import streamlit as st
import plotly.graph_objects as go
import numpy as np
import re
import subprocess
import os

# --- Seitenkonfiguration ---
st.set_page_config(page_title=" EHT Cloud-Rechner", layout="wide")

# --- Datenlade-Funktionen ---
# Die mod_time sorgt dafür, dass Streamlit den Cache clever leert, 
# sobald Rust die Dateien neu berechnet hat!
@st.cache_data
def parse_eht_output(filepath, mod_time):
    with open(filepath, "r", encoding="utf-8") as f:
        content = f.read()

    energies = {}
    if "Orbital Energies (eV):" in content:
        energy_block = content.split("Orbital Energies (eV):")[1].split("===")[0].strip()
        for line in energy_block.split("\n"):
            if "MO" in line:
                parts = line.split(":")
                mo_id = int(parts[0].replace("MO", "").strip())
                energy = float(parts[1].replace("eV", "").strip())
                energies[mo_id] = energy

    characters = {}
    if "=== All MO characters ===" in content:
        char_block = content.split("=== All MO characters ===")[1].strip()
        for line in char_block.split("\n"):
            if "MO" in line:
                mo_str, rest = line.split(":", 1)
                mo_id = int(mo_str.replace("MO", "").strip())
                char = rest.split("Character:")[1].strip()
                characters[mo_id] = char

    homo_match = re.search(r"HOMO is MO (\d+)", content)
    homo_idx = int(homo_match.group(1)) if homo_match else 0

    return energies, characters, homo_idx

@st.cache_data
def read_cube(filepath, mod_time):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    parts = lines[2].split()
    natoms = int(parts[0])
    origin = np.array([float(parts[1]), float(parts[2]), float(parts[3])])

    nx, dx, _, _ = lines[3].split()
    ny, _, dy, _ = lines[4].split()
    nz, _, _, dz = lines[5].split()
    nx, ny, nz = int(nx), int(ny), int(nz)
    dx, dy, dz = float(dx), float(dy), float(dz)

    atoms = []
    atom_start = 6
    for i in range(natoms):
        parts = lines[atom_start + i].split()
        atoms.append({
            'z': int(parts[0]),
            'pos': [float(parts[2]), float(parts[3]), float(parts[4])]
        })

    data_start = atom_start + natoms
    data = []
    for line in lines[data_start:]:
        data.extend([float(val) for val in line.split()])
    val_array = np.array(data)

    x_idx, y_idx, z_idx = np.mgrid[0:nx, 0:ny, 0:nz]
    X = origin[0] + x_idx * dx
    Y = origin[1] + y_idx * dy
    Z = origin[2] + z_idx * dz

    return atoms, X.flatten(), Y.flatten(), Z.flatten(), val_array

# --- UI Aufbau ---
st.title("☁️ Extended Hückel Cloud-Rechner & Viewer")

with st.sidebar:
    st.header("1. Molekül hochladen")
    xyz_file = st.file_uploader("Lade eine .xyz Struktur hoch", type=["xyz"])
    
    st.header("2. Darstellung")
    isovalue = st.slider("Iso-Wert (Volumen der Orbitalwolke)", min_value=0.005, max_value=0.1, value=0.03, step=0.005)
    
    st.markdown("---")
    st.markdown(" **Backend:** Rust\n\n **Frontend:** Python/Streamlit")

# --- Rust Ausführung ---
if xyz_file is not None:
    # Button zum Starten der Berechnung anzeigen
    if st.button("▶️ Quantenchemische Berechnung starten", type="primary"):
        # 1. Die hochgeladene Datei speichern
        with open("struc.xyz", "wb") as f:
            f.write(xyz_file.getvalue())
        
        st.info("Starte Rust-Backend... 🦀")
        
        # Ein leeres Feld für unser Live-Terminal in Streamlit erstellen
        log_placeholder = st.empty()
        log_text = ""

        # 2. Rust-Programm ausführen und den Output streamen
        try:
            # Popen erlaubt es uns, den Output live zu lesen
            process = subprocess.Popen(
                ["cargo", "run", "--release"],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, # Wir werfen Fehlermeldungen und normalen Output zusammen
                text=True,
                bufsize=1, # Zeilenweise puffern
                universal_newlines=True
            )

            # Zeile für Zeile auslesen und im Dashboard aktualisieren
            for line in process.stdout:
                log_text += line
                # Das log_text in einem schönen Code-Block anzeigen (Scrollbar entsteht automatisch)
                log_placeholder.code(log_text, language="bash")
            
            # Warten, bis der Prozess wirklich beendet ist
            process.wait()

            # Prüfen, ob das Programm fehlerfrei durchgelaufen ist (Exit-Code 0)
            if process.returncode == 0:
                st.success("Berechnung erfolgreich abgeschlossen!")
            else:
                st.error(f"Fehler bei der Rust-Berechnung! (Exit code: {process.returncode})")
                st.stop()

        except Exception as e:
            st.error(f"Fehler beim Starten des Backend-Prozesses: {e}")
            st.stop()
# --- Ergebnisse Anzeigen ---
# Wir prüfen, ob die Rust-Berechnung eine output-Datei hinterlassen hat
if os.path.exists("eht_output.txt"):
    st.markdown("---")
    
    # Modifikationszeit auslesen, damit Streamlit bei neuen Molekülen den Cache leert
    eht_mod_time = os.path.getmtime("eht_output.txt")
    energies, characters, homo_idx = parse_eht_output("eht_output.txt", eht_mod_time)

    if not energies:
        st.error("Konnte eht_output.txt nicht parsen. Bitte überprüfe das Format.")
        st.stop()

    col1, col2 = st.columns([1, 2])

    with col1:
        st.subheader("MO Energieschema")
        mo_list = sorted(list(energies.keys()))
        default_index = homo_idx - 1 if homo_idx > 0 else 0
        selected_mo = st.selectbox("Wähle ein MO zur 3D-Ansicht:", mo_list, index=default_index)
        
        fig2d = go.Figure()
        x_center = 0
        line_width = 0.4

        for mo_id, energy in energies.items():
            is_occupied = mo_id <= homo_idx
            
            if mo_id == selected_mo:
                color = "orange"
                width = 8
            else:
                color = "blue" if is_occupied else "red"
                width = 4

            char = characters.get(mo_id, "Unknown")
            hover_text = f"MO {mo_id}<br>Energie: {energy:.4f} eV<br>Charakter: {char}"

            fig2d.add_trace(go.Scatter(
                x=[x_center - line_width/2, x_center + line_width/2],
                y=[energy, energy], mode="lines", line=dict(color=color, width=width),
                name=f"MO {mo_id}", hoverinfo="text", hovertext=hover_text, showlegend=False
            ))

            if is_occupied:
                fig2d.add_annotation(
                    x=x_center, y=energy, text="⥮", showarrow=False, 
                    font=dict(size=16, color=color), yanchor="bottom"
                )

        fig2d.update_layout(
            yaxis_title="Energie (eV)", xaxis=dict(showticklabels=False, range=[-1, 1]),
            plot_bgcolor="white", height=600, margin=dict(l=0, r=0, t=30, b=0)
        )
        fig2d.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGray')
        
        st.plotly_chart(fig2d, width="stretch")

    with col2:
        st.subheader(f"3D Ansicht: MO {selected_mo}")
        cube_file = f"mo_{selected_mo}.cube"
        
        if os.path.exists(cube_file):
            with st.spinner('Zeichne 3D Gitter (dies kann kurz dauern)...'):
                cube_mod_time = os.path.getmtime(cube_file)
                atoms, X, Y, Z, values = read_cube(cube_file, cube_mod_time)
                
                fig3d = go.Figure()

                # Positive Phase
                fig3d.add_trace(go.Isosurface(
                    x=X, y=Y, z=Z, value=values, isomin=isovalue, isomax=isovalue,
                    surface_fill=0.7, colorscale=[[0, 'blue'], [1, 'blue']], showscale=False, name='Positive Phase'
                ))

                # Negative Phase
                fig3d.add_trace(go.Isosurface(
                    x=X, y=Y, z=Z, value=values, isomin=-isovalue, isomax=-isovalue,
                    surface_fill=0.7, colorscale=[[0, 'red'], [1, 'red']], showscale=False, name='Negative Phase'
                ))

                # Atome zeichnen
                color_map = {1: 'lightgray', 6: 'black', 7: 'blue', 8: 'red'}
                fig3d.add_trace(go.Scatter3d(
                    x=[a['pos'][0] for a in atoms], y=[a['pos'][1] for a in atoms], z=[a['pos'][2] for a in atoms],
                    mode='markers', marker=dict(size=10, color=[color_map.get(a['z'], 'green') for a in atoms]),
                    name='Atome', hoverinfo="none"
                ))

                fig3d.update_layout(
                    scene=dict(xaxis_title='X (Bohr)', yaxis_title='Y (Bohr)', zaxis_title='Z (Bohr)', aspectmode='data'),
                    height=600, margin=dict(l=0, r=0, t=0, b=0)
                )
                
                st.plotly_chart(fig3d, width="stretch")
        else:
            st.warning(f"Konnte {cube_file} nicht finden. Ist bei der Berechnung etwas schiefgegangen?")
else:
    st.info(" Bitte lade eine .xyz-Datei hoch und starte die Berechnung, um die Ergebnisse zu sehen!")
