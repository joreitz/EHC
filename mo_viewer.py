import streamlit as st
import plotly.graph_objects as go
import numpy as np
import re

st.set_page_config(page_title="EHT MO Viewer", layout="wide")

# --- Angepasste Parser-Funktionen für hochgeladene Dateien (Strings) ---
@st.cache_data
def parse_eht_output(content):
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
def read_cube(content):
    lines = content.splitlines()

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

st.title("☁️ Cloud EHT-Rechner & Viewer")

with st.sidebar:
    st.header("1. Molekül hochladen")
    xyz_file = st.file_uploader("Lade eine .xyz Datei hoch", type=["xyz"])

if xyz_file is not None:
    # 1. Die hochgeladene Datei speichern, damit Rust sie lesen kann
    with open("struc.xyz", "wb") as f:
        f.write(xyz_file.getvalue())
    
    st.success("Struktur gespeichert! Starte Berechnung...")
    
    # 2. Rust-Programm über das Terminal ausführen
    with st.spinner("Rust rechnet (Diagonalisiere Matrizen)... 🦀"):
        try:
            # Das entspricht dem Befehl "cargo run --release" im Terminal
            subprocess.run(["cargo", "run", "--release"], check=True)
        except subprocess.CalledProcessError as e:
            st.error("Fehler bei der Rust-Berechnung!")
            st.stop()
            
    st.success("Berechnung erfolgreich! Lade 3D-Modelle...")
    
    # 3. Jetzt, wo Rust fertig ist, existieren die eht_output.txt und die mo_X.cube Dateien
    # auf der Festplatte des Cloud-Servers. Wir können sie direkt einlesen!
    
    with open("eht_output.txt", "r") as f:
        eht_content = f.read()
        
    energies, characters, homo_idx = parse_eht_output(eht_content)

# --- UI Aufbau ---
st.title("🧪 Extended Hückel MO Viewer")

# NEU: Datei-Uploader in der Seitenleiste (Sidebar)
with st.sidebar:
    st.header("1. Daten hochladen")
    eht_file = st.file_uploader("Lade eht_output.txt hoch", type=["txt"])
    cube_files = st.file_uploader("Lade Cube-Dateien hoch (.cube)", type=["cube"], accept_multiple_files=True)
    
    st.header("2. Darstellung")
    isovalue = st.slider("Iso-Wert (Volumen)", min_value=0.005, max_value=0.1, value=0.03, step=0.005)

# Hauptbereich prüfen, ob die Textdatei da ist
if eht_file is None:
    st.info("👈 Bitte lade auf der linken Seite deine 'eht_output.txt' hoch, um zu starten!")
    st.stop()

# Textdatei auslesen (dekodieren aus Bytes in String)
eht_content = eht_file.getvalue().decode("utf-8")
energies, characters, homo_idx = parse_eht_output(eht_content)

# Cube Dateien in ein Dictionary sortieren (Dateiname -> Datei)
cube_dict = {f.name: f for f in cube_files} if cube_files else {}

col1, col2 = st.columns([1, 2])

with col1:
    st.subheader("MO Energieschema")
    mo_list = sorted(list(energies.keys()))
    default_index = homo_idx - 1 if homo_idx > 0 else 0
    selected_mo = st.selectbox("Wähle ein MO zur 3D-Ansicht:", mo_list, index=default_index)
    
    fig2d = go.Figure()
    x_center, line_width = 0, 0.4

    for mo_id, energy in energies.items():
        is_occupied = mo_id <= homo_idx
        if mo_id == selected_mo:
            color, width = "orange", 8
        else:
            color, width = ("blue" if is_occupied else "red"), 4

        char = characters.get(mo_id, "Unknown")
        hover_text = f"MO {mo_id}<br>Energie: {energy:.4f} eV<br>Charakter: {char}"

        fig2d.add_trace(go.Scatter(
            x=[x_center - line_width/2, x_center + line_width/2],
            y=[energy, energy], mode="lines",
            line=dict(color=color, width=width),
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
    expected_filename = f"mo_{selected_mo}.cube"
    
    if expected_filename in cube_dict:
        with st.spinner('Lade 3D Gitter...'):
            cube_content = cube_dict[expected_filename].getvalue().decode("utf-8")
            atoms, X, Y, Z, values = read_cube(cube_content)
            
            fig3d = go.Figure()

            # Positive & Negative Phasen
            fig3d.add_trace(go.Isosurface(x=X, y=Y, z=Z, value=values, isomin=isovalue, isomax=isovalue, surface_fill=0.7, colorscale=[[0, 'blue'], [1, 'blue']], showscale=False))
            fig3d.add_trace(go.Isosurface(x=X, y=Y, z=Z, value=values, isomin=-isovalue, isomax=-isovalue, surface_fill=0.7, colorscale=[[0, 'red'], [1, 'red']], showscale=False))

            # Atome zeichnen
            color_map = {1: 'lightgray', 6: 'black', 7: 'blue', 8: 'red'}
            fig3d.add_trace(go.Scatter3d(
                x=[a['pos'][0] for a in atoms], y=[a['pos'][1] for a in atoms], z=[a['pos'][2] for a in atoms],
                mode='markers', marker=dict(size=10, color=[color_map.get(a['z'], 'green') for a in atoms]), hoverinfo="none"
            ))

            fig3d.update_layout(scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z', aspectmode='data'), height=600, margin=dict(l=0, r=0, t=0, b=0))
            st.plotly_chart(fig3d, width="stretch")
    else:
        st.warning(f"Lade bitte die Datei `{expected_filename}` auf der linken Seite hoch, um das 3D-Modell zu sehen!")
