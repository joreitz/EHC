import streamlit as st
import plotly.graph_objects as go
import numpy as np
import re
import subprocess
import os

# --- Seitenkonfiguration ---
st.set_page_config(page_title="🧪 EHT Workbench", layout="wide")

# --- Funktionen für den Daten-Import (wie gehabt) ---
@st.cache_data
def parse_eht_output(filepath, mod_time):
    if not os.path.exists(filepath): return None, None, 0
    with open(filepath, "r", encoding="utf-8") as f:
        content = f.read()
    
    energies, characters = {}, {}
    if "Orbital Energies (eV):" in content:
        block = content.split("Orbital Energies (eV):")[1].split("===")[0].strip()
        for line in block.split("\n"):
            if "MO" in line:
                p = line.split(":")
                energies[int(p[0].replace("MO", "").strip())] = float(p[1].replace("eV", "").strip())
    
    if "=== All MO characters ===" in content:
        block = content.split("=== All MO characters ===")[1].strip()
        for line in block.split("\n"):
            if "MO" in line:
                m, r = line.split(":", 1)
                characters[int(m.replace("MO", "").strip())] = r.split("Character:")[1].strip()

    homo_idx = int(re.search(r"HOMO is MO (\d+)", content).group(1)) if re.search(r"HOMO is MO (\d+)", content) else 0
    return energies, characters, homo_idx

@st.cache_data
def read_cube(filepath, mod_time):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    natoms = int(lines[2].split()[0])
    origin = np.array([float(x) for x in lines[2].split()[1:4]])
    nx, dx = int(lines[3].split()[0]), float(lines[3].split()[1])
    ny, dy = int(lines[4].split()[0]), float(lines[4].split()[2])
    nz, dz = int(lines[5].split()[0]), float(lines[5].split()[3])
    
    atoms = []
    for i in range(natoms):
        p = lines[6+i].split()
        atoms.append({'z': int(p[0]), 'pos': [float(p[2]), float(p[3]), float(p[4])]})

    data = []
    for line in lines[6+natoms:]:
        data.extend([float(val) for val in line.split()])
    
    val_array = np.array(data)
    x_i, y_i, z_i = np.mgrid[0:nx, 0:ny, 0:nz]
    return atoms, (origin[0] + x_i * dx).flatten(), (origin[1] + y_i * dy).flatten(), (origin[2] + z_i * dz).flatten(), val_array

# --- UI Aufbau ---
st.title("🧪 EHT Quantenchemie Workbench")

# Beispiel-Geometrie für den Editor
default_xyz = """4
Formaldehyd Beispiel
C          0.00000        0.00000        0.00000
O          0.00000        0.00000        1.21000
H          0.94000        0.00000       -0.58000
H         -0.94000        0.00000       -0.58000"""

with st.sidebar:
    st.header("1. Geometrie Editor")
    xyz_input = st.text_area("XYZ Koordinaten hier bearbeiten:", value=default_xyz, height=300)
    
    st.header("2. Darstellung")
    isovalue = st.slider("Iso-Wert (Orbitalgröße)", 0.005, 0.1, 0.03, 0.005)
    
    if st.button("▶️ Berechnung starten", type="primary"):
        with st.spinner("Rust berechnet..."):
            with open("struc.xyz", "w") as f:
                f.write(xyz_input)
            
            log_placeholder = st.empty()
            log_text = ""
            process = subprocess.Popen(["cargo", "run", "--release"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            for line in process.stdout:
                log_text += line
                log_placeholder.code(log_text, language="bash")
            process.wait()
            if process.returncode == 0: st.success("Fertig!")
            else: st.error("Fehler in der Berechnung!")

# --- Visualisierung ---
if os.path.exists("eht_output.txt"):
    eht_time = os.path.getmtime("eht_output.txt")
    energies, characters, homo_idx = parse_eht_output("eht_output.txt", eht_time)

    col1, col2 = st.columns([1, 2])

    with col1:
        st.subheader("Energieschema")
        
        # --- DOWNLOAD BUTTON ---
        with open("eht_output.txt", "rb") as f:
            st.download_button(
                label="📄 eht_output.txt herunterladen",
                data=f,
                file_name="eht_output.txt",
                mime="text/plain"
            )
        
        # --- WICHTIG: DIESE ZEILE HAT GEFEHLT ---
        mo_list = sorted(energies.keys())
        selected_mo = st.selectbox("Wähle ein MO zur 3D-Ansicht:", mo_list, index=max(0, homo_idx-1))
        # ----------------------------------------

        # --- Verbessertes Energieniveaudiagramm ---
        fig2d = go.Figure()
        
        # Horizontale Nulllinie
        fig2d.add_hline(y=0, line_dash="dash", line_color="gray", annotation_text="0 eV", annotation_position="bottom right")

        x_center = 0
        line_width = 0.5 

        for mo_id, energy in energies.items():
            is_occupied = mo_id <= homo_idx
            
            # Jetzt ist selected_mo bekannt!
            if mo_id == selected_mo:
                color = "orange"
                width = 10
            else:
                color = "blue" if is_occupied else "red"
                width = 5

            char = characters.get(mo_id, "Unknown")
            hover_text = f"MO {mo_id}<br>Energie: {energy:.4f} eV<br>Charakter: {char}"

            fig2d.add_trace(go.Scatter(
                x=[x_center - line_width/2, x_center + line_width/2],
                y=[energy, energy],
                mode="lines",
                line=dict(color=color, width=width),
                name=f"MO {mo_id}",
                hoverinfo="text",
                hovertext=hover_text,
                showlegend=False
            ))

            if is_occupied:
                fig2d.add_annotation(
                    x=x_center, y=energy, text="⥮", showarrow=False, 
                    font=dict(size=20, color=color), yanchor="bottom", yshift=2
                )

        fig2d.update_layout(
            yaxis_title="Energie (eV)",
            xaxis=dict(showticklabels=False, range=[-0.6, 0.6]),
            plot_bgcolor="white",
            height=650,
            margin=dict(l=20, r=20, t=30, b=20)
        )
        fig2d.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightGray', zeroline=False)
        
        st.plotly_chart(fig2d, width="stretch")

    with col2:
        st.subheader(f"3D: MO {selected_mo}")
        cube_file = f"mo_{selected_mo}.cube"
        if os.path.exists(cube_file):
            atoms, X, Y, Z, vals = read_cube(cube_file, os.path.getmtime(cube_file))
            fig3d = go.Figure()
            # Orbitale
            fig3d.add_trace(go.Isosurface(x=X, y=Y, z=Z, value=vals, isomin=isovalue, isomax=isovalue, surface_fill=0.7, colorscale=[[0,'blue'],[1,'blue']], showscale=False))
            fig3d.add_trace(go.Isosurface(x=X, y=Y, z=Z, value=vals, isomin=-isovalue, isomax=-isovalue, surface_fill=0.7, colorscale=[[0,'red'],[1,'red']], showscale=False))
            
            # Atome & Bindungen
            ax, ay, az = [a['pos'][0] for a in atoms], [a['pos'][1] for a in atoms], [a['pos'][2] for a in atoms]
            bx, by, bz = [], [], []
            for i in range(len(atoms)):
                for j in range(i+1, len(atoms)):
                    d = np.linalg.norm(np.array(atoms[i]['pos']) - np.array(atoms[j]['pos']))
                    if 0.1 < d < 3.4: bx.extend([ax[i], ax[j], None]); by.extend([ay[i], ay[j], None]); bz.extend([az[i], az[j], None])
            
            fig3d.add_trace(go.Scatter3d(x=bx, y=by, z=bz, mode='lines', line=dict(color='gray', width=5), showlegend=False))
            c_map = {1: '#d3d3d3', 6: '#333333', 7: 'blue', 8: 'red', 9: 'green'}
            fig3d.add_trace(go.Scatter3d(x=ax, y=ay, z=az, mode='markers', marker=dict(size=10, color=[c_map.get(a['z'], 'magenta') for a in atoms], line=dict(color='black', width=2)), showlegend=False))
            
            fig3d.update_layout(scene=dict(aspectmode='data'), height=600, margin=dict(l=0,r=0,b=0,t=0))
            st.plotly_chart(fig3d, width="stretch")
