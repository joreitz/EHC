use nalgebra as na;
use std::collections::HashMap;
use std::fs::File;
use std::io::{Write, Result};

mod pop_analysis;
mod parser;
// use parser::{AtomData, ElementBasis, load_basis_library};

// Var's
static _BOHR_TO_ANGSTROM: f64 = 0.529177210903; // Bohr radius in Angstroms
static _ANGSTROM_TO_BOHR: f64 = 1.0 / _BOHR_TO_ANGSTROM; // Angstrom to Bohr conversion factor
static _HARTREE_TO_EV: f64 = 27.211386245988; // Hartree energy in eV
static _EV_TO_HARTREE: f64 = 1.0 / _HARTREE_TO_EV; // eV to Hartree conversion factor
static K_WH: f64 = 1.75;
#[derive(Debug)]
struct OrbitalData {
    n: u8,
    l: i8,
    _coefficient: f64,
    exponent: f64
}

#[derive(Debug)]
struct Atom {
    _label: String,
    position: [f64; 3],
    _z: u8,
    valence_orbitals: Vec<OrbitalData>,
    vsie: HashMap<(u8, u8), f64>,
}

#[derive(Debug)]
struct _Orbital {
    atom_id: usize,
    n: u8,
    l: i8,
    _m: i8,
}

#[derive(Debug)]
struct Molecule {
    atoms: Vec<Atom>,
    orbitals: Vec<_Orbital>,
    dist_matrix: Option<na::DMatrix<f64>>,
    s_mat: Option<na::DMatrix<f64>>,
    h_mat: Option<na::DMatrix<f64>>,
    energies: Option<na::DVector<f64>>,
    c_matrix: Option<na::DMatrix<f64>>,
    s_sqrt: Option<na::DMatrix<f64>>,
}

impl Atom {
    
    fn add_vsie(&mut self, n: u8, l: i8, vsie: f64) {
        self.vsie.insert((n, l.try_into().unwrap()), vsie);
    }
    
    fn get_vsie(&self, n: u8, l: i8) -> Option<f64> {
        self.vsie.get(&(n, l.try_into().unwrap())).copied()
    }

    fn new(label: &str, position: [f64; 3], z: u8, vsie: HashMap<(u8, u8), f64>) -> Self {
        Atom {
            _label: label.to_string(),
            position,
            _z: z,
            valence_orbitals: Vec::new(),
            vsie,
        }
    }

    fn add_orbital(&mut self, n: u8, l: u8, exponent: f64, coefficient: f64) {
        self.valence_orbitals.push(OrbitalData {n, l: l.try_into().unwrap(), _coefficient: coefficient, exponent});
    }

    fn find_orbital(&self, n: u8, l: u8) -> Option<&OrbitalData> {
        self.valence_orbitals.iter().find(|orbital| orbital.n == n && orbital.l == l as i8)
    }

}

impl _Orbital {
    fn _new(atom_id: usize, n: u8, l: i8, m: i8) -> Self {
        _Orbital { atom_id, n, l, _m: m }
    }

    fn _orbital_type(&self) -> String {
        match [self.l, self._m] {
            [0, 0] => "s".to_string(),
            [1, -1] => "px".to_string(),
            [1, 0] => "py".to_string(),
            [1, 1] => "pz".to_string(),
            _ => "unknown".to_string()
        }
    }

    fn get_axis(&self) -> [f64; 3] {
        match (self.l, self._m) {
            (1, -1) => [0.0, 1.0, 0.0],  // py
            (1,  0) => [0.0, 0.0, 1.0],  // pz
            (1,  1) => [1.0, 0.0, 0.0],  // px
            _ => [0.0, 0.0, 0.0],        // s hat keine Richtung
        }
    }
}

impl Molecule {

    fn calculate_distance(&self, atom_i: usize, atom_j: usize) -> f64 {
        let a1 = &self.atoms[atom_i].position;
        let a2 = &self.atoms[atom_j].position;
        let dx = a2[0] - a1[0];
        let dy = a2[1] - a1[1];
        let dz = a2[2] - a1[2];
        (dx*dx + dy*dy + dz*dz).sqrt()
    }

    fn get_exponent(&self, orbital: &_Orbital) -> f64 {
    let atom = &self.atoms[orbital.atom_id];
    atom.find_orbital(orbital.n, orbital.l.try_into().unwrap())
        .expect(&format!("Exponent not found for n={}, l={}", orbital.n, orbital.l))
        .exponent
    }

    fn get_cos_theta(&self, orbital: &_Orbital, direction: [f64; 3]) -> f64 {
        let axis = orbital.get_axis();
        axis[0]*direction[0] + axis[1]*direction[1] + axis[2]*direction[2]
    }
    
    fn direction_vector(&self, atom_i: usize, atom_j: usize) -> [f64; 3] {
        let a1 = &self.atoms[atom_i].position;
        let a2 = &self.atoms[atom_j].position;
        
        let dx = a2[0] - a1[0];
        let dy = a2[1] - a1[1];
        let dz = a2[2] - a1[2];
        
        let dist = (dx*dx + dy*dy + dz*dz).sqrt();
        [dx/dist, dy/dist, dz/dist]
    }

    fn overlap_ss(&self, n1: u8, zeta1: f64, n2: u8, zeta2: f64, r: f64) -> f64 {
        let p = (zeta1 + zeta2) / 2.0;
        let x = p * r;
        let exp_term = (-x).exp();
    
        match (n1, n2) {
            (1, 1) => exp_term * (1.0 + x + x.powi(2)/3.0),
        
            (1, 2) | (2, 1) => exp_term * (1.0 + x + 2.0*x.powi(2)/5.0 + x.powi(3)/15.0),
        
            (2, 2) => exp_term * (1.0 + x + 2.0*x.powi(2)/5.0 + 2.0*x.powi(3)/15.0 + x.powi(4)/75.0),
        
            (1, 3) | (3, 1) => exp_term * (1.0 + x + 3.0*x.powi(2)/7.0 + 2.0*x.powi(3)/21.0 + x.powi(4)/105.0),
        
            (2, 3) | (3, 2) => exp_term * (1.0 + x + 3.0*x.powi(2)/7.0 + 4.0*x.powi(3)/35.0 + 2.0*x.powi(4)/105.0 + x.powi(5)/525.0),
        
            (3, 3) => exp_term * (1.0 + x + 3.0*x.powi(2)/7.0 + 4.0*x.powi(3)/35.0 + 3.0*x.powi(4)/175.0 + 2.0*x.powi(5)/525.0 + x.powi(6)/3675.0),
        
            _ => 0.0, // Höhere n nicht implementiert
        }
    }

    fn overlap_sp(&self, n_s: u8, zeta_s: f64, n_p: u8, zeta_p: f64, r: f64, cos_theta: f64) -> f64 {
        let p = (zeta_s + zeta_p) / 2.0;
        let x = p * r;
        let exp_term = (-x).exp();
    
        let radial = match (n_s, n_p) {
            (1, 2) => x * exp_term * (1.0 + x/2.0 + x.powi(2)/10.0),
        
            (2, 2) => x * exp_term * (1.0 + x/2.0 + x.powi(2)/10.0 + x.powi(3)/60.0),
        
            (1, 3) => x * exp_term * (1.0 + x/2.0 + 3.0*x.powi(2)/28.0 + x.powi(3)/84.0),
        
            (2, 3) => x * exp_term * (1.0 + x/2.0 + 3.0*x.powi(2)/28.0 + x.powi(3)/84.0 + x.powi(4)/504.0),
        
            (3, 3) => x * exp_term * (1.0 + x/2.0 + 3.0*x.powi(2)/28.0 + 2.0*x.powi(3)/168.0 + x.powi(4)/672.0 + x.powi(5)/5040.0),
        
            // Symmetrische Fälle (p-s)
            (2, 1) => x * exp_term * (1.0 + x/2.0 + x.powi(2)/10.0),
            (3, 1) => x * exp_term * (1.0 + x/2.0 + 3.0*x.powi(2)/28.0 + x.powi(3)/84.0),
            (3, 2) => x * exp_term * (1.0 + x/2.0 + 3.0*x.powi(2)/28.0 + x.powi(3)/84.0 + x.powi(4)/504.0),
        
            _ => 0.0,
    };
    
    radial * cos_theta
    }

    fn overlap_pp(&self, atom_i: usize, atom_j: usize, s_sigma: f64, s_pi: f64) -> na::Matrix3<f64> {
        let p1 = na::Vector3::from_row_slice(&self.atoms[atom_i].position);
        let p2 = na::Vector3::from_row_slice(&self.atoms[atom_j].position);
        
        let mut z_local = p2 - p1;
        let distance = z_local.norm();
        
        if distance < 1e-10 {
            // Gleiches Atom, Überlappung ist die Einheitsmatrix
            return na::Matrix3::identity();
        }
        
        z_local /= distance; // Normieren

        // Hilfsvektor finden (meistens globale Z-Achse)
        let mut global_up = na::Vector3::new(0.0, 0.0, 1.0);
        
        // Wenn z_local fast parallel zur Z-Achse ist, nimm stattdessen X-Achse
        if z_local.z.abs() > 0.9999 {
            global_up = na::Vector3::new(1.0, 0.0, 0.0);
        }

        // y_local = z_local x global_up
        let mut y_local = z_local.cross(&global_up);
        y_local.normalize_mut();

        // x_local = y_local x z_local
        let mut x_local = y_local.cross(&z_local);
        x_local.normalize_mut();

        // Rotationsmatrix R aus den Spaltenvektoren aufbauen
        let r = na::Matrix3::from_columns(&[x_local, y_local, z_local]);

        // Lokale Überlappungsmatrix definieren
        // Beachte das Vorzeichen bei s_sigma!
        let mut s_loc = na::Matrix3::zeros();
        s_loc[(0, 0)] = s_pi;
        s_loc[(1, 1)] = s_pi;
        s_loc[(2, 2)] = -s_sigma;

        // Rücktransformation: S_global = R * S_loc * R^T
        r * s_loc * r.transpose()
    }

    fn radial_overlap_pp(&self, n1: u8, zeta1: f64, n2: u8, zeta2: f64, r: f64) -> (f64, f64) {
        let p = (zeta1 + zeta2) / 2.0;
        let x = p * r;
        let exp_term = (-x).exp();
        
        match (n1, n2) {
            (2, 2) => {
                let sigma = exp_term * (-1.0 - x - 0.2 * x.powi(2) + (2.0/15.0) * x.powi(3) + (1.0/15.0) * x.powi(4));
            
                // S_pi (Parallele Achsen -> S(0) = 1)
                let pi = exp_term * (1.0 + x + 0.4 * x.powi(2) + (1.0/15.0) * x.powi(3));
                
                (sigma, pi)
            },
            (2, 3) | (3, 2) => {
                let sigma = x.powi(2) * exp_term * (1.0 + x/3.0 + x.powi(2)/21.0 + x.powi(3)/252.0);
                let pi = x.powi(2) * exp_term * (1.0 - x/5.0 - x.powi(2)/42.0 - x.powi(3)/504.0);
                (sigma, pi)
            },
            (3, 3) => {
                let sigma = x.powi(2) * exp_term * (1.0 + x/3.0 + x.powi(2)/21.0 + x.powi(3)/168.0 + x.powi(4)/2520.0);
                let pi = x.powi(2) * exp_term * (1.0 - x/5.0 - x.powi(2)/42.0 - x.powi(3)/420.0 - x.powi(4)/5040.0);
                (sigma, pi)
            },
            _ => (0.0, 0.0),
        }
    }

    fn init(&mut self) {
        self.build_basis_set();
        self.dist_matrix = Some(self.init_dist_matrix());
        self.s_mat = Some(self.build_overlap_matrix());
        self.h_mat = Some(self.build_hamiltonian());
        self.solve();
    }
    
    fn add_atom(&mut self, atom: Atom) {
        self.atoms.push(atom);
    }

    fn build_basis_set(&mut self) {
        for (atom_id, atom) in self.atoms.iter().enumerate() {
            for orbital_data in &atom.valence_orbitals {
                for m in -orbital_data.l..=orbital_data.l {
                    self.orbitals.push(_Orbital {
                        atom_id: atom_id,
                        n: orbital_data.n,
                        l: orbital_data.l, 
                        _m: m as i8,
                    });
                }
            }
        }
    }

    fn init_dist_matrix(&self) -> na::DMatrix<f64> {
        let n= self.atoms.len();
        let mut dist_matrix = na::DMatrix::<f64>::zeros(n,n);
        
        for i in 0..n {
            for j in i..n {
                if i == j {
                    dist_matrix[(i,j)] = 0.0;
                } else {
                    let d = self.calculate_distance(i, j);
                    dist_matrix[(i,j)] = d;
                    dist_matrix[(j,i)] = d;
                }
            }
        }

        dist_matrix
    }

    fn calculate_overlap(&self, i: usize, j: usize) -> f64 {
    if i == j { return 1.0; }
    
    let o1 = &self.orbitals[i];
    let o2 = &self.orbitals[j];
    
    let r = self.dist_matrix.as_ref().unwrap()[(o1.atom_id, o2.atom_id)];
    if r < 1e-10 { return 0.0; } // Gleicher Ort
    
    let exp_i = self.get_exponent(o1);
    let exp_j = self.get_exponent(o2);
    
    let direction = self.direction_vector(o1.atom_id, o2.atom_id);
    
    match (o1.l, o2.l) {
        (0, 0) => self.overlap_ss(o1.n, exp_i, o2.n, exp_j, r),
        
        (0, 1) => {
            let cos_theta = self.get_cos_theta(o2, direction);
            -self.overlap_sp(o1.n, exp_i, o2.n, exp_j, r, cos_theta)
        },
        
        (1, 0) => {
            let cos_theta = self.get_cos_theta(o1, direction);
            self.overlap_sp(o2.n, exp_j, o1.n, exp_i, r, cos_theta)  // Vorzeichen!
        },
        
        (1, 1) => {
            let (s_sigma, s_pi) = self.radial_overlap_pp(o1.n, exp_i, o2.n, exp_j, r);
            
            // 2. Erzeuge den rotierten 3x3 Überlappungs-Block
            let p_block = self.overlap_pp(o1.atom_id, o2.atom_id, s_sigma, s_pi);
            
            // 3. Welches m entspricht welchem Matrix-Index?
            // Nach deiner Logik in `get_axis`: px (1) -> x, py (-1) -> y, pz (0) -> z
            let m_to_index = |m: i8| -> usize {
                match m {
                    1 => 0,  // px ist Index 0
                    -1 => 1, // py ist Index 1
                    0 => 2,  // pz ist Index 2
                    _ => panic!("Ungültiges m für p-Orbital"),
                }
            };
            
            let idx_1 = m_to_index(o1._m);
            let idx_2 = m_to_index(o2._m);
            
            // 4. Lies einfach den passenden Wert aus der Matrix ab!
            p_block[(idx_1, idx_2)]
        },
        
        _ => 0.0,
    }
    }

    fn build_overlap_matrix(&mut self) -> na::DMatrix<f64> {
        let n= self.orbitals.len();
        let mut s = na::DMatrix::<f64>::zeros(n,n); 
        
        for i in 0..n {
            for j in i..n {
                    s[(i,j)] = self.calculate_overlap(i, j);
                    s[(j, i)] = s[(i,j)] 
            }
        }
        println!("Checking S matrix diagonal:");
        for i in 0..s.nrows() {
            println!("S[{},{}] = {}", i, i, s[(i,i)]);
            assert!((s[(i,i)] - 1.0).abs() < 1e-6, "Diagonal not 1!");
        }

        println!("Checking S matrix off-diagonal:");
        for i in 0..s.nrows() {
            for j in 0..s.ncols() {
                if s[(i,j)].abs() > 1.0 {
                println!("ERROR: S[{},{}] = {} > 1.0", i, j, s[(i,j)]);
                }
            }
        }
        s
    }
    
    fn build_hamiltonian(&self) -> na::DMatrix<f64> {
        let n = self.orbitals.len();
        let mut hamiltonian = na::DMatrix::<f64>::zeros(n, n);

        for i in 0..n {
            for j in i..n {
                let o1 = &self.orbitals[i];
                let o2 = &self.orbitals[j];

                if i==j {
                    let atom = &self.atoms[o1.atom_id];
                hamiltonian[(i,j)] = atom.get_vsie(o1.n, o1.l).expect("VSIE NOT FOUND");
                } else {
                    let (a1, a2) = (&self.atoms[o1.atom_id], &self.atoms[o2.atom_id]);
                    let s: f64 = self.s_mat.as_ref().expect("Overlap Matrix has not been initialized.")[(i, j)];
                    hamiltonian[(i,j)] = K_WH *  s * (a1.get_vsie(o1.n, o1.l).expect("VSIE NOT FOUND")+a2.get_vsie(o2.n, o2.l).expect("VSIE NOT FOUND"))/2.0;
                    hamiltonian[(j,i)] = hamiltonian[(i,j)];
                }
            }
        }
    hamiltonian
    }

    fn solve(&mut self) {
        let h = self.h_mat.as_ref().unwrap();
        let s = self.s_mat.as_ref().unwrap();
        
        println!("Diagonalizing S matrix...");
        let s_eigen = s.clone().symmetric_eigen();
        let s_vals = &s_eigen.eigenvalues;
        let u = &s_eigen.eigenvectors;
        
        let n = s_vals.len();
        let mut s_inv_sqrt = na::DMatrix::<f64>::zeros(n, n);
        let mut s_sqrt_diag = na::DMatrix::<f64>::zeros(n, n); 
        
        for i in 0..n {
            if s_vals[i] > 1e-10 {
                s_inv_sqrt[(i,i)] = 1.0 / s_vals[i].sqrt();
                s_sqrt_diag[(i,i)] = s_vals[i].sqrt(); 
            } else {
                println!("Warning: Small/negative eigenvalue λ[{}] = {}", i, s_vals[i]);
            }
        }
        
        let x = u * s_inv_sqrt * u.transpose();
        let s_sqrt = u * s_sqrt_diag * u.transpose(); 
        
        println!("Transforming Hamiltonian...");
        let h_prime = x.transpose() * h * &x;
        
        println!("Solving eigenvalue problem...");
        let eigen = h_prime.symmetric_eigen();
        let c = x * &eigen.eigenvectors;
        
        // Ergebnisse im Struct speichern!
        self.energies = Some(eigen.eigenvalues);
        self.c_matrix = Some(c);
        self.s_sqrt = Some(s_sqrt);
    }

    fn dipole_moment(&self, q: &[f64]) -> na::Vector3<f64> {
        let mut mu = na::Vector3::zeros();
        for (i, atom) in self.atoms.iter().enumerate() {
            mu += q[i] * na::Vector3::from(atom.position); // position in Bohr, q in units of e
        }
        mu // in e·Bohr
    }
    
    fn get_mo_character(&self, mo_idx: usize) -> String {
        let c = self.c_matrix.as_ref().expect("Run solve() first!");
        let s_sqrt = self.s_sqrt.as_ref().expect("Run solve() first!");
        
        // Transform the specific MO column into the orthogonalized basis: C_perp = S^(1/2) * C
        let c_col = c.column(mo_idx);
        let c_perp = s_sqrt * c_col;

        let mut atom_pops = vec![0.0; self.atoms.len()];
        let mut p_x_pop = 0.0;
        let mut p_y_pop = 0.0;
        let mut p_z_pop = 0.0;
        let mut s_pop = 0.0;

        // Sum up the squared coefficients (which now represent exact population percentages)
        for mu in 0..self.orbitals.len() {
            let pop = c_perp[mu].powi(2);
            let orb = &self.orbitals[mu];
            
            atom_pops[orb.atom_id] += pop;
            
            if orb.l == 0 {
                s_pop += pop;
            } else if orb.l == 1 {
                match orb._m {
                    1 => p_x_pop += pop,   // px
                    -1 => p_y_pop += pop,  // py
                    0 => p_z_pop += pop,   // pz
                    _ => {}
                }
            }
        }

        // 1. Check for n (non-bonding) character
        // If more than 75% of the electron density is localized on a single atom
        let max_atom_pop = atom_pops.iter().copied().fold(f64::NAN, f64::max);
        if max_atom_pop > 0.85 {
            return "n (non-bonding)".to_string();
        }

        // 2. Check for pi character
        // If the orbital is overwhelmingly made of p-orbitals pointing in ONE direction,
        // with almost zero s-orbital mixing, it is a pi orbital. 
        // (Assuming you orient planar molecules flat along Cartesian axes).
        if p_x_pop > 0.75 && s_pop < 0.05 {
            return "pi (px)".to_string();
        } else if p_y_pop > 0.75 && s_pop < 0.05 {
            return "pi (py)".to_string();
        } else if p_z_pop > 0.75 && s_pop < 0.05 {
            return "pi (pz)".to_string();
        }

        // 3. Otherwise, it is a sigma orbital (hybridized s and p character)
        "sigma".to_string()
    }

    fn lowdin_populations(&self) {
        println!("Performing Löwdin population analysis...");
        
        let c: &nalgebra::Matrix<f64, nalgebra::Dyn, nalgebra::Dyn, nalgebra::VecStorage<f64, nalgebra::Dyn, nalgebra::Dyn>> = self.c_matrix.as_ref().expect("Run solve() first!");
        let e = self.energies.as_ref().expect("Run solve() first!");
        let s_sqrt = self.s_sqrt.as_ref().expect("Run solve() first!");
        let n_atoms = self.atoms.len();
        let n_ao = self.orbitals.len();
        let mut sorted_indices: Vec<usize> = (0..n_ao).collect();
        sorted_indices.sort_by(|&a, &b| e[a].partial_cmp(&e[b]).unwrap());
        let mut total_valectrons = 0;
        
        let mut atom_valance_electron: Vec<f64> = vec![0.0; self.atoms.len()];
        for (i, atom) in self.atoms.iter().enumerate() {
            let val_e: usize = match atom._label.as_str() {
                "H" => 1,
                "Li" => 1,
                "C" => 4,
                "Si" => 4,
                "N" => 5,
                "P" => 5,
                "O" => 6,
                "S" => 6,
                "F" => 7,
                "Cl" => 7,
                _ => panic!("Unknown element: {}", atom._label),
            };
            atom_valance_electron[i] = val_e as f64;
            total_valectrons += val_e;
        }

        let n_occ = total_valectrons / 2;
        let p_matrix = pop_analysis::orth_p_matrix(c, s_sqrt, n_occ, &sorted_indices);
        
        println!("Total valence electrons: {}, Occupied MOs: {}", total_valectrons, n_occ);
        println!("HOMO is MO {}", n_occ);
        
        let mut gap = vec![0.0; n_atoms];
        let mut q = vec![0.0; n_atoms];

        for mu in 0..n_ao {
            let atom_id = self.orbitals[mu].atom_id;
            let orbital_pop = p_matrix[(mu, mu)];
            gap[atom_id] += orbital_pop;
        }

        for atom_id in 0..n_atoms {
            q[atom_id] += atom_valance_electron[atom_id] - gap[atom_id];
        }

        println!("\nGross Atomic Populations:");
        for (i, pop) in gap.iter().enumerate() {
            println!("Atom {} ({}): {:.4} Elektronen", i, self.atoms[i]._label, pop);
        }
        println!("Partial Charges:");
        for (i, charge) in q.iter().enumerate() {
            println!("Atom {} ({}): {:.4}", i, self.atoms[i]._label, charge);
        }

        let mut wbi_matrix = na::DMatrix::<f64>::zeros(n_atoms, n_atoms);
        for mu in 0..n_ao {
            for nu in 0..n_ao {
                let a1 = self.orbitals[mu].atom_id;
                let a2 = self.orbitals[nu].atom_id;
                let p_val = p_matrix[(mu, nu)];
                wbi_matrix[(a1, a2)] += p_val * p_val;
            }
        }
        println!("\n=== Wiberg Bond Indices (Bond Orders) ===");
        println!("{:<8} | {:<8} | {:<12}", "Atom A", "Atom B", "Bond Index");
        println!("------------------------------------");
        for i in 0..n_atoms {
            for j in i+1..n_atoms {
                let bo = wbi_matrix[(i, j)];
                if bo > 0.05 {
                    println!("{:<8} | {:<8} | {:<12.4}", self.atoms[i]._label, self.atoms[j]._label, bo);
                }
            }
        }

        let mu = self.dipole_moment(&q);
        let mu_debye = mu * 2.5418;
        let mu_norm = mu_debye.norm();

        println!("\n=== Dipole moment ===");
        println!("μ = ({:.4}, {:.4}, {:.4}) Debye", mu_debye.x, mu_debye.y, mu_debye.z);
        println!("|μ| = {:.4} Debye", mu_norm);

        // --- Calculate All Allowed Transitions ---
        println!("\n=== Symmetry Allowed Transitions (TDM > 0) ===");
        println!("{:<28} | {:<28} | {:<12} | {:<30} | {:<10}", 
                 "Initial MO", "Final MO", "Delta E (eV)", "TDM Vector (x, y, z) [Debye]", "|TDM| [Debye]");
        println!("{:-<115}", "");

        for i in 0..n_occ {
            for j in n_occ..n_ao {
                let mo_i = sorted_indices[i];
                let mo_j = sorted_indices[j];
                
                let tdm = self.transition_dipole_moment(mo_i, mo_j);
                let tdm_debye = tdm * 2.5418; 
                let tdm_norm = tdm_debye.norm();
                
                if tdm_norm > 1e-4 {
                    let energy_gap = e[mo_j] - e[mo_i];
                    
                    // Fetch the character of the orbitals
                    let char_i = self.get_mo_character(mo_i);
                    let char_j = self.get_mo_character(mo_j);
                    
                    let label_i = format!("MO {:<2} ({})", i + 1, char_i);
                    let label_j = format!("MO {:<2} ({})", j + 1, char_j);

                    println!(
                        "{:<28} | {:<28} | {:<12.4} | ({:>6.3}, {:>6.3}, {:>6.3})             | {:.4}",
                        label_i, 
                        label_j, 
                        energy_gap, 
                        tdm_debye.x, tdm_debye.y, tdm_debye.z, 
                        tdm_norm
                    );
                }
            }
        }
        println!("{:-<115}", "");
    }

    fn transition_dipole_moment(&self, mo_i: usize, mo_j: usize) -> na::Vector3<f64> {
        let c = self.c_matrix.as_ref().expect("Run solve() first!");
        let s = self.s_mat.as_ref().expect("Run solve() first!");
        let n_ao = self.orbitals.len();
        
        let mut tdm = na::Vector3::zeros();
        
        for mu in 0..n_ao {
            let atom_mu = self.orbitals[mu].atom_id;
            let pos_mu = na::Vector3::from(self.atoms[atom_mu].position);
            
            for nu in 0..n_ao {
                let atom_nu = self.orbitals[nu].atom_id;
                let pos_nu = na::Vector3::from(self.atoms[atom_nu].position);
                
                // Mulliken approximation: <mu | r | nu> ≈ S_{μν} * (R_μ + R_ν) / 2
                let r_op = (pos_mu + pos_nu) / 2.0;
                
                let s_val = s[(mu, nu)];
                let c_mu_i = c[(mu, mo_i)];
                let c_nu_j = c[(nu, mo_j)];
                
                tdm += c_mu_i * c_nu_j * s_val * r_op;
            }
        }
        
        tdm 
    }   
}

fn main() {
    let basis_lib = parser::load_basis_library("basis_library.json")
        .expect("Failed to load basis");
    let atom_data = parser::read_xyz("struc.xyz")
        .expect("Failed to read XYZ");

    let mut molecule = Molecule {
        atoms: Vec::new(),
        orbitals: Vec::new(),
        dist_matrix: None,
        s_mat: None,
        h_mat: None,
        energies: None,
        c_matrix: None,
        s_sqrt: None,
    };
    
    for data in atom_data {
        let basis = basis_lib.get(&data.element)
            .expect(&format!("Element {} not found", data.element));
        
        let mut atom = Atom::new(
            &data.element,
            data.position,
            basis.atomic_number,
            HashMap::new()
        );
        
        // Orbitale hinzufügen
        for orb in &basis.orbitals {
            atom.add_orbital(orb.n, orb.l, orb.exponent, orb.coefficient);
        }
        
        // VSIE hinzufügen
        for vsie_data in &basis.vsie {
            atom.add_vsie(vsie_data[0] as u8, vsie_data[1] as i8, vsie_data[2]);
        }
        
        molecule.add_atom(atom);
    }

    molecule.init();

    fn write_output(molecule: &Molecule) -> Result<()> {
    let mut file = File::create("eht_output.txt")?;
    
    writeln!(file, "=========================================")?;
    writeln!(file, "       EXTENDED HÜCKEL CALCULATION       ")?;
    writeln!(file, "=========================================")?;
    
    // 1. Solve nur EINMAL aufrufen und Daten speichern
    let energies = molecule.energies.as_ref().unwrap();
    let c_matrix = molecule.c_matrix.as_ref().unwrap();
    let n_orbitals = energies.len();
    
    // 2. Erstelle ein Array von 0 bis n und sortiere es anhand der Energien
    let mut sorted_indices: Vec<usize> = (0..n_orbitals).collect();
    sorted_indices.sort_by(|&a, &b| energies[a].partial_cmp(&energies[b]).unwrap());

    // 3. Energien sortiert ausgeben
    writeln!(file, "\nOrbital Energies (eV):")?;
    for (rank, &idx) in sorted_indices.iter().enumerate() {
        writeln!(file, "  MO {:>2}: {:>10.4} eV", rank + 1, energies[idx])?;
    }
    
    writeln!(file, "\n=========================================\n")?;
    
    // 4. Koeffizientenmatrix sortiert ausgeben
    // Zeilen sind Atomorbitale (AO), Spalten sind Molekülorbitale (MO)
    writeln!(file, "Coefficient matrix (C):")?;
    
    // Optional: Kopfzeile für die MOs
    write!(file, "      ")?;
    for rank in 0..n_orbitals {
        write!(file, "        MO {:<2}", rank + 1)?;
    }
    writeln!(file)?;

    for i in 0..c_matrix.nrows() {
        write!(file, "AO {:>2} |", i + 1)?;
        for &idx in &sorted_indices {
            // Wir greifen auf die Spalte 'idx' zu, da das der ursprüngliche Index ist
            write!(file, "{:>12.6}", c_matrix[(i, idx)])?;
        }
        writeln!(file)?;
    }

    writeln!(file, "\n=========================================\n")?;
    
    // 5. Overlap Matrix (S) bleibt wie sie ist, da sie nur von den AOs abhängt
    writeln!(file, "Overlap Matrix (S):")?;
    let s_matrix = molecule.s_mat.as_ref().unwrap();
    for i in 0..s_matrix.nrows() {
        for j in 0..s_matrix.ncols() {
            write!(file, "{:>12.6}", s_matrix[(i,j)])?;
        }
        writeln!(file)?;
    }

    Ok(())
}
    write_output(&molecule).expect("Failed to write output");

    molecule.lowdin_populations();

}