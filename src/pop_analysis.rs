use nalgebra as na;

pub fn orth_p_matrix(
    c_matrix: &na::DMatrix<f64>,
    s_sqrt: &na::DMatrix<f64>,
    num_occ: usize,
    sorted_indices: &[usize], // <--- Das Array mit den sortierten Rängen!
) -> na::DMatrix<f64> {
    let mut p_matrix = na::DMatrix::<f64>::zeros(c_matrix.nrows(), c_matrix.nrows());
    
    for i in 0..c_matrix.nrows() {
        for j in 0..c_matrix.nrows() {
            let mut sum = 0.0;
            for k in 0..num_occ {
                // Welches ist die Spalte mit der k-tiefsten Energie?
                let mo_idx = sorted_indices[k]; 
                
                // Nutze mo_idx! Und vergiss die 2.0 für 'Closed Shell' nicht!
                sum += 2.0 * c_matrix[(i, mo_idx)] * c_matrix[(j, mo_idx)];
            }
            p_matrix[(i, j)] = sum;
        }
    }
    
    s_sqrt * p_matrix * s_sqrt
}
