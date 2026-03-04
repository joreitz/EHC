use std::fs::File;
use std::io::{BufRead, BufReader};
use anyhow::Result;
static _BOHR_TO_ANGSTROM: f64 = 0.529177210903; // Bohr radius in Angstroms
static _ANGSTROM_TO_BOHR: f64 = 1.0 / _BOHR_TO_ANGSTROM; // Angstrom to Bohr conversion factor

pub struct AtomData {
    pub element: String,
    pub position: [f64; 3],
}

pub fn read_xyz(filename: &str) -> Result<Vec<AtomData>, Box<dyn std::error::Error>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    
    // Anzahl überspringen
    lines.next();
    // Kommentar überspringen  
    lines.next();
    
    let mut atoms = Vec::new();
    
    for line in lines {
        let line = line?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        
        if parts.len() < 4 { continue; }
        
        atoms.push(AtomData {
            element: parts[0].to_string(),
            position: [
                parts[1].parse::<f64>()? * _ANGSTROM_TO_BOHR,  // ← Umrechnung!
                parts[2].parse::<f64>()? * _ANGSTROM_TO_BOHR,
                parts[3].parse::<f64>()? * _ANGSTROM_TO_BOHR,
            ],
        });
    }
    
    Ok(atoms)
}

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Deserialize, Serialize)]
pub struct OrbitalParams {
    pub n: u8,
    pub l: u8,
    pub exponent: f64,
    pub coefficient: f64,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct ElementBasis {
    pub atomic_number: u8,
    pub orbitals: Vec<OrbitalParams>,
    pub vsie: Vec<[f64; 3]>,  // [n, l, vsie_value]
}

pub fn load_basis_library(filename: &str) -> Result<HashMap<String, ElementBasis>, Box<dyn std::error::Error>> {
    let file = File::open(filename)?;
    let reader = BufReader::new(file);
    let library = serde_json::from_reader(reader)?;
    Ok(library)
}