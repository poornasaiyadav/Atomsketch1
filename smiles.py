import streamlit as st
import py3dmol
import openbabel
import cairosvg
from openbabel import pybel

# Title of the App
st.title("SMILES to 2D Structure Viewer")

# Input SMILES
smiles = st.text_input("Enter a SMILES string:", "CCO")  # Example: Ethanol

if smiles:
    # Convert SMILES to Molecule
    mol = pybel.readstring("smi", smiles)
    
    # Get Molecular Formula & Molecular Weight
    molecular_formula = mol.formula
    molecular_weight = mol.molwt
    
    # Display Molecular Information
    st.write(f"**SMILES Entered:** {smiles}")
    st.write(f"**Molecular Formula:** {molecular_formula}")
    st.write(f"**Molecular Weight:** {molecular_weight:.2f} g/mol")
    
    # Convert to 2D Structure
    mol.make3D()
    obConversion = openbabel.OBConversion()
    obConversion.SetOutFormat("svg")
    svg_str = obConversion.WriteString(mol.OBMol)
    
    # Save SVG and Render in Streamlit
    with open("mol.svg", "w") as f:
        f.write(svg_str)
    
    st.image("mol.svg", caption="2D Molecular Structure")

