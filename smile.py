import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
# NO IMPORT OF IPythonConsole
from io import BytesIO
import base64

# Set page configuration
st.set_page_config(page_title="SMILES to Structure Converter", layout="wide")

def get_molecule_image(mol, size=(300, 300)):
    """Generate image for molecule"""
    img = Draw.MolToImage(mol, size=size)
    buffered = BytesIO()
    img.save(buffered, format="PNG")
    return base64.b64encode(buffered.getvalue()).decode()

def process_smiles(smiles):
    """Process a single SMILES string and return molecule info"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    # Get molecular information
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
    mol_weight = round(Descriptors.MolWt(mol), 2)
    img_str = get_molecule_image(mol)
    
    return {
        "SMILES": smiles,
        "Molecular Formula": formula,
        "Molecular Weight": f"{mol_weight} g/mol",
        "Structure Image": img_str
    }

# App header
st.title("SMILES to Structure Converter")
st.write("Enter a SMILES string to view its chemical structure and information")

# Single SMILES input
smiles_input = st.text_input("Enter SMILES string:", value="CCO")

# Add example dropdown
examples = {
    "Ethanol": "CCO",
    "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Paracetamol": "CC(=O)NC1=CC=C(C=C1)O",
    "Ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
}

example = st.selectbox("Or select an example:", [""] + list(examples.keys()))
if example:
    smiles_input = examples[example]

# Process button
if st.button("Convert SMILES") or smiles_input:
    if smiles_input:
        # Process the SMILES
        result = process_smiles(smiles_input)
        if result:
            # Create two columns for display
            col1, col2 = st.columns([1, 1])
            
            # Display structure in first column
            with col1:
                st.subheader("Chemical Structure")
                st.image(f"data:image/png;base64,{result['Structure Image']}", use_column_width=True)
            
            # Display information in second column
            with col2:
                st.subheader("Molecule Information")
                info_df = pd.DataFrame({
                    "Property": ["SMILES", "Molecular Formula", "Molecular Weight"],
                    "Value": [
                        result["SMILES"],
                        result["Molecular Formula"],
                        result["Molecular Weight"]
                    ]
                })
                st.table(info_df)
        else:
            st.error("Invalid SMILES string. Please check your input.")

# Batch processing section
st.header("Batch Processing")
st.write("Enter multiple SMILES strings (one per line):")

batch_input = st.text_area("Batch SMILES input:", height=150)

if st.button("Process Batch"):
    if batch_input:
        # Split input by lines and process each SMILES
        smiles_list = [s.strip() for s in batch_input.split("\n") if s.strip()]
        
        if smiles_list:
            # Process all SMILES
            results = []
            mols = []
            
            for smiles in smiles_list:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    mols.append(mol)
                    
                    # Get data for table
                    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                    mol_weight = round(Descriptors.MolWt(mol), 2)
                    
                    results.append({
                        "SMILES": smiles,
                        "Molecular Formula": formula,
                        "Molecular Weight": f"{mol_weight} g/mol"
                    })
                else:
                    results.append({
                        "SMILES": smiles,
                        "Molecular Formula": "Invalid SMILES",
                        "Molecular Weight": "N/A"
                    })
            
            # Display results table
            st.subheader("Batch Results")
            results_df = pd.DataFrame(results)
            st.dataframe(results_df)
            
            # Create a grid of molecule images
            if mols:
                st.subheader("Batch Structures")
                img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(200, 200))
                buffered = BytesIO()
                img.save(buffered, format="PNG")
                img_str = base64.b64encode(buffered.getvalue()).decode()
                st.image(f"data:image/png;base64,{img_str}")
                
                # Export option
                csv = results_df.to_csv(index=False)
                b64 = base64.b64encode(csv.encode()).decode()
                st.download_button(
                    label="Download Results as CSV",
                    data=csv,
                    file_name="molecule_batch_results.csv",
                    mime="text/csv"
                )
        else:
            st.warning("No valid SMILES strings provided.")