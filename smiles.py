import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

# Streamlit App
st.title("SMILES to 2D Structure Converter 🧪")
st.write("Enter a SMILES string to generate its 2D molecular structure and properties.")

# Input for SMILES string
smiles = st.text_input("Enter SMILES:")

# Function to convert SMILES to 2D image and get properties
def smiles_to_properties(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Generate 2D structure
            img = Draw.MolToImage(mol, size=(300, 300))
            # Get molecular formula
            mol_formula = CalcMolFormula(mol)
            # Get molecular weight
            mol_weight = round(Descriptors.MolWt(mol), 2)
            return img, mol_formula, mol_weight
        else:
            return None, None, None
    except:
        return None, None, None

# Display the results
if smiles:
    img, mol_formula, mol_weight = smiles_to_properties(smiles)
    
    if img:
        st.image(img, caption="Molecular Structure", use_column_width=False)
        st.write(f"**SMILES Entered:** `{smiles}`")
        st.write(f"**Molecular Formula:** `{mol_formula}`")
        st.write(f"**Molecular Weight:** `{mol_weight} g/mol`")
    else:
        st.error("Invalid SMILES string. Please try again.")
