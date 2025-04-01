import streamlit as st
import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import Descriptors
from io import BytesIO

st.set_page_config(page_title="SMILES Descriptor Generator", layout="centered")
st.title("Molecular Descriptor Generator from SMILES")

# File uploader
uploaded_file = st.file_uploader("Upload a CSV or Excel file containing a SMILES column:", type=["csv", "xlsx"])

def is_valid_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        Chem.SanitizeMol(mol)
        return True
    except:
        return False

def compute_properties(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        return {
            "Molecular Weight (MolWt)": Descriptors.MolWt(mol),
            "LogP (Octanol-Water Partition Coefficient)": Descriptors.MolLogP(mol),
            "Topological Polar Surface Area (TPSA)": Descriptors.TPSA(mol),
            "Number of Hydrogen Bond Donors (HBD)": Descriptors.NumHDonors(mol),
            "Number of Hydrogen Bond Acceptors (HBA)": Descriptors.NumHAcceptors(mol),
            "Number of Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
            "Number of Rings": Descriptors.RingCount(mol),
            "Fraction of sp3 Carbon Atoms (FractionCSP3)": Descriptors.FractionCSP3(mol),
            "Number of Heavy Atoms": Descriptors.HeavyAtomCount(mol),
            "Number of Aromatic Rings": Descriptors.NumAromaticRings(mol),
            "Number of Aliphatic Rings": Descriptors.NumAliphaticRings(mol),
            "Number of Saturated Rings": Descriptors.NumSaturatedRings(mol),
            "Number of Heteroatoms": Descriptors.NumHeteroatoms(mol),
            "Balaban J Index": Descriptors.BalabanJ(mol),
            "Bertz Complexity": Descriptors.BertzCT(mol),
            "Chi0v (Molecular Connectivity Index)": Descriptors.Chi0v(mol),
            "Chi1v (Molecular Connectivity Index)": Descriptors.Chi1v(mol),
            "Chi2v (Molecular Connectivity Index)": Descriptors.Chi2v(mol),
            "Kappa1 (Kappa Shape Index)": Descriptors.Kappa1(mol),
            "Kappa2 (Kappa Shape Index)": Descriptors.Kappa2(mol),
            "Hall-Kier Alpha Index": Descriptors.HallKierAlpha(mol),
        }
    except:
        return {}

if uploaded_file:
    file_ext = uploaded_file.name.split(".")[-1].lower()

    if file_ext == "csv":
        df = pd.read_csv(uploaded_file)
    elif file_ext == "xlsx":
        df = pd.read_excel(uploaded_file, engine="openpyxl")
    else:
        st.error("Unsupported file type.")
        st.stop()

    if "SMILES" not in df.columns:
        st.error("Missing required 'SMILES' column.")
        st.stop()

    # Clean and validate
    df = df.copy()
    df = df[df["SMILES"].apply(is_valid_smiles)].reset_index(drop=True)

    with st.spinner("Computing molecular descriptors..."):
        descriptors = df["SMILES"].apply(compute_properties)
        df_descriptors = pd.DataFrame(list(descriptors))
        df_final = pd.concat([df, df_descriptors], axis=1)

    st.success("Descriptors generated!")
    st.dataframe(df_final.head())

    # Download button
    def convert_df(df_out):
        output = BytesIO()
        if file_ext == "csv":
            df_out.to_csv(output, index=False)
        else:
            df_out.to_excel(output, index=False, engine="openpyxl")
        output.seek(0)
        return output

    output_data = convert_df(df_final)
    st.download_button(
        label=" Download file with descriptors",
        data=output_data,
        file_name=uploaded_file.name.replace(f".{file_ext}", f"_with_descriptors.{file_ext}"),
        mime="application/octet-stream"
    )
