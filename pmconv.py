import streamlit as st
import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# from venn import venn


prot_dict = {}
prot_gene_dict = {}
ASSOC_MET_SET = None
ASSOC_PROT_SET = None
MET_SET = None
PROT_SET = None


def append_row(df, row):
    return pd.concat([
                df,
                pd.DataFrame([row], columns=row.index)]
           ).reset_index(drop=True)


st.write("# Welcome to PMconv!")

st.markdown(
    """
        ## PMconv is an open-source app framework built for proteometabolomics analyses
    """
)

st.markdown(
    """
        #### Metabolome* - metabolites associated with proteins\n
        #### Proteome* - proteins associated with metabolites 
    """
)

with st.sidebar:
    st.write(
        """
            Upload your metabolome and/or proteome profiles
            """
    )

    proteins = st.text_area("Your protein list (one protein per line)")
    st.button("Submit", disabled=not proteins, key="protein")

    if proteins is not None:
        protein_split = proteins.splitlines()
        if len(protein_split) > 0:
            PROT_SET = set(protein_split)
            st.write(len(PROT_SET))

    metabolites = st.text_area("Your metabolite list (one metabolite per line)")
    #filters = st.multiselect("Select filters", options=["Brutto", "Mol_weight", "Class", "Subclass", "Ref", "Proteins", "Origin"], default='Brutto')
    apply_filter = st.checkbox("Group by brutto formulas")
    ref_filter = st.text_input("Number of references threshold", value=0)
    nprot_filter = st.text_input("Number of associated proteins threshold", value=0)
    origin_filter = st.text_input("Origin", value="Endogenous")
    st.button("Submit", disabled=not metabolites, key="metabolite")

    if metabolites is not None:
        metabolite_split = metabolites.splitlines()
        if len(metabolite_split) > 0:
            MET_SET = set(metabolite_split)
            st.write(len(MET_SET))

if st.button("Process your profiles", type="primary"):

    if PROT_SET is not None:
        with open(
                '/home/ministreliya131/PycharmProjects/PMInteractor/protein_dictionary_hmdb.csv'
        ) as p_d:
            for line in p_d:
                line_strip = line.strip()
                key_prot = line_strip.split(",")[0]
                met_line = line_strip.split(",")[1]
                if len(met_line) > 0:
                    met_values = met_line.split("|")
                else:
                    met_values = met_line
                prot_dict.update({key_prot: met_values})
        st.write("Processing proteome profile...")
        st.write(f"{len(prot_dict)} proteins were downloaded from HMDB...")

        assoc_met = []

        with open("associated_metabolome_dictionary.csv", "w") as f2:
            for p in PROT_SET:
                if p in prot_dict:
                    assoc_met += prot_dict[p]
                    print(p, "|".join([x for x in prot_dict[p]]), sep=",", file=f2)
                else:
                    print(p, "None", sep=",", file=f2)

        ASSOC_MET_SET = assoc_met

        with open("associated_metabolome.csv", "w") as f3:
            for m in ASSOC_MET_SET:
                print(m, file=f3)

        st.write("Finished proteome processing!")

    else:
        st.write("No proteome profile was given")

    if MET_SET is not None:
        hmdb = pd.read_csv('/home/ministreliya131/PycharmProjects/PMInteractor/metabolite_dictionary_hmdb.csv', sep=";")

        met_set_df = pd.DataFrame({"Metabolite set input": list(MET_SET)})
        hmdb_met_set = pd.merge(hmdb, met_set_df,
                                left_on="Id", right_on="Metabolite set input", how="inner")

        st.write("Processing metabolome profile...")

        if apply_filter:
            hmdb_met_set["Proteins list"] = hmdb_met_set["Proteins list"].astype(str)
            met_dict_b = hmdb_met_set.groupby(by="Brutto", as_index=False)["Proteins list"].agg("|".join)
            met_dict = pd.merge(hmdb_met_set, met_dict_b, left_on="Brutto", right_on="Brutto", how="inner", suffixes=("", "_y"))
            #st.table(met_dict.head(6))
            met_dict = met_dict[met_dict["Ref"] >= int(ref_filter)]
            met_dict = met_dict[met_dict["Proteins"] >= int(nprot_filter)]
            met_dict = met_dict[met_dict["Origin"] == origin_filter]
            met_dict = met_dict[["Brutto", "Proteins list"]]
            met_dict.to_csv("associated_proteome_dictionary.csv", index=False)
        else:
            met_dict = hmdb_met_set[hmdb_met_set["Ref"] >= int(ref_filter)]
            met_dict = met_dict[met_dict["Proteins"] >= int(nprot_filter)]
            met_dict = met_dict[met_dict["Origin"] == origin_filter]
            met_dict = met_dict[["Id", "Proteins list"]]
            met_dict.to_csv("associated_proteome_dictionary.csv", index=False)


        assoc_prot = set()
        with open("associated_proteome.csv", "w") as f4:
            for line in met_dict["Proteins list"]:
                if type(line) == float:
                    #st.write(line)
                    pass
                else:
                    line_s = line.split("|")
                    for p in line_s:
                        assoc_prot.add(p)

        ASSOC_PROT_SET = assoc_prot

        # os.remove("/home/ministreliya131/PycharmProjects/PMInteractor/tmp")

        st.write("Finished metabolome processing!")

    else:
        st.write("No metabolome profile was given")

else:
    pass

col1, col2 = st.columns([1.5, 1.5])

if ASSOC_MET_SET and ASSOC_PROT_SET is not None:
    with col1:
        prot_df = pd.read_csv("associated_metabolome_dictionary.csv", sep=",", header=None)
        prot_df.rename(columns={0: "Uniprot ID",
                                1: "Associated metabolites"}, inplace=True)
        prot_df_single = prot_df[prot_df["Associated metabolites"].isna()]
        prot_df_multi = prot_df.dropna()
        prot_df_multi["Metabolite count"] = prot_df_multi["Associated metabolites"].apply(lambda x: len(x.split("|")))
        prot_df_multi_sort = prot_df_multi.sort_values(by="Metabolite count", ascending=False)

        st.write(f"##### {prot_df_single.shape[0]} proteins have no associations")
        st.write(f"##### {len(PROT_SET) - prot_df_single.shape[0]} proteins have associations with metabolites")
        st.table(prot_df_multi_sort[["Uniprot ID", "Metabolite count"]].head(6))

    with col2:
        met_df = pd.read_csv("associated_proteome_dictionary.csv", sep=",", header=None)
        met_df.rename(columns={0: "HMDB ID",
                               1: "Associated proteins"}, inplace=True)

        met_df["Protein count"] = met_df["Associated proteins"].apply(lambda x: 0 if type(x) == float else len(x.split("|")))
        met_df_sort = met_df.sort_values(by="Protein count", ascending=False)

        st.write(f"##### {len(MET_SET - set(met_df['HMDB ID'].to_list()))} metabolites have no associations")
        st.write(f"##### {len(set(met_df['HMDB ID'].to_list()))} metabolites have associations with proteins")
        st.table(met_df_sort[["HMDB ID", "Protein count"]].head(6))

else:
    pass
