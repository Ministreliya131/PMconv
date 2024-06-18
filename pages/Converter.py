import streamlit as st
import aiohttp
import asyncio
import time
import pandas as pd
import matplotlib.pyplot as plt
from venn import venn
import os
from bs4 import BeautifulSoup
import requests
import networkx as nx
import plotly.graph_objs as go
from PIL import Image
from streamlit_agraph import agraph, Node, Edge, Config

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


st.write("# Welcome to PMInteractor!")

st.markdown(
    """
        ## PMInteractor is an open-source app framework built for proteometabolomics analyses
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
    st.button("Submit", disabled=not metabolites, key="metabolite")

    if metabolites is not None:
        metabolite_split = metabolites.splitlines()
        if len(metabolite_split) > 0:
            MET_SET = set(metabolite_split)
            st.write(len(MET_SET))

if st.button("Process your profiles", type="primary"):

    if PROT_SET is not None:
        with open(
                'pages/protein_dictionary_hmdb.csv'
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

        with open("proteome_to_metabolome_dictionary.csv", "w") as f2:
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
        start_time = time.time()


        async def get_task():
            tasks = []
            async with aiohttp.ClientSession() as session:
                for p_hmdb in MET_SET:
                    task = asyncio.create_task(get_prot_and_gene(session, p_hmdb))
                    tasks.append(task)
                    await asyncio.sleep(1)
                await asyncio.gather(*tasks)


        async def get_prot_and_gene(session, p_hmdb):
            url = f"https://hmdb.ca/metabolites/{p_hmdb}/metabolite_protein_links"
            async with session.get(url) as response:
                if response.status == 200:
                    soup = BeautifulSoup(await response.text(), "html.parser")

                    mdata = []
                    table = soup.find('table', attrs={"class": "table table-condensed"
                                                               " table-striped metabolite-protein-links proteins"})
                    table_body = table.find('tbody')
                    rows = table_body.find_all('tr')
                    for row in rows:
                        cols = row.find_all('td')
                        cols = [ele.text.strip() for ele in cols]
                        mdata.append([ele for ele in cols if ele])

                    with open("pages/tmp", "a") as f:
                        for i in mdata:
                            print(p_hmdb, i[2], i[3], sep=",", file=f)


        st.write("Processing metabolome profile...")

        asyncio.run(get_task())
        end_time = time.time() - start_time
        st.write(f"Finished search with {round(end_time / 60, 2)} minutes!")

        met_dict = pd.read_csv("pages/tmp", sep=",", header=None)
        met_dict_g = met_dict.groupby(by=0, as_index=False)[1].agg("|".join)
        met_dict_g.to_csv("metabolome_to_proteome_dictionary.csv", index=False, header=None)

        assoc_prot = set()
        with open("associated_proteome.csv", "w") as f4:
            with open("pages/tmp") as f:
                for line in f:
                    line_s = line.strip().split(",")
                    assoc_prot.add(line_s[1])
                    prot_gene_dict.update({line_s[1]: line_s[2]})
                    print(line_s[1], file=f4)

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
        st.subheader("Metabolome vs Metabolome*")
        dict_venn_met = {"Metabolome": set(MET_SET), "Metabolome*": set(ASSOC_MET_SET)}
        cmap_format = ["#00BFFF", "#7B68EE"]
        venn(dict_venn_met, cmap=cmap_format, fontsize=20, figsize=(8, 8))
        st.pyplot(plt)

        prot_df = pd.read_csv("proteome_to_metabolome_dictionary.csv", sep=",", header=None)
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
        st.subheader("Proteome vs Proteome*")
        dict_venn_prot = {"Proteome": set(PROT_SET), "Proteome*": set(ASSOC_PROT_SET)}
        cmap_format = ["#DC143C", "#800000"]
        venn(dict_venn_prot, cmap=cmap_format, fontsize=20, figsize=(8, 8))
        st.pyplot(plt)

        met_df = pd.read_csv("metabolome_to_proteome_dictionary.csv", sep=",", header=None)
        met_df.rename(columns={0: "HMDB ID",
                               1: "Associated proteins"}, inplace=True)

        met_df["Protein count"] = met_df["Associated proteins"].apply(lambda x: len(x.split("|")))
        met_df_sort = met_df.sort_values(by="Protein count", ascending=False)

        st.write(f"##### {len(MET_SET - set(met_df['HMDB ID'].to_list()))} metabolites have no associations")
        st.write(f"##### {len(set(met_df['HMDB ID'].to_list()))} metabolites have associations with proteins")
        st.table(met_df_sort[["HMDB ID", "Protein count"]].head(6))

else:
    pass
