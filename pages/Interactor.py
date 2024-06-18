import streamlit as st
import pandas as pd
import os
import requests
from streamlit_agraph import agraph, Node, Edge, Config
import json

PROT_LIST = []
MET_LIST = []
ASSOC_MET_LIST = []
prot_gene_dict = {}


def append_row(df, row):
    return pd.concat([
        df,
        pd.DataFrame([row], columns=row.index)]
    ).reset_index(drop=True)


def get_key_by_value(dict, search):
    for name, val in dict.items():
        if val == search:
            return name


st.write("# Welcome to PMInteractor!")

with open("pages/gene_prot_dict.tsv") as ff:
    for line in ff:
        l_s = line.strip().split("\t")
        prot_gene_dict.update({l_s[0]: l_s[1]})

with open("pages/proteome_to_metabolome_dictionary.csv") as pm:
    for line in pm:
        ll = line.strip().split(",")
        PROT_LIST.append(ll[0])
        if ll[1] != "None":
            for j in ll[1].split("|"):
                ASSOC_MET_LIST.append(j)

with open("pages/metabolome_to_proteome_dictionary.csv") as mp:
    for line in mp:
        ll = line.strip().split(",")
        MET_LIST.append(ll[0])

gene_to_string = []
err_lst = []
for p in PROT_LIST:
    try:
        gene_to_string.append(prot_gene_dict[p])
    except KeyError:
        err_lst.append(p)

# string_api_url = "https://version-11-5.string-db.org/api"
# output_format = "json"
# method = "enrichment"
#
# request_url = "/".join([string_api_url, output_format, method])
#
# params = {
#     "identifiers": "%0d".join(gene_to_string),
#     "species": 9606,
#     "caller_identity": "PMInteractor"
# }
#
# response = requests.post(request_url, data=params)
#
# data = json.loads(response.text)
#
# func_dict = pd.DataFrame(columns=["Term", "Protein List",
#                                   "Category", "Description", "FDR"])
# for row in data:
#
#     term = row["term"]
#     preferred_names = ",".join(row["preferredNames"])
#     fdr = float(row["fdr"])
#     description = row["description"]
#     category = row["category"]
#     my_input = str(row["number_of_genes"])
#     bkg = str(row["number_of_genes_in_background"])
#
#     if fdr < 0.01:
#         # st.write("\t".join([term, category, my_input,
#         #                    bkg, preferred_names, str(fdr), description]))
#         new_row = pd.Series({"Term": term, "Protein List": preferred_names,
#                              "Category": category.upper(), "Description": description,
#                              "FDR": fdr})
#         func_dict = append_row(func_dict, new_row)
#
# options_1 = st.multiselect(
#     "Choose categories",
#     list(func_dict["Category"].unique()),
#     list(func_dict["Category"].unique())[0])
#
# opt1_df = pd.DataFrame({"Category": options_1})
# for_inter_1 = pd.merge(func_dict, opt1_df, left_on="Category", right_on="Category", how="inner")
# for_inter_1 = for_inter_1[["Term", "Protein List", "Category", "Description", "FDR"]]
#
# options_2: object = st.multiselect(
#     "Choose pathways/functions/compartments",
#     list(for_inter_1["Description"].unique()),
#     list(for_inter_1["Description"].unique())[0])
#
# opt2_df = pd.DataFrame({"Description": options_2})
# for_inter_2 = pd.merge(for_inter_1, opt2_df, left_on="Description", right_on="Description", how="inner")
# for_inter_2 = for_inter_2[["Term", "Protein List", "Category", "Description", "FDR"]]
#
# pre_prot_list = []
# for pl in for_inter_2["Protein List"]:
#     pre_prot_list += pl.split(",")
# prot_lst = list(set(pre_prot_list))
#
# plst = [get_key_by_value(prot_gene_dict, p) for p in prot_lst]

string_api_url_2 = "https://version-12-0.string-db.org/api"
output_format_2 = "tsv-no-header"
method_2 = "network"

request_url_2 = "/".join([string_api_url_2, output_format_2, method_2])

params_2 = {
    "identifiers": "%0d".join(PROT_LIST),
    "species": 9606,
    "caller_identity": "PMInteractor",
    "network_type": "physical"
}

response_2 = requests.post(request_url_2, data=params_2)

# st.write("got response")

met_prot_inter = pd.DataFrame(columns=['Source', 'Target', 'Score'])
met_prot_d = pd.read_csv("pages/proteome_to_metabolome_dictionary.csv",
                         sep=",", header=None)
for index, line in met_prot_d.iterrows():
    if line[0] in PROT_LIST:
        if str(line[1]) == 'nan':
            continue
        else:
            for s_line in line[1].split("|"):
                if s_line in MET_LIST:
                    new_row = pd.Series({'Source': prot_gene_dict[line[0]], 'Target': s_line, 'Score': 1.0})
                    met_prot_inter = append_row(met_prot_inter, new_row)

df_inter = pd.DataFrame(columns=['Source', 'Target', 'Score'])
for line in response_2.text.strip().split("\n"):
    l = line.strip().split("\t")
    p1, p2 = l[2], l[3]

    experimental_score = float(l[5])
    if experimental_score >= 0.9:
        new_row = pd.Series({'Source': p1, 'Target': p2, 'Score': experimental_score})
        df_inter = append_row(df_inter, new_row)

df_i = pd.concat([met_prot_inter, df_inter], ignore_index=True)
pre_nodes_l = set(list(df_i["Source"]) + list(df_i["Target"]))

nodes = []
edges = []
for i in pre_nodes_l:
    if "HMDB" in i:
        nodes.append(Node(id=i,
                          label=i,
                          size=25,
                          color="#0000CD",
                          shape="square"))
    else:
        nodes.append(Node(id=i,
                          label=i,
                          size=25,
                          color="#DC143C",
                          shape="dot"))
for i, lines in df_i.iterrows():
    edges.append(Edge(source=lines[0],
                      target=lines[1]))

config = Config(width=650,
                height=850,
                directed=False,
                physics=False,
                hierarchical=False,
                # **kwargs
                )

return_value = agraph(nodes=nodes,
                      edges=edges,
                      config=config)
