import requests
from streamlit_agraph import agraph, Node, Edge, Config
import pandas as pd
import streamlit as st
import os



prot_dict = {}
prot_gene_dict = {}

def append_row(d, row):
    return pd.concat([
                d,
                pd.DataFrame([row], columns=row.index)]
           ).reset_index(drop=True)

def get_key_by_value(dct, search):
    for name, val in dct.items():
        if val == search:
            return name


st.set_page_config(
    page_title="PMconv",
    layout="wide"  # 'centered' is the default, 'wide' expands to browser width
)

# Example files
xlsx_file = os.path.join(os.path.dirname(__file__), 'example', 'example.xlsx')

p_g_dict = os.path.join(os.path.dirname(__file__), 'data', 'gene_prot_dict.tsv')
with open(p_g_dict) as ff:
    for line in ff:
        l_s = line.strip().split("\t")
        prot_gene_dict.update({l_s[0]: l_s[1]})

# Main page
st.title('PMconv')
st.markdown(
    """
        #### PMconv is an open-source app built for proteometabolomics analyses
    """
)

st.markdown(
    """
        ###### Metabolome* - metabolites associated with proteins\n
        ###### Proteome* - proteins associated with metabolites 
    """
)

df = []
selection_lists = []

col1, col2 = st.columns([1.5, 1.5])

with st.sidebar:
    logo = os.path.join(os.path.dirname(__file__), 'data', 'pmconv.png')
    st.sidebar.image(logo, use_column_width=True)

    st.write(
        """
            #### Apply filters to your metabolome profile
            """
    )
    apply_filter = st.checkbox("Group by brutto formulas")
    ref_filter = st.text_input("The number of references threshold", value=0)
    nprot_filter = st.text_input("The number of associated proteins threshold", value=0)
    origin_filter = st.text_input("Origin (Endogenous, Environmental, Food or None)", value="Endogenous")
    st.button("Submit", key="metabolite")

with col1:
    # Example section
    st.subheader("📎 Example")

    st.link_button("More", "Github link")

    demo = st.checkbox("**Try example**", value=1)
    if demo:  # Demo mode

        df_prot = pd.read_excel(xlsx_file, sheet_name="Proteins")
        df_met = pd.read_excel(xlsx_file, sheet_name="Metabolites")


    else:

        with st.expander("**.xlsx template**", expanded=False):
            with open(xlsx_file, "rb") as file:  # Download .xlsx template
                st.download_button(
                    label="Download example.xlsx",
                    data=file,
                    file_name="example.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

        # Upload data section
        st.subheader("💽 Upload data")

        st.markdown(
            """
                ###### Can handle up to 200 protein identifiers
            """
        )

        uploaded_file = st.file_uploader("**Upload .xlsx file**", type="xlsx",
                                      accept_multiple_files=False)

        if uploaded_file is not None:
            try:
                df_prot = pd.read_excel(uploaded_file, sheet_name="Proteins")
            except:
                st.warning("It appears that there is only metabolome profile in your data...")
            try:
                df_met = pd.read_excel(uploaded_file, sheet_name="Metabolites")
            except:
                st.warning("It appears that there is only proteome profile in your data...")
        else:
            st.warning("Please upload a file to continue.")
            st.stop()

    if df_prot is not None:
        prot_set = set(df_prot["Uniprot IDs"].to_list())
        p_dict_hmdb = os.path.join(os.path.dirname(__file__), 'data', 'protein_dictionary_hmdb.csv')
        p_dict_loc = os.path.join(os.path.dirname(__file__), 'data', 'proteins_localization_HPA.tsv')
        with open(p_dict_hmdb) as p_d:
            for line in p_d:
                line_strip = line.strip()
                key_prot = line_strip.split(",")[0]
                met_line = line_strip.split(",")[1]
                if len(met_line) > 0:
                    met_values = met_line.split("|")
                else:
                    met_values = met_line
                prot_dict.update({key_prot: met_values})
        #st.write("Processing proteome profile...")
        #st.write(f"{len(prot_dict)} proteins were downloaded from HMDB...")

        assoc_met_dict_df = pd.DataFrame(columns=["Uniprot", "Metabolites list"])
        for p in prot_set:
            if p in prot_dict:
                new_row = pd.Series({"Uniprot": p, "Metabolites list": "|".join([x for x in prot_dict[p]])})
            else:
                new_row = pd.Series({"Uniprot": p, "Metabolites list": "None"})
            assoc_met_dict_df = append_row(assoc_met_dict_df, new_row)

        p_dict_loc_df = pd.read_csv(p_dict_loc, sep="\t")
        assoc_met_dict_df_loc = pd.merge(assoc_met_dict_df, p_dict_loc_df, 
                                         left_on="Uniprot", right_on="Entry",
                                         how="left")
        assoc_met_dict_df_loc_f = assoc_met_dict_df_loc[["Uniprot", "Secretome location",
                                                         "Subcellular main location", "Metabolites list"]]
        assoc_met_dict_csv = assoc_met_dict_df_loc_f.to_csv(sep=";", index=False).encode("utf-8")
        st.download_button(
            label="📥 Download metabolome* as CSV",
            data=assoc_met_dict_csv,
            file_name="associated_metabolome_dictionary.csv",
            mime="text/csv"
        )
        #st.write("Finished proteome processing!")

    else:
        st.write("No proteome profile was given")

    if df_met is not None:
        hmdb_data = os.path.join(os.path.dirname(__file__), 'data', 'metabolite_dictionary_hmdb.csv')
        p_dict_loc = os.path.join(os.path.dirname(__file__), 'data', 'proteins_localization_HPA.tsv')
        p_dict_loc_df = pd.read_csv(p_dict_loc, sep="\t")
        hmdb = pd.read_csv(hmdb_data, sep=";")
        hmdb_met_set = pd.merge(hmdb, df_met,
                                left_on="Id", right_on="HMDB IDs", how="inner")

        #st.write("Processing metabolome profile...")

        if apply_filter:
            hmdb_met_set["Proteins list"] = hmdb_met_set["Proteins list"].astype(str)
            met_dict_b = hmdb_met_set.groupby(by="Brutto", as_index=False)["Proteins list"].agg("|".join)
            met_dict = pd.merge(hmdb_met_set, met_dict_b, left_on="Brutto", right_on="Brutto", how="inner",
                                suffixes=("", "_y"))
            # st.table(met_dict.head(6))
            met_dict = met_dict[met_dict["Ref"] >= int(ref_filter)]
            met_dict = met_dict[met_dict["Proteins"] >= int(nprot_filter)]
            met_dict = met_dict[met_dict["Origin"] == origin_filter]
            met_dict = met_dict[["Brutto", "Proteins list"]]
            met_dict_csv = met_dict.to_csv(sep=";", index=False).encode("utf-8")

            assoc_p_l = []
            for i in met_dict["Proteins list"].astype(str):
                if len(i.split("|")) > 1:
                    for j in i.split("|"):
                        assoc_p_l.append(j)
                else:
                    assoc_p_l.append(i)

            met_dict_prot_df = pd.DataFrame({"Uniprot": assoc_p_l})
            prot_asoc_loc = pd.merge(p_dict_loc_df, met_dict_prot_df,
                                     left_on="Entry", right_on="Uniprot",
                                     how="right")
            prot_asoc_loc_f = prot_asoc_loc[["Uniprot", "Secretome location",
                                                         "Subcellular main location"]]
            prot_asoc_loc_csv = prot_asoc_loc_f.to_csv(sep=";", index=False).encode("utf-8")
            
        else:
            met_dict = hmdb_met_set[hmdb_met_set["Ref"] >= int(ref_filter)]
            met_dict = met_dict[met_dict["Proteins"] >= int(nprot_filter)]
            met_dict = met_dict[met_dict["Origin"] == origin_filter]
            met_dict = met_dict[["Id", "Proteins list"]]
            met_dict_csv = met_dict.to_csv(sep=";", index=False).encode("utf-8")

            assoc_p_l = []
            for i in met_dict["Proteins list"].astype(str):
                if len(i.split("|")) > 1:
                    for j in i.split("|"):
                        assoc_p_l.append(j)
                else:
                    assoc_p_l.append(i)

            met_dict_prot_df = pd.DataFrame({"Uniprot": assoc_p_l})
            prot_asoc_loc = pd.merge(p_dict_loc_df, met_dict_prot_df,
                                     left_on="Entry", right_on="Uniprot",
                                     how="right")
            prot_asoc_loc_f = prot_asoc_loc[["Uniprot", "Secretome location",
                                                         "Subcellular main location"]]
            prot_asoc_loc_csv = prot_asoc_loc_f.to_csv(sep=";", index=False).encode("utf-8")

        st.download_button(
            label="📥 Download proteome* as CSV",
            data=met_dict_csv,
            file_name="associated_proteome_dictionary.csv",
            mime="text/csv"
        )

        st.download_button(
            label="📥 Download proteome* proteins localization as CSV",
            data=prot_asoc_loc_csv,
            file_name="associated_proteome_dictionary_localization.csv",
            mime="text/csv"
        )

        #st.write("Finished metabolome processing!")

    else:
        st.write("No metabolome profile was given")

    with col2:
        st.subheader('Welcome to PMconv!')
        st.write('You are by default in **demo** mode.\n'
                'You can disable **Try example** on the left **📎 Example** section.')

with col2:
    st.markdown("##### Hint!")
    st.markdown("In this version of the PMconv,\n"
                "you can save an interactive graph by right-clicking\n"
                "and selecting the 'save image as' function.\n"
                "The graph can be zoomed in, zoomed out, and moved:\n"
                "the saved image will contain the fragment of interest.")

    st.markdown("### Legend")
    c1, c2 = st.columns(2)
    with c1:
        st.markdown("🔴 **Protein** (Red dot)")
    with c2:
        st.markdown("🟦 **Metabolite** (Blue square)")


if df_prot is not None and df_met is not None:
    gene_to_string = []
    err_lst = []
    protlist = list(set(df_prot["Uniprot IDs"].to_list()))
    for p in protlist:
        try:
            gene_to_string.append(prot_gene_dict[p])
        except KeyError:
            err_lst.append(p)

    string_api_url_2 = "https://version-12-0.string-db.org/api"
    output_format_2 = "tsv-no-header"
    method_2 = "network"

    request_url_2 = "/".join([string_api_url_2, output_format_2, method_2])

    params_2 = {
        "identifiers": "%0d".join(protlist),
        "species": 9606,
        "caller_identity": "PMInteractor",
        "network_type": "physical"
    }

    response_2 = requests.post(request_url_2, data=params_2)

    # st.write("got response")

    #st.write(assoc_met_dict_df.head())
    met_prot_inter = pd.DataFrame(columns=['Source', 'Target', 'Score'])
    for index, line in assoc_met_dict_df.iterrows():
        if line[0] in protlist:
            if str(line[1]) == 'nan':
                continue
            else:
                for s_line in line[1].split("|"):
                    if s_line in df_met["HMDB IDs"].to_list():
                        new_row = pd.Series({'Source': prot_gene_dict[line[0]], 'Target': s_line, 'Score': 1.0})
                        met_prot_inter = append_row(met_prot_inter, new_row)

    df_inter = pd.DataFrame(columns=['Source', 'Target', 'Score'])
    for line in response_2.text.strip().split("\n"):
        #st.write(line)
        l = line.strip().split("\t")
        p1, p2 = l[2], l[3]

        experimental_score = float(l[5])
        if experimental_score >= 0.9:
            new_row = pd.Series({'Source': p1, 'Target': p2, 'Score': experimental_score})
            df_inter = append_row(df_inter, new_row)

    df_i = pd.concat([met_prot_inter, df_inter], ignore_index=True)
    df_i_csv = df_i.to_csv(index=False).encode("utf-8")

    with col1:
        st.download_button(
            label="📥 Download protein-metabolite interaction table as CSV",
            data=df_i_csv,
            file_name="protein_metabolite_interactions.csv",
            mime="text/csv"
        )

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

else:
    with col2:
        st.warning(f"It appears that there is an error..."
                   "Please check your data and make sure that you uploaded both proteome and metabolome profiles.\n"
                   "If this does not resolve the problems, contact me by email (ministreliya13113@gmail.com) or submit an Issue on GitHub.\n\n"
                   , icon='🚨')