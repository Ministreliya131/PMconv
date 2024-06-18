import streamlit as st

def main_page():
    st.markdown("# Main page 🎈")
    st.sidebar.markdown("# Main page 🎈")

def page2():
    st.markdown("# Converter")
    st.sidebar.markdown("# Converter")

def page3():
    st.markdown("# Interactor")
    st.sidebar.markdown("# Interactor")

page_names_to_funcs = {
    "Main Page": main_page,
    "Converter": page2,
    "Interactor": page3,
}

selected_page = st.sidebar.selectbox("Select a page", page_names_to_funcs.keys())
page_names_to_funcs[selected_page]()
