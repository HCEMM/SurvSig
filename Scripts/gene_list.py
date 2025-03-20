import streamlit as st
import pandas as pd
import os
import Scripts.user_interaction as uint
import math
from st_keyup import st_keyup
import sys
import Scripts.descriptions as dsc


def organiser(t_inp, gene_name_col1):
    """
    Organize the input gene names by replacing unwanted characters with spaces,
    converting to uppercase, and removing duplicates.

    Parameters:
    t_inp (str): Input gene names as a string.
    gene_name_col1 (streamlit.column): Streamlit column for displaying gene names.

    Returns:
    pd.Series: Organized gene names.
    """
    try:
        unwanted_chars = (
            " ,\"'/()_::;[]{}&@#><%!+*?`~^|\\$€£¥₹§©®™±µ·÷×"
            "¡¿°ªºßäöüÄÖÜáéíóúàèìòùâêîôûãñõçøåœ"
            "ÁÉÍÓÚÀÈÌÒÙÂÊÎÔÛÃÑÕÇØÅŒ"
        )

        for char in unwanted_chars:
            t_inp = t_inp.replace(char, " ")
        t_inp = " ".join(t_inp.split()).upper()
        t_inp_series = pd.Series(t_inp.split()).drop_duplicates()
        return t_inp_series
    except Exception as e:
        st.error(f"Error organizing gene names: {e}")
        return pd.Series()


def name_suggestion(gene_names, df_app, name_dict, gene_name_col1, gene_name_col_exp):
    """
    Suggest correct gene names and display them in a user-friendly format.

    Parameters:
    gene_names (list): List of input gene names.
    df_app (pd.DataFrame): DataFrame of approved gene symbols.
    name_dict (dict): Dictionary mapping gene names to approved symbols.
    gene_name_col1 (streamlit.column): Streamlit column for displaying gene names.
    gene_name_col_exp (streamlit.expander): Streamlit expander for gene name suggestions.

    Returns:
    tuple: DataFrame of suggested gene names, path to saved gene file.
    """
    try:
        suggestions = []
        num_corrected = 0
        not_found_names = []

        for gene_name in gene_names:
            if gene_name not in name_dict:
                not_found_names.append(gene_name)
                continue
            elif gene_name in name_dict and gene_name in df_app['approved_symbol'].values:
                gene_name = gene_name
            elif gene_name in name_dict:
                target_names = name_dict[gene_name]
                if len(target_names) == 1:
                    if target_names[0] != gene_name:
                        suggestion = f':red[{gene_name}] :blue[(consider using ] :green[{target_names[0]}]:blue[)]'
                        options = [target_names[0], 'Remove gene']

                        option = gene_name_col_exp.selectbox(suggestion, options, key="option_copy"+gene_name,
                                                             help=f"https://www.genecards.org/Search/Keyword?queryString={gene_name}")
                        if option == 'Remove gene':
                            continue
                        gene_name = target_names[0]
                        num_corrected += 1
                else:
                    suggestion = f':red[{gene_name}] :blue[(consider using] :green[{", ".join(target_names)}]:blue[)]'
                    options = target_names + ['Remove gene']
                    option = gene_name_col_exp.selectbox(suggestion, options,
                                                         key="option_copy_multi"+gene_name,
                                                         help=f"https://www.genecards.org/Search/Keyword?queryString={gene_name}")
                    if option == 'Remove gene':
                        continue
                    gene_name = option
                    num_corrected += 1

            suggestions.append(gene_name)

        name_list = pd.DataFrame(suggestions, columns=['Gene Names'])

        if len(name_list) != 0:
            gene_name_col1.subheader(':blue[Gene Names Actually List]')

            num_cols = min(max(math.ceil(len(name_list) / 10), 1), 4)
            gene_list_cols = gene_name_col1.columns(num_cols)
            n = math.ceil(len(name_list) / num_cols)

            for i in range(num_cols):
                start = i * n
                end = (i + 1) * n if i < num_cols - 1 else None
                df = pd.DataFrame(name_list[start:end])
                gene_list_cols[i].dataframe(df, use_container_width=True, hide_index=True)

        if num_corrected != 0:
            gene_name_col1.warning(f":blue[Number of corrected genes: ] :red[{num_corrected}]")

        if num_corrected == 0 and len(name_list) >= 1:
            gene_name_col1.success(":blue[All of your genes are ] :red[correct!]")

        if len(not_found_names) != 0:
            gene_name_col1.warning(
                f":blue[Number of ] :red[NOT ] :blue[found genes: ] :red[{str(len(not_found_names))}]")
            gene_name_col1.warning(f":blue[The following names/terms are ] :red[NOT ] :blue[in ] "
                                   f":blue[the gene dictionary: ] :red[{not_found_names}]")

        user_id = uint.get_session_id()
        gene_file = uint.naming(user_id)[38]
        name_list.to_csv(gene_file, index=False, header=False, encoding="utf")
        not_dict_genes = uint.naming(user_id)[61]
        not_found_names = pd.Series(not_found_names)
        not_found_names.to_csv(not_dict_genes, index=False, header=False, encoding="utf")

        return name_list, gene_file
    except Exception as e:
        st.error(f"Error in name suggestion: {e}")
        return pd.DataFrame(), ""


def name_suggestion_cluster(gene_names_with_cluster, df_app, name_dict, gene_name_col1, gene_name_col_exp):
    """
    Suggest correct gene names and display them along with cluster information.

    Parameters:
    gene_names_with_cluster (pd.DataFrame): DataFrame of input gene names with cluster information.
    df_app (pd.DataFrame): DataFrame of approved gene symbols.
    name_dict (dict): Dictionary mapping gene names to approved symbols.
    gene_name_col1 (streamlit.column): Streamlit column for displaying gene names.
    gene_name_col_exp (streamlit.expander): Streamlit expander for gene name suggestions.

    Returns:
    tuple: DataFrame of suggested gene names with clusters, path to saved gene file.
    """
    try:
        gene_names_with_cluster.drop_duplicates(subset=gene_names_with_cluster.columns[0], inplace=True)
        gene_names_with_cluster.reset_index(drop=True, inplace=True)
        suggestions = []
        num_corrected = 0
        not_found_names = []
        i = 0

        for index, row in gene_names_with_cluster.iterrows():
            gene_name = row.iloc[0]
            cluster_info = row.iloc[1]
            i = i + 1

            if gene_name not in name_dict:
                not_found_names.append(gene_name)
                continue
            else:
                if gene_name in df_app['approved_symbol'].values:
                    suggestions.append([gene_name, cluster_info])
                elif gene_name in name_dict:
                    target_names = name_dict[gene_name]
                    if len(target_names) == 1 and target_names[0] != gene_name:
                        suggestion = f':red[{gene_name}] :blue[(consider using ] :green[{target_names[0]}]:blue[)]'
                        options = [target_names[0], 'Remove gene']
                    else:
                        suggestion = f':red[{gene_name}] :blue[(consider using] :green[{", ".join(target_names)}]:blue[)]'
                        options = target_names + ['Remove gene']

                    option = gene_name_col_exp.selectbox(suggestion, options, key="option_copy_cluster"+gene_name+str(i),
                                                         help=f"https://www.genecards.org/Search/Keyword?queryString={gene_name}")
                    if option == 'Remove gene':
                        continue
                    else:
                        gene_name = option
                        num_corrected += 1
                        suggestions.append([gene_name, cluster_info])
                else:
                    suggestions.append([gene_name, cluster_info])

        name_list = pd.DataFrame(suggestions, columns=['Samples', 'Cluster'])
        gene_name_col1.subheader(':blue[Updated Gene Names and Clusters]')
        selected_clusters = gene_name_col1.multiselect(":blue[Choose the ] :red[Desired Cluster(s)]",
                                                       options=name_list["Cluster"].unique(), key="cluster_multi",
                                                       placeholder="Choose Cluster(s)",
                                                       default=name_list["Cluster"].unique())

        if len(selected_clusters) == 0:
            st.error(":blue[Please Select ] :red[Minimum One Cluster]")
        else:
            name_list = name_list[name_list['Cluster'].isin(selected_clusters)]
            name_list.drop_duplicates(keep='first', inplace=True, subset=['Samples'])
            gene_name_col1.dataframe(name_list, use_container_width=True)

            if num_corrected > 0:
                gene_name_col1.warning(f":blue[Number of corrected genes: ] :red[{num_corrected}]")

            if len(not_found_names) > 0:
                gene_name_col1.warning(f":blue[Number of NOT found genes: ] :red[{len(not_found_names)}]")
                gene_name_col1.warning(
                    f":blue[The following names/terms are NOT in the gene dictionary: ] :red[{', '.join(not_found_names)}]")
            user_id = uint.get_session_id()
            gene_file = uint.naming(user_id)[38]
            name_list.to_csv(gene_file, index=False, encoding="utf")
            not_dict_genes = uint.naming(user_id)[61]
            not_found_names = pd.Series(not_found_names)
            not_found_names.to_csv(not_dict_genes, index=False, header=False, encoding="utf")

            return name_list, gene_file
    except Exception as e:
        st.error(f"Error in name suggestion with clusters: {e}")
        return pd.DataFrame(), ""


@st.cache_resource(show_spinner='Gene Dictionary & List Maker are Loading')
def load_gene_name_dict():
    """
    Load the gene name dictionary into a DataFrame.

    Returns:
    tuple: DataFrame of approved gene symbols, dictionary mapping gene names to approved symbols.
    """
    try:
        df = pd.read_csv('source_data/approved_symbols_sep.txt', dtype=str)
        headers = df.columns
        df = df.map(lambda x: x.upper() if isinstance(x, str) else x)
        df.columns = headers

        name_dict = {}
        for i in range(len(df)):
            target_name = df.loc[i, 'approved_symbol']
            for col_name in df.columns:
                if col_name.startswith('approved_symbol') or col_name.startswith(
                        'previous_symbol_') or col_name.startswith(
                        'alias_symbol_'):
                    old_name = df.loc[i, col_name]
                    if pd.notnull(old_name):
                        if old_name not in name_dict:
                            name_dict[old_name] = []
                        name_dict[old_name].append(target_name)

        return df, name_dict
    except Exception as e:
        st.error(f"Error loading gene name dictionary: {e}")
        return pd.DataFrame(), {}


def app(user_session_id, df, option_dataset, uploaded_file):
    """
    Main function for the gene list application.

    Parameters:
    user_session_id (str): User session ID.
    df (pd.DataFrame): DataFrame containing gene expression data.
    option_dataset (str): Option for dataset.
    uploaded_file (UploadedFile): Uploaded file containing gene names.

    Returns:
    pd.DataFrame: Selected genes DataFrame.
    """
    try:
        st.warning(':blue[Avoid duplicates gene names, because the application ] :red[keeps only the first match ] :blue[, deletes others]')
        if uploaded_file is not None:
            uploaded_file.seek(0)

        df_app, name_dict = load_gene_name_dict()
        gene_name_col1, gene_name_col2 = st.columns(2)
        gene_name_col2.caption("", help=dsc.gl_text())
        gene_name_col_exp = gene_name_col2.expander(":blue[Hide Options for Mapping]", expanded=True)
        selected_genes = None

        gene_names = []

        if uploaded_file is None:
            up_e = st.empty()
            up_e.warning(":blue[Please ] :red[upload your gene list ] :blue[or ] :red[you can create one]")

            st.warning(
                ":blue[If you upload a file, the provided genes will be ] :red[REMOVED]! :blue[Please upload your file first, and then ] :red[ADD ] :blue[genes manually.]")

            gene_input = gene_name_col1.text_input(
                ":blue[Add ] :red[More Genes]:blue[. You Can Copy Gene List Here!]",  key="gl_maker_feature")

            t_inp = organiser(gene_input, gene_name_col1)
            gene_names = t_inp.tolist()
            gene_names = list(set(gene_names))


            if len(gene_names) <= 2500:
                name_list, gene_file = name_suggestion(gene_names, df_app, name_dict, gene_name_col1,
                                                       gene_name_col_exp)

                if len(name_list) != 0:
                    up_e.empty()

                with open(gene_file, "r") as f:
                    gene_ids = f.read().splitlines()

                selected_genes = df[df.index.isin(gene_ids)]
                selected_genes = selected_genes.dropna(axis=1)

                difference = list(set(name_list["Gene Names"]) - set(selected_genes.index.tolist()))
                corrected_genes = list(set(gene_names) - set(selected_genes.index.tolist()))
                corrected_genes = set(corrected_genes) - set(difference)

                not_database_genes = uint.naming(user_session_id)[39]
                corrected_genes_list = uint.naming(user_session_id)[40]

                if len(difference) > 0:
                    st.warning(
                        f':blue[The following gene(s) is/are ] :red[NOT ] :blue[it this dataset: ] :red[{difference}]')
                    difference = pd.Series(difference)
                    difference.to_csv(not_database_genes, index=False, header=None)
                elif selected_genes.shape[1] != 0 and len(difference) > 0:
                    gene_name_col1.success("All genes what provide, are in the dataset!")
                elif selected_genes.shape[1] == 0:
                    gene_name_col1.error(
                        ":blue[Please ] :red[add gene(s) ] :blue[or ] :red[try other gene(s) name]")

                if len(corrected_genes) > 0:
                    corrected_genes = pd.Series(list(corrected_genes))
                    corrected_genes.to_csv(corrected_genes_list, index=False, header=None)

                gene_file_dl = pd.Series(selected_genes.index)
                gene_list_file_dl = gene_file_dl.to_csv(index=False, sep="\t", header=None)

                if len(gene_file_dl) < 1:
                    dl_btn = True
                else:
                    dl_btn = False
                gene_name_col1.header(" ", divider="blue")
                down_load_add_pdf = st.download_button(
                    label="Download Gene List as CSV",
                    key='download_add_btn1',
                    data=gene_list_file_dl,
                    file_name='gene_list.csv',
                    mime='text/csv',
                    disabled=dl_btn
                )

                os.remove(gene_file)
            else:
                st.sidebar.warning(
                    ":blue[Maximum number of genes is ] :red[2500. ] :blue[Please reduce number of genes in your list]")
                sys.exit(1)

            if selected_genes.shape[0] == 0:
                st.error(":blue[You did ] :red[NOT ] :blue[provide gene(s) ] :red[OR ]"
                         " :blue[gene(s) not in the database. Try other gene(s)!]")

            else:
                return selected_genes


        if uploaded_file is not None:
            try:
                if uploaded_file.name.endswith(".tsv"):
                    gene_names = pd.read_csv(uploaded_file, sep='\t')
                elif uploaded_file.name.endswith(".txt"):
                    gene_names = pd.read_csv(uploaded_file, header=None, sep='delimiter')
                elif uploaded_file.name.endswith(".csv"):
                    gene_names = pd.read_csv(uploaded_file, header=None, sep=',')

                if len(gene_names.index) <= 2500:
                    gene_names = gene_names.iloc[:, 0]
                    gene_names = gene_names.drop_duplicates(keep='first')
                    gene_names.replace(r'^\s+', '', regex=True, inplace=True)
                    gene_names.replace(r'\s+$', '', regex=True, inplace=True)
                    gene_names = gene_names.str.cat(sep=' ')
                    gene_names = organiser(gene_names, gene_name_col1)
                    gene_input = gene_name_col1.text_input(
                        ":blue[Add ] :red[More Genes]:blue[. You Can Copy Gene List Here!]",
                        key="gl_maker")
                    gene_name_col1.header(" ", divider="blue")

                    if len(gene_input) > 2500:
                        st.warning(
                            ":blue[Maximum number of genes is ] :red[2500. ] :blue[Please reduce number of genes in your list]")
                    else:
                        t_inp = organiser(gene_input, gene_name_col1)
                        t_inp = t_inp.tolist()

                        gene_names = pd.concat([gene_names, pd.Series(t_inp)])
                        gene_names.drop_duplicates(keep='first', inplace=True)
                        name_list, gene_file = name_suggestion(gene_names, df_app, name_dict, gene_name_col1,
                                                               gene_name_col_exp)

                        with open(gene_file, "r") as f:
                            gene_ids = f.read().splitlines()

                        selected_genes = df[df.index.isin(gene_ids)]
                        selected_genes.dropna(axis=0, inplace=True, how='all')
                        selected_genes.dropna(axis=1, inplace=True, how='all')

                        difference = list(set(name_list["Gene Names"]) - set(selected_genes.index.tolist()))
                        corrected_genes = list(set(gene_names) - set(selected_genes.index.tolist()))
                        corrected_genes = set(corrected_genes) - set(difference)

                        not_database_genes = uint.naming(user_session_id)[39]
                        corrected_genes_list = uint.naming(user_session_id)[40]

                        if len(difference) > 0:
                            st.warning(
                                f':blue[The following gene(s) is/are ] :red[NOT ] :blue[it this dataset: ] :red[{difference}]')
                            difference = pd.Series(difference)
                            difference.to_csv(not_database_genes, index=False, header=None)
                        elif selected_genes.shape[1] != 0 and len(difference) > 0:
                            gene_name_col1.success("All genes what provide, are in the dataset!")
                        elif selected_genes.shape[1] == 0:
                            gene_name_col1.error(":blue[Something went ] :red[wrong! ] :blue[Please add ] "
                                                 ":red[gene(s) or try other gene(s) name]")

                        if len(corrected_genes) > 0:
                            corrected_genes = pd.Series(list(corrected_genes))
                            corrected_genes.to_csv(corrected_genes_list, index=False, header=None)

                        gene_file_dl = pd.Series(selected_genes.index)
                        gene_list_file_dl = gene_file_dl.to_csv(index=False, sep="\t", header=None)

                        if len(gene_file_dl) < 1:
                            dl_btn = True
                        else:
                            dl_btn = False

                        down_load_add_pdf = st.download_button(
                            label="Download Gene List as CSV",
                            key='download_add_btn',
                            data=gene_list_file_dl,
                            file_name='gene_list.csv',
                            mime='text/csv',
                            disabled=dl_btn
                        )
                        os.remove(gene_file)

                    if len(selected_genes.index) == 0:
                        st.error(
                            ":red[Wrong file format, or no find your genes. ] :blue[Please upload another file and take care of the format.]")

                    else:
                        return selected_genes

                else:
                    st.sidebar.warning(
                        ":blue[Maximum number of genes is ] :red[2500. ] :blue[Please reduce number of genes in your list]")
                    sys.exit(1)
            except Exception as e:
                st.error(f"Error processing uploaded file: {e}")

    except Exception as e:
        st.error(f"Error in gene list app: {e}")


def gene_finder():
    """group
    Function to find a gene by name.

    Returns:
    None
    """
    try:
        gene_name = st_keyup(key="gene_name", placeholder="Type a Gene Name and Find it", label=None)
        if gene_name is not None:
            st.link_button("Search", url=f"https://www.genecards.org/Search/Keyword?queryString={gene_name}")
            st.header(" ", divider="blue")
    except Exception as e:
        st.error(f"Error in gene finder: {e}")


def app_maker(user_session_id, df, option_dataset):
    """
    Function to create an app for gene list management.

    Parameters:
    user_session_id (str): User session ID.
    df (pd.DataFrame): DataFrame containing gene expression data.
    option_dataset (str): Option for dataset.

    Returns:
    pd.DataFrame: Selected genes DataFrame.
    """
    try:
        st.warning(':blue[Avoid duplicates gene names, because the application ] :red[keeps only the first match ] :blue[, deletes others]')
        df_app, name_dict = load_gene_name_dict()

        gene_name_col1, gene_name_col2 = st.columns(2)
        gene_name_col_exp = gene_name_col2.expander(":blue[Hide Options for Mapping]", expanded=True)
        selected_genes = None

        uploaded_file = st.sidebar.file_uploader('', label_visibility="visible",
                                                   type=('txt', 'csv', 'tsv'), key="maker_uploader")
        st.sidebar.subheader(" ", divider="blue")
        gene_names = []

        if uploaded_file is not None:
            try:
                if uploaded_file.name.endswith(".tsv"):
                    gene_names = pd.read_csv(uploaded_file, sep='\t')
                elif uploaded_file.name.endswith(".txt"):
                    gene_names = pd.read_csv(uploaded_file, header=None, sep='delimiter')
                elif uploaded_file.name.endswith(".csv"):
                    gene_names = pd.read_csv(uploaded_file, header=None, sep=',')

                if len(gene_names.index) <= 2500:
                    try:
                        gene_names = gene_names.iloc[:, 0]
                        gene_names = gene_names.drop_duplicates(keep='first')
                        gene_names.replace(r'^\s+', '', regex=True, inplace=True)
                        gene_names.replace(r'\s+$', '', regex=True, inplace=True)
                        gene_names = gene_names.str.cat(sep=' ')
                        gene_names = organiser(gene_names, gene_name_col1)

                        gene_names.drop_duplicates(keep='first', inplace=True)

                        name_list, gene_file = name_suggestion(gene_names, df_app, name_dict, gene_name_col1,
                                                               gene_name_col_exp)

                        with open(gene_file, "r") as f:
                            gene_ids = f.read().splitlines()

                        selected_genes = df[df.index.isin(gene_ids)]
                        selected_genes.dropna(axis=0, inplace=True, how='all')
                        selected_genes.dropna(axis=1, inplace=True, how='all')

                        difference = list(set(name_list["Gene Names"]) - set(selected_genes.index.tolist()))
                        corrected_genes = list(set(gene_names) - set(selected_genes.index.tolist()))
                        corrected_genes = set(corrected_genes) - set(difference)

                        not_database_genes = uint.naming(user_session_id)[39]
                        corrected_genes_list = uint.naming(user_session_id)[40]

                        if len(difference) > 0:
                            st.warning(
                                f':blue[The following gene(s) is/are ] :red[NOT ] :blue[it this dataset: ] :red[{difference}]')
                            difference = pd.Series(difference)
                            difference.to_csv(not_database_genes, index=False, header=None)
                        elif selected_genes.shape[1] != 0 and len(difference) > 0:
                            gene_name_col1.success("All genes what provide, are in the dataset!")
                        elif selected_genes.shape[1] == 0:
                            gene_name_col1.error(":blue[Something went ] :red[wrong! ] :blue[Please add ] "
                                                 ":red[gene(s) or try other gene(s) name]")

                        if len(corrected_genes) > 0:
                            corrected_genes = pd.Series(list(corrected_genes))
                            corrected_genes.to_csv(corrected_genes_list, index=False, header=None)

                        gene_file_dl = pd.Series(selected_genes.index)
                        gene_list_file_dl = gene_file_dl.to_csv(index=False, sep="\t", header=None)

                        if len(gene_file_dl) < 1:
                            dl_btn = True
                        else:
                            dl_btn = False

                        down_load_add_pdf = st.download_button(
                            label="Download Gene List as CSV",
                            key='download_add_btn_maker',
                            data=gene_list_file_dl,
                            file_name='gene_list.csv',
                            mime='text/csv',
                            disabled=dl_btn
                        )
                        os.remove(gene_file)

                        if len(selected_genes.index) == 0:
                            st.error(
                                ":red[Wrong file format, or no find your genes. ] :blue[Please upload another file and take care of the format.]")
                        else:
                            return selected_genes

                    except Exception as e:
                        st.error(f"Error processing gene names: {e}")
                else:
                    st.sidebar.warning(
                        ":blue[Maximum number of genes is ] :red[2500. ] :blue[Please reduce number of genes in your list]")
                    sys.exit(1)

            except Exception as e:
                st.error(f"Error reading uploaded file: {e}")

    except Exception as e:
        st.error(f"Error in app maker: {e}")


def cluster_app(user_session_id, user_gene_cluster_path, df, option_dataset, uploaded_file):
    """
    Function to handle gene clustering with uploaded file.

    Parameters:
    user_session_id (str): User session ID.
    user_gene_cluster_path (str): Path to save user gene cluster file.
    df (pd.DataFrame): DataFrame containing gene expression data.
    option_dataset (str): Option for dataset.
    uploaded_file (UploadedFile): Uploaded file containing gene names with clusters.

    Returns:
    pd.DataFrame: Selected genes DataFrame.
    """
    try:
        uploaded_file.seek(0)
        st.warning(
            ':blue[If you want to ] :red[specify your own gene clusters] :blue[, please pay attention to the file ]'
            ':blue[structure. The ] :red[first column contains the gene names] :blue[and ] :red[the other columns ] '
            ':red[the clusters. ] :blue[Cluster column names ] :red[MUST ] :blue[contain one of the words ] '
            ':red["cluster", "group" or "class". ] :blue[For example, Gene_Clusters.] :red[IMPORTANT: ]'
            ':blue[Avoid duplicates gene names, because the application ] :red[keeps only the first match ] :blue[, deletes others]')
        final_df = None
        selected_cluster_column = None
        df_app, name_dict = load_gene_name_dict()
        gene_name_col1, gene_name_col2 = st.columns(2)

        gene_name_col2.caption("", help=dsc.gl_text())
        gene_name_col_exp = gene_name_col2.expander(":blue[Hide Options for Mapping]", expanded=True)
        selected_genes = None

        gene_names = []

        if uploaded_file is None:
            st.warning(":blue[Please ] :red[upload your gene list with clusters! ]")
        else:
            try:
                if uploaded_file.name.endswith(".tsv"):
                    gene_names = pd.read_csv(uploaded_file, sep='\t', header=None)
                elif uploaded_file.name.endswith(".txt"):
                    gene_names = pd.read_csv(uploaded_file, header=None,
                                             sep='delimiter')
                elif uploaded_file.name.endswith(".csv"):
                    gene_names = pd.read_csv(uploaded_file, header=None, sep=',')

                gene_names.columns = gene_names.iloc[0]
                gene_names.columns = gene_names.columns.fillna("Unknown")
                gene_names = gene_names[1:].map(lambda x: x.upper() if isinstance(x, str) else x)
                gene_names.reset_index(drop=True, inplace=True)

                cluster_columns = [col for col in gene_names.columns if
                                   any(sub in col.lower() for sub in ["class", "group", "cluster"])]
                if cluster_columns:
                    selected_cluster_column = gene_name_col1.selectbox(
                        ":red[Choose a Cluster-Related Column ] :blue[to Use:]",
                        options=cluster_columns, key="user_select_clust")

                    index_values = gene_names.iloc[:, 0]
                    selected_data = gene_names[[selected_cluster_column]]
                    merge_gene_names = selected_data.set_index(index_values)

                    if len(gene_names) <= 2500:
                        gene_names_series = gene_names.iloc[:, 0].drop_duplicates().str.strip()
                        concatenated_gene_names = gene_names_series.str.cat(sep=' ')
                        organised_gene_names = pd.Series(organiser(concatenated_gene_names, gene_name_col1))
                        final_df = merge_gene_names[merge_gene_names.index.isin(organised_gene_names)]
                        final_df.dropna(how='any', inplace=True)
                        final_df.reset_index(inplace=True)
                    else:
                        st.sidebar.warning(
                            ":blue[Maximum number of genes is ] :red[2500. ] :blue[Please reduce number of genes in your ] "
                            ":blue[list]")
                        sys.exit(1)

                    updated_df, gene_file = name_suggestion_cluster(final_df, df_app, name_dict, gene_name_col1,
                                                                    gene_name_col_exp)

                    updated_df.drop_duplicates(subset=['Samples'], keep='first', inplace=True)
                    selected_genes = df[df.index.isin(updated_df['Samples'])]
                    selected_genes = selected_genes.dropna(axis=1)
                    selected_genes.dropna(axis=0, inplace=True, how='all')
                    selected_genes.dropna(axis=1, inplace=True, how='all')

                    difference = list(set(updated_df["Samples"]) - set(selected_genes.index.tolist()))
                    corrected_genes = list(set(gene_names.iloc[:, 0]) - set(selected_genes.index.tolist()))
                    corrected_genes = set(corrected_genes) - set(difference)

                    not_database_genes = uint.naming(user_session_id)[39]
                    corrected_genes_list = uint.naming(user_session_id)[40]

                    if len(difference) > 0:
                        st.warning(
                            f':blue[The following gene(s) is/are ] :red[NOT ] :blue[it this dataset: ] :red[{difference}]')
                        difference = pd.Series(difference)
                        difference.to_csv(not_database_genes, index=False, header=None)
                    elif selected_genes.shape[1] == 0:
                        gene_name_col1.error(
                            ":blue[Something went ] :red[wrong! ] :blue[Please add ] :red[gene(s) or try other gene(s) name]")

                    if len(corrected_genes) > 0:
                        corrected_genes = pd.Series(list(corrected_genes))
                        corrected_genes.to_csv(corrected_genes_list, index=False, header=None)

                    gene_file_dl = pd.Series(selected_genes.index)
                    gene_list_file_dl = gene_file_dl.to_csv(index=False, sep="\t", header=None)

                    dl_btn = len(gene_file_dl) < 1

                    down_load_add_pdf = st.download_button(
                        label="Download Gene List as CSV",
                        key='download_add_btn_maker',
                        data=gene_list_file_dl,
                        file_name='gene_list.csv',
                        mime='text/csv',
                        disabled=dl_btn
                    )

                    os.remove(gene_file)

                else:
                    st.warning(":red[NO ] :blue[cluster-related column found in your gene list.]")
                    sys.exit(1)

                if len(selected_genes.index) == 0:
                    st.error(
                        ":red[Wrong file format, or no find your genes. ] :blue[Please upload another file and take care of the format.]")
                else:
                    updated_df = updated_df[updated_df['Samples'].isin(selected_genes.index)]
                    updated_df.set_index("Samples", inplace=True)
                    updated_df.reset_index(inplace=True)
                    updated_df.to_csv(user_gene_cluster_path, index=False)

                    return selected_genes

            except Exception as e:
                st.error(
                    f":blue[Something Went] :red[wrong! ] :blue[Try Another File, or Choose Minimum One Cluster, If ]"
                    f":blue[ You Use Your Own Clustering Results!]")
    except Exception as e:
        st.error(f"Error in cluster app: {e}")


def cluster_checker(uploaded_file):
    """
    Check if the uploaded file contains cluster-related columns.

    Parameters:
    uploaded_file (UploadedFile): Uploaded file containing gene names with clusters.

    Returns:
    bool: True if cluster-related columns are found, False otherwise.
    """
    try:
        global gene_names
        cluster_columns = False
        user_cluster = False

        if uploaded_file is not None:
            try:
                uploaded_file.seek(0)

                if uploaded_file.size == 0:
                    st.error("Uploaded file is empty.")
                    return False

                if uploaded_file.name.endswith(".tsv"):
                    gene_names = pd.read_csv(uploaded_file, header=None, sep='\t')
                elif uploaded_file.name.endswith(".txt"):
                    gene_names = pd.read_csv(uploaded_file, header=None, sep='\t')
                elif uploaded_file.name.endswith(".csv"):
                    gene_names = pd.read_csv(uploaded_file, header=None, sep=',')

                if gene_names.empty:
                    st.error("No data found in the file.")
                    return False
                gene_names.columns = gene_names.iloc[0]
                gene_names.columns = gene_names.columns.fillna("Unknown")
                gene_names = gene_names[1:].map(lambda x: x.upper() if isinstance(x, str) else x)
                gene_names.reset_index(drop=True, inplace=True)

                cluster_columns = [col for col in gene_names.columns if
                                   any(sub in col.lower() for sub in ["class", "group", "cluster"])]

                user_cluster = bool(cluster_columns)

            except pd.errors.EmptyDataError:
                st.error("No columns to parse from file.")
                return False
            except Exception as e:
                st.error(f"An error occurred: {e}")
                return False
        return user_cluster
    except Exception as e:
        st.error(f"Error in cluster checker: {e}")
        return False


def col_maker(selected_genes):
    # Print out the selected genes
    name_list = pd.DataFrame(selected_genes.index)

    if len(name_list) != 0:
        st.subheader(":blue[The Chosen Genes]")

        # Determine the number of columns
        num_cols = min(max(math.ceil(len(name_list) / 10), 1), 4)
        gene_list_cols = st.columns(num_cols)

        # Determine the length of each part
        n = math.ceil(len(name_list) / num_cols)

        # Add each part to a column
        for i in range(num_cols):
            start = i * n
            end = (i + 1) * n if i < num_cols - 1 else None
            df = pd.DataFrame(name_list[start:end])
            gene_list_cols[i].dataframe(df, use_container_width=True, hide_index=True)
