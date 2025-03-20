import pandas as pd
from Scripts.r_code_settings import file_copy
from streamlit_lottie import st_lottie_spinner
import json
import streamlit as st
from multiprocessing import Process, Queue
import Scripts.descriptions as dsc
from Scripts.rpy2_heatmap_plots import gene_ont_all_cluster_kegg, gene_ont_all_gl_kegg, gene_ont_all_cluster_go, \
    gene_ont_all_gl_go


def run_analysis(target_func, args, queue):
    try:
        target_func(*args)
        queue.put(True)
    except Exception as e:
        queue.put(False)
        print(e)


@st.fragment
def gene_ont(user_id, dataset_option, ht_row_cluster_path, df_index, gene_list, gene_clust_chb):
    """
    Perform Gene Ontology analysis.

    Parameters:
        user_id (str): ID of the user.
        dataset_option (str): Selected dataset option.
        ht_row_cluster_path (str): Path to the row cluster data.
        df_index (list): List of data frame indices.
        gene_list (list): List of genes.
        gene_clust_chb (bool): toggle for gene clustering.

    Returns:
        None
    """
    try:
        sub_ont = 'None'
        go_form = st.form("go_form")
        go_col = go_form.columns(2)
        gene_list_options = ["Gene List"]

        if gene_clust_chb:
            gene_list_options = ["Gene List", "Cluster"]

        gs_opt = go_col[0].selectbox(
            ":blue[Choose the ] :red[Gene Set]",
            options=gene_list_options,
            key="gene_list",
            help=dsc.go_geneset()
        )

        go_fun = go_col[0].selectbox(
            ":blue[Choose the Statistical Method for ] :red[Enrichment Analysis]",
            options=("enrichGO", "enrichKEGG"),
            key="go_fun",
            help=dsc.go_fun(),
            index=0
        )

        if go_fun == 'enrichGO':
            sub_ont = go_col[0].selectbox(
                ":blue[Choose the ] :red[Sub Ontology]",
                options=("BP", "MF", "CC", "ALL"),
                key="SubOnt",
                help=dsc.sub_ont(),
                index=0
            )

        p_value_cutoff = go_col[1].select_slider(
            ":red[P Value ] :blue[Cutoff]",
            options=(0.001, 0.01, 0.05),
            key='p_value_cutoff',
            help=dsc.p_value_cutoff(),
            value=0.05
        )

        q_value_cutoff = go_col[1].slider(
            ":red[Q Value ] :blue[Cutoff]",
            min_value=0.1,
            max_value=0.5,
            value=0.2,
            step=0.1,
            key='q_value_cutoff',
            help=dsc.q_value_cutoff()
        )

        corr_meth = go_col[1].selectbox(
            ":blue[Select the ] :red[Correlation Method]",
            key='corr_meth_go',
            help=dsc.corr_meth_go(),
            options=("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
            index=6
        )

        max_gene = go_col[1].number_input(
            ":blue[Select the ] :red[Maximum Number of Gene]",
            key='max_gene',
            help=dsc.max_gene(),
            min_value=100,
            max_value=1000,
            value=500,
            step=1
        )

        max_category = go_col[1].slider(
            ":blue[Select the ] :red[Maximum Number of Pathways ] :blue[for Plot]",
            key='max_category',
            help=dsc.max_category_description(),
            min_value=1,
            max_value=25,
            value=10,
            step=1
        )

        order = go_col[1].selectbox(
            ":blue[Select the ] :red[Base of Ordering] :blue[(only for Gene List)]",
            key='order',
            help=dsc.order_by_description(),
            options=["GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "Count"],
            index=0
        )

        title =  go_col[1].text_input(":blue[Choose the ] :red[Title]", key="title", value="Enriched Pathway Terms",
                                      help=dsc.gene_ontology_title())

        st.warning(":blue[This process can take up to a ] :red[few minutes. ] :blue[Please be patient.]")

        if gs_opt == "Gene List":
            start_btn = go_form.form_submit_button("Start Analysing")
            if start_btn:
                lottie_url = "style_items/surv_sig_lottie.json"
                with open(lottie_url, "r") as f:
                    animation = json.load(f)
                try:
                    # Use a Queue to check if the analysis was successful
                    queue = Queue()
                    analysis_success = False

                    with st_lottie_spinner(animation, height=200, width=200, key="plot_spinner"):
                        if go_fun == 'enrichGO':
                            r_process = Process(
                                target=run_analysis,
                                args=(gene_ont_all_gl_go,
                                      (user_id, dataset_option, ht_row_cluster_path, df_index, gene_list, sub_ont, go_fun,
                                       max_gene, corr_meth, q_value_cutoff, p_value_cutoff, title, max_category,
                                       order), queue)
                            )
                        else:
                            r_process = Process(
                                target=run_analysis,
                                args=(gene_ont_all_gl_kegg,
                                      (user_id, dataset_option, ht_row_cluster_path, df_index, gene_list, sub_ont, go_fun,
                                       max_gene, corr_meth, q_value_cutoff, p_value_cutoff, title, max_category,
                                       order), queue)
                            )
                        r_process.start()
                        r_process.join()
                        analysis_success = queue.get()

                    if analysis_success:
                        result_path_go = f"result_data/result_data_{user_id}/gene_ont_result_df_{user_id}.csv"
                        with open(f"result_data/result_data_{user_id}/gene_ont_{user_id}.pdf", "rb") as pdf_file:
                            pdfbyte = pdf_file.read()

                        st.download_button(
                            label="Download DotPlot",
                            data=pdfbyte,
                            file_name="gene_ont.pdf",
                            mime='application/octet-stream'
                        )

                        st.header(":blue[Gene Ontology Plot]", divider="blue")
                        img_cols = st.columns(3)
                        img_cols[1].image(f"result_data/result_data_{user_id}/gene_ont_{user_id}.png")

                        st.header(":blue[Gene Ontology Result Table]", divider="blue",
                                  help='If the value is too small, it will appear as 0 in the table.')
                        result_go = pd.read_csv(result_path_go, index_col=0)
                        st.dataframe(result_go, use_container_width=True)
                    else:
                        st.warning(":red[No pathways found! Try other enrichment analysis methods or settings.]")

                except Exception as e:
                    st.error(f"Error during analysis: {e}")

        elif gs_opt == "Cluster":
            start_btn = go_form.form_submit_button("Start Analysing")
            if start_btn:
                lottie_url = "style_items/surv_sig_lottie.json"
                with open(lottie_url, "r") as f:
                    animation = json.load(f)
                try:
                    # Use a Queue to check if the analysis was successful
                    queue = Queue()
                    analysis_success = False

                    with st_lottie_spinner(animation, height=200, width=200, key="plot_spinner"):
                        if go_fun == 'enrichGO':
                            r_process = Process(
                                target=run_analysis,
                                args=(gene_ont_all_cluster_go,
                                      (user_id, dataset_option, ht_row_cluster_path, df_index, gene_list, sub_ont, go_fun,
                                       max_gene, corr_meth, q_value_cutoff, p_value_cutoff, title, max_category), queue)
                            )
                        else:
                            r_process = Process(
                                target=run_analysis,
                                args=(gene_ont_all_cluster_kegg,
                                      (user_id, dataset_option, ht_row_cluster_path, df_index, gene_list, sub_ont, go_fun,
                                       max_gene, corr_meth, q_value_cutoff, p_value_cutoff, title, max_category), queue)
                            )
                        r_process.start()
                        r_process.join()
                        analysis_success = queue.get()
                        print(analysis_success)

                    if analysis_success:
                        result_path_go = f"result_data/result_data_{user_id}/gene_ont_result_df_{user_id}.csv"
                        with open(f"result_data/result_data_{user_id}/gene_ont_{user_id}.pdf", "rb") as pdf_file:
                            pdfbyte = pdf_file.read()

                        st.download_button(
                            label="Download Result",
                            data=pdfbyte,
                            file_name="gene_ont.pdf",
                            mime='application/octet-stream'
                        )

                        st.header(":blue[Gene Ontology Plot]", divider="blue")
                        img_cols = st.columns(3)
                        img_cols[1].image(f"result_data/result_data_{user_id}/gene_ont_{user_id}.png")

                        st.header(":blue[Gene Ontology Result Table]", divider="blue",
                                  help='If the value is too small, it will appear as 0 in the table.')
                        result_go = pd.read_csv(result_path_go, index_col=0)
                        st.dataframe(result_go, use_container_width=True)
                    else:
                        st.warning(":red[No pathways found! Try other enrichment analysis methods or settings.]")


                except Exception as e:
                    st.error(f"Error during analysis: {e}")
        file_copy(user_id)

    except Exception as e:
        st.error(f"Error in Gene Ontology analysis: {e}")
