import streamlit as st

# Page configuration
st.set_page_config(
    page_title="SurvSig",
    page_icon="style_items/favicon.svg",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        'Report a bug': "mailto:kolos.nemes@hcemm.eu",
        'About': "https://www.hcemm.eu/teams/genomic-instability-and-cancer/cancer-genomics-and-epigenetics-core-group/"
    }
)
from Scripts.style import style, warning_tab, landing_animation, display_free_use_message

landing_animation()

import time
time.sleep(2.5) # Sleep for 3 seconds


import pandas as pd
import numpy as np
import Scripts.descriptions as dscript
import Scripts.data_loader as dl
import Scripts.user_interaction as uint
import Scripts.gene_list as gl
from Scripts.data_preprocessing import data_preprocess, post_process
from Scripts.dimension_reduction import dim_reductor, dim_reductor_gene
from Scripts.clustering import clustering, clustering_genes
from Scripts.multivariate_analysis import mv_db
from Scripts.plotly_plots import plotly_heatmap
from Scripts.subtype_selector import subtype_select, subtypes_selector_tcga, subtypes_selector_tcga_sg, subtype_select_sg
from Scripts.chi_sqrt_test import chi_sqrt_test
from Scripts.r_code_settings import downloader, file_copy
from Scripts.ssgsea import ssgsea_analysis, ssgsea_compare
from Scripts.gene_ontology import gene_ont
from Scripts.single_gene import sg_analysis_surv
from Scripts.gene_compare import gene_comp
from Scripts.gene_set_finder import gene_finder
from Scripts.help_tab import help_tab
from Scripts.survival_info import survival_info
import copy
import io
import time
from Scripts.rpy2_heatmap_plots import r_libraries
from Scripts.cookies import  display_cookie_banner


style()
display_cookie_banner()

user_session_id = uint.get_session_id()
file_name_list = uint.naming(user_session_id)


def main_app():
    try:
            # ########## Main part of app ########## #

            uint.del_files_py(24000)

            ht_expression_path = file_name_list[0]
            selected_genes_name = file_name_list[0]
            selected_genes_name_nmf = file_name_list[1]
            selected_genes_dm_name = file_name_list[3]
            column_cluster_path = file_name_list[4]
            ht_row_cluster_path = file_name_list[6]
            ht_top_annotation_path = file_name_list[2]
            ht_png = file_name_list[7]
            ht_pdf = file_name_list[9]
            surv_dataframe_path = file_name_list[5]
            surv_png = file_name_list[8]
            surv_pdf = file_name_list[10]
            surv_df = surv_dataframe_path
            multi_db_name = file_name_list[11]
            cont_name = file_name_list[12]
            comp_name = file_name_list[13]
            scatter_2d_name = file_name_list[14]
            scatter_3d_name = file_name_list[15]
            nmf_cons_name = file_name_list[16]
            sub_hm_name = file_name_list[17]
            chi_hm_name = file_name_list[18]
            chi_bp_name = file_name_list[19]
            scatter_2d_name_gene = file_name_list[20]
            scatter_3d_name_gene = file_name_list[21]
            multi_info_name = file_name_list[22]
            forest_pdf = file_name_list[23]
            forest_png = file_name_list[24]
            gene_sig_name = file_name_list[25]
            ssgsea_violin_box = file_name_list[26]
            ssgsea_result_path = file_name_list[28]
            selected_genes_dm_name_gene = file_name_list[31]
            user_gene_cluster_path = file_name_list[32]
            violoin_comapre_violin_path = file_name_list[34]


            def delete_files(del_time):
                # Placeholder functions for deletion
                uint.del_files_py_user(del_time)
                uint.del_file_py_zip(del_time)

            # Main app logic
            dataset_type = st.sidebar.selectbox(
                ':blue[Select the ] :red[Collection]', index=0,
                options=("Neuroendocrine Lung Cancer", "TCGA"), key="type_data", help=dscript.exp_type())

            if dataset_type == 'Neuroendocrine Lung Cancer':
                option_dataset = st.sidebar.selectbox(
                    ':blue[Select a ] :red[Dataset]', index=0,
                    options=("George-SCLC", "Anish-SCLC", "Jiang-SCLC", "George-LCNEC", "Fernandez-Carcinoid",
                             "Alcala-Carcinoid", 'Liu-SCLC', 'Rousseaux-Mixed(non-NE)'),
                    key="opt_data", help=dscript.cohort_description())
                st.sidebar.header(" ", divider="blue")
            else:
                option_dataset = dataset_type


            if 'last_option_dataset' not in st.session_state or not np.array_equal(
                    st.session_state['last_option_dataset'], dataset_type):
                delete_files(1)
            st.session_state['last_option_dataset'] = copy.deepcopy(dataset_type)

            selected_genes = None
            selected_genes1 = None
            option_dm = None
            clusters_db = None
            clust_nmf = None
            selected_genes_dim_red = None
            clusters_db_gene = None
            gene_clust_chb = False
            uploaded_file = None

            r_libraries() # Call after streamlit, to go to the main thread


            if option_dataset == "TCGA":
                start = time.time()

                tab1, tab2, tab4, tab6, tab3, tab9, tab8, tab7, tab5, tab11, tab10 = st.tabs(
                    ["Clustering :robot_face:", "Genes Handling 🧬",
                     "Subtypes 🔍", "ssGSEA Survival  🔬",
                     "ssGSEA Compare  :chart_with_upwards_trend:",
                     "Gene Compare ⚗️",
                     "Gene Survival :petri_dish:",
                     "Ontology :test_tube:",
                     "Multivariate & Chi" + '$$^2$$ ' + ":bar_chart:",
                     "Gene Set Finder 🧩",
                     "Help 	:question:"])

                with tab10:
                    help_tab()

                with tab1:
                    tcga_df, sample_type, cancer_types_list, tcga_clinics, tcga_subtypes, tcga_surv = dl.data_read()

                    target_cancer_exp_df, cancer_opt = dl.cancer_type_separator(tcga_df, sample_type, cancer_types_list,
                                                                                gene_sig_exp=None, key="main")

                    tcga_surv, tcga_clinics = dl.tcga_surv_time_corrector(tcga_surv, target_cancer_exp_df, tcga_clinics)

                    tcga_surv, tcga_clinics = dl.tcga_survival_selector(tcga_surv, tcga_clinics.copy())

                    if len(target_cancer_exp_df.columns) > 250:
                        st.sidebar.warning(
                            ":blue[You work with a ] :red[high number ] :blue[of samples, which increases the running time!]")

                    del tcga_df

                with tab2:

                    uploader_container = st.sidebar.empty()
                    gene_set_finder_container = st.sidebar.empty()
                    demo_conatainer_container = st.sidebar.empty()

                    feature_chb = gene_set_finder_container.toggle(":red[Find a ] :blue[Gene Set]",
                                                                   key="feature_1", help=dscript.dataset_ana_intro())

                    uploaded_file = uploader_container.file_uploader(':blue[Upload your ] :red[gene list]',
                                                                     label_visibility="visible",
                                                                     type=('txt', 'csv', 'tsv'),
                                                                     key="gene_uploader")

                    demo_toggle = demo_conatainer_container.toggle(":blue[Start ] :red[Demo]", key="demo_toggle1",
                                                                   help=dscript.demo_gene_list_help())

                    # Handle file input
                    if demo_toggle:
                        # Simulate an uploaded file-like object for the demo file
                        with open("source_data/NE50-signatures.csv", "rb") as demo_file:
                            demo_content = demo_file.read()
                        # Create a BytesIO object to simulate file upload
                        uploaded_file = io.BytesIO(demo_content)
                        # Mimic the UploadedFile attributes
                        uploaded_file.name = "NE50-signatures.csv"
                        uploaded_file.size = len(demo_content)

                    st.sidebar.subheader(" ", divider="blue")

                    if feature_chb:
                        uploader_container.empty()
                        demo_conatainer_container.empty()

                    if demo_toggle:
                        uploader_container.empty()
                        gene_set_finder_container.empty()
                        st.toast("🖱️ :blue[Click to the ] :red[Draw Heatmap and Survival Plot ] :blue[Button!] 🖱️")

                    with tab11:
                        if not feature_chb:
                            st.warning(
                                ':blue[Please Swich On the: ] :red["Find a Gene Set" ] :blue[Toggle on the Sidebar]')

                    if feature_chb:
                        with tab11:
                            selected_genes = gene_finder(option_dataset, target_cancer_exp_df.copy(), gene_sig_name,
                                                         user_session_id, tcga_subtypes.copy(), tcga_clinics.copy())

                            if selected_genes is not None:
                                pass
                            else:
                                with tab1:
                                    st.warning(
                                        ":red[No ] :blue[ any significant genes!] :blue[Please ] :red[change the settings of analisys! ]"
                                        ":blue[ If you use gene list please ] :red[upload one!]")
                                with tab2:
                                    st.warning(
                                        ":red[No ] :blue[ any significant genes!] :blue[Please ] :red[change the settings of analisys! ]"
                                        ":blue[ If you use gene list please ] :red[upload one!]")
                                with tab3:
                                    st.warning(
                                        ":red[No ] :blue[ any significant genes!] :blue[Please ] :red[change the settings of analisys! ]"
                                        ":blue[ If you use gene list please ] :red[upload one!]")
                                with tab4:
                                    st.warning(
                                        ":red[No ] :blue[ any significant genes!] :blue[Please ] :red[change the settings of analisys! ]"
                                        ":blue[ If you use gene list please ] :red[upload one!]")
                                with tab5:
                                    st.warning(
                                        ":red[No ] :blue[ any significant genes!] :blue[Please ] :red[change the settings of analisys! ]"
                                        ":blue[ If you use gene list please ] :red[upload one!]")
                                with tab6:
                                    st.warning(
                                        ":red[No ] :blue[ any significant genes!] :blue[Please ] :red[change the settings of analisys! ]"
                                        ":blue[ If you use gene list please ] :red[upload one!]")
                                with tab7:
                                    st.warning(
                                        ":red[No ] :blue[ any significant genes!] :blue[Please ] :red[change the settings of analisys! ]"
                                        ":blue[ If you use gene list please ] :red[upload one!]")
                                with tab8:
                                    st.warning(
                                        ":red[No ] :blue[ any significant genes!] :blue[Please ] :red[change the settings of analisys! ]"
                                        ":blue[ If you use gene list please ] :red[upload one!]")
                                with tab9:
                                    st.warning(
                                        ":red[No ] :blue[ any significant genes!] :blue[Please ] :red[change the settings of analisys! ]"
                                        ":blue[ If you use gene list please ] :red[upload one!]")
                                with tab10:
                                    st.warning(
                                        ":red[No ] :blue[ any significant genes!] :blue[Please ] :red[change the settings of analisys! ]"
                                        ":blue[ If you use gene list please ] :red[upload one!]")
                                with tab11:
                                    st.warning(
                                        ":red[No ] :blue[ any significant genes!] :blue[Please ] :red[change the settings of analisys! ]"
                                        ":blue[ If you use gene list please ] :red[upload one!]")
                                st.stop()

                        gl.col_maker(selected_genes)

                        cluster_gene_list_toggle = False

                        if uploaded_file:
                            feature_chb = False

                    else:
                        cluster_gene_list_toggle = False

                        cluster_check = gl.cluster_checker(uploaded_file)

                        if demo_toggle:
                            demo_cluster = True
                            cluster_check = True
                        else:
                            demo_cluster = False

                        if uploaded_file is not None and cluster_check:
                            cluster_gene_list_toggle = st.toggle(":blue[Use Gene List with ] :red[Clusters]",
                                                                 key="cluster_toggle_tcga", value=demo_cluster)
                            if cluster_gene_list_toggle:
                                selected_genes = gl.cluster_app(user_session_id, user_gene_cluster_path,
                                                                option_dataset=option_dataset, df=target_cancer_exp_df,
                                                                uploaded_file=uploaded_file)
                            else:
                                selected_genes = gl.app(user_session_id, option_dataset=option_dataset,
                                                        df=target_cancer_exp_df,
                                                        uploaded_file=uploaded_file)

                        else:
                            selected_genes = gl.app(user_session_id, option_dataset=option_dataset,
                                                    df=target_cancer_exp_df,
                                                    uploaded_file=uploaded_file)

                    if selected_genes is not None:
                        # Check if min and max are the same
                        all_same = selected_genes.values.min() == selected_genes.values.max()

                        if all_same:
                            st.warning(":red[All genes expression values are same!]")
                            with tab1:
                                if all_same:
                                    st.warning(":red[All genes expression values are same!]")
                            st.stop()
                        else:
                            pass

                with tab1:
                    if selected_genes is None:
                        dscript.gene_list_warning()
                        display_free_use_message()
                        warning_tab()
                        img_container = st.container(key='main_img_NE')
                        img_container.image("style_items/main_fig.svg", use_container_width=False, width=750)

                with tab4:
                    if selected_genes is not None and selected_genes.shape[0] >= 1:
                        if selected_genes.shape[0] > 2500:
                            st.sidebar.error(":blue[Too many genes. Please ] :red[reduce ] :blue[the number of genes!]")
                            exit()  # Halt the script execution without causing a Streamlit error
                        elif 500 <= selected_genes.shape[0] <= 2500:
                            st.sidebar.warning(":blue[Too many genes] :red[slow down] :blue[application runtimes]")
                            filtered_target_cancer, tcga_anno = subtypes_selector_tcga(selected_genes,
                                                                                       tcga_clinics.copy(),
                                                                                       tcga_subtypes,
                                                                                       ht_top_annotation_path,
                                                                                       key='main')
                        else:
                            filtered_target_cancer, tcga_anno = subtypes_selector_tcga(selected_genes,
                                                                                       tcga_clinics.copy(),
                                                                                       tcga_subtypes,
                                                                                       ht_top_annotation_path,
                                                                                       key='main')
                    else:
                        filtered_target_cancer = pd.DataFrame()

                with tab5:
                    if selected_genes is None:
                        dscript.gene_list_warning()
                    if len(filtered_target_cancer.columns) < 8 and selected_genes is not None:
                        st.warning(":blue[The number of samples should be at ] :red[least 8.]")
                    warning_tab()
                with tab3:
                    if selected_genes is None:
                        dscript.gene_list_warning()
                    if len(filtered_target_cancer.columns) < 8 and selected_genes is not None:
                        st.warning(":blue[The number of samples should be at ] :red[least 8.]")
                    warning_tab()
                with tab6:
                    if selected_genes is None:
                        dscript.gene_list_warning()
                    if len(filtered_target_cancer.columns) < 8 and selected_genes is not None:
                        st.warning(":blue[The number of samples should be at ] :red[least 8.]")
                    warning_tab()
                with tab7:
                    if selected_genes is None:
                        dscript.gene_list_warning()
                    if len(filtered_target_cancer.columns) < 8 and selected_genes is not None:
                        st.warning(":blue[The number of samples should be at ] :red[least 8.]")
                    warning_tab()
                with tab1:
                    if selected_genes is not None:
                        selected_genes1, selected_genes_nmf, df_names, random_state_samples, df_names_gene = \
                            post_process(filtered_target_cancer, selected_genes_name, selected_genes_name_nmf,
                                         user_session_id)
                        if selected_genes1 is None and filtered_target_cancer.shape[1] > 8:
                            st.error(
                                ":blue[Gene(s) ] :red[not in the database ] :blue[or you did] :red[not select genes. ] "
                                ":blue[You can also try another ] "
                                ":blue[database or switch off Z-scoring]")
                        else:
                            pass

                        if selected_genes1 is not None:
                            gene_clust_exp = st.sidebar.expander(':blue[Gene Clustering]')
                            sample_clust_exp = st.sidebar.expander(':blue[Sample Clustering]', expanded=True)
                            gene_cont = gene_clust_exp.container()
                            gene_clust_col1, gene_clust_col2 = st.sidebar.columns(2)
                            adv_exp = st.expander(":blue[Advanced Settings]")
                            adv_exp_gene = st.expander(":blue[Advanced Settings of Gene Clustering]")
                            adv_col1, adv_col2 = adv_exp.columns(2)
                            adv_col1_gene, adv_col2_gene = adv_exp_gene.columns(2)
                            e = st.columns(2, gap="large")

                            gene_fake_cluster = pd.DataFrame()
                            gene_fake_cluster["Names"] = df_names_gene
                            gene_fake_cluster["Cluster"] = "Genes"
                            gene_fake_cluster.to_csv(file_name_list[6], index=False)

                            if selected_genes1.shape[1] >= 8:

                                if cluster_gene_list_toggle:
                                    gene_clust_use = True
                                else:
                                    gene_clust_use = False

                                gene_clust_chb = gene_cont.toggle(":red[Clustering ] :blue[of Genes]",
                                                                    key="gene_clust_chb", value=gene_clust_use)

                                if gene_clust_chb:
                                    ht_row_cluster_path = file_name_list[6]

                                    gene_selected_genes1 = selected_genes1.copy()
                                    gene_selected_genes1 = gene_selected_genes1.transpose()

                                    gene_selected_genes_nmf = selected_genes_nmf.copy()
                                    gene_selected_genes_nmf = gene_selected_genes_nmf.transpose()

                                    selected_genes_dim_red_gene, option_dm_gene, lat_vec_gene = dim_reductor_gene(
                                        gene_selected_genes1,
                                        gene_selected_genes_nmf,
                                        random_state_samples,
                                        df_names_gene,
                                        selected_genes_dm_name_gene,
                                        user_session_id,
                                        selected_genes_name_nmf,
                                        ht_row_cluster_path,
                                        option_dataset,
                                        adv_col1_gene,
                                        adv_col2_gene,
                                        ht_png,
                                        ht_pdf, surv_png,
                                        surv_pdf,
                                        gene_clust_col1,
                                        clinic_df=tcga_surv,
                                        surv_dataframe_path=surv_df,
                                        gene_exp=gene_cont)

                                    clusters_db_gene = clustering_genes(selected_genes_dim_red_gene, selected_genes1,
                                                                        df_names_gene,
                                                                        ht_row_cluster_path, lat_vec_gene,
                                                                        random_state_samples, adv_col1_gene,
                                                                        ht_expression_path, gene_clust_col2,
                                                                        scatter_2d_name_gene,
                                                                        scatter_3d_name_gene, gene_exp=gene_cont,
                                                                        adv_col2_gene=adv_col2_gene,
                                                                        option_dm=option_dm,
                                                                        user_session_id=user_session_id,
                                                                        column_cluster_path=ht_row_cluster_path,
                                                                        selected_genes_dm_path=selected_genes_dm_name_gene,
                                                                        user_gene_cluster_path=user_gene_cluster_path,
                                                                        cluster_toggle=cluster_gene_list_toggle,
                                                                        chb_gene_clust=gene_clust_chb, e=e)

                                    if clusters_db_gene is not None:
                                        with tab2:
                                            st.subheader(":blue[Gene Clusters]", divider="blue")
                                            st.dataframe(clusters_db_gene, use_container_width=True, hide_index=True)

                            else:
                                    st.sidebar.warning(
                                        ":blue[Not enough genes for clustering them! You also can not use ] :red[ssGSEA (4)] "
                                        ":blue[and] "
                                        ":red[ Gene Ontology (10)]")


                            selected_genes_dim_red, option_dm, lat_vec, clust_nmf = dim_reductor(selected_genes1,
                                                                                                 selected_genes_nmf,
                                                                                                 random_state_samples,
                                                                                                 df_names,
                                                                                                 selected_genes_dm_name,
                                                                                                 user_session_id,
                                                                                                 selected_genes_name_nmf,
                                                                                                 column_cluster_path,
                                                                                                 option_dataset,
                                                                                                 adv_col1,
                                                                                                 adv_col2,
                                                                                                 ht_expression_path,
                                                                                                 ht_top_annotation_path,
                                                                                                 ht_png, ht_pdf,
                                                                                                 surv_png,
                                                                                                 surv_pdf,
                                                                                                 ht_row_cluster_path,
                                                                                                 tcga_anno, cancer_opt,
                                                                                                 sidebar_expander=sample_clust_exp,
                                                                                                 clinic_df=tcga_surv,
                                                                                                 surv_dataframe_path=surv_df,
                                                                                                 nmf_cons_name=nmf_cons_name,
                                                                                                 e=e,
                                                                                                 gene_clust_chb = gene_clust_chb)

                            if option_dm != 'NMF Clustering':
                                clusters_db = clustering(selected_genes_dim_red, option_dm, selected_genes1,
                                                         option_dataset,
                                                         df_names, column_cluster_path, lat_vec, user_session_id,
                                                         random_state_samples, adv_col1, adv_col2, ht_expression_path,
                                                         column_cluster_path, ht_top_annotation_path,
                                                         ht_png, ht_pdf, surv_png, surv_pdf, ht_row_cluster_path,
                                                         tcga_anno, cancer_opt, scatter_2d_name=scatter_2d_name,
                                                         scatter_3d_name=scatter_3d_name, clinic_df=tcga_surv,
                                                         surv_dataframe_path=surv_df, sidebar_expander=sample_clust_exp,
                                                         selected_genes_dm_path=selected_genes_dm_name,
                                                         cluster_toggle=gene_clust_chb, e=e)

                            end = time.time()
                            wait_time = end - start

                            new_results_surv = pd.read_csv(surv_df, nrows=25)
                            if 'last_results_surv' not in st.session_state or not np.array_equal(
                                    st.session_state['last_results_surv'], new_results_surv):
                                delete_files(wait_time + 2)
                            st.session_state['last_results_surv'] = copy.deepcopy(new_results_surv)

                            new_results_gene = clusters_db_gene
                            if 'last_results_gene' not in st.session_state or not np.array_equal(
                                    st.session_state['last_results_gene'], new_results_gene):
                                delete_files(wait_time + 2)
                            st.session_state['last_results_gene'] = copy.deepcopy(new_results_gene)

                            new_results_hm_zscore = pd.read_csv(ht_expression_path, nrows=25)
                            if 'last_results_hm_z' not in st.session_state or not np.array_equal(
                                    st.session_state['last_results_hm_z'], new_results_hm_zscore):
                                delete_files(wait_time + 2)
                            st.session_state['last_results_hm_z'] = copy.deepcopy(new_results_hm_zscore)

                            new_results = clusters_db
                            if 'last_results' not in st.session_state or not np.array_equal(
                                    st.session_state['last_results'], new_results):
                                delete_files(wait_time + 2)
                            st.session_state['last_results'] = copy.deepcopy(new_results)

                            new_results_nmf = clust_nmf
                            if 'last_results_nmf' not in st.session_state or not np.array_equal(
                                    st.session_state['last_results_nmf'], new_results_nmf):
                                delete_files(wait_time + 2)
                            st.session_state['last_results_nmf'] = copy.deepcopy(new_results_nmf)

                            new_results_dr = selected_genes_dim_red[:25]
                            if 'last_results_dr' not in st.session_state or not np.array_equal(
                                    st.session_state['last_results_dr'], new_results_dr):
                                delete_files(wait_time + 2)
                            st.session_state['last_results_dr'] = copy.deepcopy(new_results_dr)

                            try:
                                plot_col1, plot_col2 = st.columns(2)
                                plot_col1.subheader(":blue[Clustering Results and Heatmap]", divider="blue")
                                plot_col2.subheader(":blue[Survival Results between Clusters]", divider="blue")
                                st.subheader("", divider="blue")
                                plot_col1.image(ht_png, use_container_width=True, output_format='PNG')
                                data_surv = pd.read_csv(surv_dataframe_path)
                                try:
                                    plot_col2.image(surv_png, use_container_width=True, output_format='PNG')
                                    survival_info(surv_dataframe_path, user_session_id, surv_type="cluster")
                                except:
                                    plot_col2.image("style_items/error.svg", use_container_width=True)
                                    plot_col2.warning(":blue[Insufficient data available for] :red[Survival Analysis].")
                            except:
                                st.info(
                                    ":blue[Let's click ]" ':red[ "Draw Heatmap and Survival Plot" ] :blue[button on the sidebar]')

                            file_copy(user_session_id)

                            with tab6:
                                ssgsea_chb = st.toggle(":blue[Use ] :red[ssGSEA Survival]", key='ssgsea_cl',
                                                         help=dscript.ssgsea_surv_intro())
                                st.header(" ", divider="blue")
                                if selected_genes1.shape[1] >= 1:
                                    if ssgsea_chb:
                                        if gene_clust_chb:
                                            cluster_gene_db = pd.read_csv(ht_row_cluster_path, index_col=0)
                                        else:
                                            cluster_gene_db = pd.DataFrame()
                                            cluster_gene_db["Samples"] = df_names_gene
                                            cluster_gene_db["Cluster"] = "All Genes"
                                            cluster_gene_db.set_index("Samples", inplace=True)

                                        ssgsea_data = target_cancer_exp_df.loc[:, target_cancer_exp_df.columns.isin(
                                            filtered_target_cancer.columns)]
                                        ssgsea_analysis(random_state_samples, ssgsea_data, cluster_gene_db,
                                                        option_dataset, tcga_clinics.copy(), ssgsea_violin_box,
                                                        ssgsea_result_path, user_session_id, tcga_anno,
                                                        gene_cluster=gene_clust_chb, subset_gene_exp_df=selected_genes,
                                                        st_anno=tcga_anno, sample_cluster=clusters_db)

                                else:
                                    st.warning(":blue[Too few gene(s) for ssGSEA! At least ] :red[4 genes]!")

                            with tab3:
                                ssgsea_comp_chb = st.toggle(":blue[Use ] :red[ssGSEA Compare]", key='ssgsea_comp',
                                                              help=dscript.ssgsea_comp_intro())
                                st.header(" ", divider="blue")

                                if selected_genes1.shape[1] >= 1:
                                    if ssgsea_comp_chb:
                                        if gene_clust_chb:
                                            cluster_gene_db = pd.read_csv(ht_row_cluster_path, index_col=0)
                                        else:
                                            cluster_gene_db = pd.DataFrame()
                                            cluster_gene_db["Samples"] = df_names_gene
                                            cluster_gene_db["Cluster"] = "All Genes"
                                            cluster_gene_db.set_index("Samples", inplace=True)

                                        ssgsea_data = target_cancer_exp_df.loc[:, target_cancer_exp_df.columns.isin(
                                            filtered_target_cancer.columns)]
                                        ssgsea_compare(target_cancer_exp_df, violoin_comapre_violin_path,
                                                       random_state_samples, ssgsea_data,
                                                       cluster_gene_db, option_dataset, tcga_anno.copy(),
                                                       ssgsea_result_path, user_session_id, tcga_anno,
                                                       gene_cluster=gene_clust_chb,
                                                       st_anno=tcga_anno,
                                                       sample_cluster=clusters_db )

                                else:
                                    st.warning(":blue[Too few gene(s) for ssGSEA! At least ] :red[4 genes]!")


                            with tab5:
                                chi_ch = st.toggle(":blue[Use ] :red[Multivariate ] :blue[and ] "
                                                     ":red[Chi Square Test]", key="chi_ch",
                                                     help=dscript.multi_intro())
                                st.header(" ", divider="blue")

                                if chi_ch:
                                    if clusters_db is not None and len(clusters_db) != 0:
                                        mv_db(clusters_db, selected_genes1, multi_db_name,
                                              user_session_id, option_dataset, multi_info_name,
                                              forest_pdf, forest_png, tcga_anno.copy(), tcga_clinics.copy())

                                    elif clusters_db is None and option_dm == "NMF Clustering":
                                        mv_db(clust_nmf, selected_genes1, multi_db_name,
                                              user_session_id, option_dataset, multi_info_name,
                                              forest_pdf, forest_png, tcga_anno.copy(), tcga_clinics.copy())

                                    elif clusters_db is None:
                                        pass

                                if chi_ch:
                                    if option_dm != 'NMF Clustering':
                                        chi_sqrt_test(clusters_db, tcga_anno.copy(),
                                                      cont_name, comp_name, chi_hm_name,
                                                      chi_bp_name)
                                    else:
                                        chi_sqrt_test(clust_nmf, tcga_anno.copy(),
                                                      cont_name, comp_name, chi_hm_name,
                                                      chi_bp_name)

                            with tab7:
                                gen_ont = st.toggle(":blue[Use  ] :red[Gene Ontology]", key="gen_ont_chb",
                                                      help=dscript.go_intro())
                                st.header(" ", divider="blue")

                                if selected_genes1.shape[1] >= 10:
                                    if gen_ont:
                                        gene_ont(user_session_id, option_dataset, ht_row_cluster_path,
                                                 target_cancer_exp_df.index, selected_genes.index, gene_clust_chb)

                                else:
                                    st.warning(":blue[Increase the number of genes: ] :red[minimum 10 genes!]")

                        else:
                            pass
                    else:
                        pass

                    with tab8:
                        warning_tab()
                        single_gene_chb = st.toggle(":blue[Use ] :red[Gene Survival]", key="sg_chb",
                                                      help=dscript.sg_intro())
                        st.subheader(" ", divider="blue")
                        if single_gene_chb:

                            sg_exp = st.expander(":blue[Subtype Selection for Gene Survival]")
                            filtered_target_cancer_sg, tcga_anno_sg = subtypes_selector_tcga_sg(target_cancer_exp_df,
                                                                                                tcga_clinics,
                                                                                                tcga_subtypes,
                                                                                                ht_top_annotation_path,
                                                                                                sg_exp,
                                                                                                key='sg')
                            if len(filtered_target_cancer_sg.columns) >= 8:
                                sg_analysis_surv(user_session_id,
                                                 filtered_target_cancer_sg,
                                                 tcga_clinics.copy(), option_dataset,
                                                 tcga_anno_sg.copy())

                            else:
                                st.warning(":blue[The number of samples should be at ] :red[least 8.]")

                    with tab9:
                        warning_tab()
                        comp_gene_chb = st.toggle(":blue[Use ] :red[Gene Compare]", key="gc_chb",
                                                    help=dscript.gene_comp_into())
                        st.subheader(" ", divider="blue")
                        if comp_gene_chb:
                            cg_exp = st.expander(":blue[Subtype Selection for Gene Compare]")
                            filtered_target_cancer_sg, tcga_anno_sg = subtypes_selector_tcga_sg(target_cancer_exp_df,
                                                                                                tcga_clinics,
                                                                                                tcga_subtypes,
                                                                                                ht_top_annotation_path,
                                                                                                cg_exp,
                                                                                                key='gc')

                            columns_to_drop = [col for col in tcga_anno_sg.columns if tcga_anno_sg[col].nunique() < 2]
                            tcga_anno_sg = tcga_anno_sg.drop(columns=columns_to_drop)

                            if len(filtered_target_cancer_sg.columns) >= 8:
                                gene_comp(user_session_id, filtered_target_cancer_sg, tcga_anno_sg, option_dataset,
                                          tcga_anno_sg, clusters_db)



                            else:
                                st.warning(":blue[The number of samples should be at ] :red[least 8.]")

            elif (option_dataset == "Anish-SCLC" or option_dataset == "Jiang-SCLC" or option_dataset == "George-LCNEC"
                  or option_dataset == "Fernandez-Carcinoid" or option_dataset == "Alcala-Carcinoid"
                  or option_dataset == 'Liu-SCLC' or option_dataset == 'George-SCLC'
                  or option_dataset == "Rousseaux-Mixed(non-NE)"):

                start = time.time()

                tab1, tab2, tab4, tab6, tab3, tab9, tab8, tab7, tab5, tab11, tab10 = st.tabs(
                    ["Clustering :robot_face:", "Genes Handling 🧬",
                     "Subtypes 🔍", "ssGSEA Survival  🔬",
                     "ssGSEA Compare  :chart_with_upwards_trend:",
                     "Gene Compare ⚗️",
                     "Gene Survival :petri_dish:",
                     "Ontology :test_tube:",
                     "Multivariate & Chi" + '$$^2$$ ' + ":bar_chart:",
                     "Gene Set Finder 🧩",
                     "Help 	:question:"])

                with tab10:
                    help_tab()

                with tab1:
                    if option_dataset == 'Anish-SCLC':
                        df, clinic_df, clinic_st = dl.anish_sclc_data_load()
                        name_list = []
                    elif option_dataset == 'George-SCLC':
                        df, clinic_df, clinic_st = dl.george_sclc_data_load()
                        name_list = []
                    elif option_dataset == 'George-LCNEC':
                        df, clinic_df, clinic_st = dl.george_lcnec_data_load()
                        name_list = []
                    elif option_dataset == 'Fernandez-Carcinoid':
                        df, clinic_df, clinic_st = dl.carcinoid_data_load()
                        name_list = []
                    elif option_dataset == 'Liu-SCLC':
                        df, clinic_df, clinic_st = dl.liu_data_load()
                        name_list = []
                    elif option_dataset == 'Alcala-Carcinoid':
                        surv_diff_ch = st.empty()
                        surv_diff_ch.warning(':blue[LC1/2/3  samples lack survival data. It is possible that ] '
                                             ':red[an entire cluster is excluded from the analysis. ] '
                                             ':blue[You can use on the sidebar, ] :red["LC Exclusion" ] '
                                             ':blue[checkbox to ] :red[remove ] '
                                             ':blue[these samples from the analysis.]')
                        exclusion_chb = st.sidebar.checkbox(":blue[LC ] :red[Exclusion]", key="ex_chb")

                        if exclusion_chb:
                            surv_diff_ch.empty()

                        df, clinic_df, clinic_st = dl.alcala_data_load(exclusion_chb)
                        name_list = []
                    elif option_dataset == 'Rousseaux-Mixed(non-NE)':
                        surv_diff_ch = st.empty()
                        surv_diff_ch.warning(
                            ':blue[Non-Tumoral-Lung (NTL)  samples lack survival data. It is possible that ] '
                            ':red[an entire cluster is excluded from the analysis. ] '
                            ':blue[You can use on the sidebar, ] :red["NTL Exclusion" ] '
                            ':blue[checkbox to ] :red[remove ] '
                            ':blue[these samples from the analysis.]')
                        exclusion_chb = st.sidebar.checkbox(":blue[NTL ] :red[Exclusion]", key="ex_chb_ntl")

                        if exclusion_chb:
                            surv_diff_ch.empty()

                        surv_type = dl.rousseaux_surv_type_selector()
                        new_results_rousseax_surv_type = surv_type
                        df, clinic_df, clinic_st = dl.rousseaux_data_load(exclusion_chb, surv_type)
                        name_list = []
                    else:
                        df, clinic_df, clinic_st = dl.jiang_sclc_data_load()
                        name_list = []

                with tab2:
                    uploader_container = st.sidebar.empty()
                    gene_set_finder_container = st.sidebar.empty()
                    demo_conatainer_container = st.sidebar.empty()

                    feature_chb = gene_set_finder_container.toggle(":red[Find a ] :blue[Gene Set]",
                                                    key="feature_2", help=dscript.dataset_ana_intro())

                    uploaded_file = uploader_container.file_uploader(':blue[Upload your ] :red[gene list]',
                                                         label_visibility="visible",
                                                         type=('txt', 'csv', 'tsv'),
                                                         key="gene_uploader")

                    demo_toggle = demo_conatainer_container.toggle(":blue[Start ] :red[Demo]", key="demo_toggle",
                                                                   help=dscript.demo_gene_list_help())

                    # Handle file input
                    if demo_toggle:
                        # Simulate an uploaded file-like object for the demo file
                        with open("source_data/NE50-signatures.csv", "rb") as demo_file:
                            demo_content = demo_file.read()
                        # Create a BytesIO object to simulate file upload
                        uploaded_file = io.BytesIO(demo_content)
                        # Mimic the UploadedFile attributes
                        uploaded_file.name = "NE50-signatures.csv"
                        uploaded_file.size = len(demo_content)

                    if feature_chb:
                        uploader_container.empty()
                        demo_conatainer_container.empty()

                    if demo_toggle:
                        uploader_container.empty()
                        gene_set_finder_container.empty()
                        st.toast("🖱️ :blue[Click to the ] :red[Draw Heatmap and Survival Plot ] :blue[Button!] 🖱️")

                    st.sidebar.subheader(" ", divider="blue")

                    with tab11:
                        if not feature_chb:
                            st.warning(
                                ':blue[Please Swich On the: ] :red["Find a Gene Set" ] :blue[Toggle on the Sidebar]')

                    if feature_chb:
                        with tab11:
                            selected_genes = gene_finder(option_dataset, df, gene_sig_name, user_session_id, clinic_st,
                                                         clinic_df)

                            if selected_genes is not None:
                                pass
                            else:
                                with tab1:
                                    st.warning(
                                        ":red[No ] :blue[ any significant genes!] :blue[Please ] :red[change the settings of analisys! ]"
                                        ":blue[ If you use gene list please ] :red[upload one!]")
                                with tab2:
                                    st.warning(
                                        ":red[No ] :blue[ any significant genes!] :blue[Please ] :red[change the settings of analisys! ]"
                                        ":blue[ If you use gene list please ] :red[upload one!]")
                                with tab3:
                                    st.warning(
                                        ":red[No ] :blue[ any significant genes!] :blue[Please ] :red[change the settings of analisys! ]"
                                        ":blue[ If you use gene list please ] :red[upload one!]")
                                with tab4:
                                    st.warning(
                                        ":red[No ] :blue[ any significant genes!] :blue[Please ] :red[change the settings of analisys! ]"
                                        ":blue[ If you use gene list please ] :red[upload one!]")
                                with tab5:
                                    st.warning(
                                        ":red[No ] :blue[ any significant genes!] :blue[Please ] :red[change the settings of analisys! ]"
                                        ":blue[ If you use gene list please ] :red[upload one!]")
                                with tab6:
                                    st.warning(
                                        ":red[No ] :blue[ any significant genes!] :blue[Please ] :red[change the settings of analisys! ]"
                                        ":blue[ If you use gene list please ] :red[upload one!]")
                                with tab7:
                                    st.warning(
                                        ":red[No ] :blue[ any significant genes!] :blue[Please ] :red[change the settings of analisys! ]"
                                        ":blue[ If you use gene list please ] :red[upload one!]")
                                with tab8:
                                    st.warning(
                                        ":red[No ] :blue[ any significant genes!] :blue[Please ] :red[change the settings of analisys! ]"
                                        ":blue[ If you use gene list please ] :red[upload one!]")
                                with tab9:
                                    st.warning(
                                        ":red[No ] :blue[ any significant genes!] :blue[Please ] :red[change the settings of analisys! ]"
                                        ":blue[ If you use gene list please ] :red[upload one!]")
                                with tab10:
                                    st.warning(
                                        ":red[No ] :blue[ any significant genes!] :blue[Please ] :red[change the settings of analisys! ]"
                                        ":blue[ If you use gene list please ] :red[upload one!]")
                                with tab11:
                                    st.warning(
                                        ":red[No ] :blue[ any significant genes!] :blue[Please ] :red[change the settings of analisys! ]"
                                        ":blue[ If you use gene list please ] :red[upload one!]")
                                st.stop()

                        gl.col_maker(selected_genes)

                        cluster_gene_list_toggle = False

                        if uploaded_file:
                            feature_chb = False

                    else:
                        cluster_gene_list_toggle = False

                        cluster_check = gl.cluster_checker(uploaded_file)

                        if demo_toggle:
                            demo_cluster = True
                            cluster_check = True
                        else:
                            demo_cluster = False

                        if uploaded_file is not None and cluster_check:
                            cluster_gene_list_toggle = st.toggle(":blue[Use Gene List with ] :red[Clusters]",
                                                                 key="cluster_toggle_ne", value=demo_cluster)
                            if cluster_gene_list_toggle:
                                selected_genes = gl.cluster_app(user_session_id, user_gene_cluster_path,
                                                                option_dataset=option_dataset, df=df,
                                                                uploaded_file=uploaded_file)
                            else:
                                selected_genes = gl.app(user_session_id, option_dataset=option_dataset, df=df,
                                                        uploaded_file=uploaded_file)

                        else:
                            selected_genes = gl.app(user_session_id, option_dataset=option_dataset, df=df,
                                                    uploaded_file=uploaded_file)
                    if selected_genes is not None:
                        # Check if min and max are the same
                        all_same = selected_genes.values.min() == selected_genes.values.max()

                        if all_same:
                            st.warning(":red[All genes expression values are same!]")
                            with tab1:
                                if all_same:
                                    st.warning(":red[All genes expression values are same!]")
                            st.stop()
                        else:
                            pass
                    if selected_genes is not None and selected_genes.shape[0] >= 1:

                        if selected_genes.shape[0] > 2500:
                            st.sidebar.error(":blue[Too many genes. Please ] :red[reduce ] :blue[the number of genes!]")
                            exit()  # Halt the script execution without causing a Streamlit error
                        elif 500 <= selected_genes.shape[0] <= 2500:
                            st.sidebar.warning(":blue[Too many genes] :red[slow down] :blue[application runtimes]")
                            selected_genes_anno = selected_genes.copy()
                            selected_genes_anno = selected_genes_anno.transpose()
                            if option_dataset == 'Anish-SCLC':
                                st_annotation = dl.data_load_anno_anish(selected_genes_anno)
                            elif option_dataset == 'George-SCLC':
                                st_annotation = dl.data_load_anno_george(selected_genes_anno)
                            elif option_dataset == 'George-LCNEC':
                                st_annotation = dl.data_load_anno_george_lcnec(selected_genes_anno)
                            elif option_dataset == 'Fernandez-Carcinoid':
                                st_annotation = dl.data_load_anno_carcinoid(selected_genes_anno)
                            elif option_dataset == 'Liu-SCLC':
                                st_annotation = dl.data_load_anno_liu(selected_genes_anno)
                            elif option_dataset == 'Alcala-Carcinoid':
                                st_annotation = dl.data_load_anno_alcala(selected_genes_anno)
                            elif option_dataset == 'Rousseaux-Mixed(non-NE)':
                                st_annotation = dl.data_load_anno_rousseaux(selected_genes_anno)
                            else:
                                st_annotation = dl.data_load_anno_jiang(selected_genes_anno)
                        else:
                            selected_genes_anno = selected_genes.copy()
                            selected_genes_anno = selected_genes_anno.transpose()
                            if option_dataset == 'Anish-SCLC':
                                st_annotation = dl.data_load_anno_anish(selected_genes_anno)
                            elif option_dataset == 'George-SCLC':
                                st_annotation = dl.data_load_anno_george(selected_genes_anno)
                            elif option_dataset == 'George-LCNEC':
                                st_annotation = dl.data_load_anno_george(selected_genes_anno)
                            elif option_dataset == 'Fernandez-Carcinoid':
                                st_annotation = dl.data_load_anno_carcinoid(selected_genes_anno)
                            elif option_dataset == 'Liu-SCLC':
                                st_annotation = dl.data_load_anno_liu(selected_genes_anno)
                            elif option_dataset == 'Alcala-Carcinoid':
                                st_annotation = dl.data_load_anno_alcala(selected_genes_anno)
                            elif option_dataset == 'Rousseaux-Mixed(non-NE)':
                                st_annotation = dl.data_load_anno_rousseaux(selected_genes_anno)
                            else:
                                st_annotation = dl.data_load_anno_jiang(selected_genes_anno)

                with tab1:
                    if selected_genes is None:
                        dscript.gene_list_warning()
                        display_free_use_message()
                        warning_tab()
                        img_container = st.container(key='main_img_NE')
                        img_container.image("style_items/main_fig.svg", use_container_width=False, width=750)

                with tab4:
                    if selected_genes is None:
                        dscript.gene_list_warning()
                    if selected_genes is not None and selected_genes.shape[0] >= 1:

                        selected_genes_filtered, name_list, merged_df = subtype_select(selected_genes, st_annotation,
                                                                                       selected_genes_anno,
                                                                                       ht_top_annotation_path,
                                                                                       option_dataset,
                                                                                       clinic_st)

                        hm_chb = st.toggle(':red[Display Annotations]', key='hm_preview', value=True)
                        st.subheader(" ", divider="blue")
                        if len(selected_genes_filtered.columns) > 0:
                            if hm_chb:
                                plotly_heatmap(merged_df, name_list, option_dataset, sub_hm_name)
                        else:
                            st.warning(":blue[The number of samples should be at ] :red[least 8.]")

                with tab5:
                    if selected_genes is None:
                        dscript.gene_list_warning()
                    if len(name_list) < 8 and selected_genes is not None:
                        st.warning(":blue[The number of samples should be at ] :red[least 8.]")
                    warning_tab()
                with tab6:
                    if selected_genes is None:
                        dscript.gene_list_warning()
                    if len(name_list) < 8 and selected_genes is not None:
                        st.warning(":blue[The number of samples should be at ] :red[least 8.]")
                    warning_tab()
                with tab3:
                    if selected_genes is None:
                        dscript.gene_list_warning()
                    if len(name_list) < 8 and selected_genes is not None:
                        st.warning(":blue[The number of samples should be at ] :red[least 8.]")
                    warning_tab()
                with tab7:
                    if selected_genes is None:
                        dscript.gene_list_warning()
                    if len(name_list) < 8 and selected_genes is not None:
                        st.warning(":blue[The number of samples should be at ] :red[least 8.]")
                    warning_tab()
                with tab1:
                    if selected_genes is not None and selected_genes.shape[0] >= 1:
                        selected_genes1, selected_genes_nmf, df_names, random_state_samples, df_names_gene = \
                            data_preprocess(selected_genes_name, selected_genes_name_nmf, selected_genes_filtered,
                                            name_list, user_session_id)
                        if selected_genes1 is None and selected_genes_filtered.shape[1] > 8:
                            st.error(
                                ":blue[Gene(s) ] :red[not in the database ] :blue[or you did] :red[not select genes. ] "
                                ":blue[You can also try another ] "
                                ":blue[database or switch off Z-scoring]")
                        else:
                            pass

                        if selected_genes.shape[1] >= 8 and selected_genes1 is not None:
                            gene_clust_exp = st.sidebar.expander(':blue[Gene Clustering]')
                            sample_clust_exp = st.sidebar.expander(':blue[Sample Clustering]', expanded=True)
                            gene_cont = gene_clust_exp.container()
                            gene_clust_col1, gene_clust_col2 = st.sidebar.columns(2)
                            adv_exp = st.expander(":blue[Advanced Settings]")
                            adv_col1, adv_col2 = adv_exp.columns(2)
                            adv_exp_gene = st.expander(":blue[Advanced Settings of Gene Clustering]")
                            adv_col1_gene, adv_col2_gene = adv_exp_gene.columns(2)
                            e = st.columns(2, gap="large")

                            gene_fake_cluster = pd.DataFrame()
                            gene_fake_cluster["Names"] = df_names_gene
                            gene_fake_cluster["Cluster"] = "Genes"
                            gene_fake_cluster.to_csv(file_name_list[6], index=False)

                            if selected_genes1.shape[1] >= 8:
                                if cluster_gene_list_toggle:
                                    gene_clust_use = True
                                else:
                                    gene_clust_use = False

                                gene_clust_chb = gene_cont.toggle(":red[Clustering ] :blue[of Genes]",
                                                                    key="gene_clust_chb", value=gene_clust_use)

                                if gene_clust_chb:
                                    ht_row_cluster_path = file_name_list[6]

                                    gene_selected_genes1 = selected_genes1.copy()
                                    gene_selected_genes1 = gene_selected_genes1.transpose()

                                    gene_selected_genes_nmf = selected_genes_nmf.copy()
                                    gene_selected_genes_nmf = gene_selected_genes_nmf.transpose()

                                    # Call dim_reductor_gene and check if outputs have changed
                                    selected_genes_dim_red_gene, option_dm_gene, lat_vec_gene = dim_reductor_gene(
                                        gene_selected_genes1,
                                        gene_selected_genes_nmf,
                                        random_state_samples,
                                        df_names_gene,
                                        selected_genes_dm_name_gene,
                                        user_session_id,
                                        selected_genes_name_nmf,
                                        ht_row_cluster_path,
                                        option_dataset,
                                        adv_col1_gene,
                                        adv_col2_gene,
                                        ht_png, ht_pdf,
                                        surv_png,
                                        surv_pdf,
                                        gene_clust_col1,
                                        clinic_df=clinic_df,
                                        surv_dataframe_path=surv_df,
                                        gene_exp=gene_cont)

                                    # Call clustering_genes and check if output has changed
                                    clusters_db_gene = clustering_genes(selected_genes_dim_red_gene,
                                                                        selected_genes1,
                                                                        df_names_gene,
                                                                        ht_row_cluster_path,
                                                                        lat_vec_gene,
                                                                        random_state_samples,
                                                                        adv_col1_gene, ht_expression_path,
                                                                        gene_clust_col2,
                                                                        scatter_2d_name_gene,
                                                                        scatter_3d_name_gene,
                                                                        gene_exp=gene_cont,
                                                                        adv_col2_gene=adv_col2_gene,
                                                                        option_dm=option_dm,
                                                                        user_session_id=user_session_id,
                                                                        column_cluster_path=ht_row_cluster_path,
                                                                        selected_genes_dm_path=selected_genes_dm_name_gene,
                                                                        user_gene_cluster_path=user_gene_cluster_path,
                                                                        cluster_toggle=cluster_gene_list_toggle,
                                                                        chb_gene_clust=gene_clust_chb, e = e)

                                    if clusters_db_gene is not None:
                                        with tab2:
                                            st.subheader(":blue[Gene Clusters]", divider="blue")
                                            st.dataframe(clusters_db_gene, use_container_width=True, hide_index=True)

                            else:
                                st.sidebar.warning(
                                    ":blue[Not enough genes for clustering them! You also can not use ] :red[ssGSEA (4)] "
                                    ":blue[and ] "
                                    ":red[ Gene Ontology (10)]")



                            selected_genes_dim_red, option_dm, lat_vec, clust_nmf = dim_reductor(selected_genes1,
                                                                                                 selected_genes_nmf,
                                                                                                 random_state_samples,
                                                                                                 df_names,
                                                                                                 selected_genes_dm_name,
                                                                                                 user_session_id,
                                                                                                 selected_genes_name_nmf,
                                                                                                 column_cluster_path,
                                                                                                 option_dataset,
                                                                                                 adv_col1,
                                                                                                 adv_col2,
                                                                                                 ht_expression_path,
                                                                                                 ht_top_annotation_path,
                                                                                                 ht_png, ht_pdf,
                                                                                                 surv_png,
                                                                                                 surv_pdf,
                                                                                                 ht_row_cluster_path,
                                                                                                 sidebar_expander=sample_clust_exp,
                                                                                                 selected_subtypes=None,
                                                                                                 cancer_opt=None,
                                                                                                 clinic_df=clinic_df,
                                                                                                 surv_dataframe_path=surv_df,
                                                                                                 nmf_cons_name=nmf_cons_name,
                                                                                                 e=e,
                                                                                                 gene_clust_chb = gene_clust_chb)

                            if option_dm != 'NMF Clustering':
                                clusters_db = clustering(selected_genes_dim_red, option_dm,
                                                         selected_genes1, option_dataset,
                                                         df_names, column_cluster_path, lat_vec,
                                                         user_session_id,
                                                         random_state_samples, adv_col1, adv_col2, ht_expression_path,
                                                         column_cluster_path, ht_top_annotation_path,
                                                         ht_png, ht_pdf, surv_png, surv_pdf, ht_row_cluster_path,
                                                         scatter_2d_name=scatter_2d_name,
                                                         scatter_3d_name=scatter_3d_name,
                                                         sidebar_expander=sample_clust_exp,
                                                         selected_subtypes=None, cancer_opt=None,
                                                         clinic_df=clinic_df, surv_dataframe_path=surv_df,
                                                         selected_genes_dm_path=selected_genes_dm_name,
                                                         cluster_toggle=gene_clust_chb, e=e)

                            end = time.time()
                            wait_time = end - start

                            if option_dataset == 'Rousseaux-Mixed(non-NE)':
                                new_results_surv = pd.read_csv(surv_df, nrows=25)
                                if 'last_results_surv' not in st.session_state or not np.array_equal(
                                        st.session_state['last_results_surv'], new_results_surv):
                                    delete_files(wait_time + 2)
                                st.session_state['last_results_surv'] = copy.deepcopy(new_results_surv)

                            new_results_ne_gene = clusters_db_gene
                            if 'last_results_ne_gene' not in st.session_state or not np.array_equal(
                                    st.session_state['last_results_ne_gene'], new_results_ne_gene):
                                delete_files(wait_time + 2)
                                st.session_state['last_results_ne_gene'] = new_results_ne_gene

                            new_results_hm_zscore_ne = pd.read_csv(ht_expression_path, nrows=25)
                            if 'last_results_hm_z_ne' not in st.session_state or not np.array_equal(
                                    st.session_state['last_results_hm_z_ne'], new_results_hm_zscore_ne):
                                delete_files(wait_time + 2)
                            st.session_state['last_results_hm_z_ne'] = copy.deepcopy(new_results_hm_zscore_ne)

                            new_results_ne = clusters_db
                            if 'last_results_ne' not in st.session_state or not np.array_equal(
                                    st.session_state['last_results_ne'], new_results_ne):
                                delete_files(wait_time + 2)
                                st.session_state['last_results_ne'] = new_results_ne

                            new_results_ne_nmf = clust_nmf
                            if 'last_results_ne_nmf' not in st.session_state or not np.array_equal(
                                    st.session_state['last_results_ne_nmf'], new_results_ne_nmf):
                                delete_files(wait_time + 2)
                                st.session_state['last_results_ne_nmf'] = copy.deepcopy(new_results_ne_nmf)

                            new_results_ne_dr = selected_genes_dim_red[:25]
                            if 'last_results_ne_dr' not in st.session_state or not np.array_equal(
                                    st.session_state['last_results_ne_dr'], new_results_ne_dr):
                                delete_files(wait_time + 2)
                                st.session_state['last_results_ne_dr'] = copy.deepcopy(new_results_ne_dr)

                            try:
                                plot_col1, plot_col2 = st.columns(2)
                                plot_col1.subheader(":blue[Clustering Results and Heatmap]", divider="blue")
                                plot_col2.subheader(":blue[Survival Results between Clusters]", divider="blue")
                                st.subheader("", divider="blue")
                                plot_col1.image(ht_png, use_container_width=True, output_format='PNG')
                                try:
                                    plot_col2.image(surv_png, use_container_width=True, output_format='PNG')
                                    survival_info(surv_dataframe_path, user_session_id, surv_type="cluster")
                                except:
                                    plot_col2.image("style_items/error.svg", use_container_width=True)
                                    plot_col2.warning(":blue[Insufficient data available for] :red[Survival Analysis].")

                            except:
                                st.info(
                                    ":blue[Let's click ]" ':red[ "Draw Heatmap and Survival Plot" ] :blue[button on the sidebar]')

                            file_copy(user_session_id)

                            with tab6:
                                ssgsea_chb = st.toggle(":blue[Use ] :red[ssGSEA Survival]", key='ssgsea_cl',
                                                         help=dscript.ssgsea_surv_intro())
                                st.header(" ", divider="blue")
                                if selected_genes1.shape[1] >= 1:
                                    if ssgsea_chb:
                                        if gene_clust_chb:
                                            cluster_gene_db = pd.read_csv(ht_row_cluster_path, index_col=0)
                                        else:
                                            cluster_gene_db = pd.DataFrame()
                                            cluster_gene_db["Samples"] = df_names_gene
                                            cluster_gene_db["Cluster"] = "All Genes"
                                            cluster_gene_db.set_index("Samples", inplace=True)

                                        ssgsea_data = df.loc[:, df.columns.isin(selected_genes_filtered.columns)]
                                        if option_dataset != 'Rousseaux-Mixed(non-NE)':
                                            ssgsea_analysis(random_state_samples, ssgsea_data, cluster_gene_db,
                                                            option_dataset, clinic_df, ssgsea_violin_box,
                                                            ssgsea_result_path, user_session_id, anno_file=clinic_st,
                                                            gene_cluster=gene_clust_chb,
                                                            subset_gene_exp_df=selected_genes,
                                                            st_anno=st_annotation, sample_cluster=clusters_db)
                                        else:
                                            ssgsea_analysis(random_state_samples, ssgsea_data, cluster_gene_db,
                                                            option_dataset, clinic_df, ssgsea_violin_box,
                                                            ssgsea_result_path, user_session_id, anno_file=clinic_st,
                                                            gene_cluster=gene_clust_chb,
                                                            subset_gene_exp_df=selected_genes,
                                                            st_anno=st_annotation, sample_cluster=clusters_db)

                                else:
                                    st.warning(":blue[Too few gene(s) for ssGSEA! At least ] :red[4 genes]!")

                            with tab3:
                                ssgsea_comp_chb = st.toggle(":blue[Use ] :red[ssGSEA Compare]", key='ssgsea_comp',
                                                              help=dscript.ssgsea_comp_intro())
                                st.header(" ", divider="blue")

                                if selected_genes1.shape[1] >= 1:
                                    if ssgsea_comp_chb:
                                        if gene_clust_chb:
                                            cluster_gene_db = pd.read_csv(ht_row_cluster_path, index_col=0)
                                        else:
                                            cluster_gene_db = pd.DataFrame()
                                            cluster_gene_db["Samples"] = df_names_gene
                                            cluster_gene_db["Cluster"] = "All Genes"
                                            cluster_gene_db.set_index("Samples", inplace=True)

                                        ssgsea_data = df.loc[:, df.columns.isin(selected_genes_filtered.columns)]
                                        if option_dataset != 'Rousseaux-Mixed(non-NE)':
                                            ssgsea_compare(df, violoin_comapre_violin_path,
                                                           random_state_samples, ssgsea_data, cluster_gene_db,
                                                           option_dataset,
                                                           st_annotation, ssgsea_result_path, user_session_id,
                                                           anno_file=clinic_st,
                                                           gene_cluster=gene_clust_chb,
                                                           st_anno=st_annotation,
                                                           sample_cluster=clusters_db)
                                        else:
                                            ssgsea_compare(df, violoin_comapre_violin_path,
                                                           random_state_samples, ssgsea_data,
                                                           cluster_gene_db,
                                                           option_dataset,
                                                           st_annotation, ssgsea_result_path, user_session_id,
                                                           anno_file=clinic_st,
                                                           gene_cluster=gene_clust_chb,
                                                           st_anno=st_annotation,
                                                           sample_cluster=clusters_db)
                                else:
                                    st.warning(":blue[Too few gene(s) for ssGSEA! At least ] :red[4 genes]!")

                            with tab5:
                                chi_ch = st.toggle(":blue[Use ] :red[Multivariate ] :blue[and ] "
                                                     ":red[Chi Square Test]", key="chi_ch",
                                                     help=dscript.multi_intro())
                                st.header(" ", divider="blue")

                                if chi_ch:
                                    if clusters_db is not None and len(clusters_db) != 0:
                                        if option_dataset != 'Rousseaux-Mixed(non-NE)':
                                            mv_db(clusters_db, selected_genes1, multi_db_name,
                                                  user_session_id, option_dataset, multi_info_name,
                                                  forest_pdf, forest_png, anno_file=st_annotation,
                                                  clinic_df=None)
                                        else:
                                            mv_db(clusters_db, selected_genes1, multi_db_name,
                                                  user_session_id, option_dataset, multi_info_name,
                                                  forest_pdf, forest_png, anno_file=clinic_st,
                                                  clinic_df=clinic_df)

                                    elif clusters_db is None and option_dm == "NMF Clustering":
                                        if option_dataset != 'Rousseaux-Mixed(non-NE)':
                                            mv_db(clust_nmf, selected_genes1, multi_db_name,
                                                  user_session_id, option_dataset, multi_info_name,
                                                  forest_pdf, forest_png, anno_file=None,
                                                  clinic_df=None)
                                        else:
                                            mv_db(clust_nmf, selected_genes1, multi_db_name,
                                                  user_session_id, option_dataset, multi_info_name,
                                                  forest_pdf, forest_png, anno_file=clinic_st,
                                                  clinic_df=clinic_df)

                                    elif clusters_db is None:
                                        pass

                                if chi_ch:
                                    if option_dm != 'NMF Clustering':
                                        chi_sqrt_test(clusters_db, merged_df, cont_name,
                                                      comp_name, chi_hm_name,
                                                      chi_bp_name)
                                    else:
                                        chi_sqrt_test(clust_nmf, merged_df, cont_name,
                                                      comp_name, chi_hm_name,
                                                      chi_bp_name)

                            with tab7:
                                def start_go(user_session_id, option_dataset, ht_row_cluster_path, df, selected_genes,
                                             gene_clust_chb):
                                    if selected_genes1.shape[1] >= 10:
                                        gene_ont(user_session_id, option_dataset, ht_row_cluster_path, df.index,
                                                 selected_genes.index, gene_clust_chb)
                                    else:
                                        st.warning(":blue[Increase the number of genes: ] :red[minimum 10 genes!]")

                                gen_ont = st.toggle(":blue[Use  ] :red[Gene Ontology]", key="gen_ont_chb",
                                                      help=dscript.go_intro())
                                st.header(" ", divider="blue")
                                if gen_ont:
                                    start_go(user_session_id, option_dataset, ht_row_cluster_path, df,
                                             selected_genes, gene_clust_chb)
                        else:
                            pass
                    else:
                        pass

                    with tab8:
                        warning_tab()
                        single_gene_chb = st.toggle(":blue[Use ] :red[Gene Survival]", key="sg_chb",
                                                      help=dscript.sg_intro())
                        st.subheader(" ", divider="blue")
                        if single_gene_chb:
                            if option_dataset == 'Anish-SCLC':
                                sg_anno = dl.data_load_anno_anish(df.transpose())
                            elif option_dataset == 'George-SCLC':
                                sg_anno = dl.data_load_anno_george(df.transpose())
                            elif option_dataset == 'George-LCNEC':
                                sg_anno = dl.data_load_anno_george_lcnec(df.transpose())
                            elif option_dataset == 'Fernandez-Carcinoid':
                                sg_anno = dl.data_load_anno_carcinoid(df.transpose())
                            elif option_dataset == 'Liu-SCLC':
                                sg_anno = dl.data_load_anno_liu(df.transpose())
                            elif option_dataset == 'Alcala-Carcinoid':
                                sg_anno = dl.data_load_anno_alcala(df.transpose())
                            elif option_dataset == 'Rousseaux-Mixed(non-NE)':
                                sg_anno = dl.data_load_anno_rousseaux(df.transpose())
                            else:
                                sg_anno = dl.data_load_anno_jiang(df.transpose())
                            sg_exp = st.expander(":blue[Subtype Selection for Gene Survival]")

                            selected_genes_filtered, name_list, merged_df = subtype_select_sg(df.transpose(),
                                                                                              sg_anno.copy(),
                                                                                              df.copy().transpose(),
                                                                                              option_dataset,
                                                                                              clinic_st.copy(), sg_exp,
                                                                                              key='sg')
                            if len(name_list) >= 8:
                                if option_dataset != 'Rousseaux-Mixed(non-NE)':
                                    sg_analysis_surv(user_session_id, selected_genes_filtered, clinic_df.copy(),
                                                     option_dataset, anno_file=None)
                                else:
                                    sg_analysis_surv(user_session_id, selected_genes_filtered, clinic_df.copy(),
                                                     option_dataset,  anno_file=clinic_st)
                            else:
                                st.warning(":blue[The number of samples should be at ] :red[least 8.]")

                    with tab9:
                        warning_tab()
                        comp_gene_chb = st.toggle(":blue[Use ] :red[Gene Compare]", key="gc_chb",
                                                    help=dscript.gene_comp_into())
                        st.subheader(" ", divider="blue")
                        if comp_gene_chb:
                            if option_dataset == 'Anish-SCLC':
                                gc_anno = dl.data_load_anno_anish(df.transpose())
                            elif option_dataset == 'George-SCLC':
                                gc_anno = dl.data_load_anno_george(df.transpose())
                            elif option_dataset == 'George-LCNEC':
                                gc_anno = dl.data_load_anno_george_lcnec(df.transpose())
                            elif option_dataset == 'Fernandez-Carcinoid':
                                gc_anno = dl.data_load_anno_carcinoid(df.transpose())
                            elif option_dataset == 'Liu-SCLC':
                                gc_anno = dl.data_load_anno_liu(df.transpose())
                            elif option_dataset == 'Alcala-Carcinoid':
                                gc_anno = dl.data_load_anno_alcala(df.transpose())
                            elif option_dataset == 'Rousseaux-Mixed(non-NE)':
                                gc_anno = dl.data_load_anno_rousseaux(df.transpose())
                            else:
                                gc_anno = dl.data_load_anno_jiang(df.transpose())

                            cg_exp = st.expander(":blue[Subtype Selection for Gene Compare]")

                            selected_genes_filtered, name_list, merged_df = subtype_select_sg(df.transpose(),
                                                                                              gc_anno.copy(),
                                                                                              df.copy().transpose(),
                                                                                              option_dataset,
                                                                                              clinic_st.copy(), cg_exp,
                                                                                              key='gene_comp')
                            if len(name_list) >= 8:
                                gene_comp(user_session_id, selected_genes_filtered, gc_anno, option_dataset, clinic_st,
                                          clusters_db)
                            else:
                                st.warning(":blue[The number of samples should be at ] :red[least 8.]")

            download_exp = st.sidebar.expander(":blue[Download Result]")

            dl_name = download_exp.text_input(
                ":blue[Give a Name for Your ] :red[Results]",
                key="ne_download",
                value="My_Results",
                max_chars=100,
                help=dscript.download()
            )

            # Create a placeholder for the info message
            info_placeholder = download_exp.empty()

            # Display the info message in the placeholder
            info_placeholder.info(
                "Please enter a name and click 'Accept and Prepare Download' to generate your results."
            )

            accept_dl = download_exp.button("Accept and Prepare Download")

            if accept_dl:
                # Clear the info message when the button is clicked
                info_placeholder.empty()
                # Proceed with the download
                downloader(user_session_id, dl_name, download_exp)

    except:
        st.warning(':blue[Something went ] :red[wrong]:red[! ] :blue[Please ] :red[refresh ] '
                   ':blue[the the webpage or ] :red[ try other settings. ] '
                   ':red[Thank you for your patience!]')


if __name__ == "__main__":
    main_app()
