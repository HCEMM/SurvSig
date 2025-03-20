import streamlit as st

@st.cache_data()
def help_tab():
    img_cols2 = st.columns(5)
    img_cols = st.columns(3)
    img_cols[1].image("style_items/main_fig.svg", use_container_width=True)
    img_cols[1].video("videos/WelcometoS.mp4")

    #Shorts
    st.header(":blue[Tutorial Short Videos] :movie_camera:", divider=True)

    #Expanders
    vid_epx1 = st.expander(":blue[**DATASETS**    ] "
                           ":red[    *Click Here to Watch Videos*]")
    vid_epx1.header(":red[DATASETS]")
    vid_cols1 = vid_epx1.columns(3)
    vid_cols1[0].subheader(':blue[Dataset Selection]')
    vid_cols1[0].video("videos/dataset.mp4")

    vid_epx2 = st.expander(":blue[**GENE SIGNATURE**    ] "
                           ":red[    *Click Here to Watch Videos*]")
    vid_cols2 = vid_epx2.columns(3)
    vid_cols2[0].subheader(':blue[Intro]')
    vid_cols2[0].video("videos/gene_sig_intro.mp4")
    vid_epx2.header(":red[USE OWN LIST]")
    vid_cols2 = vid_epx2.columns(3)
    vid_cols2[0].subheader(':blue[File Upload]')
    vid_cols2[0].video("videos/file_uplaod.mp4")
    vid_cols2[1].subheader(':blue[Gene Names Copy]')
    vid_cols2[1].video("videos/gene_names_copy.mp4")
    vid_cols2[2].subheader(':blue[Gene Handling]')
    vid_cols2[2].video("videos/gene_handling.mp4")
    vid_epx2.header(":red[CREATE GENE LIST]")
    vid_cols2 = vid_epx2.columns(3)
    vid_cols2[0].subheader(':blue[Dataset Analysing Intro]')
    vid_cols2[0].video("videos/dataset_analsing_intro.mp4")
    vid_cols2[1].subheader(':blue[Analyse the Dataset]')
    vid_cols2[1].video("videos/analyse_dataset.mp4")
    vid_cols2[2].subheader(':blue[Clinical Detail Based Analysing]')
    vid_cols2[2].video("videos/clin_based.mp4")

    vid_epx3 = st.expander(":blue[**SUBTYPE SELECTION**    ] "
                           ":red[    *Click Here to Watch Videos*]")
    vid_epx3.subheader(":red[SUBTYPE SELECTION]")
    vid_cols3 = vid_epx3.columns(3)
    vid_cols3[0].subheader(':blue[How to Select Subtypes]')
    vid_cols3[0].video("videos/st_selector.mp4")

    vid_epx4 = st.expander(":blue[**CLUSTERING**    ] "
                           ":red[    *Click Here to Watch Videos*]")
    vid_epx4.subheader(":red[CLUSTERING DATA]")
    vid_cols4 = vid_epx4.columns(3)
    vid_cols4[0].subheader(':blue[How to Clustering Data]')
    vid_cols4[0].video("videos/clust_dim.mp4")

    vid_epx5 = st.expander(":blue[**HEATMAP/SURVIVAL ANALYSIS**    ] "
                           ":red[    *Click Here to Watch Videos*]")
    vid_epx5.subheader(":red[SURVIVAL ANALYSIS & HEATMAP]")
    vid_cols5 = vid_epx5.columns(3)
    vid_cols5[0].subheader(':blue[How to Draw Plots]')
    vid_cols5[0].video("videos/heatmap_survival.mp4")

    vid_epx6 = st.expander(":blue[**ssGSEA ANALYSIS**    ] "
                           ":red[    *Click Here to Watch Videos*]")
    vid_epx6.subheader(":red[ssGSEA ANALYSIS]")
    vid_cols6 = vid_epx6.columns(3)
    vid_cols6[0].subheader(':blue[ssGSEA Survival Analysis]')
    vid_cols6[0].video("videos/ssgsea_surv.mp4")
    vid_cols6[1].subheader(':blue[ssGSEA Compare]')
    vid_cols6[1].video("videos/ssgsea_compare.mp4")

    vid_epx6 = st.expander(":blue[**GENE ANALYSIS**    ] "
                           ":red[    *Click Here to Watch Videos*]")
    vid_epx6.subheader(":red[SINGLE GENE ANALYSIS]")
    vid_cols6 = vid_epx6.columns(3)
    vid_cols6[0].subheader(':blue[Single Gene Survival Analysis]')
    vid_cols6[0].video("videos/gene_surv.mp4")
    vid_cols6[1].subheader(':blue[Gene Compare]')
    vid_cols6[1].video("videos/gene_compare.mp4")

    vid_epx7 = st.expander(":blue[**PATHWAY ANALYSIS**    ] "
                           ":red[    *Click Here to Watch Videos*]")
    vid_epx7.subheader(":red[GO & KEGG]")
    vid_cols7 = vid_epx7.columns(3)
    vid_cols7[0].subheader(':blue[How to Use Pathway Analysis]')
    vid_cols7[0].video("videos/pathway_gene.mp4")

    vid_epx8 = st.expander(":blue[**MULTIVARIATE ANALYSIS**    ] "
                           ":red[    *Click Here to Watch Videos*]")
    vid_epx8.subheader(":red[MULTIVARIATE]")
    vid_cols8 = vid_epx8.columns(3)
    vid_cols8[0].subheader(':blue[How to Use Multivariate Analysis and Forest Plot]')
    vid_cols8[0].video("videos/multivariate.mp4")

    vid_epx9 = st.expander(":blue[**CHI2 TEST**    ] "
                           ":red[    *Click Here to Watch Videos*]")
    vid_epx9.subheader(":red[CHI2 TEST]")
    vid_cols9 = vid_epx9.columns(3)
    vid_cols9[0].subheader(':blue[How to Use Chi2 Test and Setting the Plots]')
    vid_cols9[0].video("videos/chi2.mp4")
