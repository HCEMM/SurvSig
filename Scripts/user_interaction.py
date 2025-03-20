import streamlit as st
import string
import os
import datetime
import random
import uuid
from datetime import datetime
import shutil


# Define function to generate a random string
def generate_random_string(length):
    """
    Generate a random string of specified length.

    Args:
        length (int): The length of the random string to generate.

    Returns:
        str: A random string of specified length.
    """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for _ in range(length))


# Define function to get the user session ID
def get_session_id():
    """
    Get or create a user session ID.

    Returns:
        str: The user session ID.
    """
    if "session_id" not in st.session_state:
        st.session_state.session_id = uuid.uuid4()
    user_id = str(st.session_state.session_id)
    return user_id


# Define function to delete files older than a certain time in a directory
def del_files_py(remove_time):
    """
    Delete files and directories from the "./result_data" directory
    that are older than the specified time in seconds.

    Args:
        remove_time (int): Maximum age of files and directories to keep in seconds.
    """
    now = datetime.now()
    del_dir = "result_data"

    if not os.path.isdir(del_dir):
        print(f"Directory {del_dir} does not exist")
        return

    for filename in os.listdir(del_dir):
        filepath = os.path.join(del_dir, filename)
        try:
            creation_time = datetime.fromtimestamp(os.path.getctime(filepath))
            time_difference = now - creation_time
            if time_difference.total_seconds() > remove_time:
                if os.path.isfile(filepath):
                    os.remove(filepath)
                elif os.path.isdir(filepath):
                    shutil.rmtree(filepath)
        except OSError as e:
            print(f"Error: {e}")


# Define function to delete files for a specific user older than a certain time
def del_files_py_user(remove_time):
    """
    Delete files and directories from the "./result_data" directory
    that are older than the specified time in seconds for a specific user.

    Args:
        remove_time (int): Maximum age of files and directories to keep in seconds.
    """
    now = datetime.now()
    user_id = get_session_id()
    del_dir = f"result_data/result_data_{user_id}"

    if not os.path.isdir(del_dir):
        print(f"Directory {del_dir} does not exist")
        return

    for filename in os.listdir(del_dir):
        filepath = os.path.join(del_dir, filename)
        try:
            creation_time = datetime.fromtimestamp(os.path.getctime(filepath))
            time_difference = now - creation_time
            if time_difference.total_seconds() > remove_time:
                if os.path.isfile(filepath):
                    os.remove(filepath)
                elif os.path.isdir(filepath):
                    shutil.rmtree(filepath)
        except OSError as e:
            print(f"Error: {e}")


# Define function to delete a specific .zip file older than a certain time
def del_file_py_zip(remove_time):
    """
    Delete a specific .zip file from the "./result_data" directory
    that is older than the specified time in seconds.

    Args:
        remove_time (int): Maximum age of the .zip file to keep in seconds.
    """
    now = datetime.now()
    user_id = get_session_id()
    zip_filepath = os.path.join(os.getcwd(), "result_data", f"result_data_{user_id}.zip")

    try:
        creation_time = datetime.fromtimestamp(os.path.getctime(zip_filepath))
        time_difference = now - creation_time
        if time_difference.total_seconds() > remove_time:
            os.remove(zip_filepath)
    except:
        pass


# Define function to make unique names for files for users
def naming(user_session_id):
    """
    Generate a list of unique file paths for a user session.

    Args:
        user_session_id (str): The user session ID.

    Returns:
        list: A list of file paths.
    """
    user_dir = f"/result_data_{user_session_id}"
    try:
        os.makedirs(f"result_data{user_dir}", exist_ok=True)
    except OSError as error:
        print(f"Directory '{user_dir}' already exists")

    selected_genes_name = f"result_data{user_dir}/selected_genes_{user_session_id}.csv"  # 0
    selected_genes_name_nmf = f"result_data{user_dir}/selected_genes_nmf_{user_session_id}.csv"  # 1
    anno_file_name = f"result_data{user_dir}/anno_{user_session_id}.csv"  # 2
    selected_genes_dm_name = f"result_data{user_dir}/selected_genes_dm_{user_session_id}.csv"  # 3
    cluster_name = f"result_data{user_dir}/cluster_{user_session_id}.csv"  # 4
    surv_df_name = f"result_data{user_dir}/surv_{user_session_id}.csv"  # 5
    gene_cluster = f"result_data{user_dir}/cluster_gene_{user_session_id}.csv"  # 6
    heatmap_png_name = f"result_data{user_dir}/heatmap_{user_session_id}.png"  # 7
    surv_plot_png_name = f"result_data{user_dir}/surv_plot_{user_session_id}.png"  # 8
    heatmap_pdf_name = f"result_data{user_dir}/heatmap_{user_session_id}.pdf"  # 9
    surv_plot_pdf_name = f"result_data{user_dir}/surv_plot_{user_session_id}.pdf"  # 10
    multi_db_name = f"result_data{user_dir}/multi_db_{user_session_id}.csv"  # 11
    cont_name = f"result_data{user_dir}/cont_mtx_{user_session_id}.csv"  # 12
    comp_name = f"result_data{user_dir}/comp_df_{user_session_id}.csv"  # 13
    scatter_2d_name = f"result_data{user_dir}/scatter_2d_{user_session_id}.pdf"  # 14
    scatter_3d_name = f"result_data{user_dir}/scatter_3d_{user_session_id}.pdf"  # 15
    nmf_cons_name = f"result_data{user_dir}/nmf_cons_{user_session_id}.pdf"  # 16
    sub_hm_name = f"result_data{user_dir}/sub_hm_{user_session_id}.pdf"  # 17
    chi_hm_name = f"result_data{user_dir}/chi2_hm_{user_session_id}.pdf"  # 18
    chi_bp_name = f"result_data{user_dir}/chi2_bp_{user_session_id}.pdf"  # 19
    scatter_2d_name_gene = f"result_data{user_dir}/scatter_2d_gene_{user_session_id}.pdf"  # 20
    scatter_3d_name_gene = f"result_data{user_dir}/scatter_3d_gene_{user_session_id}.png"  # 21
    multi_info_name = f"result_data{user_dir}/multi_var_info_{user_session_id}.tsv"  # 22
    forest_pdf = f"result_data{user_dir}/forest_{user_session_id}.pdf"  # 23
    forest_png = f"result_data{user_dir}/forest_{user_session_id}.png"  # 24
    gene_sig_name = f"result_data{user_dir}/gene_list_maker_info_{user_session_id}.txt"  # 25
    strip_ssgsea = f"result_data{user_dir}/violin_box_{user_session_id}.pdf"  # 26
    corr_ssgsea = f"result_data{user_dir}/corr_{user_session_id}.pdf"  # 27
    ssgsea_result_path = f"result_data{user_dir}/ssgsea_result_{user_session_id}.csv"  # 28
    anno_color_corr = f"result_data{user_dir}/color_corr_{user_session_id}.csv"  # 29
    drug_corr_path = f"result_data{user_dir}/drug_corr_{user_session_id}.tsv"  # 30
    selected_genes_dm_name_gene = f"result_data{user_dir}/selected_genes_dm_gene_{user_session_id}.csv"  # 31
    hm_ssgsea_groups_path = f"result_data{user_dir}/hm_ssgsea_groups_{user_session_id}.csv"  # 32
    user_gene_cluster_path = f"result_data{user_dir}/user_gene_cluster_{user_session_id}.csv"  # 33
    compare_violin_ssgsea = f"result_data{user_dir}/violin_box_compare_{user_session_id}.pdf"  # 34
    anno_color_corr_gene = f"result_data{user_dir}/color_corr_gene_{user_session_id}.csv"  # 35
    final_gene_list_path = f"result_data{user_dir}/final_gene_list_{user_session_id}.csv"  # 36
    file_path_ssgsea_p = f"result_data{user_dir}/1000_lowest_gene_{user_session_id}.csv"  # 37
    gene_file = f"result_data{user_dir}/uploaded_gene_list_{user_session_id}.csv"  # 38
    not_database_genes = f"result_data{user_dir}/not_dataset_genes_{user_session_id}.csv"  # 39
    corrected_genes_list = f"result_data{user_dir}/corrected_genes_{user_session_id}.csv"  # 40
    cluster_color_path_streamlit = f"result_data{user_dir}/cluster_colors_{user_session_id}.csv"  # 41
    anno_colors_path = f"result_data{user_dir}/ht_annotation_colors_{user_session_id}.csv"  # 42
    multi_db_name_sg = f"result_data{user_dir}/multi_db_single_gene_{user_session_id}.csv"  # 43
    multi_info_name_sg = f"result_data{user_dir}/multi_var_info_single_gene_{user_session_id}.tsv"  # 44
    forest_pdf_sg = f"result_data{user_dir}/forest_single_gene_{user_session_id}.pdf"  # 45
    forest_png_sg = f"result_data{user_dir}/forest_single_gene_{user_session_id}.png"  # 46
    surv_colors_path = f"result_data{user_dir}/surv_single_gene_colors_{user_session_id}.csv"  # 47
    surv_dataframe_path = f"result_data{user_dir}/surv_single_gene_{user_session_id}.csv"  # 48
    surv_pdf_sg = f"result_data{user_dir}/surv_plot_single_gene_{user_session_id}.pdf"  # 49
    surv_png_sg = f"result_data{user_dir}/surv_plot_single_gene_{user_session_id}.png"  # 50
    sg_dot_name = f"result_data{user_dir}/sg_dot_name_{user_session_id}.pdf"  # 51
    pv_dot_name = f"result_data{user_dir}/p_value_sg_{user_session_id}.pdf"  # 52
    ssgsea_dot_path = f"result_data{user_dir}/ssgsea_dot_{user_session_id}.pdf"  # 53
    ht_expression_path = f"result_data{user_dir}/hm_exp_data_ssgsea_{user_session_id}.csv"  # 54
    anno_colors_path2 = f"result_data{user_dir}/ht_annotation_colors2_{user_session_id}.csv"  # 55
    ht_expression_path_ssgsea = f"result_data{user_dir}/hm_exp_data_ssgsea2_{user_session_id}.csv"  # 56
    ht_top_annotation_path2 = f"result_data{user_dir}/anno2_{user_session_id}.csv"  # 57
    ssgsea_row_cluster = f"result_data{user_dir}/ssgsea_row_{user_session_id}.csv"  # 58
    ht_png = f"result_data{user_dir}/ht_png_ssgsea_{user_session_id}.png"  # 59
    ht_pdf = f"result_data{user_dir}/ht_pdf_ssgsea_{user_session_id}.pdf"  # 60
    not_dict_genes = f"result_data{user_dir}/not_dict_genes_{user_session_id}.csv" # 61
    survival_info_cluster = f"result_data{user_dir}/survival_info_{user_session_id}.csv" # 62
    survival_info_ssgsea = f"result_data{user_dir}/survival_info_ssgsea_{user_session_id}.csv" # 63
    survival_info_sg = f"result_data{user_dir}/survival_info_sg_{user_session_id}.csv" # 64
    survival_info_cluster_pw = f"result_data{user_dir}/survival_info_pw_{user_session_id}.csv" # 65
    survival_info_ssgsea_pw = f"result_data{user_dir}/survival_info_ssgsea_pw_{user_session_id}.csv" # 66
    survival_info_sg_pw = f"result_data{user_dir}/survival_info_sg_pw_{user_session_id}.csv" # 67

    names_list = [
        selected_genes_name, selected_genes_name_nmf, anno_file_name, selected_genes_dm_name, cluster_name,
        surv_df_name, gene_cluster, heatmap_png_name, surv_plot_png_name, heatmap_pdf_name,
        surv_plot_pdf_name, multi_db_name, cont_name, comp_name, scatter_2d_name, scatter_3d_name,
        nmf_cons_name, sub_hm_name, chi_hm_name, chi_bp_name, scatter_2d_name_gene, scatter_3d_name_gene,
        multi_info_name, forest_pdf, forest_png, gene_sig_name, strip_ssgsea, corr_ssgsea, ssgsea_result_path,
        anno_color_corr, drug_corr_path, selected_genes_dm_name_gene, hm_ssgsea_groups_path,
        user_gene_cluster_path, compare_violin_ssgsea, anno_color_corr_gene, final_gene_list_path,
        file_path_ssgsea_p, gene_file, not_database_genes, corrected_genes_list, cluster_color_path_streamlit,
        anno_colors_path, multi_db_name_sg, multi_info_name_sg, forest_pdf_sg, forest_png_sg, surv_colors_path,
        surv_dataframe_path, surv_pdf_sg, surv_png_sg, sg_dot_name, pv_dot_name, ssgsea_dot_path,
        ht_expression_path, anno_colors_path2, ht_expression_path_ssgsea, ht_top_annotation_path2,
        ssgsea_row_cluster, ht_png, ht_pdf, not_dict_genes, survival_info_cluster, survival_info_ssgsea,
        survival_info_sg, survival_info_cluster_pw, survival_info_ssgsea_pw, survival_info_sg_pw]

    return names_list
