import pandas as pd
import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import numpy as np
import kaleido
import random
import string
import math
from itertools import cycle, combinations
from scipy.stats import mannwhitneyu


def george_plotly(annotation, df_names):
    """
    Generate data for heatmap for George-SCLC dataset.

    Args:
        annotation (pd.DataFrame): Annotation DataFrame with metadata.
        df_names (list): List of sample names to include.

    Returns:
        tuple: Transformed annotation DataFrame, heatmap data, text data, and custom colorscale.
    """
    annotation = annotation[annotation.index.isin(df_names)]
    annotation = annotation[["NAPY", "NE", "Smoker", "Sex", "UICC"]]
    annotation.fillna('nan', inplace=True, axis=1)
    annotation_num = annotation.copy()
    annotation_num = annotation_num.replace(['ASCL1', 'NEUROD1', 'POU2F3', 'YAP1'], ['0', '1', '2', '3'])
    annotation_num = annotation_num.replace(['non-NE', 'NE'], ['4', '5'])
    annotation_num = annotation_num.replace(['Current', 'Former', 'Never', 'NA_Smoker'], ['6', '7', '8', '9'])
    annotation_num = annotation_num.replace(['Male', 'Female'], ['10', '11'])
    annotation_num = annotation_num.replace(['I', 'IA', 'IB', 'II', 'IIA', 'IIB', 'IIIA', 'IIIB', 'IV', 'NA_UICC'],
                                            ['12', '13', '14', '15', '16', '17', '18', '19', '20', '21'])
    annotation_num = annotation_num.astype('int32')
    annotation_num[["NAPY_STR", "NE_SCORE_STR", "SMOKER_STR", "SEX_STR", 'UICC_STR']] = annotation[["NAPY", "NE",
                                                                                                    "Smoker", "Sex",
                                                                                                    "UICC"]]
    annotation_num = annotation_num.replace(['ASCL1', 'NEUROD1', 'POU2F3', 'YAP1'], ['A', 'N', 'P', 'Y'])
    annotation_num = annotation_num.replace(['non-NE'], ['n-NE'])
    annotation_num = annotation_num.replace(['Male', 'Female'], ['M', 'F'])
    annotation_num = annotation_num.replace(['Current', 'Former', 'Never'], ['C', 'F', 'N'])
    annotation_num = annotation_num.sort_values(by=["NAPY", "NE"])

    unique_values = sorted(pd.unique(annotation_num[["NAPY", "NE", "Smoker", "Sex", "UICC"]].values.ravel()))
    mapping_dict = {old_value: new_value for new_value, old_value in enumerate(unique_values)}
    for column in ["NAPY", "NE", "Smoker", "Sex", "UICC"]:
        annotation_num[column] = annotation_num[column].map(mapping_dict)

    z = annotation_num[["NAPY", "NE", "Smoker", "Sex", "UICC"]].transpose().values
    z_text = annotation_num[["NAPY_STR", "NE_SCORE_STR", "SMOKER_STR", "SEX_STR", "UICC_STR"]].transpose().values

    # This is your color map:
    color_dict = {0: "#FD8A8A", 1: "#8FEBC8", 2: "#8E60DB", 3: "#EFFA7A", 4: "#84CBD9", 5: "#F9C65F",
                  6: "#F2A359", 7: "#64A5DE", 8: "#98D19C", 9: "#a9a9ab", 10: "#7df8ff", 11: "#c79558",
                  12: "#7AB1E7", 13: "#8FBCEC", 14: "#A4C7F1", 15: "#B9D2F6", 16: "#CEDDFB",
                  17: "#E3E8FF", 18: "#D1C2D2", 19: "#BF9CAB", 20: "#AD769E", 21: "#a9a9ab"}

    custom_color_dict = {mapping_dict[k]: color_dict[k] for k in color_dict.keys() if k in unique_values}
    custom_colorscale = list(custom_color_dict.values())

    return annotation_num, z, z_text, custom_colorscale


def liu_plotly(annotation, df_names):
    """
    Generate data for heatmap for Liu-SCLC dataset.

    Args:
        annotation (pd.DataFrame): Annotation DataFrame with metadata.
        df_names (list): List of sample names to include.

    Returns:
        tuple: Transformed annotation DataFrame, heatmap data, text data, and custom colorscale.
    """
    annotation = annotation[annotation.index.isin(df_names)]
    annotation = annotation[["NAPY", "NE", "Smoke", "Gender", "TNM_Stage", "EMT"]]
    annotation.fillna('nan', inplace=True, axis=1)
    annotation_num = annotation.copy()
    annotation_num = annotation_num.replace(['ASCL1', 'NEUROD1', 'POU2F3', 'YAP1'], ['0', '1', '2', '3'])
    annotation_num = annotation_num.replace(['non-NE', 'NE'], ['4', '5'])
    annotation_num = annotation_num.replace(['yes', 'no', 'unknown'], ['6', '8', '9'])
    annotation_num = annotation_num.replace(['Male', 'Female'], ['10', '11'])
    annotation_num = annotation_num.replace(['Ia2', 'Ia3', 'Ib', 'IIa', 'IIb', 'IIb', 'IIIa', 'IIIb', 'IV', 'nan'],
                                            ['12', '13', '14', '15', '16', '17', '18', '19', '20', '21'])
    annotation_num = annotation_num.replace(['Epithelial', 'Mesenchymal'], ['22', '23'])
    annotation_num = annotation_num.astype('int32')
    annotation_num[["NAPY_STR", "NE_SCORE_STR", "SMOKE_STR", 'GENDER_STR', "TNM_STR", "EMT_STR"]] = annotation[["NAPY",
                                                                                                                "NE",
                                                                                                                "Smoke",
                                                                                                                "Gender",
                                                                                                                "TNM_Stage",
                                                                                                                "EMT"]]
    annotation_num = annotation_num.replace(['ASCL1', 'NEUROD1', 'POU2F3', 'YAP1'], ['A', 'N', 'P', 'Y'])
    annotation_num = annotation_num.replace(['non-NE'], ['n-NE'])
    annotation_num = annotation_num.replace(['Male', 'Female'], ['M', 'F'])
    annotation_num = annotation_num.replace(['Current', 'Former', 'Never'], ['C', 'F', 'N'])
    annotation_num = annotation_num.replace(['Epithelial', 'Mesenchymal'], ['E', 'M'])
    annotation_num = annotation_num.sort_values(by=["NAPY", "NE"])

    unique_values = sorted(pd.unique(annotation_num[["NAPY", "NE", "Smoke", "Gender", "TNM_Stage",
                                                     "EMT"]].values.ravel()))
    mapping_dict = {old_value: new_value for new_value, old_value in enumerate(unique_values)}
    for column in ["NAPY", "NE", "Smoke", "Gender", "TNM_Stage", "EMT"]:
        annotation_num[column] = annotation_num[column].map(mapping_dict)

    z = annotation_num[["NAPY", "NE", "Smoke", "Gender", "TNM_Stage", "EMT"]].transpose().values
    z_text = annotation_num[
        ["NAPY_STR", "NE_SCORE_STR", "SMOKE_STR", 'GENDER_STR', "TNM_STR", "EMT_STR"]].transpose().values

    # This is your color map:
    color_dict = {0: "#FD8A8A", 1: "#8FEBC8", 2: "#8E60DB", 3: "#EFFA7A", 4: "#84CBD9", 5: "#F9C65F",
                  6: "#F2A359", 7: "#64A5DE", 8: "#98D19C", 9: "#a9a9ab", 10: "#7df8ff", 11: "#c79558",
                  12: "#7AB1E7", 13: "#8FBCEC", 14: "#A4C7F1", 15: "#B9D2F6", 16: "#CEDDFB",
                  17: "#E3E8FF", 18: "#D1C2D2", 19: "#BF9CAB", 20: "#AD769E", 21: "#a9a9ab", 22: "#9DE481",
                  23: "#674E62"}

    custom_color_dict = {mapping_dict[k]: color_dict[k] for k in color_dict.keys() if k in unique_values}
    custom_colorscale = list(custom_color_dict.values())

    return annotation_num, z, z_text, custom_colorscale


def jiang_plotly(annotation, df_names):
    """Generate data for heatmap for Jiang-SCLC dataset.

    Args:
        annotation (pd.DataFrame): Annotation DataFrame with metadata.
        df_names (list): List of sample names to include.

    Returns:
        tuple: Transformed annotation DataFrame, heatmap data, text data, and custom colorscale."""

    annotation = annotation[annotation.index.isin(df_names)]
    annotation = annotation[["NAPY", "NE", 'Gender', 'Smoke', 'UICC stage', 'EMT', 'Chemotherapy']]
    annotation.fillna('nan', inplace=True, axis=1)
    annotation_num = annotation.copy()
    annotation_num = annotation_num.replace(['ASCL1', 'NEUROD1', 'POU2F3', 'YAP1'], ['0', '1', '2', '3'])
    annotation_num = annotation_num.replace(['non-NE', 'NE'], ['4', '5'])
    annotation_num = annotation_num.replace(['Female', 'Male'], ['6', '7'])
    annotation_num = annotation_num.replace(['No', 'Yes'], ['9', '10'])
    annotation_num = annotation_num.replace(['I', 'IA', 'IB', 'II', 'IIA', 'IIB', 'IIIA', 'IIIB', 'IV'],
                                            ['12', '13', '14', '15', '16', '17', '18', '19', '20'])
    annotation_num = annotation_num.replace(['Epithelial', 'Mesenchymal'], ['22', '23'])
    annotation_num = annotation_num.replace(['Naive', 'Chemo'], ['25', '26'])
    annotation_num = annotation_num.replace(['IIIA/IIIB'], ['28'])

    annotation_num = annotation_num.astype('int')
    annotation_num[["NAPY_STR", "NE_SCORE_STR", "SEX_STR", "SMOKER_STR", 'UICC_STR', 'EMT_STR', 'CHEMO_STR']] = \
        annotation[["NAPY", "NE", 'Gender', 'Smoke', 'UICC stage', 'EMT', 'Chemotherapy']]
    annotation_num = annotation_num.replace(['ASCL1', 'NEUROD1', 'POU2F3', 'YAP1'], ['A', 'N', 'P', 'Y'])
    annotation_num = annotation_num.replace(['non-NE'], ['n-NE'])
    annotation_num = annotation_num.replace(['Male', 'Female'], ['M', 'F'])
    annotation_num = annotation_num.replace(['Yes', 'No'], ['Y', 'N'])
    annotation_num = annotation_num.sort_values(by=["NAPY", "NE"])
    annotation_num = annotation_num.replace(['Epithelial', 'Mesenchymal'], ['E', 'M'])
    annotation_num = annotation_num.replace(['Naive', 'Chemo'], ['N', 'C'])

    unique_values = sorted(
        pd.unique(
            annotation_num[["NAPY", "NE", 'Gender', 'Smoke', 'UICC stage', 'EMT', 'Chemotherapy']].values.ravel()))
    mapping_dict = {old_value: new_value for new_value, old_value in enumerate(unique_values)}
    for column in ["NAPY", "NE", 'Gender', 'Smoke', 'UICC stage', 'EMT', 'Chemotherapy']:
        annotation_num[column] = annotation_num[column].map(mapping_dict)

    z = annotation_num[["NAPY", "NE", 'Gender', 'Smoke', 'UICC stage', 'EMT', 'Chemotherapy']].transpose().values
    z_text = annotation_num[
        ["NAPY_STR", "NE_SCORE_STR", "SEX_STR", "SMOKER_STR", 'UICC_STR', 'EMT_STR', 'CHEMO_STR']].transpose().values
    # This is your color map:
    color_dict = {0: "#FD8A8A", 1: "#8FEBC8", 2: "#8E60DB", 3: "#EFFA7A", 4: "#84CBD9", 5: "#F9C65F", 6: "#F2A359",
                  7: "#64A5DE", 9: "#98D19C", 10: "#F2A359", 12: "#7AB1E7", 13: "#8FBCEC",
                  14: "#A4C7F1", 15: "#B9D2F6", 16: "#CEDDFB", 17: "#E3E8FF", 18: "#D1C2D2", 19: "#BF9CAB",
                  20: "#AD769E", 22: "#9DE481", 23: "#674E62", 25: "#98D19C",
                  26: "#F2A359", 28: "#c0a7c2"}

    custom_color_dict = {mapping_dict[k]: color_dict[k] for k in color_dict.keys() if k in unique_values}
    custom_colorscale = list(custom_color_dict.values())

    return annotation_num, z, z_text, custom_colorscale


def lcnec_plotly(annotation, df_names):
    """
    Generate data for heatmap for George-LCNEC dataset.

    Args:
        annotation (pd.DataFrame): Annotation DataFrame with metadata.
        df_names (list): List of sample names to include.

    Returns:
        tuple: Transformed annotation DataFrame, heatmap data, text data, and custom colorscale.
    """
    annotation = annotation[annotation.index.isin(df_names)]
    annotation = annotation[["NAPY", "NE", 'Gender', 'UICC Stage', 'EMT']]
    annotation.fillna('nan', inplace=True, axis=1)
    annotation_num = annotation.copy()
    annotation_num = annotation_num.replace(['ASCL1', 'NEUROD1', 'POU2F3', 'YAP1'], ['0', '1', '2', '3'])
    annotation_num = annotation_num.replace(['non-NE', 'NE'], ['4', '5'])
    annotation_num = annotation_num.replace(['Female', 'Male'], ['6', '7'])
    annotation_num = annotation_num.replace(['Ia', 'Ib', 'IIa', 'IIb', 'III', 'IIIa', 'IIIb', 'IV', 'NA_UICC'],
                                            ['12', '13', '14', '15', '17', '18', '19', '20', '21'])
    annotation_num = annotation_num.replace(['Epithelial', 'Mesenchymal'], ['22', '23'])

    annotation_num = annotation_num.astype('float').astype('int')
    annotation_num[["NAPY_STR", "NE_SCORE_STR", "SEX_STR", 'UICC_STR', 'EMT_STR']] = \
        annotation[["NAPY", "NE", 'Gender', 'UICC Stage', 'EMT']]
    annotation_num = annotation_num.replace(['ASCL1', 'NEUROD1', 'POU2F3', 'YAP1'], ['A', 'N', 'P', 'Y'])
    annotation_num = annotation_num.replace(['non-NE'], ['n-NE'])
    annotation_num = annotation_num.replace(['Male', 'Female'], ['M', 'F'])
    annotation_num = annotation_num.sort_values(by=["NAPY", "NE"])
    annotation_num = annotation_num.replace(['Epithelial', 'Mesenchymal'], ['E', 'M'])

    unique_values = sorted(
        pd.unique(
            annotation_num[["NAPY", "NE", 'Gender', 'UICC Stage', 'EMT']].values.ravel()))
    mapping_dict = {old_value: new_value for new_value, old_value in enumerate(unique_values)}
    for column in ["NAPY", "NE", 'Gender', 'UICC Stage', 'EMT']:
        annotation_num[column] = annotation_num[column].map(mapping_dict)

    z = annotation_num[["NAPY", "NE", 'Gender', 'UICC Stage', 'EMT']].transpose().values
    z_text = annotation_num[["NAPY_STR", "NE_SCORE_STR", "SEX_STR", 'UICC_STR', 'EMT_STR']].transpose().values
    # This is your color map:
    color_dict = {0: "#FD8A8A", 1: "#8FEBC8", 2: "#8E60DB", 3: "#EFFA7A", 4: "#84CBD9", 5: "#F9C65F", 6: "#F2A359",
                  7: "#64A5DE", 9: "#98D19C", 10: "#F2A359", 12: "#7AB1E7", 13: "#8FBCEC",
                  14: "#A4C7F1", 15: "#B9D2F6", 16: "#CEDDFB", 17: "#E3E8FF", 18: "#D1C2D2", 19: "#BF9CAB",
                  20: "#AD769E", 21: "#a9a9ab", 22: "#9DE481", 23: "#674E62", 25: "#98D19C",
                  26: "#F2A359", 28: "#c0a7c2"}

    custom_color_dict = {mapping_dict[k]: color_dict[k] for k in color_dict.keys() if k in unique_values}
    custom_colorscale = list(custom_color_dict.values())

    return annotation_num, z, z_text, custom_colorscale


def carcinoid_plotly(annotation, df_names):
    """
    Generate data for heatmap for Fernandez-Carcinoid dataset.

    Args:
        annotation (pd.DataFrame): Annotation DataFrame with metadata.
        df_names (list): List of sample names to include.

    Returns:
        tuple: Transformed annotation DataFrame, heatmap data, text data, and custom colorscale.
    """

    annotation = annotation[annotation.index.isin(df_names)]
    annotation = annotation[["NAPY", "NE", 'Gender', 'UICC Stage', 'EMT']]
    annotation.fillna('nan', inplace=True, axis=1)
    annotation_num = annotation.copy()
    annotation_num = annotation_num.replace(['ASCL1', 'NEUROD1', 'POU2F3', 'YAP1'], ['0', '1', '2', '3'])
    annotation_num = annotation_num.replace(['non-NE', 'NE'], ['4', '5'])
    annotation_num = annotation_num.replace(['Female', 'Male'], ['6', '7'])
    annotation_num = annotation_num.replace(['I', 'IA', 'IA1', 'IA2', 'IA3', 'IB', 'IIA', 'IIB', 'IIIA', 'NA_UICC'],
                                            ['12', '13', '14', '15', '17', '18', '19', '24', '20', '21'])
    annotation_num = annotation_num.replace(['Epithelial', 'Mesenchymal'], ['22', '23'])

    annotation_num = annotation_num.astype('float').astype('int')
    annotation_num[["NAPY_STR", "NE_SCORE_STR", "SEX_STR", 'UICC_STR', 'EMT_STR']] = \
        annotation[["NAPY", "NE", 'Gender', 'UICC Stage', 'EMT']]
    annotation_num = annotation_num.replace(['ASCL1', 'NEUROD1', 'POU2F3', 'YAP1'], ['A', 'N', 'P', 'Y'])
    annotation_num = annotation_num.replace(['non-NE'], ['n-NE'])
    annotation_num = annotation_num.replace(['Male', 'Female'], ['M', 'F'])
    annotation_num = annotation_num.sort_values(by=["NAPY", "NE"])
    annotation_num = annotation_num.replace(['Epithelial', 'Mesenchymal'], ['E', 'M'])

    unique_values = sorted(
        pd.unique(
            annotation_num[["NAPY", "NE", 'Gender', 'UICC Stage', 'EMT']].values.ravel()))
    mapping_dict = {old_value: new_value for new_value, old_value in enumerate(unique_values)}
    for column in ["NAPY", "NE", 'Gender', 'UICC Stage', 'EMT']:
        annotation_num[column] = annotation_num[column].map(mapping_dict)

    z = annotation_num[["NAPY", "NE", 'Gender', 'UICC Stage', 'EMT']].transpose().values
    z_text = annotation_num[["NAPY_STR", "NE_SCORE_STR", "SEX_STR", 'UICC_STR', 'EMT_STR']].transpose().values
    # This is your color map:
    color_dict = {0: "#FD8A8A", 1: "#8FEBC8", 2: "#8E60DB", 3: "#EFFA7A", 4: "#84CBD9", 5: "#F9C65F", 6: "#F2A359",
                  7: "#64A5DE", 9: "#98D19C", 10: "#F2A359", 12: "#7AB1E7", 13: "#8FBCEC",
                  14: "#A4C7F1", 15: "#B9D2F6", 16: "#B9D2F6", 17: "#CEDDFB", 18: "#E3E8FF", 19: "#D1C2D2",
                  20: "#AD769E", 21: "#a9a9ab", 22: "#9DE481", 23: "#674E62", 24: "#BF9CAB", 25: "#98D19C",
                  26: "#F2A359", 28: "#c0a7c2"}

    custom_color_dict = {mapping_dict[k]: color_dict[k] for k in color_dict.keys() if k in unique_values}
    custom_colorscale = list(custom_color_dict.values())

    return annotation_num, z, z_text, custom_colorscale


def alcala_plotly(annotation, df_names):
    """
    Generate data for heatmap for Alcala-Carcinoid dataset.

    Args:
        annotation (pd.DataFrame): Annotation DataFrame with metadata.
        df_names (list): List of sample names to include.

    Returns:
        tuple: Transformed annotation DataFrame, heatmap data, text data, and custom colorscale.
    """
    annotation = annotation[annotation.index.isin(df_names)]
    annotation = annotation[["NAPY", "NE", 'Gender', 'EMT']]
    annotation.fillna('nan', inplace=True, axis=1)
    annotation_num = annotation.copy()
    annotation_num = annotation_num.replace(['ASCL1', 'NEUROD1', 'POU2F3', 'YAP1'], ['0', '1', '2', '3'])
    annotation_num = annotation_num.replace(['non-NE', 'NE'], ['4', '5'])
    annotation_num = annotation_num.replace(['Female', 'Male'], ['6', '7'])
    annotation_num = annotation_num.replace(['Epithelial', 'Mesenchymal'], ['22', '23'])

    annotation_num = annotation_num.astype('float').astype('int')
    annotation_num[["NAPY_STR", "NE_SCORE_STR", "SEX_STR", 'EMT_STR']] = \
        annotation[["NAPY", "NE", 'Gender', 'EMT']]
    annotation_num = annotation_num.replace(['ASCL1', 'NEUROD1', 'POU2F3', 'YAP1'], ['A', 'N', 'P', 'Y'])
    annotation_num = annotation_num.replace(['non-NE'], ['n-NE'])
    annotation_num = annotation_num.replace(['Male', 'Female'], ['M', 'F'])
    annotation_num = annotation_num.sort_values(by=["NAPY", "NE"])
    annotation_num = annotation_num.replace(['Epithelial', 'Mesenchymal'], ['E', 'M'])

    unique_values = sorted(
        pd.unique(
            annotation_num[["NAPY", "NE", 'Gender', 'EMT']].values.ravel()))
    mapping_dict = {old_value: new_value for new_value, old_value in enumerate(unique_values)}
    for column in ["NAPY", "NE", 'Gender', 'EMT']:
        annotation_num[column] = annotation_num[column].map(mapping_dict)

    z = annotation_num[["NAPY", "NE", 'Gender', 'EMT']].transpose().values
    z_text = annotation_num[["NAPY_STR", "NE_SCORE_STR", "SEX_STR", 'EMT_STR']].transpose().values
    # This is your color map:
    color_dict = {0: "#FD8A8A", 1: "#8FEBC8", 2: "#8E60DB", 3: "#EFFA7A", 4: "#84CBD9", 5: "#F9C65F", 6: "#F2A359",
                  7: "#64A5DE", 22: "#9DE481", 23: "#674E62"}

    custom_color_dict = {mapping_dict[k]: color_dict[k] for k in color_dict.keys() if k in unique_values}
    custom_colorscale = list(custom_color_dict.values())

    return annotation_num, z, z_text, custom_colorscale


def rousseaux_plotly(annotation, df_names):
    """
    Generate data for heatmap for Rousseaux-Mixed(non-NE) dataset.

    Args:
        annotation (pd.DataFrame): Annotation DataFrame with metadata.
        df_names (list): List of sample names to include.

    Returns:
        tuple: Transformed annotation DataFrame, heatmap data, text data, and custom colorscale.
    """
    annotation = annotation[annotation.index.isin(df_names)]
    annotation = annotation[["NAPY", "NE", 'Gender', 'EMT']]
    annotation[["NAPY", "NE", 'Gender', 'EMT']] = annotation[["NAPY", "NE", 'Gender', 'EMT']].replace('NTL', np.nan)
    annotation.fillna('nan', inplace=True, axis=1)

    annotation_num = annotation.copy()
    annotation_num = annotation_num.replace(['ASCL1', 'NEUROD1', 'POU2F3', 'YAP1'], ['0', '1', '2', '3'])
    annotation_num = annotation_num.replace(['non-NE', 'NE'], ['4', '5'])
    annotation_num = annotation_num.replace(['F', 'M'], ['6', '7'])
    annotation_num = annotation_num.replace(['Epithelial', 'Mesenchymal'], ['22', '23'])
    annotation_num = annotation_num.replace(['NTL_Gender'], ['24'])

    annotation_num = annotation_num.astype('float').astype('int')
    annotation_num[["NAPY_STR", "NE_SCORE_STR", "SEX_STR", 'EMT_STR']] = \
        annotation[["NAPY", "NE", 'Gender', 'EMT']]
    annotation_num = annotation_num.replace(['ASCL1', 'NEUROD1', 'POU2F3', 'YAP1'], ['A', 'N', 'P', 'Y'])
    annotation_num = annotation_num.replace(['non-NE'], ['n-NE'])
    annotation_num = annotation_num.sort_values(by=["NAPY", "NE"])
    annotation_num = annotation_num.replace(['Epithelial', 'Mesenchymal'], ['E', 'M'])

    unique_values = sorted(
        pd.unique(
            annotation_num[["NAPY", "NE", 'Gender', 'EMT']].values.ravel()))
    mapping_dict = {old_value: new_value for new_value, old_value in enumerate(unique_values)}
    for column in ["NAPY", "NE", 'Gender', 'EMT']:
        annotation_num[column] = annotation_num[column].map(mapping_dict)

    z = annotation_num[["NAPY", "NE", 'Gender', 'EMT']].transpose().values
    z_text = annotation_num[["NAPY_STR", "NE_SCORE_STR", "SEX_STR", 'EMT_STR']].transpose().values
    # This is your color map:
    color_dict = {0: "#FD8A8A", 1: "#8FEBC8", 2: "#8E60DB", 3: "#EFFA7A", 4: "#84CBD9", 5: "#F9C65F", 6: "#F2A359",
                  7: "#64A5DE", 22: "#9DE481", 23: "#674E62", 24: "#a9a9ab"}

    custom_color_dict = {mapping_dict[k]: color_dict[k] for k in color_dict.keys() if k in unique_values}
    custom_colorscale = list(custom_color_dict.values())

    return annotation_num, z, z_text, custom_colorscale


def anish_plotly(annotation, df_names):
    """
    Generate data for heatmap for Anish-SCLC dataset.

    Args:
        annotation (pd.DataFrame): Annotation DataFrame with metadata.
        df_names (list): List of sample names to include.

    Returns:
        tuple: Transformed annotation DataFrame, heatmap data, text data, and custom colorscale.
    """
    annotation = annotation[annotation.index.isin(df_names)]
    annotation = annotation[["NAPY", "NE", "Sex", "Smoke", "Tumor_Stage"]]
    annotation.fillna('nan', inplace=True, axis=1)
    annotation_num = annotation.copy()
    annotation_num = annotation_num.replace(['ASCL1', 'NEUROD1', 'POU2F3', 'YAP1'], ['0', '1', '2', '3'])
    annotation_num = annotation_num.replace(['non-NE', 'NE'], ['4', '5'])
    annotation_num = annotation_num.replace(['F', 'M', 'NA_Sex'], ['6', '7', '8'])
    annotation_num = annotation_num.replace(['No', 'Yes', 'NA_Smoking'], ['9', '10', '11'])
    annotation_num = annotation_num.replace(['Extensive stage', 'Limited stage'], ['12', '13'])

    annotation_num = annotation_num.astype('int32')

    annotation_num[["NAPY_STR", "NE_SCORE_STR", "SEX_STR", "SMOKE_STR", 'TS_STR']] = \
        annotation[["NAPY", "NE", "Sex", "Smoke", "Tumor_Stage"]]
    annotation_num = annotation_num.replace(['ASCL1', 'NEUROD1', 'POU2F3', 'YAP1'], ['A', 'N', 'P', 'Y'])
    annotation_num = annotation_num.replace(['non-NE'], ['n-NE'])
    annotation_num = annotation_num.replace(['Extensive stage', 'Limited stage'], ['Ext', 'Lim'])
    annotation_num = annotation_num.replace(['Yes', 'No'], ['Y', 'N'])
    annotation_num = annotation_num.sort_values(by=["NAPY", "NE"])

    unique_values = sorted(pd.unique(annotation_num[["NAPY", "NE", "Sex", "Smoke", "Tumor_Stage"]].values.ravel()))
    mapping_dict = {old_value: new_value for new_value, old_value in enumerate(unique_values)}
    for column in ["NAPY", "NE", "Sex", "Smoke", "Tumor_Stage"]:
        annotation_num[column] = annotation_num[column].map(mapping_dict)

    z = annotation_num[["NAPY", "NE", "Sex", "Smoke", "Tumor_Stage"]].transpose().values
    z_text = annotation_num[["NAPY_STR", "NE_SCORE_STR", "SEX_STR", "SMOKE_STR", 'TS_STR']].transpose().values
    # This is your color map:
    color_dict = {0: "#FD8A8A", 1: "#8FEBC8", 2: "#8E60DB", 3: "#EFFA7A", 4: "#84CBD9", 5: "#F9C65F", 6: "#c79558",
                  7: "#7df8ff", 8: "#a9a9ab", 9: "#98D19C", 10: "#F2A359", 11: "#a9a9ab", 12: "#f74e05", 13: "#fcacf0"}

    custom_color_dict = {mapping_dict[k]: color_dict[k] for k in color_dict.keys() if k in unique_values}
    custom_colorscale = list(custom_color_dict.values())

    return annotation_num, z, z_text, custom_colorscale


def plotly_heatmap(annotation, df_names, df_option, sub_hm_name):
    """
    Create and display a heatmap using Plotly based on the dataset option.

    Args:
        annotation (pd.DataFrame): Annotation DataFrame with metadata.
        df_names (list): List of sample names to include.
        df_option (str): Dataset option identifier.
        sub_hm_name (str): Path to save the heatmap image.

    Returns:
        None
    """
    annotation_num = None
    z = None
    z_text = None
    y_text = None
    custom_colorscale = None

    if df_option == "George-SCLC":
        annotation_num, z, z_text, custom_colorscale = george_plotly(annotation, df_names)
        y_text = annotation_num[["NAPY", "NE", "Smoker", "Sex", "UICC"]].columns
    elif df_option == "Anish-SCLC":
        annotation_num, z, z_text, custom_colorscale = anish_plotly(annotation, df_names)
        y_text = annotation_num[["NAPY", "NE", "Sex", "Smoke", "Tumor_Stage"]].columns
    elif df_option == "Liu-SCLC":
        annotation_num, z, z_text, custom_colorscale = liu_plotly(annotation, df_names)
        y_text = annotation_num[["NAPY", "NE", "Smoke", "Gender", "TNM_Stage", "EMT"]].columns
    elif df_option == "George-LCNEC":
        annotation_num, z, z_text, custom_colorscale = lcnec_plotly(annotation, df_names)
        y_text = annotation_num[["NAPY", "NE", 'Gender', 'UICC Stage', 'EMT']].columns
    elif df_option == "Fernandez-Carcinoid":
        annotation_num, z, z_text, custom_colorscale = carcinoid_plotly(annotation, df_names)
        y_text = annotation_num[["NAPY", "NE", 'Gender', 'UICC Stage', 'EMT']].columns
    elif df_option == "Alcala-Carcinoid":
        annotation_num, z, z_text, custom_colorscale = alcala_plotly(annotation, df_names)
        y_text = annotation_num[["NAPY", "NE", 'Gender', 'EMT']].columns
    elif df_option == "Jiang-SCLC":
        annotation_num, z, z_text, custom_colorscale = jiang_plotly(annotation, df_names)
        y_text = annotation_num[["NAPY", "NE", 'Gender', 'Smoke', 'UICC stage', 'EMT', 'Chemotherapy']].columns
    elif df_option == 'Rousseaux-Mixed(non-NE)':
        annotation_num, z, z_text, custom_colorscale = rousseaux_plotly(annotation, df_names)
        y_text = annotation_num[["NAPY", "NE", 'Gender', 'EMT']].columns
    elif df_option == "TCGA":
        y_text = None
        return st.warning(":blue[Annotation preview is ] :red[NOT ] :blue[available at the moment]")

    hm = [go.Heatmap(z=z, text=z_text, y=y_text,
                     hoverinfo='text', colorscale=custom_colorscale, showscale=False,
                     xgap=3, ygap=3, showlegend=False, textfont=dict(family="sans serif", size=9, color="White"))]

    fig = go.Figure(data=hm)
    fig.update_traces(text=z_text, texttemplate="%{text}")

    fig.update_layout(
        autosize=True,
        width=500,
        height=500,
        plot_bgcolor='#FFFFFF',
        showlegend=False,
        xaxis=dict(zeroline=False, showgrid=False, showticklabels=False),
        yaxis=dict(zeroline=False, showgrid=False, showticklabels=True),
    )

    # Update axis properties
    fig.update_xaxes(
        tickfont=dict(color='black'),  # Set tick font color
        title_font=dict(color='black')  # Set title font color
    )
    fig.update_yaxes(
        tickfont=dict(color='black'),  # Set tick font color
        title_font=dict(color='black')  # Set title font color
    )

    fig.write_image(sub_hm_name, width=1920, height=1080, scale=4)
    try:
        st.plotly_chart(fig, use_container_width=True)
    except:
        st.image("style_items/error.svg", use_container_width=True)
        st.warning(":blue[There is problem with the ] :red[chart! ] :blue[Please try another ] :red[settings ] "
                   ":blue[or ] :red[refresh ] :blue[page]")


def clust_scatter(dim_red, clust_labels, sample_names_scatter_2d, scatter_2d_name, sample_type, cluster_toggle, e):
    """
    Create a 2D scatter plot of clustering results using Plotly.

    Args:
        dim_red (np.array): 2D array of dimensionality reduced data.
        clust_labels (np.array): Array of cluster labels.
        sample_names_scatter_2d (list): List of sample names.
        scatter_2d_name (str): Path to save the scatter plot image.
        sample_type (str): Type of samples ('gene' or other).
        cluster_toggle (bool): Whether to toggle cluster labels.
        :param e: st.columns

    Returns:
        None

    """

    if dim_red.shape[0] > 500:
        sct_mark_size = 10
    elif 500 > dim_red.shape[0] > 150:
        sct_mark_size = 12.5
    else:
        sct_mark_size = 15

    # Create a dataframe with the dimensionality reduced results and cluster labels
    if sample_type == "gene":
        title = "Clustering Result of Gene"
    else:
        title = "Clustering Result of Samples"

    df_scutter = pd.DataFrame(
        {'x': dim_red[:, 0], 'y': dim_red[:, 1], 'Clusters': clust_labels, 'Samples': sample_names_scatter_2d})

    # Sort the dataframe by cluster labels and convert them to strings
    df_scutter = df_scutter.sort_values(by=['Clusters'], ascending=True)
    df_scutter["Clusters"] = df_scutter["Clusters"].astype(str)

    if len(df_scutter["Clusters"].unique()) > 9:
        cds = px.colors.qualitative.Light24
    else:
        cds = px.colors.qualitative.Set1

    # Create the scatter plot using Plotly
    fig = px.scatter(df_scutter, x='x', y='y', color="Clusters", render_mode="svg", hover_name='Samples',
                     color_discrete_sequence=cds, width=500, height=500)

    # Customize the plot with marker size, line width and color, plot background color, and title
    fig.update_traces(marker_size=sct_mark_size, marker_line_width=2, marker_line_color="black", marker_opacity=0.85)
    fig.update_layout(plot_bgcolor='#FFFFFF', title=title, scattergap=1.0)
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', tickcolor='black', title_standoff=True,
                     gridcolor='black')
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', tickcolor='black', title_standoff=True,
                     gridcolor='black')

    # Update axis properties
    fig.update_xaxes(
        tickfont=dict(color='black', family='Arial'),  # Set tick font color and font style
        title_font=dict(color='black', family='Arial')  # Set title font color and font style
    )
    fig.update_yaxes(
        tickfont=dict(color='black', family='Arial'),  # Set tick font color and font style
        title_font=dict(color='black', family='Arial')  # Set title font color and font style
    )

    # Update the figure layout with the equal axis ranges and legend font properties
    fig.update_layout(
        xaxis=dict(
            autorange=True
        ),
        yaxis=dict(
            autorange=True
        ),
        autosize=False,
        width=550,  # You can adjust the size to fit your needs
        height=550,
        font=dict(family="Arial"),  # Set Arial font style for the entire plot
        legend_title=dict(
            font=dict(size=20, color="black", family="Arial")  # Set legend title properties
        ),
        legend=dict(
            font=dict(size=16, color="black", family="Arial")  # Set legend font properties
        )
    )

    fig.update_layout(
        legend=dict(
            x=-0.1,  # Position the legend to the right
            y=1,  # Center the legend vertically
            xanchor="right",  # Align the left edge of the legend box to the x position
            yanchor="top",  # Align the center of the legend box to the y position
            font=dict(size=16, color="black", family="Arial"),
            bgcolor="rgba(255,255,255,0.8)",  # Optional: Add a semi-transparent background
            bordercolor="black",  # Optional: Add a border
            borderwidth=1  # Optional: Set border width
        )
    )


    fig.write_image(scatter_2d_name, width=600, height=600, scale=8)
    # Display the plot in a column of the app

    if not cluster_toggle:
        solo = st.columns(3)
        try:
            solo[1].plotly_chart(fig, use_container_width=True)
        except:
            st.image("style_items/error.svg", use_container_width=True)
            st.warning(":blue[There is problem with the ] :red[chart! ] :blue[Please try another ] :red[settings ] "
                       ":blue[or ] :red[refresh ] :blue[page]")
    else:
        if sample_type == "gene":
            try:
                e[1].plotly_chart(fig, use_container_width=True)
            except:
                e[1].image("style_items/error.svg", use_container_width=True)
                e[1].warning(":blue[There is problem with the ] :red[chart! ] :blue[Please try another ] :red[settings ] "
                           ":blue[or ] :red[refresh ] :blue[page]")
        else:
            try:
                e[0].plotly_chart(fig, use_container_width=True)
            except:
                e[0].image("style_items/error.svg", use_container_width=True)
                e[0].warning(
                    ":blue[There is problem with the ] :red[chart! ] :blue[Please try another ] :red[settings ] "
                    ":blue[or ] :red[refresh ] :blue[page]")


# Define function to create a 3D scatter plot with cluster labels
def clust_scatter3d(dim_red, clust_labels, sample_names_scatter_3d, sample_type, cluster_toggle, e):
    """
    Create a 3D scatter plot of clustering results using Plotly.

    Args:
        dim_red (np.array): 3D array of dimensionality reduced data.
        clust_labels (np.array): Array of cluster labels.
        sample_names_scatter_3d (list): List of sample names.
        sample_type (str): Type of samples ('gene' or other).
        cluster_toggle (bool): Whether to toggle cluster labels.

    Returns:
        None
    """
    # Determine scatter marker size based on the number of samples
    if dim_red.shape[0] > 500:
        sct_mark_size = 9
    elif 500 > dim_red.shape[0] > 150:
        sct_mark_size = 10
    else:
        sct_mark_size = 12.5

    # Set the title based on sample type
    title = "Clustering Result of Gene" if sample_type == "gene" else "Clustering Result of Samples"

    # Create DataFrame for scatter plot
    df_scutter = pd.DataFrame(
        {'x': dim_red[:, 0], 'y': dim_red[:, 1], 'z': dim_red[:, 2], 'Clusters': clust_labels,
         'Samples': sample_names_scatter_3d})

    # Configure image export options
    config = {
        'toImageButtonOptions': {
            'format': 'svg',
            'filename': 'custom_image',
            'height': 500,
            'width': 700,
            'scale': 1
        }
    }

    df_scutter = df_scutter.sort_values(by=['Clusters'], ascending=True)
    df_scutter["Clusters"] = df_scutter["Clusters"].astype(str)

    cds = px.colors.qualitative.Light24 if len(df_scutter["Clusters"].unique()) > 9 else px.colors.qualitative.Set1

    # Create the 3D scatter plot
    fig = px.scatter_3d(df_scutter, x='x', y='y', z='z', color="Clusters", hover_name='Samples',
                        color_discrete_sequence=cds, width=1000, height=650)
    fig.update_layout(
        title=title,
        font=dict(family="Arial"),
        legend_title=dict(font=dict(size=20, color="black", family="Arial")),
        legend=dict(font=dict(size=16, color="black", family="Arial")),
        plot_bgcolor='#FFFFFF',
        scene_camera=dict(
        eye=dict(x=3, y=3, z=1)  # Zoomed-out perspective
        )
    )
    fig.update_traces(marker_size=sct_mark_size, marker_line_width=2, marker_line_color="black", marker_opacity=1)

    # Animation for rotation
    frames = []
    duration_seconds = 5
    total_frames = duration_seconds * 600
    total_rotations = 10

    for t in range(total_frames):
        angle = (t / total_frames) * 360 * total_rotations
        camera = dict(
            up=dict(x=0, y=0, z=1),
            center=dict(x=0, y=0, z=0),
            eye=dict(x=2.5 * math.cos(math.radians(angle)), y=2.5 * math.sin(math.radians(angle)), z=0.5)
        )
        frames.append(dict(layout=dict(scene_camera=camera)))

    fig.frames = frames
    fig.update_layout(
        updatemenus=[
            {
                "buttons": [
                    {
                        "args": [None,
                                 {"frame": {"duration": 50, "redraw": True, "autoplay": True}, "fromcurrent": True,
                                  "transition": {"duration": 50}}],
                        "label": "Play",
                        "method": "animate"
                    },
                    {
                        "args": [[None], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate",
                                          "transition": {"duration": 0}}],
                        "label": "Pause",
                        "method": "animate"
                    }
                ],
                "direction": "left",
                "pad": {"r": 10, "t": 87},
                "showactive": True,
                "type": "buttons",
                "x": 0.1,
                "xanchor": "left",
                "y": 0,
                "yanchor": "bottom"
            }
        ]
    )

    fig.update_layout(
        legend=dict(
            x=0,  # Position the legend to the right
            y=1,  # Center the legend vertically
            xanchor="left",  # Align the left edge of the legend box to the x position
            yanchor="top",  # Align the center of the legend box to the y position
            font=dict(size=16, color="black", family="Arial"),
            bgcolor="rgba(255,255,255,0.8)",  # Optional: Add a semi-transparent background
            bordercolor="black",  # Optional: Add a border
            borderwidth=1  # Optional: Set border width
        )
    )

    if not cluster_toggle:
        # Display the plot
        solo = st.columns(3)
        try:
            solo[1].plotly_chart(fig, use_container_width=True, vertical_alignment="bottom", config=config)
        except Exception as ex:
            solo[1].image("style_items/error.svg", use_container_width=True)
            solo[1].warning(":blue[There is problem with the ] :red[chart! ] :blue[Please try another ] :red[settings ] "
                      ":blue[or ] :red[refresh ] :blue[page]")
            print(f"Error: {ex}")
    else:
        if sample_type == "gene":
            try:
                e[1].plotly_chart(fig, use_container_width=True, vertical_alignment="bottom", config=config)
            except Exception as ex:
                e[1].image("style_items/error.svg", use_container_width=True)
                e[1].warning(
                    ":blue[There is problem with the ] :red[chart! ] :blue[Please try another ] :red[settings ] "
                    ":blue[or ] :red[refresh ] :blue[page]")
                print(f"Error: {ex}")
        if sample_type == "sample":
            try:
                e[0].plotly_chart(fig, use_container_width=True, vertical_alignment="bottom", config=config)
            except Exception as ex:
                e[0].image("style_items/error.svg", use_container_width=True)
                e[0].warning(
                    ":blue[There is problem with the ] :red[chart! ] :blue[Please try another ] :red[settings ] "
                    ":blue[or ] :red[refresh ] :blue[page]")
                print(f"Error: {ex}")


def plot_heatmap(p_value_df, chi_hm_name, chi_plot_col):
    """
    Create and display a heatmap of p-values using Plotly.

    Args:
        p_value_df (pd.DataFrame): DataFrame of p-values.
        chi_hm_name (str): Path to save the heatmap image.
        chi_plot_col (list): List of Streamlit columns for plotting.

    Returns:
        None
    """
    p_val_cols = st.columns(3)

    significance_level = p_val_cols[1].slider(":blue[Set the ] :red[Level of Significance]", min_value=0.01,
                                              max_value=0.05, value=0.05, key='sig_lvl')

    # Set the diagonal and upper triangular elements to NaN
    for i in range(p_value_df.shape[0]):
        for j in range(i, p_value_df.shape[1]):
            p_value_df.iloc[i, j] = np.nan
            if p_value_df.iloc[i, j] == 0:
                p_value_df.iloc[i, j] = np.nan

    # Reverse the order of the rows
    p_value_df = p_value_df.iloc[::-1]

    # Define custom colorscale
    custom_color_scale = [
        [0.0, 'red'],
        [0.25, 'orange'],
        [0.5, 'yellow'],
        [0.75, 'blue'],
        [1.0, 'white']
    ]

    # Create a heatmap using Plotly Express
    fig = px.imshow(p_value_df,
                    origin='lower',
                    range_color=[0, significance_level],
                    color_continuous_scale=custom_color_scale,
                    title='Heatmap of Contingency Matrix')

    # Annotate non-NaN values
    for i in range(p_value_df.shape[0]):
        for j in range(p_value_df.shape[1]):
            value = p_value_df.iloc[i, j]
            if not pd.isna(value):
                text = "NS" if value >= significance_level else "{:.2e}".format(value)
                fig.add_annotation(
                    x=j, y=i,
                    text=text,
                    showarrow=False,
                    font=dict(size=16, color='black', family='Arial'),
                    align="center",
                    xref="x",
                    yref="y",
                    ax=0,
                    ay=0
                )

    # Update axis properties
    fig.update_xaxes(
        title_text='Cluster',
        nticks=len(p_value_df.columns) + 1,
        tickfont=dict(color='black', family='Arial', size=14),
        title_font=dict(color='black', family='Arial', size=16)
    )
    fig.update_yaxes(
        title_text='Cluster',
        nticks=len(p_value_df.index) + 1,
        tickfont=dict(color='black', family='Arial', size=14),
        title_font=dict(color='black', family='Arial', size=16)
    )

    fig.update_layout(
        font=dict(family="Arial"),
        coloraxis_colorbar=dict(
            tickvals=[0, significance_level / 2, significance_level],
            ticktext=[str(0), str(significance_level / 2), str(significance_level)],
            tickfont=dict(color='black', family='Arial', size=14),
            titlefont=dict(color='black', family='Arial', size=16)
        ),
        plot_bgcolor='lightgrey',
        xaxis=dict(showgrid=False),
        yaxis=dict(showgrid=False),
        title=dict(font=dict(size=20, color='black', family='Arial'))
    )

    # Show the plot
    try:
        chi_plot_col[0].plotly_chart(fig, use_container_width=True)
    except Exception as ex:
        chi_plot_col[0].image("style_items/error.svg", use_container_width=True)
        chi_plot_col[0].warning(
            ":blue[There is problem with the ] :red[chart! ] :blue[Please try another ] :red[settings ] "
            ":blue[or ] :red[refresh ] :blue[page]")
        print(f"Error: {ex}")

    # Save the plot
    try:
        fig.write_image(chi_hm_name, width=1200, height=1000, scale=1)
    except Exception as ex:
        print(f"Error saving image: {ex}")


def plot_gene_correlation(x_gene, y_gene, x_gene_name, y_gene_name, corr_coef, p_value, color_palette, mono_color,
                          merged_df, comp_type, user_id):
    """
    Generates a scatter plot to visualize the correlation between the expression levels of two genes,
    optionally coloring the dots based on a selected palette.

    Parameters:
    - x_gene: Expression data for the x-axis gene.
    - y_gene: Expression data for the y-axis gene.
    - x_gene_name: Name of the x-axis gene.
    - y_gene_name: Name of the y-axis gene.
    - corr_coef: Correlation coefficient for annotation.
    - p_value: P-value for annotation.
    - color_palette: Selected coloring palette or 'Monochrome'.
    - mono_color: Selected color for monochrome plots.
    - merged_df: DataFrame with sample annotations, including the 'NAPY' column.
    - comp_type: Comparison type (e.g., "ssgsea" or "gc_ssgsea").
    - user_id: User identifier for session-specific operations.

    Returns:
    - fig: Plotly figure object for the scatter plot.
    """

    if x_gene.shape[0] > 500:
        sct_mark_size = 10
    elif 500 > x_gene.shape[0] > 150:
        sct_mark_size = 12.5
    else:
        sct_mark_size = 15

    # Set plot titles and labels based on comparison type
    if comp_type == "ssgsea":
        plot_title = " Cluster"
        plot_text = "ssGSEA Score"
        plot_text_2 = " Cluster"
        plot_title_gene = ""
        plot_text_gene = ""
        plot_text_2_gene = ""
    elif comp_type == "gc_ssgsea":
        plot_title = " Cluster"
        plot_text = "ssGSEA Score"
        plot_text_2 = " Cluster"
        plot_title_gene = ""
        plot_text_gene = "Expression"
        plot_text_2_gene = ""
    else:
        plot_title = ""
        plot_text = "Expression"
        plot_text_2 = ""
        plot_title_gene = ""
        plot_text_gene = ""
        plot_text_2_gene = ""

    # Prepare DataFrame for plotting
    df_plot = pd.DataFrame({x_gene_name: x_gene, y_gene_name: y_gene}).reset_index()

    for column in df_plot.columns:
        if column in merged_df.columns:
            if column == x_gene_name and column != y_gene_name:
                x_gene_name = x_gene_name + "_genes"
                df_plot = df_plot.rename(columns={column: column + "_genes"})
            if column != x_gene_name and column == y_gene_name:
                y_gene_name = y_gene_name + "_genes"
                df_plot = df_plot.rename(columns={column: column + "_genes"})
            else:
                df_plot = df_plot.rename(columns={column: column + "_genes"})
                y_gene_name = y_gene_name + "_genes"
                x_gene_name = x_gene_name + "_genes"

    if color_palette == "Sample_Cluster":
        if 'S0' in merged_df["Sample_Cluster"].unique():
            # Shift Set1 palette (up to 9 colors) for 'S0'
            if (len(merged_df["Sample_Cluster"].unique()) - 1) <= 9:
                # Shift Set1 colors
                color_map = {
                    'S0': '#E41A1C', 'S1': '#377EB8', 'S2': '#4DAF4A', 'S3': '#984EA3',
                    'S4': '#FF7F00', 'S5': '#FFFF33', 'S6': '#A65628', 'S7': '#F781BF', 'S8': '#999999'
                }
            else:
                # Shift Light24 palette for 'S0' (more than 9 colors)
                color_map = {
                    'S0': '#1F77B4', 'S1': '#AEC7E8', 'S2': '#FF7F0E', 'S3': '#FFBB78',
                    'S4': '#2CA02C', 'S5': '#98DF8A', 'S6': '#D62728', 'S7': '#FF9896',
                    'S8': '#9467BD', 'S9': '#C5B0D5', 'S10': '#8C564B', 'S11': '#C49C94',
                    'S12': '#E377C2', 'S13': '#F7B6D2', 'S14': '#7F7F7F', 'S15': '#C7C7C7'
                }
        else:
            if (len(merged_df["Sample_Cluster"].unique()) - 1) <= 9:
                # Set1 palette (up to 9 colors)
                color_map = {
                    'S1': '#E41A1C', 'S2': '#377EB8', 'S3': '#4DAF4A', 'S4': '#984EA3',
                    'S5': '#FF7F00', 'S6': '#FFFF33', 'S7': '#A65628', 'S8': '#F781BF', 'S9': '#999999'
                }
            else:
                # Light24 palette (more than 9 colors)
                color_map = {
                    'S1': '#1F77B4', 'S2': '#AEC7E8', 'S3': '#FF7F0E', 'S4': '#FFBB78',
                    'S5': '#2CA02C', 'S6': '#98DF8A', 'S7': '#D62728', 'S8': '#FF9896',
                    'S9': '#9467BD', 'S10': '#C5B0D5', 'S11': '#8C564B', 'S12': '#C49C94',
                    'S13': '#E377C2', 'S14': '#F7B6D2', 'S15': '#7F7F7F', 'S16': '#C7C7C7'
                }


    else:
    # Define color map
        color_map = {
            "ADC": "#1f77b4", "LCNEC": "#9beaf2", "CARCI": "#ff7f0e", "SQC": "#ffb5e9", "BAS": "#2ca02c", "NTL": "#000000",
            "LCC": "#7a37b8", "Other": "#87888a", "SCLC": "#d62728",
            'ASCL1': '#FD8A8A', 'YAP1': '#EFFA7A', 'NEUROD1': '#8FEBC8', 'POU2F3': '#8E60DB', 'NE': '#FCCA53',
            'non-NE': '#84CBD9', 'Female': '#C79558', 'Male': '#7DF8FF', 'F': '#C79558', 'M': '#7DF8FF',
            'nan': '#a9a9ab', 'Epithelial': '#9DE481', 'Epithelial-Mesenchymal': '#D3B7D9', 'Mesenchymal': '#674E62',
            'Current': '#F2A359', 'Former': '#64A5DE', 'Never': '#98D19C', 'No': '#98D19C', 'Yes': '#F2A359',
            'Passive': '#ae95b8', 'no': '#98D19C', 'yes': '#F2A359', 'unknown': '#a9a9ab', 'Extensive stage': '#f74e05',
            'Limited stage': '#fcacf0', 'I': '#7AB1E7', 'IA': '#8FBCEC', 'IB': '#A4C7F1', 'II': '#B9D2F6',
            'IIIA/IIIB': "#c0a7c2", 'Ia': '#8FBCEC', 'Ib': '#A4C7F1', 'IIa': '#CEDDFB', 'IIb': '#E3E8FF',
            'III': "#D1C2D2", 'IIIa': '#c0a7c2', 'IIIb': '#BF9CAB', 'IIIB': '#BF9CAB', 'IV': '#AD769E', 'IA1': '#A4C7F1',
            'IA3': '#CEDDFB', 'IIA': '#D1C2D2', 'IIB': '#BF9CAB', 'IIIA': '#AD769E', 'MALE': '#7DF8FF',
            'FEMALE': '#C79558', '[Not Available]': '#a9a9ab', '[Not Evaluated]': '#a9a9ab', '[Unknown]': '#a9a9ab',
            'AMER. IND. / ALASKA NATIVE': '#B3FFD9', 'ASIAN': '#CCE6FF', 'BLACK / AFR. AMER.': '#D9B3FF',
            'NAT. HAWAIIAN / PACIFIC ISLANDER': '#FFCCE6', 'WHITE': '#FFD9B3', 'Asian': '#CCE6FF', 'IA2': '#B9D2F6',
            'Black or African American': '#D9B3FF', 'White': '#FFD9B3', 'Naive': '#98D19C', 'Chemo': '#F2A359',
            'Ia2': '#8FBCEC', 'Ia3': '#A4C7F1', "Typical": "#FF9900", "Atypical": "#0099FF",
            'Supra_carcinoid': '#ff8cd9', 'Supra_Carcinoid': '#ff8cd9', 'NA_Smoker': '#a9a9ab', 'NA_UICC': '#a9a9ab',
            'NA_EC': '#a9a9ab', 'NA_Sex': '#a9a9ab', 'NA_Race': '#a9a9ab', 'NA_Smoking': '#a9a9ab',
            'NA_Smoke': '#a9a9ab', 'NA_Mol_ST': '#a9a9ab', 'Caucasian:': '#FFD9B3', 'Type 1 LCNEC': '#00A20B',
            'Type 2 LCNEC': '#AF5600', 'Supra_carcinoid_HS': '#ff8cd9', 'Supra_carcinoid_MC': '#ff8cd9',
            'Supra_carcinoid_mol_clust': '#CC00CC', 'LC1': '#FF0009', 'LC2': '#45029C', 'LC3': '#046220',
            'Carcinoid-A1': '#0099FF', 'Carcinoid-A2': '#FF9900', 'Carcinoid-B': '#66CC00', 'Carcinoid': '#CC00CC',
            'NTL_TNM': '#a9a9ab', 'NTL_M': '#a9a9ab', 'NTL_N': '#a9a9ab', 'NTL_T_Stage': '#a9a9ab',
            'NTL_Gender': '#a9a9ab', 'non_CARCI': '#a9a9ab'
        }

    color_pal = ['#0099FF', '#FF9900', '#66CC00', '#FF66CC', '#00CC66',
                 '#FF3300', '#CC00CC', '#33CCCC', '#FFCC00', '#99FF33',
                 '#FF0066', '#33CC00', '#FF99CC', '#00CCFF', '#CC3300',
                 '#66CCFF', "#CDEEAA", "#EFB960", "#C3AED6", "#874C62",
                 "#EA907A", "#2192FF", "#557571", "#F9FD50", "#B2EC5D"]

    if color_palette != 'Monochrome':
        merged_df.fillna('nan', inplace=True)
        df_plot = df_plot.merge(merged_df[[color_palette]], left_on='index', right_index=True, how='left')

        # Sort categories alphabetically or by another criterion here
        df_plot = df_plot.sort_values(by=color_palette)

        all_categories = df_plot[color_palette].dropna().unique()
        color_cycle = cycle(color_pal)  # Cycle through your dynamic palette

        full_color_map = {cat: color_map.get(cat, next(color_cycle)) for cat in all_categories}

        df_plot['color'] = df_plot[color_palette].map(full_color_map)

        if comp_type != 'gc_ssgsea':
            fig = px.scatter(df_plot, x=x_gene_name, y=y_gene_name, color=color_palette,
                             color_discrete_map=full_color_map,
                             title=f'Correlation between {x_gene_name} and {y_gene_name}' + plot_title,
                             labels={x_gene_name: plot_text + f' of {plot_text_2} {x_gene_name}',
                                     y_gene_name: plot_text + f' of {plot_text_2} {y_gene_name}'})
        else:
            fig = px.scatter(df_plot, x=x_gene_name, y=y_gene_name, color=color_palette,
                             color_discrete_map=full_color_map,
                             title=f'Correlation between {x_gene_name} and {y_gene_name}' + plot_title,
                             labels={x_gene_name: plot_text_gene + f' of {plot_text_2_gene} {x_gene_name}',
                                     y_gene_name: plot_text + f' of {plot_text_2} {y_gene_name}'})
    else:
        if comp_type != 'gc_ssgsea':
            fig = px.scatter(df_plot, x=x_gene_name, y=y_gene_name,
                             title=f'Correlation between {x_gene_name} and {y_gene_name}' + plot_title,
                             labels={x_gene_name: plot_text + f' of {plot_text_2} {x_gene_name}',
                                     y_gene_name: plot_text + f' of {plot_text_2} {y_gene_name}'})
        else:
            fig = px.scatter(df_plot, x=x_gene_name, y=y_gene_name,
                             title=f'Correlation between {x_gene_name} and {y_gene_name}' + plot_title,
                             labels={x_gene_name: plot_text_gene + f' of {plot_text_2_gene} {x_gene_name}',
                                     y_gene_name: plot_text + f' of {plot_text_2} {y_gene_name}'})
        fig.update_traces(marker=dict(color=mono_color))

    fig.update_layout(
        font=dict(family="Arial", size=16),
        legend_title=dict(font=dict(size=16, color='black', family="Arial")),
        legend=dict(font=dict(size=14, color='black', family="Arial"))
    )
    fig.update_traces(marker_size=sct_mark_size, marker_line_width=1, marker_line_color="black", marker_opacity=0.75)

    fig.update_xaxes(
        tickfont=dict(color='black'),
        title_font=dict(color='black')
    )
    fig.update_yaxes(
        tickfont=dict(color='black'),
        title_font=dict(color='black')
    )

    # Annotate correlation coefficient and p-value
    fig.add_annotation(
        text=f"Corr. Coef.: {corr_coef:.2f}  P-value: {p_value:.2e}",
        xref="paper", yref="paper",
        x=0.05, y=1.1, showarrow=False,
        font=dict(size=14, color='black', family="Arial", )
    )

    m, b = np.polyfit(df_plot[x_gene_name], df_plot[y_gene_name], 1)
    x = np.linspace(df_plot[x_gene_name].min(), df_plot[x_gene_name].max(), 100)
    fig.add_trace(go.Scatter(x=x, y=m * x + b, mode='lines', name='Fitted Line', line=dict(color='red')))

    # Save the plot as PDF
    user_dir = "/result_data_" + user_id
    if comp_type == "ssgsea":
        file_path = "result_data" + user_dir + "/ssgsea_corr_" + user_id + ".pdf"
    elif comp_type == "gc_ssgsea":
        file_path = "result_data" + user_dir + "/gene_ssgsea_corr_" + user_id + ".pdf"
    else:
        file_path = "result_data" + user_dir + "/gene_compare_corr_" + user_id + ".pdf"
    fig.write_image(file_path, width=600, height=600, scale=1)

    return fig


def plot_component_counts(component_counts_df, chi_bp_name, chi_plot_col):
    """
    Create and display a stacked bar chart of component counts using Plotly.

    Args:
        component_counts_df (pd.DataFrame): DataFrame with component counts.
        chi_bp_name (str): Path to save the bar plot image.
        chi_plot_col (list): List of Streamlit columns for plotting.

    Returns:
        None
    """
    # Reset index to make 'Cluster' a column in the DataFrame
    color_palettes = component_counts_df.columns
    component_counts_df = component_counts_df.reset_index()
    component_counts_df.rename(columns={'index': 'Cluster'}, inplace=True)

    # Predefined color map
    color_map = {
        "ADC": "#1f77b4", "LCNEC": "#9beaf2", "CARCI": "#ff7f0e", "SQC": "#ffb5e9", "BAS": "#2ca02c", "NTL": "#000000",
        "LCC": "#7a37b8", "Other": "#87888a", "SCLC": "#d62728", 'NA': '#e6e6e8', 'NaN': '#e6e6e8', 'Na': '#e6e6e8',
        'ASCL1': '#FD8A8A', 'YAP1': '#EFFA7A', 'NEUROD1': '#8FEBC8', 'POU2F3': '#8E60DB', 'NE': '#FCCA53',
        'non-NE': '#84CBD9', 'Female': '#c79558', 'Male': '#7df8ff', 'F': '#c79558', 'M': '#7df8ff',
        'nan': '#e6e6e8', 'Epithelial': '#9DE481', 'Epithelial-Mesenchymal': '#D3B7D9', 'Mesenchymal': '#674E62',
        'Current': '#F2A359', 'Former': '#64A5DE', 'Never': '#98D19C', 'No': '#98D19C', 'Yes': '#F2A359',
        'Passive': '#ae95b8', 'no': '#98D19C', 'yes': '#F2A359', 'unknown': '#a9a9ab', 'Extensive stage': '#f74e05',
        'Limited stage': '#fcacf0', 'I': '#7AB1E7', 'IA': '#8FBCEC', 'IB': '#A4C7F1', 'II': '#B9D2F6',
        'IIIA/IIIB': "#c0a7c2", 'Ia': '#8FBCEC', 'Ib': '#A4C7F1', 'IIa': '#CEDDFB', 'IIb': '#E3E8FF',
        'III': "#D1C2D2", 'IIIa': '#c0a7c2', 'IIIb': '#BF9CAB', 'IIIB': '#BF9CAB', 'IV': '#AD769E', 'IA1': '#A4C7F1',
        'IA3': '#CEDDFB', 'IIA': '#D1C2D2', 'IIB': '#BF9CAB', 'IIIA': '#AD769E', 'MALE': '#7df8ff',
        'FEMALE': '#c79558', '[Not Available]': '#a9a9ab', '[Not Evaluated]': '#a9a9ab', '[Unknown]': '#a9a9ab',
        'AMER. IND. / ALASKA NATIVE': '#B3FFD9', 'ASIAN': '#CCE6FF', 'BLACK / AFR. AMER.': '#D9B3FF',
        'NAT. HAWAIIAN / PACIFIC ISLANDER': '#FFCCE6', 'WHITE': '#FFD9B3', 'Asian': '#CCE6FF', 'IA2': '#B9D2F6',
        'Black or African American': '#D9B3FF', 'White': '#FFD9B3', 'Naive': '#98D19C', 'Chemo': '#F2A359',
        'Ia2': '#8FBCEC', 'Ia3': '#A4C7F1', "Typical": "#FF9900", "Atypical": "#0099FF", 'Supra_carcinoid': '#ff8cd9',
        'Supra_Carcinoid': '#ff8cd9', 'NA_Smoker': '#a9a9ab', 'NA_UICC': '#a9a9ab', 'NA_EC': '#a9a9ab',
        'NA_Sex': '#a9a9ab', 'NA_Race': '#a9a9ab', 'NA_Smoking': '#a9a9ab', 'NA_Smoke': '#a9a9ab',
        'NA_Mol_ST': '#a9a9ab', 'Caucasian:': '#FFD9B3', 'Type 1 LCNEC': '#00A20B', 'Type 2 LCNEC': '#AF5600',
        'Supra_carcinoid_HS': '#ff8cd9', 'Supra_carcinoid_MC': '#ff8cd9', 'Supra_carcinoid_mol_clust': '#CC00CC',
        'LC1': '#FF0009', 'LC2': '#45029C', 'LC3': '#046220', 'Carcinoid-A1': '#0099FF', 'Carcinoid-A2': '#FF9900',
        'Carcinoid-B': '#66CC00', 'Carcinoid': '#CC00CC', 'NA_TNM': '#a9a9ab', 'NTL_TNM': '#a9a9ab',
        'NTL_M': '#a9a9ab', 'NTL_N': '#a9a9ab', 'NTL_T_Stage': '#a9a9ab', 'NTL_Gender': '#a9a9ab', 'non_CARCI': '#a9a9ab'
    }

    # Define fallback colors
    color_pal = ['#0099FF', '#FF9900', '#66CC00', '#FF66CC', '#00CC66',
                 '#FF3300', '#CC00CC', '#33CCCC', '#FFCC00', '#99FF33',
                 '#FF0066', '#33CC00', '#FF99CC', '#00CCFF', '#CC3300',
                 '#66CCFF', "#CDEEAA", "#EFB960", "#C3AED6", "#874C62",
                 "#EA907A", "#2192FF", "#557571", "#F9FD50", "#B2EC5D"]

    color_cycle = cycle(color_pal)

    # Create a mapping of categories to their colors
    all_categories = color_palettes
    full_color_map = {cat: color_map.get(cat, next(color_cycle)) for cat in all_categories}

    # Map the colors to a new column in the DataFrame
    color_list = pd.Series(color_palettes.map(full_color_map))

    # Use the color column for the color in bar plots
    color_palette_bar = color_list.unique().tolist()

    # Melt the DataFrame to make it suitable for a stacked bar chart
    component_counts_melted = component_counts_df.melt(id_vars='Cluster',
                                                       value_vars=component_counts_df.columns[1:],
                                                       var_name='Component',
                                                       value_name='Count')

    # Create a stacked bar chart
    fig = px.bar(component_counts_melted,
                 x='Cluster',
                 y='Count',
                 color='Component',
                 title='Cluster Component Counts',
                 labels={'Cluster': 'Cluster'},
                 color_discrete_sequence=color_palette_bar,
                 height=400)

    # Update axis properties
    fig.update_xaxes(
        tickfont=dict(color='black', family="Arial", size=14),  # Set tick font color, family, and size
        title_font=dict(color='black', family="Arial", size=16)  # Set title font color, family, and size
    )
    fig.update_yaxes(
        tickfont=dict(color='black', family="Arial", size=14),  # Set tick font color, family, and size
        title_font=dict(color='black', family="Arial", size=16)  # Set title font color, family, and size
    )

    # Update layout for legend
    fig.update_layout(
        font=dict(family="Arial"),
        legend=dict(
            title=dict(font=dict(size=16, color='black', family="Arial")),  # Set legend title font
            font=dict(size=14, color='black', family="Arial")  # Set legend font
        )
    )

    # Save the figure as an image
    fig.write_image(chi_bp_name, width=600, height=600, scale=1)

    try:
        # Show the plot
        chi_plot_col[1].plotly_chart(fig, use_container_width=True)
    except Exception as ex:
        chi_plot_col[1].image("style_items/error.svg", use_container_width=True)
        chi_plot_col[1].warning(
            ":blue[There is a problem with the ] :red[chart! ] :blue[Please try another ] :red[setting ] "
            ":blue[or ] :red[refresh ] :blue[the page]")
        print(f"Error: {ex}")


def plot_custom_strip(df_melt, strip_path, anno_file, st_anno, ssgsea_col, key, sample_cluster):
    """
    Create and save a custom strip plot using ssGSEA results.

    Args:
        df_melt (pd.DataFrame): Melted DataFrame containing ssGSEA results.
        strip_path (str): Path to save the strip plot image.
        anno_file (pd.DataFrame): Annotation file DataFrame.
        st_anno (pd.DataFrame): Streamlit annotation DataFrame.
        ssgsea_col (st.delta_generator.DeltaGenerator): Streamlit column for ssGSEA.
        key (str): Unique key for the Streamlit components.

    Returns:
        go.Figure: Generated Plotly figure.
    """
    # Determine scatter marker size based on the number of samples
    if df_melt.shape[0] > 500:
        sct_mark_size = 10
    elif 500 > df_melt.shape[0] > 150:
        sct_mark_size = 12.5
    else:
        sct_mark_size = 15

    # Drop columns with only one unique value
    cols_to_drop_anno = anno_file.columns[anno_file.nunique() < 2]
    cols_to_drop_st = st_anno.columns[st_anno.nunique() < 2]

    anno_file = anno_file.drop(columns=cols_to_drop_anno)
    st_anno = st_anno.drop(columns=cols_to_drop_st)

    overlap_columns = anno_file.columns.intersection(st_anno.columns)

    # Rename overlapping columns in clini_df before merging
    st_anno_renamed = st_anno.rename(columns={col: f"{col}_from_clinic_df" for col in overlap_columns})

    # Merge anno_file and clini_df_renamed
    merged_df = pd.merge(anno_file, st_anno_renamed, left_index=True, right_index=True, how='outer')

    # Handle overlapping columns
    for col in overlap_columns:
        merged_df.drop(col, axis=1, inplace=True)  # Drop anno_file's original overlapping columns
        merged_df.rename(columns={f"{col}_from_clinic_df": col}, inplace=True)  # Rename clini_df's columns back

    merged_df = merged_df[merged_df.index.isin(df_melt["sample"].unique())]

    # Reset the index of df_single to turn the index into a column
    merged_df.index.rename('sample', inplace=True)
    df_single_reset = merged_df.reset_index()

    # Merge df_multi with df_single_reset on the 'sample' column, keeping all rows from df_multi
    merged_df = pd.merge(df_melt, df_single_reset, on='sample', how='left')
    merged_df.set_index("sample", inplace=True)
    merged_df['Sample_Cluster'] = sample_cluster['Cluster'].reindex(merged_df.index)
    merged_df.fillna('NA', inplace=True)
    anno_groups = list(merged_df.columns)
    anno_groups.remove('ssGSEA')

    # Select coloring palette
    ssgsea_color_pal = ssgsea_col.selectbox(":blue[Select the ] :red[Coloring Palette]",
                                            options=anno_groups, key="ssgsea_color_pal" + key, index=0)

    if ssgsea_color_pal == 'Cluster':
        jitter_color = 'Cluster'
    else:
        jitter_color = merged_df[ssgsea_color_pal]

    clusters_unique = df_melt['Cluster'].unique()
    num_clusters = len(clusters_unique)
    if num_clusters <= 9:
        color_palette = px.colors.qualitative.Set1
    else:
        color_palette = px.colors.qualitative.Light24

    if ssgsea_color_pal == "Cluster":
        color_palette_ssgsea = color_palette
    else:
        if ssgsea_color_pal == "Sample_Cluster":
            if 'S0' in merged_df["Sample_Cluster"].unique():
                # Shift Set1 palette (up to 9 colors) for 'S0'
                if (len(merged_df["Sample_Cluster"].unique()) - 1) <= 9:
                    # Shift Set1 colors
                    color_map = {
                        'S0': '#E41A1C', 'S1': '#377EB8', 'S2': '#4DAF4A', 'S3': '#984EA3',
                        'S4': '#FF7F00', 'S5': '#FFFF33', 'S6': '#A65628', 'S7': '#F781BF', 'S8': '#999999'
                    }
                else:
                    # Shift Light24 palette for 'S0' (more than 9 colors)
                    color_map = {
                        'S0': '#1F77B4', 'S1': '#AEC7E8', 'S2': '#FF7F0E', 'S3': '#FFBB78',
                        'S4': '#2CA02C', 'S5': '#98DF8A', 'S6': '#D62728', 'S7': '#FF9896',
                        'S8': '#9467BD', 'S9': '#C5B0D5', 'S10': '#8C564B', 'S11': '#C49C94',
                        'S12': '#E377C2', 'S13': '#F7B6D2', 'S14': '#7F7F7F', 'S15': '#C7C7C7'
                    }
            else:
                if (len(merged_df["Sample_Cluster"].unique()) - 1) <= 9:
                    # Set1 palette (up to 9 colors)
                    color_map = {
                        'S1': '#E41A1C', 'S2': '#377EB8', 'S3': '#4DAF4A', 'S4': '#984EA3',
                        'S5': '#FF7F00', 'S6': '#FFFF33', 'S7': '#A65628', 'S8': '#F781BF', 'S9': '#999999'
                    }
                else:
                    # Light24 palette (more than 9 colors)
                    color_map = {
                        'S1': '#1F77B4', 'S2': '#AEC7E8', 'S3': '#FF7F0E', 'S4': '#FFBB78',
                        'S5': '#2CA02C', 'S6': '#98DF8A', 'S7': '#D62728', 'S8': '#FF9896',
                        'S9': '#9467BD', 'S10': '#C5B0D5', 'S11': '#8C564B', 'S12': '#C49C94',
                        'S13': '#E377C2', 'S14': '#F7B6D2', 'S15': '#7F7F7F', 'S16': '#C7C7C7'
                    }

        else:
            # Predefined color map
            color_map = {
                "ADC": "#1f77b4", "LCNEC": "#9beaf2", "CARCI": "#ff7f0e", "SQC": "#ffb5e9", "BAS": "#2ca02c",
                "NTL": "#000000",
                "LCC": "#7a37b8", "Other": "#87888a", "SCLC": "#d62728", 'NA': '#e6e6e8', 'NaN': '#e6e6e8',
                'Na': '#e6e6e8',
                'ASCL1': '#FD8A8A', 'YAP1': '#EFFA7A', 'NEUROD1': '#8FEBC8', 'POU2F3': '#8E60DB', 'NE': '#FCCA53',
                'non-NE': '#84CBD9', 'Female': '#c79558', 'Male': '#7df8ff', 'F': '#c79558', 'M': '#7df8ff',
                'nan': '#e6e6e8', 'Epithelial': '#9DE481', 'Epithelial-Mesenchymal': '#D3B7D9',
                'Mesenchymal': '#674E62',
                'Current': '#F2A359', 'Former': '#64A5DE', 'Never': '#98D19C', 'No': '#98D19C', 'Yes': '#F2A359',
                'Passive': '#ae95b8', 'no': '#98D19C', 'yes': '#F2A359', 'unknown': '#a9a9ab',
                'Extensive stage': '#f74e05',
                'Limited stage': '#fcacf0', 'I': '#7AB1E7', 'IA': '#8FBCEC', 'IB': '#A4C7F1', 'II': '#B9D2F6',
                'IIIA/IIIB': "#c0a7c2", 'Ia': '#8FBCEC', 'Ib': '#A4C7F1', 'IIa': '#CEDDFB', 'IIb': '#E3E8FF',
                'III': "#D1C2D2", 'IIIa': '#c0a7c2', 'IIIb': '#BF9CAB', 'IIIB': '#BF9CAB', 'IV': '#AD769E',
                'IA1': '#A4C7F1',
                'IA3': '#CEDDFB', 'IIA': '#D1C2D2', 'IIB': '#BF9CAB', 'IIIA': '#AD769E', 'MALE': '#7df8ff',
                'FEMALE': '#c79558', '[Not Available]': '#a9a9ab', '[Not Evaluated]': '#a9a9ab', '[Unknown]': '#a9a9ab',
                'AMER. IND. / ALASKA NATIVE': '#B3FFD9', 'ASIAN': '#CCE6FF', 'BLACK / AFR. AMER.': '#D9B3FF',
                'NAT. HAWAIIAN / PACIFIC ISLANDER': '#FFCCE6', 'WHITE': '#FFD9B3', 'Asian': '#CCE6FF', 'IA2': '#B9D2F6',
                'Black or African American': '#D9B3FF', 'White': '#FFD9B3', 'Naive': '#98D19C', 'Chemo': '#F2A359',
                'Ia2': '#8FBCEC', 'Ia3': '#A4C7F1', "Typical": "#FF9900", "Atypical": "#0099FF",
                'Supra_carcinoid': '#ff8cd9', 'Supra_Carcinoid': '#ff8cd9', 'NA_Smoker': '#a9a9ab',
                'NA_UICC': '#a9a9ab', 'NA_EC': '#a9a9ab', 'NA_Sex': '#a9a9ab', 'NA_Race': '#a9a9ab',
                'NA_Smoking': '#a9a9ab', 'NA_Smoke': '#a9a9ab', 'NA_Mol_ST': '#a9a9ab', 'Caucasian:': '#FFD9B3',
                'Type 1 LCNEC': '#00A20B', 'Type 2 LCNEC': '#AF5600', 'Supra_carcinoid_HS': '#ff8cd9',
                'Supra_carcinoid_MC': '#ff8cd9', 'Supra_carcinoid_mol_clust': '#CC00CC', 'LC1': '#FF0009',
                'LC2': '#45029C', 'LC3': '#046220', 'Carcinoid-A1': '#0099FF', 'Carcinoid-A2': '#FF9900',
                'Carcinoid-B': '#66CC00', 'Carcinoid': '#CC00CC', 'NA_TNM': '#a9a9ab', 'NTL_TNM': '#a9a9ab',
                'NTL_M': '#a9a9ab', 'NTL_N': '#a9a9ab', 'NTL_T_Stage': '#a9a9ab', 'NTL_Gender': '#a9a9ab',
                'non_CARCI': '#a9a9ab'
            }

        # Define fallback colors
        color_pal = ['#0099FF', '#FF9900', '#66CC00', '#FF66CC', '#00CC66',
                     '#FF3300', '#CC00CC', '#33CCCC', '#FFCC00', '#99FF33',
                     '#FF0066', '#33CC00', '#FF99CC', '#00CCFF', '#CC3300',
                     '#66CCFF', "#CDEEAA", "#EFB960", "#C3AED6", "#874C62",
                     "#EA907A", "#2192FF", "#557571", "#F9FD50", "#B2EC5D"]

        # Prepare a cycle of fallback colors from color_pal
        color_cycle = cycle(color_pal)

        # Identify all unique categories within the selected column
        all_categories = merged_df[ssgsea_color_pal].dropna().unique()

        # Create a mapping of categories to their colors
        full_color_map = {cat: color_map.get(cat, next(color_cycle)) for cat in all_categories}

        # Map the colors to a new column in the DataFrame
        color_list = pd.Series(merged_df[ssgsea_color_pal].map(full_color_map))

        # Use the color column for the color in strip plots
        color_palette_ssgsea = color_list.unique().tolist()

    # Define the range for the violin plot based on the min and max values of the data
    min_val = df_melt['ssGSEA'].min()
    max_val = df_melt['ssGSEA'].max()

    if len(df_melt["Cluster"].unique()) > 1:

        # Let users select the order of violin plots
        violin_order = st.radio(
            ":blue[Select the ] :red[Order ] :blue[of Violin Plots]",
            options=["Median", "Standard Deviation", "Custom"],
            key="violin_order" + key,
            horizontal=True
        )

        if violin_order in ["Median", "Standard Deviation"]:
            # Let users select the sort order (ascending or descending)
            sort_order = st.radio(
                ":blue[Select the ] :red[Sort Order]",
                options=["Ascending", "Descending"],
                key="sort_order" + key,
                horizontal=True
            )
            # Determine if ascending or descending
            ascending = True if sort_order == "Ascending" else False
        else:
            ascending = False  # Default value when sort order is not applicable

        if violin_order == "Median":
            cluster_order_method = df_melt.groupby('Cluster')['ssGSEA'].median()
            sorted_clusters = cluster_order_method.sort_values(ascending=ascending).index.tolist()
        elif violin_order == "Standard Deviation":
            cluster_order_method = df_melt.groupby('Cluster')['ssGSEA'].std()
            sorted_clusters = cluster_order_method.sort_values(ascending=ascending).index.tolist()
        else:
            clusters_unique = df_melt['Cluster'].unique().tolist()
            st.write("### :blue[Select Clusters in Desired Order:]")
            sorted_clusters = st.multiselect(
                ":blue[Select Clusters in the desired order ] "
                ":red[Select all Clusters from 'Select Cluster' Selection'])]",
                options=clusters_unique,
                default=[],
                key="custom_order_multiselect" + key
            )
            # Check if all clusters are selected
            missing_clusters = set(clusters_unique) - set(sorted_clusters)
            if missing_clusters:
                st.warning(f":red[You have not selected all clusters. ] :blue[Missing clusters: ] "
                           f":red[{', '.join(missing_clusters)}]")
                st.stop()  # Optionally stop execution until all clusters are selected

        # Set 'Cluster' as categorical with the desired order
        df_melt['Cluster'] = pd.Categorical(
            df_melt['Cluster'],
            categories=sorted_clusters,
            ordered=True
        )

    else:
        sorted_clusters = df_melt['Cluster'].unique().tolist()

    # Generate violin plots with sorted clusters
    fig = go.Figure()
    for idx, cluster in enumerate(sorted_clusters):
        df_subset = df_melt[df_melt['Cluster'] == cluster]
        fig.add_trace(go.Violin(
            y=df_subset['ssGSEA'],
            x=df_subset['Cluster'],
            fillcolor=color_palette[idx % len(color_palette)],
            opacity=0.5,
            line_color='black',
            line_width=3,
            box_visible=False,
            points=False,
            meanline_visible=False,
            showlegend=False,
            legendgroup=cluster,
            width=0,
            spanmode='hard',
            span=[min_val * 1.5, max_val * 1.5]
        ))

    # Generate strip plots with category order
    temp_fig = px.strip(
        df_melt,
        x='Cluster',
        y='ssGSEA',
        color=jitter_color,
        title="ssGSEA Results",
        stripmode='overlay',
        color_discrete_sequence=color_palette_ssgsea,
        category_orders={'Cluster': sorted_clusters}
    )
    for trace in temp_fig.data:
        fig.add_trace(trace)

    # Generate box plots with sorted clusters
    for idx, cluster in enumerate(sorted_clusters):
        df_subset = df_melt[df_melt['Cluster'] == cluster]
        fig.add_trace(go.Box(
            y=df_subset['ssGSEA'],
            x=df_subset['Cluster'],
            fillcolor="white",
            line=dict(color='black', width=3),
            opacity=0.5,
            showlegend=False,
            boxpoints=False,
            legendgroup=cluster,
            width=0.05
        ))

    # Update x-axis category order
    fig.update_xaxes(categoryorder='array', categoryarray=sorted_clusters)
    fig.update_layout(title="GSEA Results")
    fig.update_traces(
        marker_size=sct_mark_size,
        marker_line_width=1,
        marker_line_color="black",
        marker_opacity=0.75
    )
    fig.update_yaxes(title_text="ssGSEA")
    fig.update_xaxes(title_text="Cluster", type='category')

    # Update axis properties
    fig.update_xaxes(
        tickfont=dict(color='black'),
        title_font=dict(color='black')
    )
    fig.update_yaxes(
        tickfont=dict(color='black'),
        title_font=dict(color='black')
    )

    # Save the plot as an image
    fig.write_image(strip_path, width=800, height=800)

    return fig


def plot_custom_fusion(gene, gene_name, color_palette, mono_color, merged_df, comp_type, user_id, axis_id, key,
                       multi_col):
    """
    Creates a violin plot to visualize the distribution of gene expression levels across different
    categories, optionally coloring the dots based on a selected palette.

    Parameters:
    - gene: Expression data for the gene.
    - gene_name: Name of the gene.
    - color_palette: Selected coloring palette or 'Monochrome'.
    - mono_color: Selected color for monochrome plots.
    - merged_df: DataFrame with sample annotations, including the 'NAPY' column.
    - comp_type: Comparison type (e.g., "ssgsea").
    - user_id: User identifier for session-specific operations.
    - axis_id: Identifier for the axis (e.g., "x_gene" or "y_gene").

    Returns:
    - fig: Plotly figure object for the violin plot.
    - p_value_df: DataFrame containing p-values for category comparisons.
    """
    # ------------------------- Preprocessing -------------------------

    # Create df_plot DataFrame with gene expression
    df_plot = pd.DataFrame({gene_name: gene})

    # Rename overlapping columns in df_plot if they exist in merged_df
    for column in df_plot.columns:
        if column in merged_df.columns:
            df_plot = df_plot.rename(columns={column: f"{column}_genes"})
            gene_name = f"{gene_name}_genes"

    # Merge df_plot with merged_df
    df_plot = df_plot.merge(merged_df, left_index=True, right_index=True)

    # Define plot title and text based on comp_type
    if comp_type == "ssgsea":
        plot_title = "ssGSEA Score"
        plot_text = "Cluster"
    else:
        plot_title = "Expression"
        plot_text = ""

    # ------------------------- Color Mapping -------------------------

    if color_palette == "Sample_Cluster":
        if 'S0' in merged_df["Sample_Cluster"].unique():
            # Shift Set1 palette (up to 9 colors) for 'S0'
            if (len(merged_df["Sample_Cluster"].unique()) - 1) <= 9:
                color_map = {
                    'S0': '#E41A1C', 'S1': '#377EB8', 'S2': '#4DAF4A', 'S3': '#984EA3',
                    'S4': '#FF7F00', 'S5': '#FFFF33', 'S6': '#A65628', 'S7': '#F781BF', 'S8': '#999999'
                }
            else:
                # Shift Light24 palette for 'S0' (more than 9 colors)
                color_map = {
                    'S0': '#1F77B4', 'S1': '#AEC7E8', 'S2': '#FF7F0E', 'S3': '#FFBB78',
                    'S4': '#2CA02C', 'S5': '#98DF8A', 'S6': '#D62728', 'S7': '#FF9896',
                    'S8': '#9467BD', 'S9': '#C5B0D5', 'S10': '#8C564B', 'S11': '#C49C94',
                    'S12': '#E377C2', 'S13': '#F7B6D2', 'S14': '#7F7F7F', 'S15': '#C7C7C7'
                }
        else:
            if (len(merged_df["Sample_Cluster"].unique()) - 1) <= 9:
                # Set1 palette (up to 9 colors)
                color_map = {
                    'S1': '#E41A1C', 'S2': '#377EB8', 'S3': '#4DAF4A', 'S4': '#984EA3',
                    'S5': '#FF7F00', 'S6': '#FFFF33', 'S7': '#A65628', 'S8': '#F781BF', 'S9': '#999999'
                }
            else:
                # Light24 palette (more than 9 colors)
                color_map = {
                    'S1': '#1F77B4', 'S2': '#AEC7E8', 'S3': '#FF7F0E', 'S4': '#FFBB78',
                    'S5': '#2CA02C', 'S6': '#98DF8A', 'S7': '#D62728', 'S8': '#FF9896',
                    'S9': '#9467BD', 'S10': '#C5B0D5', 'S11': '#8C564B', 'S12': '#C49C94',
                    'S13': '#E377C2', 'S14': '#F7B6D2', 'S15': '#7F7F7F', 'S16': '#C7C7C7'
                }
    else:
        # Predefined color map
        color_map = {
            "ADC": "#1f77b4", "LCNEC": "#9beaf2", "CARCI": "#ff7f0e", "SQC": "#ffb5e9", "BAS": "#2ca02c", "NTL": "#000000",
            "LCC": "#7a37b8", "Other": "#87888a", "SCLC": "#d62728",
            'ASCL1': '#FD8A8A', 'YAP1': '#EFFA7A', 'NEUROD1': '#8FEBC8', 'POU2F3': '#8E60DB', 'NE': '#FCCA53',
            'non-NE': '#84CBD9', 'Female': '#c79558', 'Male': '#7df8ff', 'F': '#c79558', 'M': '#7df8ff',
            'nan': '#a9a9ab', 'Epithelial': '#9DE481', 'Epithelial-Mesenchymal': '#D3B7D9', 'Mesenchymal': '#674E62',
            'Current': '#F2A359', 'Former': '#64A5DE', 'Never': '#98D19C', 'No': '#98D19C', 'Yes': '#F2A359',
            'Passive': '#ae95b8', 'no': '#98D19C', 'yes': '#F2A359', 'unknown': '#a9a9ab', 'Extensive stage': '#f74e05',
            'Limited stage': '#fcacf0', 'I': '#7AB1E7', 'IA': '#8FBCEC', 'IB': '#A4C7F1', 'II': '#B9D2F6',
            'IIIA/IIIB': "#c0a7c2", 'Ia': '#8FBCEC', 'Ib': '#A4C7F1', 'IIa': '#CEDDFB', 'IIb': '#E3E8FF',
            'III': "#D1C2D2", 'IIIa': '#c0a7c2', 'IIIb': '#BF9CAB', 'IIIB': '#BF9CAB', 'IV': '#AD769E', 'IA1': '#A4C7F1',
            'IA3': '#CEDDFB', 'IIA': '#D1C2D2', 'IIB': '#BF9CAB', 'IIIA': '#AD769E', 'MALE': '#7df8ff',
            'FEMALE': '#c79558', '[Not Available]': '#a9a9ab', '[Not Evaluated]': '#a9a9ab', '[Unknown]': '#a9a9ab',
            'AMER. IND. / ALASKA NATIVE': '#B3FFD9', 'ASIAN': '#CCE6FF', 'BLACK / AFR. AMER.': '#D9B3FF',
            'NAT. HAWAIIAN / PACIFIC ISLANDER': '#FFCCE6', 'WHITE': '#FFD9B3', 'Asian': '#CCE6FF', 'IA2': '#B9D2F6',
            'Black or African American': '#D9B3FF', 'White': '#FFD9B3', 'Naive': '#98D19C', 'Chemo': '#F2A359',
            'Ia2': '#8FBCEC', 'Ia3': '#A4C7F1', "Typical": "#FF9900", "Atypical": "#0099FF",
            'Supra_carcinoid': '#ff8cd9', 'Supra_Carcinoid': '#ff8cd9', 'NA_Smoker': '#a9a9ab', 'NA_UICC': '#a9a9ab',
            'NA_EC': '#a9a9ab', 'NA_Sex': '#a9a9ab', 'NA_Race': '#a9a9ab', 'NA_Smoking': '#a9a9ab',
            'NA_Smoke': '#a9a9ab', 'NA_Mol_ST': '#a9a9ab', 'Caucasian:': '#FFD9B3', 'Type 1 LCNEC': '#00A20B',
            'Type 2 LCNEC': '#AF5600', 'Supra_carcinoid_HS': '#ff8cd9', 'Supra_carcinoid_MC': '#ff8cd9',
            'Supra_carcinoid_mol_clust': '#CC00CC', 'LC1': '#FF0009', 'LC2': '#45029C', 'LC3': '#046220',
            'Carcinoid-A1': '#0099FF', 'Carcinoid-A2': '#FF9900', 'Carcinoid-B': '#66CC00', 'Carcinoid': '#CC00CC',
            'NTL_TNM': '#a9a9ab', 'NTL_M': '#a9a9ab', 'NTL_N': '#a9a9ab', 'NTL_T_Stage': '#a9a9ab',
            'NTL_Gender': '#a9a9ab', 'non_CARCI': '#a9a9ab'
        }

    # Define fallback colors
    color_pal = [
        '#0099FF', '#FF9900', '#66CC00', '#FF66CC', '#00CC66',
        '#FF3300', '#CC00CC', '#33CCCC', '#FFCC00', '#99FF33',
        '#FF0066', '#33CC00', '#FF99CC', '#00CCFF', '#CC3300',
        '#66CCFF', "#CDEEAA", "#EFB960", "#C3AED6", "#874C62",
        "#EA907A", "#2192FF", "#557571", "#F9FD50", "#B2EC5D"
    ]

    if color_palette == 'Monochrome':
        categories_unique = ['All']
        df_plot['Category'] = 'All'
        colors = { 'All': mono_color }
        color_palette_ssgsea = [mono_color]
    else:
        df_plot['Category'] = df_plot[color_palette]
        categories_unique = sorted(df_plot['Category'].dropna().unique())
        color_cycle = cycle(color_pal)  # Prepare a cycle of fallback colors
        colors = {category: color_map.get(category, next(color_cycle)) for category in categories_unique}

        # Identify all unique categories within the selected column
        all_categories = df_plot[color_palette].dropna().unique()

        # Create a mapping of categories to their colors
        full_color_map = {cat: color_map.get(cat, next(color_cycle)) for cat in all_categories}

        # Map the colors to a new column in the DataFrame
        color_list = pd.Series(df_plot[color_palette].map(full_color_map))

        # Use the color column for the color in strip plots
        color_palette_ssgsea = color_list.unique().tolist()

    # Define the range for the violin plot based on the min and max values of the data
    min_val = df_plot[gene_name].min()
    max_val = df_plot[gene_name].max()

    # ------------------------- Custom Ordering -------------------------

    # Check the number of unique categories
    unique_categories = df_plot['Category'].unique()
    num_unique_categories = len(unique_categories)

    if num_unique_categories > 1:
        # Let users select the order of categories using multiselect
        sorted_categories = multi_col.multiselect(
                            ":blue[Select Categories in the desired order on Violin Plots] "
                            ":red[Select all Clusters from 'Select Cluster' Selection']",
            options=sorted(categories_unique),  # Providing a sorted list for initial selection
            default=sorted(categories_unique),  # Pre-selecting all categories in sorted order
            key="custom_order_multiselect_fusion_" + key  # Ensure unique key per user
        )

        # Check if all categories are selected
        if len(sorted_categories) != num_unique_categories:
            missing_categories = set(categories_unique) - set(sorted_categories)
            multi_col.warning(f":red[You have not selected all clusters. ] :blue[Missing clusters: ] "
                       f":red[{', '.join(missing_categories)}]")
            st.stop()  # Stop execution until all categories are selected
    else:
        # Only one unique category exists
        sorted_categories = unique_categories.tolist()

    # Set 'Category' as categorical with the desired order
    df_plot['Category'] = pd.Categorical(
        df_plot['Category'],
        categories=sorted_categories,
        ordered=True
    )

    # ------------------------- Plotting -------------------------

    fig = go.Figure()

    for category in sorted_categories:
        df_subset = df_plot[df_plot['Category'] == category]
        color = colors[category] if color_palette != 'Monochrome' else mono_color
        current_min = df_subset[gene_name].min()
        current_max = df_subset[gene_name].max()

        # Add Box trace (for jitter/strip plot)
        fig.add_trace(go.Box(
            x=df_subset['Category'],
            y=df_subset[gene_name],
            width=0.15,
            line=dict(color='black', width=3),
            fillcolor="white",
            opacity=0.4,
            boxpoints=False,  # Disable box points
            jitter=0.4,
            showlegend=False,
            name=f"{category} Box"
        ))

        # Add Violin trace
        fig.add_trace(go.Violin(
            x=df_subset['Category'],
            y=df_subset[gene_name],
            line_width=2,
            name=category,
            line_color='black',
            fillcolor=color,
            opacity=0.75,
            box_visible=False,
            meanline_visible=False,
            showlegend=False,
            legendgroup=category,
            width=0,
            points='all',
            pointpos=0,
            spanmode='hard',
            span=[current_min, current_max],
            marker=dict(color=color, size=10, opacity=1, line=dict(width=1.5, color='black'))
        ))

    # Update layout for bold axis titles and Arial font style
    fig.update_layout(
        title=f"{gene_name} {plot_text} {plot_title} Distribution",
        xaxis_title=dict(text=color_palette if color_palette != "Monochrome" else "Category",
                        font=dict(family="Arial", size=16, color='black')),
        yaxis_title=dict(text=gene_name, font=dict(family="Arial", size=16, color='black')),
        font=dict(family="Arial"),
        violinmode='overlay',
        boxmode='overlay'
    )

    fig.update_xaxes(
        tickfont=dict(color='black'),  # Set tick font color
        title_font=dict(color='black'),  # Set title font color
        categoryorder='array',
        categoryarray=sorted_categories  # Ensure x-axis follows the sorted order
    )
    fig.update_yaxes(
        tickfont=dict(color='black'),  # Set tick font color
        title_font=dict(color='black')  # Set title font color
    )

    fig.update_traces(marker_line_width=1)

    # ------------------------- Statistical Comparisons -------------------------

    # Calculate the p-values for all unique pairs of categories
    p_value_records = []
    p_values = {}
    for (cat1, cat2) in combinations(sorted_categories, 2):
        data1 = df_plot[df_plot['Category'] == cat1][gene_name]
        data2 = df_plot[df_plot['Category'] == cat2][gene_name]
        stat, p = mannwhitneyu(data1, data2, alternative='two-sided')
        p_values[(cat1, cat2)] = p

        # Store the results as a tuple in the list
        p_value_records.append((cat1, cat2, p))

    p_value_df = pd.DataFrame(p_value_records, columns=['Category1', 'Category2', 'P-Value'])
    p_value_df.sort_values(by=['P-Value'], inplace=True)

    p_value_df['P-Value'] = p_value_df['P-Value'].apply(lambda x: f'{x:.2e}')

    # For annotation's position
    y_values = df_plot[gene_name].values
    y_min, y_max = y_values.min(), y_values.max()

    # Set the first annotation 10% above the maximum y value
    first_annotation_y = y_max + 0.1 * (y_max - y_min)

    # Calculate the range for subsequent annotations
    significant_pvals = [p for p in p_values.values() if p < 0.05]
    num_annotations = len(significant_pvals)
    if num_annotations > 1:
        annotation_spacing = (y_max + 0.45 * (y_max - y_min) - first_annotation_y) / (num_annotations - 1)
    else:
        annotation_spacing = 0

    # Add annotations on the figure based on the p-values
    current_y = first_annotation_y
    for (cat1, cat2), p in p_values.items():
        if p < 0.05:  # If p-value is less than 0.05, it is considered significant
            x1 = sorted_categories.index(cat1)
            x2 = sorted_categories.index(cat2)

            # Add the significance bar
            fig.add_shape(type='line',
                          x0=x1,
                          x1=x2,
                          y0=current_y,
                          y1=current_y,
                          line=dict(color='black', width=1))

            # Determine the significance text
            if 0.01 < p < 0.05:
                text_anno = '*'
            elif 0.001 < p < 0.01:
                text_anno = '**'
            else:
                text_anno = '***'

            # Add the significance text
            fig.add_annotation(x=(x1 + x2) / 2,
                               y=current_y,
                               text=text_anno,
                               showarrow=False,
                               font=dict(size=16, color='black'))

            # Increment y position for next annotation
            current_y += annotation_spacing

    # ------------------------- Saving the Plot -------------------------

    user_dir = f"/result_data_{user_id}"
    if comp_type == "ssgsea":
        file_path = f"result_data{user_dir}/ssgsea_comp_violin_{user_id}{axis_id}.pdf"
    elif comp_type == "gc_ssgsea":
        file_path = f"result_data{user_dir}/gene_ssgsea_corr_{user_id}{axis_id}.pdf"
    else:
        file_path = f"result_data{user_dir}/gene_compare_corr_{user_id}{axis_id}.pdf"

    fig.write_image(file_path, width=600, height=600, scale=1)

    return fig, p_value_df


def dot_scatter_ssgsea(merged, color_list_plotly, th_num, th_num1, th_num2, sep_opt, ssgsea_dot_path):
    """
    Create a scatter plot of ssGSEA scores with additional formatting.

    Parameters:
    merged (DataFrame): The input data containing ssGSEA scores and group information.
    color_list_plotly (list): List of colors for the plot.
    th_num (float): Threshold number for a horizontal line.
    th_num1 (float): First threshold number for a horizontal line.
    th_num2 (float): Second threshold number for a horizontal line.
    sep_opt (str): Option to determine the type of horizontal lines to add.
    ssgsea_dot_path (str): Path to save the scatter plot image.

    Returns:
    fig (Figure): A Plotly figure object.
    """

    # Determine marker size based on the number of rows in the merged DataFrame
    if merged.shape[0] > 500:
        sct_mark_size = 10
    elif 500 > merged.shape[0] > 150:
        sct_mark_size = 12.5
    else:
        sct_mark_size = 15

    fig = px.strip(merged, y="ssGSEA", color='Group', color_discrete_sequence=color_list_plotly, orientation='v',
                   stripmode='overlay')
    fig.update_layout(legend=dict(traceorder='reversed'))

    # Add horizontal lines based on the sep_opt parameter
    if sep_opt != 'Percentage decomposition (lower and upper limit)':
        fig.add_hline(y=th_num, line_color="red", line_width=2, line_dash="dash")
    else:
        fig.add_hline(y=th_num1, line_color="green", line_width=2, line_dash="dash")
        fig.add_hline(y=th_num2, line_color="red", line_width=2, line_dash="dash")

    fig.update_traces(marker_size=sct_mark_size, marker_line_width=1, marker_line_color="black")

    # Update axis properties
    fig.update_xaxes(
        tickfont=dict(color='black', family="Arial"),  # Set tick font color and family
        title_font=dict(color='black', family="Arial", size=16)  # Set title font color, family, and size
    )
    fig.update_yaxes(
        tickfont=dict(color='black', family="Arial"),  # Set tick font color and family
        title_font=dict(color='black', family="Arial", size=16)  # Set title font color, family, and size
    )

    # Update layout for legend
    fig.update_layout(
        font=dict(family="Arial"),
        legend=dict(
            title=dict(font=dict(size=16, color='black', family="Arial"))  # Set legend title font
        )
    )

    # Save the plot as an image
    fig.write_image(ssgsea_dot_path, width=1200, height=1000, scale=1)

    return fig


def fdr_p_plot(results_df, lower_limit, upper_limit, pv_dot_name):
    """
    Create a scatter plot of P_Value and Corrected_P_Value with additional formatting.

    Parameters:
    results_df (DataFrame): The input data containing P_Value and Corrected_P_Value.
    lower_limit (float): Lower limit for threshold percentage.
    upper_limit (float): Upper limit for threshold percentage.
    pv_dot_name (str): The name of the output image file.

    Returns:
    fig_p_value (Figure): A Plotly figure object.
    """

    # Create the scatter plot for P_Value with lines connecting the points
    trace1 = go.Scatter(
        x=results_df.index,  # Use Threshold_Percentage for x-axis values
        y=results_df['P_Value'],
        mode='lines+markers',  # This will draw both lines and markers
        marker=dict(color='blue'),
        name='P_Value',
        hoverinfo='text',
        text=results_df['hover_text_pv']  # Updated hover text
    )

    # Create the scatter plot for Corrected_P_Value with lines connecting the points
    trace2 = go.Scatter(
        x=results_df.index,  # Use Threshold_Percentage for x-axis values
        y=results_df['P_Value_FDR'],
        mode='lines+markers',  # This will draw both lines and markers
        marker=dict(color='red'),
        name='P_Value_FDR',
        hoverinfo='text',
        text=results_df['hover_text_pvfdr']  # Updated hover text
    )

    layout = go.Layout(
        title='P Value and Corrected P Value by Threshold Percentage',
        xaxis=dict(
            title='Threshold Percentage',
            range=[0, 100],  # Set x-axis to range from 0 to 100
            tick0=0,  # Starting tick
            dtick=5  # Interval between ticks
        ),
        yaxis=dict(title='P Value', showgrid=False, type="log", tickvals=[0, 0.01, 0.05, 0.1, 1],
                   ticktext=['0', '0.01', '0.05', '0.1', '1'], autorange=True),
        hovermode='closest',
        shapes=[
            # Rectangle for the area below the lower limit
            {'type': 'rect', 'x0': 0, 'y0': 0, 'x1': lower_limit, 'y1': 1, 'xref': 'x', 'yref': 'paper',
             'fillcolor': 'lightgrey', 'opacity': 0.5, 'layer': 'below', 'line_width': 0},
            # Rectangle for the area above the upper limit
            {'type': 'rect', 'x0': upper_limit, 'y0': 0, 'x1': 100, 'y1': 1, 'xref': 'x', 'yref': 'paper',
             'fillcolor': 'lightgrey', 'opacity': 0.5, 'layer': 'below', 'line_width': 0},
            # Line Vertical for lower limit
            {'type': 'line', 'x0': lower_limit, 'y0': 0, 'x1': lower_limit, 'y1': 1, 'xref': 'x', 'yref': 'paper',
             'line': {'color': 'darkgrey', 'width': 1}},
            # Line Vertical for upper limit
            {'type': 'line', 'x0': upper_limit, 'y0': 0, 'x1': upper_limit, 'y1': 1, 'xref': 'x', 'yref': 'paper',
             'line': {'color': 'darkgrey', 'width': 1}},
        ]
    )

    # Combine the traces and layout in a Figure
    fig_p_value = go.Figure(data=[trace1, trace2], layout=layout)

    # Determine marker size based on the number of rows in the results DataFrame
    point_sizer = len(results_df.index)
    if point_sizer >= 50:
        msize = 5
        min_size_p = 7.5
    else:
        msize = 10
        min_size_p = 15

    fig_p_value.update_traces(marker_size=msize, marker_line_width=1, marker_line_color="black")

    # Update axis properties
    fig_p_value.update_xaxes(
        tickfont=dict(color='black'),  # Set tick font color
        title_font=dict(color='black')  # Set title font color
    )
    fig_p_value.update_yaxes(
        tickfont=dict(color='black'),  # Set tick font color
        title_font=dict(color='black')  # Set title font color
    )

    fig_p_value.add_hline(y=0.05, line_color="red", line_width=1, line_dash="dash")

    # Highlight the minimum P_Value point
    fig_p_value.add_trace(
        go.Scatter(
            mode='markers',
            x=[results_df['P_Value'].idxmin()],
            y=[results_df['P_Value'].min()],
            marker=dict(
                symbol="circle-dot",
                color='yellow',
                size=min_size_p,
                line=dict(
                    color='black',
                    width=2
                )
            ),
            showlegend=False
        )
    )

    # Save the plot as an image
    fig_p_value.write_image(pv_dot_name, width=1200, height=1000, scale=1)

    return fig_p_value


def dot_scatter_sg(merged_df, selected_gene, cluster_colors, th_num, th_num1, th_num2, sep_opt):
    """
    Create a scatter plot for the expression profile of a selected gene with additional formatting.

    Parameters:
    merged_df (DataFrame): The input data containing gene expression values and group information.
    selected_gene (str): The gene for which the expression profile is to be plotted.
    cluster_colors (list): List of colors for the plot.
    th_num (float): Threshold number for a horizontal line.
    th_num1 (float): First threshold number for a horizontal line.
    th_num2 (float): Second threshold number for a horizontal line.
    sep_opt (str): Option to determine the type of horizontal lines to add.

    Returns:
    fig (Figure): A Plotly figure object.
    """

    # Determine marker size based on the number of rows in the merged DataFrame
    if merged_df.shape[0] > 500:
        sct_mark_size = 10
    elif 500 > merged_df.shape[0] > 150:
        sct_mark_size = 12.5
    else:
        sct_mark_size = 15

    # Create strip plot
    fig = px.strip(merged_df, y=selected_gene, title=f"Expression profile of {selected_gene}",
                   color='Group', color_discrete_sequence=cluster_colors, orientation='v',
                   stripmode='overlay', hover_data=merged_df[[selected_gene]])

    fig.update_traces(marker_size=sct_mark_size, marker_line_width=1, marker_line_color="black")

    # Update axis properties
    fig.update_xaxes(
        tickfont=dict(color='black', family="Arial"),  # Set tick font color and family
        title_font=dict(color='black', family="Arial", size=16)  # Set title font color, family, and size
    )
    fig.update_yaxes(
        tickfont=dict(color='black', family="Arial"),  # Set tick font color and family
        title_font=dict(color='black', family="Arial", size=16)  # Set title font color, family, and size
    )

    # Update layout for legend
    fig.update_layout(
        font=dict(family="Arial"),
        legend=dict(
            title=dict(font=dict(size=16, color='black', family="Arial"))  # Set legend title font
        )
    )

    if merged_df['Group'].iloc[0] == 'Low' or merged_df['Group'].iloc[0] == 'Intermediate':
        fig.update_layout(legend=dict(traceorder='reversed'))

    # Add horizontal lines based on the sep_opt parameter
    if sep_opt != 'Percentage decomposition (lower and upper limit)':
        fig.add_hline(y=th_num, line_color="red", line_width=2, line_dash="dash")
    else:
        fig.add_hline(y=th_num1, line_color="green", line_width=2, line_dash="dash")
        fig.add_hline(y=th_num2, line_color="red", line_width=2, line_dash="dash")

    return fig


def con_hm(_model, cmap, sample_names1, nmf_cons_name, e, gene_clust_chb):
    """
    Visualize the consensus matrix heatmap.

    Parameters:
    _model: Trained NMF model.
    cmap: Color map for the heatmap.
    sample_names1: Sample names for the heatmap.
    nmf_cons_name: Filename for saving the heatmap.

    Returns:
    None
    """
    con_mat_h = pd.DataFrame(_model.consensus_matrix_h['sim1'])
    fig = px.imshow(con_mat_h, color_continuous_scale=cmap,
                    labels=dict(x=con_mat_h.empty, color="Similarity"),
                    x=sample_names1, y=sample_names1,
                    title="Samples")
    # Update axis properties
    fig.update_xaxes(tickfont=dict(color='black'), title_font=dict(color='black'))
    fig.update_yaxes(tickfont=dict(color='black'), title_font=dict(color='black'))

    # Save the heatmap as an image
    fig.write_image(nmf_cons_name, width=800, height=800, scale=2)

    if gene_clust_chb:
        try:
            # Display the plot in Streamlit
            e[0].plotly_chart(fig, use_container_width=True, theme=None)
        except Exception as e:
            e[0].image("style_items/error.svg", use_column_width=True)
            e[0].warning(":blue[There is a problem with the ] :red[chart! ] :blue[Please try another ] :red[setting ] "
                            ":blue[or ] :red[refresh ] :blue[the page]")
            print(f"Error: {e}")
    else:
        con_hom_col = st.columns(3)
        try:
            # Display the plot in Streamlit
            con_hom_col[1].plotly_chart(fig, use_container_width=True, theme=None)
        except Exception as e:
            con_hom_col[1].image("style_items/error.svg", use_column_width=True)
            con_hom_col[1].warning(":blue[There is a problem with the ] :red[chart! ] :blue[Please try another ] :red[setting ] "
                            ":blue[or ] :red[refresh ] :blue[the page]")
            print(f"Error: {e}")
