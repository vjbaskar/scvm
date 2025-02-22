#!/usr/bin/env python3
import pandas as df


def clean_summary(df):
    """
    This function cleans the summary file generated by cellranger count
    :param df: pandas dataframe
    """
    df['Reads Mapped to Genome'] = df['Reads Mapped to Genome'].str.replace("%","").astype(float)
    df['Reads Mapped Confidently to Genome'] = df['Reads Mapped Confidently to Genome'].str.replace("%","").astype(float)
    df['Reads Mapped Confidently to Intronic Regions'] = df['Reads Mapped Confidently to Intronic Regions'].str.replace("%","").astype(float)
    df['Fraction Reads in Cells'] = df['Fraction Reads in Cells'].str.replace("%","").astype(float)
    df['Total Genes Detected'] = df['Total Genes Detected'].str.replace(",","").astype(int)
    df['Median UMI Counts per Cell'] = df['Median UMI Counts per Cell'].str.replace(",","").astype(int)
    df['Estimated Number of Cells'] = df['Estimated Number of Cells'].str.replace(",","").astype(int)
    df['Mean Reads per Cell'] = df['Mean Reads per Cell'].str.replace(",","").astype(int)
    df['Median Genes per Cell'] = df['Median Genes per Cell'].str.replace(",","").astype(int)
    df['Number of Reads'] = df['Number of Reads'].str.replace(",","").astype(int)
    df['Reads Mapped Confidently to Transcriptome'] = df['Reads Mapped Confidently to Transcriptome'].str.replace("%","").astype(float)
    df['Sequencing Saturation'] = df['Sequencing Saturation'].str.replace("%","").astype(float)
    df['Reads Mapped Confidently to Intergenic Regions'] = df['Reads Mapped Confidently to Intergenic Regions'].str.replace("%","").astype(float)
    df['Reads Mapped Confidently to Intronic Regions'] = df['Reads Mapped Confidently to Intronic Regions'].str.replace("%","").astype(float)
    df['Reads Mapped Confidently to Intronic Regions'] = df['Reads Mapped Confidently to Intronic Regions'].str.replace("%","").astype(float)

    return df

