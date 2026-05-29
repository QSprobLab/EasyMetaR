#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Merge positive and negative ion mode metabolomics tables.

The script expects the input tables to be in *wide* format: each column after the first represents a sample.
The columns are:
- ``compound`` (metabolite name)
- ``mz`` (theoretical m/z)
- ``rt`` (theoretical retention time)
- ``adduct`` (adduct)
- sample columns with intensities.

The merged output feature table has rows as metabolites and columns as samples.  The feature info table contains the metadata.
"""

import argparse
import os
import pandas as pd


def load_table(path):
    """Load a wide-format table and return a DataFrame.
    """
    df = pd.read_csv(path)
    if 'compound' not in df.columns:
        raise ValueError(f"Missing 'compound' column in {path}")
    return df


def merge_tables(pos, neg):
    # Load tables
    pos_df = load_table(pos)
    neg_df = load_table(neg) if neg else None
    # Identify sample columns (all except meta columns)
    meta = {'compound', 'mz', 'rt', 'adduct'}
    pos_samples = [c for c in pos_df.columns if c not in meta]
    if neg_df is not None:
        neg_samples = [c for c in neg_df.columns if c not in meta]
    else:
        neg_samples = []
    # Melt into long format
    def melt(df, samples, mode):
        long = df.melt(id_vars=['compound','mz','rt','adduct'], value_vars=samples,
                        var_name='sample', value_name='intensity')
        long['mode'] = mode
        return long
    pos_long = melt(pos_df, pos_samples, 'pos')
    neg_long = melt(neg_df, neg_samples, 'neg') if neg_df is not None else None
    long_df = pd.concat([pos_long, neg_long], ignore_index=True) if neg_df is not None else pos_long
    return long_df


def write_feature_table(long_df, out_tab, out_info):
    # Feature table: rows=compound, columns=samples
    pivot = long_df.pivot_table(index='compound', columns='sample', values='intensity', aggfunc='first').reset_index()
    pivot.to_csv(out_tab, index=False)
    # Feature info: compound, mz, rt, adduct, mode
    info = long_df.drop_duplicates(subset=['compound','mz','rt','adduct','mode'])
    info.to_csv(out_info, index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pos', required=True, help='Positive mode peak table CSV')
    parser.add_argument('--neg', required=False, help='Negative mode peak table CSV')
    parser.add_argument('--output-tab', required=True, help='Output feature table CSV')
    parser.add_argument('--output-info', required=True, help='Output feature info CSV')
    parser.add_argument('--db', required=False, help='Metabolite database CSV (unused, kept for compatibility)')
    args = parser.parse_args()

    long_df = merge_tables(args.pos, args.neg)
    write_feature_table(long_df, args.output_tab, args.output_info)

if __name__ == "__main__":
    main()
