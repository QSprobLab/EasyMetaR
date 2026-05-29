#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Process the raw feature table: treat zeros as missing, impute with 1/10 of the minimum non‑zero value per row, and normalise each column to a total intensity of 1000.
"""

import pandas as pd
import numpy as np


def process_feature_table(input_file='feature_table_raw.csv', output_file='feature_table.csv'):
    df = pd.read_csv(input_file)
    # Preserve first column (compound)
    first_col = df.iloc[:, 0]
    numeric = df.iloc[:, 1:].copy()
    # Replace zeros with NaN
    numeric = numeric.replace(0, np.nan)
    # Impute: 1/10 of row minimum
    row_mins = numeric.min(axis=1, skipna=True)
    fill_values = row_mins / 10
    numeric = numeric.fillna(fill_values)
    # Normalise columns to 1000
    col_sums = numeric.sum(axis=0)
    norm_df = numeric.div(col_sums, axis=1) * 1000
    # Concatenate first column
    out_df = pd.concat([first_col, norm_df], axis=1)
    out_df.columns = df.columns
    out_df.to_csv(output_file, index=False)
    print(f"Processed {input_file} → {output_file}")

if __name__ == "__main__":
    process_feature_table()
