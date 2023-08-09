import streamlit as st
import pandas as pd
import numpy as np


def custom_format(n):
    if n > 0.001:
        return f"{n:.3f}"
    else:
        return f"{n:.3e}"


def main():
    data = {'p-value': [0.01, 0.0001, 0.1, 0.00001, 0.5],
            'fdr': [0.02, 0.0002, 0.2, 0.00002, 0.6],
            'other_column': [1, 2, 3, 4, 5]}
    df = pd.DataFrame(data)

    st.dataframe(df.style.format({'p-value': custom_format, 'fdr': custom_format}))

    return


if __name__ == '__main__':
    main()
