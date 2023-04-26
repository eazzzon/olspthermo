import pandas as pd
import numpy as np
import sys, os

if not os.path.exists("OSAT") and os.path.exists(
    "../OSAT"
):  # hack to allow scripts to be placed in subdirectories next to AlOlivThermo:
    sys.path.insert(1, os.path.abspath(".."))

from OSAT import thermometers

df_ol = pd.read_excel("input_ol.xlsx")
df_sp = pd.read_excel("input_sp.xlsx")

results= thermometers.models(df_ol, df_sp).compute()

results.to_excel('results_pi082.xlsx', index=False)
