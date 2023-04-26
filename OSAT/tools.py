import numpy as np
import pandas as pd

def vec_mc(comp, std, i, mc):
    """
    vectorized MC calculation model

    Parameters:
    --------
    comp: :class: `pandas.Dataframe`
    std: :class: `pandas.Dataframe`
    i: :class: `numpy.array`
        index for iteration
    mc: :class: `numpy.array`
        number of monte carlo simulation

    Returns:
    ---------
    df: :class: `pandas.Dataframe`
        expanded dataframe after monte carlo simulation
    """
    comp = comp.copy()
    std = std.copy()
    index = i
    index = np.ones(mc) * index

    comp_mc =  comp.iloc[i,:]
    std_mc = std.iloc[i,:]
    comp_mc = np.ones((mc,1)) * comp_mc.to_numpy()
    std_mc = np.random.normal(0,1,(mc, len(std_mc))) * std_mc.to_numpy()
    comp_mc = comp_mc + std_mc
    df = pd.DataFrame(comp_mc)
    df.index = index
    return df




