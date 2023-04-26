import numpy as np
import pandas as pd
import periodictable as pt


def normalize(df, comps):
    """
    Normalize oxides data (weight percentage or moles)

    Parameters
    -----------

    df : :class:`pandas.DataFrame`
        dataframe to be normalized

    comps : :class:`list`
        list of optional oxides to be normalized
    Returns
    -------
        dfmc : :class:`pandas.DataFrame`
        normalized dataframe
    """
    dfmc = df.copy(deep = True)
    if comps:
        cmpnts = [c for c in comps if c in dfmc.columns]
        dfmc.loc[:, cmpnts] = 100 * dfmc.loc[:, cmpnts].divide(
            dfmc.loc[:, cmpnts].sum(axis=1).replace(0, np.nan), axis=0
            )
    else:
        dfmc = 100 * dfmc.divide(dfmc.sum(axis=1).replace(0, 100.0), axis=0)
    return dfmc

def weight2mole(df, comps=None, norm_fractions = True):
    """
    Transform weight percentage to mole (fractions)

    Parameters
    -----------

        df : :class:`pandas.DataFrame`
            dataframe in weight percentage to convert to mole fractions

        comps : :class:`list`
            list of optional oxides to be converted

        norm_fractions : :class:`boolean`
            default is True
            if norm_fractions:
                sum and convert moles to 100
            else:
                moles equal division of oxides and total.
    Returns
    -------
        :class:`pandas.DataFrame`
        converted dataframe
    """
    if comps == None:
            comps = ['SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'P2O5', 'Na2O', 'K2O']
    dfmc = df.copy()
    dfmc = dfmc[comps]
    owms = [pt.formula(c).mass for c in dfmc.columns]
    if norm_fractions:
        return normalize(dfmc.div(owms), comps)
    else:
        return dfmc.div(owms)

def get_oxygen_n(df, comps):
    """
    calculate number of oxygens in oxides

    Parameters
    -----------

    df : :class:`pandas.DataFrame`
        dataframe in weight percentage to calculate oxygens in each oxide

    comps : :class:`list`
        list of optional oxides to be calculated

    Returns
    -------
        dfmc : :class:`pandas.DataFrame`
        dataframe of number of oxygens

        nocs : :class:`list`
        oxygen numbers of oxides
    """
    dfmc = weight2mole(df, comps, norm_fractions = False)
    nocs = [pt.formula(c).atoms.get(pt.elements.symbol('O')) for c in dfmc.columns]
    dfmc = dfmc.mul(nocs)
    return dfmc, nocs

def get_cation_n(df, comps):
    """
    calculate number of oxygens in oxides

    Parameters
    -----------

    df : :class:`pandas.DataFrame`
        dataframe in weight percentage to calculate cations in each oxide

    comps : :class:`list`
        list of optional oxides to be calculated

    Returns
    -------
        dfmc : :class:`pandas.DataFrame`
        dataframe of number of cations
        ncas : :class:`list`
        cation numbers of oxides
    """
    dfmc = weight2mole(df, comps, norm_fractions = False)
    ncas = [list(pt.formula(c).atoms.values())[0] for c in dfmc.columns]
    dfmc = dfmc.mul(ncas)
    return dfmc, ncas









