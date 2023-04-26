import numpy as np
import pandas as pd
from OSAT.cocvert import weight2mole, get_oxygen_n, get_cation_n
import scipy.stats as st


def _thermo(df, a, b, c, d, e, f, g, h, i, j, k, l, m, n, p, r):
    """
    fitting equation for the thermodynamic model
    """
    return (
        a
        * (df["lnkdal"] + r)
        / (
            +b * df["reg_X2(1+X4-X2)"]
            + c * df["reg_(1-X2)X3"]
            + d * df["reg_(1-X2)X4"]
            + e * df["reg_(1-X2)X5"]
            + f * df["reg_X3(X3+X4+X5)"]
            + g * df["reg_X4(X3+X4+X5)"]
            + h * df["reg_X5(X3+X4+X5)"]
            + i * df["reg_X3X4"]
            + j * df["reg_X3X5"]
            + k * df["reg_X4X5"]
            + l * df["reg_X2pow0.5"]
            + m * df["reg_X2pow2"]
            + n * df["reg_X2"]
        )
        + p
    )


_coeff = np.array(
    [
        -1.68118521e-01,
        3.90197643e-01,
        9.19774257e-03,
        -2.49266643e00,
        6.46891489e-02,
        -3.07131662e-02,
        -4.14143648e00,
        -4.27844629e-01,
        4.63710203e00,
        5.45307660e-02,
        1.08027441e01,
        -6.30188015e-01,
        -5.93120349e-01,
        1.48669890e00,
        2.71379230e00,
        6.54194424e-01,
    ]
)  # fitted coefficients for thermodynamic model


def _conSpn(df):
    """
    calcualte spinel components used in the thermometer

    Parameters
    -----------------
    df : : class: `pandas.DataFrame`
        spinel dataframe
    Return
    ---------------------------
    dfmc : : class: `pandas.DataFrame`
        recalculated X terms and regression terms of spinel

    """

    comps = ["SiO2", "TiO2", "Al2O3", "Cr2O3", "FeO", "MnO", "MgO"]

    dfmc = df.copy(deep=True)
    dfmc = dfmc[comps]
    dfmc = get_cation_n(dfmc, comps)[0].copy()  # dfmc has changed to moles
    # calculate Fe3+, Fe2+, spinel components following MELTS (Gualda et al. 2012, Sack and Ghiorso 1991)
    dfmc["sumcat"] = (
        dfmc["MgO"]
        + dfmc["Al2O3"]
        + dfmc["TiO2"]
        + dfmc["Cr2O3"]
        + dfmc["MnO"]
        + dfmc["FeO"]
    )
    dfmc["sumchg"] = (
        2 * (dfmc["MgO"] + dfmc["MnO"])
        + 3 * (dfmc["Al2O3"] + dfmc["Cr2O3"])
        + 4 * dfmc["TiO2"]
    )
    dfmc["fe3"] = 0
    dfmc.loc[(8 * dfmc["sumcat"] / 3 - dfmc["sumchg"] - 2 * dfmc["FeO"]) > 0, "fe3"] = (
        8 * dfmc["sumcat"] / 3 - dfmc["sumchg"] - 2 * dfmc["FeO"]
    )
    dfmc["fe2"] = dfmc["FeO"] - dfmc["fe3"]
    dfmc.loc[dfmc["fe2"] < 0, "fe3"] = dfmc[
        "FeO"
    ]  # calculate for super oxidized spinel (in air)
    dfmc.loc[dfmc["fe2"] < 0, "fe2"] = 0
    dfmc["fTet"] = 1
    dfmc["fOct"] = 1
    dfmc.loc[dfmc["fe3"] == 0, "fTet"] = (
        (1.0 / 3.0) * dfmc["sumcat"] + dfmc["TiO2"]
    ) / (dfmc["fe2"] + dfmc["MgO"])
    dfmc.loc[dfmc["fe3"] == 0, "fOct"] = (
        (2.0 / 3.0)
        * dfmc["sumcat"]
        / (2.0 * dfmc["TiO2"] + dfmc["Al2O3"] + dfmc["Cr2O3"])
    )
    dfmc["proj"] = (dfmc["fe2"] + dfmc["MnO"] + dfmc["MgO"]) / (
        dfmc["fe2"] + dfmc["MgO"] + dfmc["MnO"]
    )  # proj < 1 if NiO and CoO is in the compositions
    # calculate spinel components following Sack and Ghiorso 1991
    dfmc["mFeCr2O4"] = dfmc["Cr2O3"] * dfmc["fOct"] / 2
    dfmc["mFe/MnAl2O4"] = (
        dfmc["fe2"] * dfmc["fTet"]
        + dfmc["MnO"] * dfmc["fTet"]
        - dfmc["Cr2O3"] * dfmc["fOct"] / 2.0
        - dfmc["fe3"] / 2.0
        - 2.0 * dfmc["TiO2"] * dfmc["fOct"]
    ) * dfmc["proj"]

    dfmc["mFe3O4"] = dfmc["fe3"] / 2
    dfmc["mMgAl2O4"] = dfmc["MgO"] * dfmc["proj"] * dfmc["fTet"]
    dfmc["mFe2TiO4"] = dfmc["TiO2"] * dfmc["fOct"]
    dfmc["sum_comps"] = (
        dfmc["mFeCr2O4"]
        + dfmc["mFe/MnAl2O4"]
        + dfmc["mFe3O4"]
        + dfmc["mMgAl2O4"]
        + +dfmc["mFe2TiO4"]
    )  # note that the code calcualte the moles of the components of the spinel but not mole fractions, \
    # one should use components / sum_comps to get identical results as the web calculator
    dfmc["X2"] = dfmc["mMgAl2O4"] / dfmc["sum_comps"]  
    dfmc["X3"] = dfmc["mFeCr2O4"] / dfmc["sum_comps"]  
    dfmc["X4"] = dfmc["mFe2TiO4"] / dfmc["sum_comps"]  
    dfmc["X5"] = dfmc["mFe3O4"] / dfmc["sum_comps"]  

    dfmc.columns = [x + "_xca" for x in dfmc.columns]

    # organised terms
    dfmc["reg_X2(1+X4-X2)"] = (1 - dfmc["X2_xca"]) * (
        1 + dfmc["X4_xca"] - dfmc["X2_xca"]
    )
    dfmc["reg_(1-X2)X3"] = (1 - dfmc["X2_xca"]) * dfmc["X3_xca"]
    dfmc["reg_(1-X2)X4"] = (1 - dfmc["X2_xca"]) * dfmc["X4_xca"]
    dfmc["reg_(1-X2)X5"] = (1 - dfmc["X2_xca"]) * dfmc["X5_xca"]
    dfmc["reg_X3(X3+X4+X5)"] = dfmc["X3_xca"] * (
        dfmc["X3_xca"] + dfmc["X4_xca"] + dfmc["X5_xca"]
    )
    dfmc["reg_X4(X3+X4+X5)"] = dfmc["X4_xca"] * (
        dfmc["X3_xca"] + dfmc["X4_xca"] + dfmc["X5_xca"]
    )
    dfmc["reg_X5(X3+X4+X5)"] = dfmc["X5_xca"] * (
        dfmc["X3_xca"] + dfmc["X4_xca"] + dfmc["X5_xca"]
    )
    dfmc["reg_X3X4"] = dfmc["X3_xca"] * dfmc["X4_xca"]
    dfmc["reg_X3X5"] = dfmc["X3_xca"] * dfmc["X5_xca"]
    dfmc["reg_X4X5"] = dfmc["X4_xca"] * dfmc["X5_xca"]
    dfmc["reg_X2pow0.5"] = dfmc["X2_xca"] ** 0.5
    dfmc["reg_X2pow2"] = dfmc["X2_xca"] ** 2
    dfmc["reg_X2"] = dfmc["X2_xca"]
    return dfmc


def _ztest(df):
    """
    Z-test protocol to select the favored temperature estimation from different models

    The protocol follows as:
    1) if T_thermo > T_empirical, T_chosen = T_thermo

    2) if not, Z-score = |T_thermo - T_empirical| / (sig_empirical^2 + sig_thermo^2)

    A width of 0.6744897 * 2 standard deviation will be required that 50% mass of a normal distribution lie between 
    mean -  |T_empirical - T_thermo|/2 and mean +  |T_empirical - T_thermo|/2

    The protocol will select the comparison with higher Z-score to support a larger possibility that the empirical model is different from thermodynamic model.

    Parameter
    -----------------------------
    df : : class: `pandas.DataFrame`
        ol-sp dataframe
    Return
    ---------------------------
    t_chosen: : class: `list`
        chosen temperature
    err_chosen : class: `list`
        chosen temperature error
    z_thermo_cr: class: `pandas.Series`
        z-score of thermodynamic model and KdCr model
    z_thermo_al: class: `pandas.Series`
        z-score of thermodynamic model and KdAl model
    """

    dfmc = df.copy()
    t_thermo = dfmc["t_thermo"]
    t_kdcr = dfmc["t_kdcr"]
    t_kdal = dfmc["t_kdal"]
    err_thermo = 23.9
    err_kdcr = 34.2
    err_kdal = 43.3
    z_thermo_cr = abs(t_thermo - t_kdcr) / (np.sqrt(err_thermo ** 2 + err_kdcr ** 2))
    z_thermo_al = abs(t_thermo - t_kdal) / (np.sqrt(err_thermo ** 2 + err_kdal ** 2))
    t_chosen_3 = []
    err_chosen_3 = []

    t_chosen_2 = []
    err_chosen_2 = []

    for idx, row in dfmc[["t_thermo", "t_kdcr", "t_kdal"]].iterrows():
        if (row["t_thermo"] > row["t_kdcr"]) and (row["t_thermo"] > row["t_kdal"]):
            t = row["t_thermo"]
            err = 23.9
        elif z_thermo_al[idx] > z_thermo_cr[idx] > 0.6744897 * 2:  # 0.6744897 comes from p-value of 50% mass of a (half) normal distribution
            t = row["t_kdal"]
            err = 43.3
        elif z_thermo_cr[idx] > z_thermo_al[idx] > 0.6744897 * 2:
            t = row["t_kdcr"]
            err = 34.2
        elif (z_thermo_cr[idx] < 0.6744897 * 2) and (z_thermo_al[idx] < 0.6744897 * 2):
            t = row["t_thermo"]
            err = 23.9
        elif z_thermo_al[idx] > 0.6744897 * 2 > z_thermo_cr[idx]:
            t = row["t_kdal"]
            err = 43.3
        elif z_thermo_cr[idx] > 0.6744897 * 2 > z_thermo_al[idx]:
            t = row["t_kdcr"]
            err = 34.2
        t_chosen_3.append(t)
        err_chosen_3.append(err)

    for idx, row in dfmc[["t_thermo", "t_kdal"]].iterrows():
        if (row["t_thermo"] > row["t_kdal"]):
            t = row["t_thermo"]
            err = 23.9
        elif z_thermo_al[idx] > 0.6744897 * 2:
            t = row["t_kdal"]
            err = 43.3
        else:
            t = row["t_thermo"]
            err = 23.9
        t_chosen_2.append(t)
        err_chosen_2.append(err)

    return t_chosen_3, err_chosen_3,t_chosen_2, err_chosen_2, z_thermo_cr, z_thermo_al


class models:
    """
    class of Al-in-olivine thermometer, intergrated the three newly proposed thermometers and Coogan et al. (2014)
    """

    def __init__(self, df_ol, df_sp):
        self.df_sp = df_sp.fillna(0).copy()
        # self.df_sp = pd.read_excel('input_sp.xlsx')
        self.df_sp_conv = _conSpn(self.df_sp)
        self.df_sp = self.df_sp.join(self.df_sp_conv)
        self.df_sp["Cr#"] = (
            self.df_sp["Cr2O3"]
            / 151.99
            / (self.df_sp["Cr2O3"] / 151.99 + self.df_sp["Al2O3"] / 101.96)
        )
        # self.df_sp.columns = [x.split(".")[0] for x in self.df_sp.columns]
        self.df_ol = df_ol.fillna(0).copy()
        self.df_ol.columns = [x + "_ol" for x in self.df_ol.columns]
        self.df_ol.rename(columns={"sample_no_ol": "sample_no"}, inplace=True)
        self.df_ol_sp = self.df_ol.merge(self.df_sp, on="sample_no", how="left")
        self.df_ol_sp["lnkdal"] = np.log(
            self.df_ol_sp["Al2O3_ol"] / self.df_ol_sp["Al2O3"]
        )

        self.df_ol_sp["lnkdcr"] = np.log(
            self.df_ol_sp["Cr2O3_ol"] / self.df_ol_sp["Cr2O3"]
        )

        self.df_ol_sp["Fo"] = (
            self.df_ol_sp["MgO_ol"]
            / 40.3044
            / (self.df_ol_sp["MgO_ol"] / 40.3044 + self.df_ol_sp["FeO_ol"] / 71.844)
        )

    def compute(self):

        self.df_ol_sp["t_thermo"] = 10000 / _thermo(self.df_ol_sp, *_coeff) - 273.15
        self.df_ol_sp["t_kdcr"] = (
            10000
            / (
                0.0488
                - 0.6572 * self.df_ol_sp["lnkdal"]
                - 0.3886 * self.df_ol_sp["lnkdcr"]
                + 0.5427 * self.df_ol_sp["Cr#"]
            )
            - 273.15
        )
        self.df_ol_sp["t_kdal"] = (
            10000
            / (
                0.7395
                - 0.8654 * self.df_ol_sp["lnkdal"]
                + 1.1438 * self.df_ol_sp["Cr#"]
            )
            - 273.15
        )
        self.df_ol_sp["t_coogan"] = (
            10000
            / (0.575 + 0.884 * self.df_ol_sp["Cr#"] - 0.897 * self.df_ol_sp["lnkdal"])
            - 273.15
        )

        t_chosen_3, err_chosen_3,t_chosen_2, err_chosen_2, z_thermo_cr, z_thermo_al = _ztest(self.df_ol_sp)
        # z_test results
        self.df_ol_sp[
            ["t_zThermoAlCr", "err_zThermoAlCr","t_zThermoAl", "err_zThermoAl", "z_thermo_cr", "z_thermo_al"]
        ] = pd.DataFrame(
            {
                "t_zThermoAlCr": t_chosen_3,
                "err_zThermoAlCr": err_chosen_3,
                "t_zThermoAl": t_chosen_2,
                "err_zThermoAl": err_chosen_2,
                "z_thermo_cr": z_thermo_cr,
                "z_thermo_al": z_thermo_al,
            }
        )
        self.df_ol_sp["p_thermo_cr"] = st.norm.sf(abs(self.df_ol_sp["z_thermo_cr"]))
        self.df_ol_sp["p_thermo_al"] = st.norm.sf(abs(self.df_ol_sp["z_thermo_al"]))
        return self.df_ol_sp
