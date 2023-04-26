import numpy as np
import pandas as pd
from .cocvert import weight2mole, get_oxygen_n, get_cation_n


def conSpn(self, df, comps=None):
        """
          double *e,      /* comp of spinel in moles of elements                      */
          double *m,      /* comp of spinel in moles of endmember components          */
          double *r,      /* comp of spinel in terms of the independent comp var      */
          double *x,      /* comp of spinel in mole fractions of endmember comp       */
          double **dm,    /* Jacobian matrix: dm[i][j] = dr[i]/dm[j]                  */
          double ***d2m,  /* vector of matrices: d2m[i][j][k] = d2r[i]/dm[j]dm[k]     */
          double **dr,    /* Jacobian matrix: dr[i][j] = dx[i]/dr[j]                  */
          double ****d3m) /* 3rd deriv matrix: d3m[i][j][k][l]=d3r[i]/dm[j]dm[k]dm[l] */

        """

        if comps == None:
            comps = ['SiO2', 'TiO2', 'Al2O3', 'Cr2O3', 'FeO', 'MnO', 'MgO']

        dfmc = df.copy(deep = True)
        dfmc = dfmc[comps]
        dfmc = get_cation_n(dfmc, comps)[0].copy()  # dfmc has changed to moles

        dfmc['sumcat'] = dfmc['MgO'] + dfmc['Al2O3'] + dfmc['TiO2'] + dfmc['Cr2O3'] + dfmc['MnO'] + dfmc['FeO']
        dfmc['sumchg'] = 2 * (dfmc['MgO'] + dfmc['MnO']) + 3 * (dfmc['Al2O3'] + dfmc['Cr2O3']) + 4 * dfmc['TiO2']
        dfmc['fe3'] = 0
        dfmc.loc[(8 * dfmc['sumcat'] / 3 - dfmc['sumchg'] - 2 * dfmc['FeO']) > 0, 'fe3'] = 8 * dfmc['sumcat'] / 3 - dfmc['sumchg'] - 2 * dfmc['FeO']
        dfmc['fe2'] = dfmc['FeO'] - dfmc['fe3']
        dfmc['fTet'] = 1
        dfmc['fOct'] = 1
        dfmc.loc[dfmc['fe3']==0, 'fTet'] = ((1.0/3.0) * dfmc['sumcat'] + dfmc['TiO2']) / (dfmc['fe2'] + dfmc['MgO'])
        dfmc.loc[dfmc['fe3']==0, 'fOct'] = (2.0/3.0) * dfmc['sumcat'] / (2.0 * dfmc['TiO2'] + dfmc['Al2O3'] + dfmc['Cr2O3'])
        dfmc['proj'] = (dfmc['fe2'] + dfmc['MnO'] + dfmc['MgO'])/(dfmc['fe2'] + dfmc['MgO'] + dfmc['MnO'])  # proj < 1 if NiO and CoO is in the compositions
        dfmc['mFeCr2O4'] = dfmc['Cr2O3'] * dfmc['fOct'] / 2
        dfmc['mFe/MnAl2O4'] = (dfmc['fe2'] * dfmc['fTet'] + dfmc['MnO'] * dfmc['fTet'] - dfmc['Cr2O3'] * dfmc['fOct'] / 2.0 - dfmc['fe3'] / 2.0 - 2.0 * dfmc['TiO2'] * dfmc['fOct']) * dfmc['proj']

        dfmc['mFe3O4'] = dfmc['fe3'] / 2
        dfmc['mMgAl2O4'] = dfmc['MgO'] * dfmc['proj'] * dfmc['fTet']
        dfmc['mFe2TiO4'] = dfmc['TiO2'] * dfmc['fOct']
        dfmc['sum_comps'] = dfmc['mFeCr2O4'] + dfmc['mFe/MnAl2O4'] + dfmc['mFe3O4'] + dfmc['mMgAl2O4'] + + dfmc['mFe2TiO4']  # note that the code calcualte the moles of the components of the spinel but not mole fractions, \
                                                                                                                             # one should use components / sum_comps to get identical results as the web calculator
        dfmc['X2'] = dfmc['mMgAl2O4'] / dfmc['sum_comps']  #r[0]
        dfmc['X3'] = dfmc['mFeCr2O4'] / dfmc['sum_comps']  #r[1]
        dfmc['X4'] = dfmc['mFe2TiO4'] / dfmc['sum_comps']  #r[2]
        dfmc['X5'] = dfmc['mFe3O4'] / dfmc['sum_comps']    #r[3]
        dfmc['totAl'] = 2 * (1 - dfmc['X3'] - dfmc['X4'] - dfmc['X5'])
        dfmc['totCr'] = 2 * dfmc['X3']
        dfmc['totFe2'] = 1.0 - dfmc['X2'] + dfmc['X4']
        dfmc['totFe3'] = 2.0 * dfmc['X5']
        dfmc['totMg']  = dfmc['X2']
        dfmc['totTi']  = dfmc['X4']
        dfmc['ratio'] = 2.0 - dfmc['totCr'] - dfmc['totTi']      # available oct / available tet sites

        dfmc['xmg2oct'] = dfmc['totMg']  * dfmc['ratio']/(1.0 + dfmc['ratio'])
        dfmc['xfe2oct'] = dfmc['totFe2'] * dfmc['ratio']/(1.0 + dfmc['ratio'])
        dfmc['xal3oct'] = dfmc['totAl']  * dfmc['ratio']/(1.0 + dfmc['ratio'])
        dfmc['xfe3oct'] = dfmc['totFe3'] * dfmc['ratio'] / (1.0 + dfmc['ratio'])

        dfmc['xmg2tet'] = dfmc['totMg']  - dfmc['xmg2oct']
        dfmc['xfe2tet'] = dfmc['totFe2'] - dfmc['xfe2oct']
        dfmc['xal3tet'] = dfmc['totAl']  - dfmc['xal3oct']
        dfmc['xfe3tet'] = dfmc['totFe3'] - dfmc['xfe3oct']

        dfmc['xmg2oct'] /= 2.0
        dfmc['xfe2oct'] /= 2.0
        dfmc['xal3oct'] /= 2.0
        dfmc['xfe3oct'] /= 2.0


        dfmc['s1'] = dfmc['xmg2tet'] - 2 * dfmc['xmg2oct']
        dfmc['s2'] = dfmc['xal3oct'] - dfmc['xal3tet'] / 2
        dfmc['s4'] = dfmc['xfe3oct'] - dfmc['xfe3tet'] / 2

        dfmc['xmg2tet'] = (dfmc['X2'] + dfmc['s1'])/2.0
        dfmc['xfe2tet'] = dfmc['X4'] - 0.5*dfmc['X2'] - 0.5*dfmc['s1'] + dfmc['s2'] + dfmc['X3'] + dfmc['s4']
        dfmc['xal3tet'] = 1.0 - dfmc['X3'] - dfmc['X4'] - dfmc['X5'] - dfmc['s2']
        dfmc['xfe3tet'] = dfmc['X5'] - dfmc['s4']
        dfmc['xmg2oct'] = (dfmc['X2'] - dfmc['s1'])/4.0
        dfmc['xfe2oct'] = (2.0 - dfmc['X2'] + dfmc['s1'] - 2.0*dfmc['s2'] - 2.0*dfmc['X3'] - 2.0*dfmc['s4'])/4.0
        dfmc['xal3oct'] = (1.0 - dfmc['X3'] - dfmc['X4'] - dfmc['X5'] + dfmc['s2'])/2.0
        dfmc['xfe3oct'] = (dfmc['X5'] + dfmc['s4'])/2.0
        dfmc['xcr3oct'] = dfmc['X3']
        dfmc['xti4oct'] = dfmc['X4']/2.0
        dfmc.columns = [x + '_xca' for x in dfmc.columns]  # leave like this for the regression use

        # organised terms
        dfmc['reg_X2'] = dfmc['X2_xca']
        dfmc['reg_X4'] = dfmc['X4_xca']
        dfmc['reg_X5'] = dfmc['X5_xca']
        dfmc['reg_X2X4'] = dfmc['X4_xca'] * dfmc['X2_xca']
        dfmc['reg_X3X4'] = dfmc['X3_xca'] * dfmc['X4_xca']
        dfmc['reg_X2sqrt'] = dfmc['X2_xca']**0.5
        return dfmc

class AlOlivThermo:

    def eq17_zhang2022(self, df):
        df_sp = df.iloc[:, 18::].copy()
        df_sp.columns = [x.split('.')[0] for x in df_sp.columns]
        df_ol = df.iloc[:, 2:16].copy()
        df_sp = df_sp.join(conSpn(df_sp))
        df_sp['Cr#'] = df_sp['Cr2O3'] / 151.99 / (df_sp['Cr2O3'] / 151.99 + df_sp['Al2O3'] / 101.96)
        df_sp['lnkd'] = np.log(df_ol['Al2O3'] / df_sp['Al2O3'])
        df_sp['phi_sp'] = 8.948 * df_sp['reg_X2'] + 18.707 * df_sp['reg_X4'] + 2.065 * df_spdfmc['reg_X5'] \\
        -16.853 * df_sp['reg_X2X4'] - 17.570 * df_sp['reg_X3X4'] - 16.074 * df_sp['reg_X2sqrt'] + 0.904 * df_sp['Cr#']

        return 10000 / (8.578 - 0.756 * df_sp['lnkd'] + df_sp['phi_sp']) -273.15

    def eq18_zhang2022(self, df):

        df_sp = df.iloc[:, 18::].copy()
        df_sp.columns = [x.split('.')[0] for x in df_sp.columns]
        df_ol = df.iloc[:, 2:16].copy()
        df_sp = df_sp.join(conSpn(df_sp))
        df_sp['Cr#'] = df_sp['Cr2O3'] / 151.99 / (df_sp['Cr2O3'] / 151.99 + df_sp['Al2O3'] / 101.96)
        df_sp['lnkd'] = np.log(df_ol['Al2O3'] / df_sp['Al2O3'])

        return 10000 / (0.832 + 0.938 * df_sp['Cr#'] - 0.878 * df_sp['lnkd']) - 273.15

    def coogan2014(self, df):

        df_sp = df.iloc[:, 18::].copy()
        df_sp.columns = [x.split('.')[0] for x in df_sp.columns]
        df_ol = df.iloc[:, 2:16].copy()
        df_sp = df_sp.join(conSpn(df_sp))
        df_sp['Cr#'] = df_sp['Cr2O3'] / 151.99 / (df_sp['Cr2O3'] / 151.99 + df_sp['Al2O3'] / 101.96)
        df_sp['lnkd'] = np.log(df_ol['Al2O3'] / df_sp['Al2O3'])

        return 10000 / (0.575 + 0.884 * df_sp['Cr#'] - 0.897 * df_sp['lnkd']) - 273.15



