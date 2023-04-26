import pandas as pd
import numpy as np
import olivthermo

AlOlivThermo = olivthermo.AlOlivThermo()

df_temp_cry = pd.read_excel('AlOlvNatural.xlsx')

df_temp_cry['T_eq17'] = AlOlivThermo.eq17_zhang2022(df_temp_cry)['T_cry']
df_temp_cry['T_eq18'] = AlOlivThermo.eq18_zhang2022(df_temp_cry)['T_cry']
df_temp_cry['T_coogan'] = AlOlivThermo.coogan2014(df_temp_cry)['T_cry']

df_temp_cry.to_excel('t_cry.xlsx', index=False)