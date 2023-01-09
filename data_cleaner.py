#%%

import os 
import pandas as pd
import datetime
import seaborn as sns

sns.set()


#%%
T2Month = {'T1':'03', 'T2':'06', 'T3':'09', 'T4':'12'}

def clean_data():
    for file in os.listdir('data'):
        if file.endswith('.csv'):
            if file.startswith('insee'):
                name = file.split('.')[0][6:]
                df = pd.read_csv('data/' + file, sep=';', header=3)
                df.drop(df.columns[2], axis=1, inplace=True)
                df.rename( columns={'Unnamed: 1':name}, inplace=True )
                if 'T' in df.Period.iloc[0]:
                    df.drop(df.tail(1).index, axis=0, inplace=True)
                    df.Period = df.Period.apply(lambda x: f"{x.split('-')[0]}-{T2Month[x.split('-')[1]]}" )
                df.Period = df.Period.astype('datetime64[ns]')
                df = df.set_index('Period')
            if file.startswith('wb'):
                name = file.split('.')[0][3:]
                df = pd.read_csv('data/' + file, sep=',', header=2)
                df = df[df['Country Name'] == 'France'].iloc[:,4:-1].T
                df.index = pd.to_datetime(df.index)
                df.rename( columns={df.columns[0]:name}, inplace=True )
            df.asfreq('Q', method='bfill', how='start')
            print(df.columns)
            df = (df - df.mean())/df.std()
            try :
                main_df = pd.merge(main_df, df, how='outer', left_index=True, right_index=True)
            except NameError:
                main_df = df
    main_df.to_csv('data/main_df.csv')

#%%
clean_data()
#%%
df = pd.read_csv("data/insee_Foreign_trade_index_Exports_Food_products.csv", sep=';', header=3)
                
# %%
df.Period = df.Period.astype('datetime64[ns]')

# %%
df = df.set_index('Period')

# %%
df = df.asfreq('Q', method='bfill', how='start')
# %%
df.plot()
# %%
df
# %%
df1 = pd.read_csv('data/insee_domestic_demand_on_GDP.csv', sep=';', header=3)
# %%

# %%

df1.Period.iloc[0].split("-")
# %%
for el in df1.Period.iloc[:-1]:
    a,b = el.split("-")
    print(f"{a}-{T2Month[b]}")
# %%
df1.plot()
# %%
df = pd.read_csv('data/wb_CentralGovDebt.csv', sep=',', header=2)
# %%

# %%
main_df = pd.read_csv('data/main_df.csv', index_col=0)
# %%
main_df.plot()
# %%
main_df['GDP(dollars)']
# %%
