import pandas as pd

df = pd.read_csv('aggregated_data.csv')

df2 = df.rename(columns = {'Unnamed: 0' :'date',
                           'Produit interieur brut (PIB)':'PIB',
                           'Importations':'MTR',
                           'Depense de consommation des menages':'PCR',
                           'Formation brute de capital fixe':'capital_fixe',
                           'Exportations':'XTR',
       'Monthly economic survey in services - Past trend of activity - All services activities - SA series':'past_service_activities',
                           'Consumer price index - Base 2015 - All households - France - Food':'price_food',
                           'Consumer price index - Base 2015 - All households - France - Energy':'price_energy',
                           'Monthly economic survey in services - Expected trend of demand - All services activities - SA series':'expected_service_activities',
                           'Monthly economic outlook survey in the building industry - Utilisation rate of production capacities in % - Overall - SA series':'rate_production',
                           'Consumer price index - Base 2015 - All households - France - Manufactured products':'price_products',
       'Monthly economic outlook survey in the building industry - Past trend in the workforce - Overall - SA series':'past_trend_workforce',
                           'Consumer price index - Base 2015 - All households - France - Services':'price_servies',
                           'Monthly economic outlook survey in the building industry - Trend of past activity - Overall - SA series':'building_industry_past',
                           'Monthly economic outlook survey in the building industry - Trend of expected activity - Overall - SA series':'building_industry_expected',
                           'Rent reference index (RRI)':'rent',
                           'SA-WDA industrial production index (base 100 in 2015) - Manufacturing (NAF rev. 2, level A10, item CZ)':'manufacturing_index',
                           'Monthly economic outlook survey in the building industry - Judgment on the order book level - Overall - SA series': 'building_industry_overall',
                           'Monthly economic survey in services - Expected trend of activity - All services activities - SA series':'trend_services',
'ILO unemployment rate - Total - Metropolitan France - SA data': 'unemployment_rate',
                           'Business climate summary indicator - All sectors - Metropolitan France': 'buisness_climate'
                           }, inplace = False)

df2.to_csv('nice_names_data.csv', index=False)