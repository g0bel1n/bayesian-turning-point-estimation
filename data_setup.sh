#!/bin/bash

printf 'Data Aggregation script \n'

DIR="data"
if [ -d "$DIR" ]; then
  # Take action if $DIR exists. #
  echo "Directory $DIR exists."
else
  #  Control will jump here if $DIR does NOT exists #
  echo "Directory $DIR does not exists."
  mkdir $DIR
fi

cd data
today=$(date +'%d%m%Y')
#declare -a keys = ("001625988","001626407")
series_insee=('001625988;insee_Foreign_trade_index_Imports_agriculture' 
        '001626407;insee_Foreign_trade_index_Exports_Food_products'
        '010565680;insee_domestic_demand_on_GDP'
        )

series_worldbank=('wb_GDP(dollars);https://api.worldbank.org/v2/en/indicator/NY.GDP.MKTP.CD?downloadformat=csv'
                    'wb_CentralGovDebt;https://api.worldbank.org/v2/en/indicator/GC.DOD.TOTL.GD.ZS?downloadformat=csv'
                    'wb_GNI(dollars);https://api.worldbank.org/v2/en/indicator/NY.GNP.MKTP.CD?downloadformat=csv'
                    'wb_ConsumerPrices;https://api.worldbank.org/v2/en/indicator/FP.CPI.TOTL.ZG?downloadformat=csv'
                    'wb_LaborParticipation;https://api.worldbank.org/v2/en/indicator/SL.TLF.CACT.ZS?downloadformat=csv'
                    'wb_Unemployment;https://api.worldbank.org/v2/en/indicator/SL.UEM.TOTL.ZS?downloadformat=csv')


for val in "${series_insee[@]}"; do

    IFS=";" read -r -a arr <<< "${val}"
    key="${arr[0]}"
    name="${arr[1]}"
    if [ -f "$name.csv" ]; then
      # Take action if $DIR exists. #
      echo "Serie $name exists."
    else
      #  Control will jump here if $DIR does NOT exists #
      echo "Serie $name does not exists."
      curl -JLO "https://www.insee.fr/en/statistiques/serie/telecharger/csv/${key}?ordre=chronologique&transposition=donneescolonne&periodeDebut=1&anneeDebut=1999&periodeFin=5&anneeFin=2022&revision=sansrevisions"  -P data
      unzip serie_${key}_${today}.zip  && rm serie_${key}_${today}.zip
      mv *_values.csv ${name}.csv
      rm characteristics.csv
    fi
done

for val in "${series_worldbank[@]}"; do

    IFS=";" read -r -a arr <<< "${val}"
    key="${arr[1]}"
    name="${arr[0]}"
    if [ -f "$name.csv" ]; then
      # Take action if $DIR exists. #
      echo "Serie $name exists."
    else
      #  Control will jump here if $DIR does NOT exists #
      echo "Serie $name does not exists."
      curl -o ${name}.zip -JL $key --create-dirs
      unzip ${name}.zip && rm ${name}.zip
      rm -rf Metadata_*
      mv API_*.csv ${name}.csv
    fi
done

cd ..
  
