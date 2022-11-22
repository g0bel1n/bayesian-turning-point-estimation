#!/bin/bash

printf 'Data Aggregation script'

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
series=('001625988;Foreign_trade_index_Imports_agriculture' 
        '001626407;Foreign_trade_index_Exports_Food_products'
        '010565680;domestic_demand_on_GDP'
        )


for val in "${series[@]}"; do

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

cd ..
  
