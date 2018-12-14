Pink and Chum Salmon Spawner-Recruit Database
=============================================


Summary 
------- 

This database contains stock-recruitment data for 46 pink and 53 chum salmon
stocks throughout their North American ranges. 

The database is composed of six comma separated value files (.csv):
  1. chum_info.csv
  2. chum_data.csv
  3. chum_sources.csv
  4. pink_info.csv
  5. pink_data.csv
  6. pink_sources.csv

Data were compiled and processed by Michael Malick as part of his PhD research.
Please see below for a description of each file and documentation for each
column/variable contained in the files. Several researchers are currently
conducting research using these data, therefore, before starting any research it
is advised that you contact the maintainer of the data (Michael Malick
<malickmj@gmail.com>) to make sure similar research is not already being pursued
by other users of the data. 



chum_info.csv and pink_info.csv
-------------------------------

These files contain detailed information about each salmon stock in the database
(e.g., ocean entry locations, region, source, comments, etc.). The `stock.id`
variable links each stock across data files, similar to a relational database
key. 
    
#### Columns
   
 - `stock.id`:     Unique stock identifier. The stock.id links each stock
                   across all files

 - `species`:      Species of salmon

 - `stock`:        Name of the salmon stock

 - `region`:       Geographic region the stock is located in. There are a total
                   of 15 regions (see below)

 - `sub.region`:   Geographic sub-region the stock is located in. Sub-regions
                   are unique for each species

 - `jurisdiction`: The jurisdiction that is primarily responsible for the stock
                   (i.e., state or province)

 - `lat`:          Latitude of ocean entry location

 - `lon`:          Longitude of ocean entry location

 - `source.id`:    Data source identifier. Links to the 'chum_sources.csv' and
                   'pink_sources.csv' files

 - `comment`:      Additional comments regarding the data for the stock



chum_data.csv and pink_data.csv
-------------------------------

The spawner and recruitment data series for each stock. All data are
aligned by brood year (year of spawning) and are in number of salmon. 


#### Columns
 
 - `stock.id`:   See chum_info.csv

 - `species`:    See chum_info.csv

 - `stock`:      See chum_info.csv

 - `region`:     See chum_info.csv

 - `sub.region`: See chum_info.csv

 - `brood.yr`:   Brood year, that is, the year in which the number of spawners
                 were enumerated on their way to the spawning grounds

 - `spawners`:   Total spawners for the brood year

 - `recruits`:   Total recruitment (i.e., offspring) resulting from the
                 spawners for the brood year. This is the sum of the recruit.x
                 columns for chum salmon

 - `use`:        Should this data point be used for analysis? For completeness,
                 all available brood years were included. However, some years
                 of data were deemed unreliable, were missing either spawner
                 or recruitment data, or were anomaously high or low values
                 that were suspicious (whenever possible these data points were
                 clarified with the data provider). Therefore, this column is
                 used to identify the reliable data points that should be
                 included in analyses. A 1 indicates that the data are useable
                 and a blank cell indicates the data should not be included in
                 the analyses.

 - `recruits.x`  The number of recruits returning at a given age (ages 2-6).

 - `age`         The age assumption used to calculate total recruits for chum
                 salmon. There are three values for this variable: 'data'
                 indicates that age data were used for the particular year;
                 'avg' indicates that the average age proportions for the given
                 stock were used for the particular year; and 'age4' indicates
                 that all recruits were assumed to be age-4. The average method
                 was only used if there were at least data for half of the
                 brood years for that stock. This follows the methods outlined
                 in Pyper (2002). [Pyper, B.J., Mueter, F.J., Peterman, R.M.,
                 Blackbourn, D.J., and Wood, C.C. 2002. Spatial covariation in
                 survival rates of Northeast Pacific chum salmon. Transactions
                 of the American Fisheries Society 131: 343–363.]



chum_sources.csv and pink_sources.csv
-------------------------------------

Information about each data source. These sources correspond to the `source.id`
column of the 'chum_info.csv' and 'pink_info.csv' files. Each source represents
the primary individual who provided the data.

#### Columns
 
 - `source.id`: See chum_info.csv

 - `source`:    Gives the name and work location of the person who provided the
                data or the citation of the publication the data were obtained
                from.



Regions
-------

Each stock is included in one of 15 regions. The regions are consistent across
species (note: sub-regions are species specific):

1.  Outside Washington (WA)
2.  Inside Washington (WA)
3.  Fraser River
4.  BC South
5.  BC Central
6.  BC North
7.  Southeast Alaska (SEAK)
8.  Yakutat
9.  Prince William Sound (PWS)
10. Cook Inlet
11. Kodiak
12. Chignik
13. AK Peninsula
14. Bristol Bay
15. Arctic-Yukon-Kuskokwim (AYK)


