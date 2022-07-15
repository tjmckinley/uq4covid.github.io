The file `dataProcess.R` cleans up the data and formats correctly for modelling.

## LTLA data request

Deaths within 28 days of a positive test

https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=cumDeaths28DaysByDeathDate&format=csv

## Region

Deaths by age

https://api.coronavirus.data.gov.uk/v2/data?areaType=region&metric=newDeaths28DaysByDeathDateAgeDemographics&format=csv

## LAD to Region lookup

https://geoportal.statistics.gov.uk/datasets/ons::local-authority-district-to-region-april-2019-lookup-in-england/explore

## NHS Region

Cumulative admissions by age, and total hospital cases. Age data is not complete, so total over
age-classes is not the same as total hospitalisations and will need to be corrected for.

https://api.coronavirus.data.gov.uk/v2/data?areaType=nhsRegion&metric=cumAdmissionsByAge&format=csv

https://api.coronavirus.data.gov.uk/v2/data?areaType=nhsRegion&metric=hospitalCases&metric=cumAdmissions&format=csv

