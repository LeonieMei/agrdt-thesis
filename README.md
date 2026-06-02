# Code for SARS-CoV-2 rapid antigen test analysis as part of my PhD thesis

This repository contains code for analyzing antigen rapid diagnostic test (Ag-RDT) data
from Charité employees. The results are presented in my PhD thesis titled *Temporal analysis of SARS-CoV-2 
rapid antigen test and PCR data*. The data were collected under the supervision of and provided by J. Horn and J. Seybold (Charité Berlin). The raw data are not 
included but a description is provided. The main findings of the project were published in
[SARS-CoV-2 rapid antigen test sensitivity and viral load in newly symptomatic 
hospital employees, December 2020 to February 2022 - an observational study](https://www.thelancet.com/journals/lanmic/article/PIIS2666-5247(23)00412-3/fulltext).

## Data

There are seven data files in the `data` directory.

### agrdtDataThesisT1.tsv (N gene PCR target) and agrdtDataThesisT2.tsv (E gene PCR target)

These files contain Ag-RDT, RT-PCR, symptom and other data from the Charité cohort with 
the 
following columns (and more):

* `age`: The person's age.
* `agrdt`: Specifies the result of the Ag-RDT (0: negative, 1: positive).
* `agrdtYN`: Specifies whether an Ag-RDT was performed.
* `binDaysPostOnset`: The time period in which the test was performed relative to 
     symptom onset (0-2 or 2-7 days post symptom onset).
* `daysPostOnset`: The number of days between the onset of symptoms and presentation 
     at the test centre (with a granularity of 12h during the first 24h of symptoms 
     onset, and a granularity of 24h thereafter).
* `gender`: The person's gender (male or female).
* `hasTyping`: Specifies whether a typing RT-PCR was performed.
* `immunYN`: 
     Specifies whether the person was immunized at least once (through 
    vaccination and/or prior infection) at the time of testing.
* `immun2YN`: Specifies whether the person was immunized at least twice (through 
     vaccination and/or prior infection) at least 2 weeks prior to the time of 
    testing. Only people with no prior immunization will get a value of 0.
* `immunN`: A person's count of prior immunizations (vaccination or prior infection).
* `immunNatLeast`: A person's minimum count of prior immunizations (for some people, 
     current immunization status wasn't available at the end of the study).
* `infectionKey`: A unique ID, given to every infection (composed of the first 
    positive PCR's ID, the date of the PCR and the personHash).
* `isFirstPosPcr`: Specifies whether the RT-PCR was the first positive one in an 
    infection.
* `pcrDate`: The date the RT-PCR was performed.
* `pcrId`: A unique ID, given to every RT-PCR test.
* `pcrPositive`: Specifies whether the RT-PCR was positive.
* `personHash`: A unique ID for each person.
* `reasonPres`: The reason for presentation at the test centre.
* `recovered`: Specifies 
     whether a person had a prior infection at the time of testing.
* `samplingMonth`: The time of sampling given in the format YYYY-MM.
* `samplingMonth2`: The time of sampling given in two-months intervals (starting in Dec. 2020).
* `surveyData`: Specifies whether data from the questionnaire are available.
* `symptoms`: Specifies whether a person experienced Covid-19 symptoms.
* `testDevice`: Specifies the test device that was used (Abbott Panbio COVID-19 Rapid Test Device or Roche SARS-CoV-2 Rapid Antigen Test).
* `vaccN`: A person's count of prior vaccinations.
* `vaccNatLeast`: A person's minimum count of prior vaccination (for some people, 
     current vaccination status wasn't available at the end of the study).
* `variant`: The SARS-CoV-2 variant of the sample, based on typing RT-PCR or 
    epidemiological assignment.
* `vl`: The logarithm base 10 of the RT-PCR estimated number of viral RNA molecules 
    / mL.


... and information on the presence of specific symptoms, their duration and severity;
on prior contact to potential Covid-19 cases; and other infection-specific data.

### abbottVsRocheWildtype.tsv

Data from testing SARS-CoV-2 Wildtype samples using the Abbott Panbio COVID-19 Rapid Test Device or the Roche SARS-CoV-2 Rapid Antigen Test. Contains
 the following columns:
 
* `agrdt`: Specifies the result of the Ag-RDT (0: negative, 1: positive).
* `test`: The test device used, either the Abbott Panbio COVID-19 Rapid Test Device or the Roche SARS-CoV-2 Rapid Antigen Test (0: Abbott, 1: Roche).
* `variant`: The SARS-CoV-2 variant of the sample, based on typing RT-PCR.
* `vl`: The logarithm base 10 of the RT-PCR estimated number of viral RNA molecules / mL.
* `zVl`: Log<sub>10</sub> viral load z-scores.
 
 
### abbottVsRocheOmicron.tsv

Data from testing SARS-CoV-2 Omicron samples using the Abbott Panbio COVID
-19 Rapid Test Device or the Roche SARS-CoV-2 Rapid Antigen Test. Contains
 the following columns:
 
 * `agrdt`: Specifies the result of the Ag-RDT (0: negative, 1: positive).
 * `batch`: Specifies the batch (tests were performed in two batches).
 * `test`: The test device used, either the Abbott Panbio COVID-19 Rapid Test Device or the Roche SARS-CoV-2 Rapid Antigen Test (0: Abbott, 1: Roche).
 * `variant`: The SARS-CoV-2 variant of the sample, based on typing RT-PCR.
* `vl`: The logarithm base 10 of the RT-PCR estimated number of viral RNA molecules / mL.
 * `zVl`: Log<sub>10</sub> viral load z-scores.

### nantigenViralload.tsv

Data from the Charité Institute of Virology cohort. Contains temporal viral 
load and N-antigen measurements. The columns are:

* `Day (symptom onset)`: Specifies the number of days between the day of sampling and symptom onset.
* `Infection ID`: Identifier for each infection, composed of the person ID and the infection number.
* `Infection number`: Specifies the number of a person's infection (first, second, ...).
* `Person ID`: Identifier for each participating individual.
* `Target`: Specifies the molecular source of the measurement (N-antigen or viral RNA).
* `Value`: The measurement (log<sub>10</sub> viral load or log<sub>10</sub> COI).


### rapidTestVariantsSummary.tsv

Data from evaluating rapid test performance on supernatants from cells infected with B.1, B.1617.2 and Omicron BA.1 and BA.2 SARS-CoV-2. The file contains the following columns:

* `Dilution`: The sequential dilution step, with a 1:10 ratio.
* `Experiment no.`: Indicates the repetition number. I.e. the same (diluted) virus 
   stock sample was applied to three different rapid test devices.
* `PFU`: Number of plaque-forming units for a given (diluted) sample of virus stock.
* `Result`: The Ag-RDT testline intensity observed for a given sample, from 0 to 5.
* `RNA`: Viral load measurement for a given sample.
* `Variant`: The SARS-CoV-2 variant of a given sample.



### Nowcast_R_aktuell.csv (accessed May 27, 2025) 

Data containing estimated 7-day R<sub>t</sub> and case numbers in Germany, provided by the Robert Koch-Institut (RKI).
Source: https://github.com/robert-koch-institut/SARS-CoV-2-Nowcasting_und_-R-Schaetzung


## Code

### Directory `agrdt`

Functions for data manipulation (`data.py`), plotting (`plotting.py`), statistical analyses (`regression.py`) and table generation (`tables.py`); and files containing data manipulation parameters and plotting parameter specifications (`dataParams.py`, `plotParams.py`).
 
### Directory `notebooks`

#### Analyses on Charité cohort data
* `agrdtStats.py`: Counts, simple statistics, Ag-RDT specificity computation.
* `comprehensiveModels.py`: Comprehensive statistical analyses, estimating 
  associations between Ag-RDT sensitivity and multiple potential influencing factors.
* `immunization.py`: Analysis on temporal changes in immunization throughout the study.
* `sensitivity.py`: Simple statistical analyses on Ag-RDT sensitivity differences 
  between different infection subgroups and over time.
* `symptoms.py`: Analyses on reported symptoms.
* `viralLoads.py`: Analyses on changes in log<sub>10</sub> viral load over time and 
  between infection subgroups; median and skewness analyses.

#### Other
* `abbottVsRoche.py`: Analysis comparing the performance of the Abbott Panbio COVID-19 Rapid Test Device and the Roche SARS-CoV-2 Rapid Antigen Test.
* `rapidTestVariantExperiments.py`: Analysis of Ag-RDT performance on cell culture supernatants of cells infected with different viral variants.
* `RNA-N-ratios.py`: Analyses on temporal changes in the ratio of N-antigen to viral RNA throughout infection.
 


## License
This project is licensed under the terms of the MIT license.
 


