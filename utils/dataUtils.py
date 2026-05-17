import numpy as np
import pandas as pd
from math import floor, ceil
from datetime import date, datetime
from scipy.spatial.distance import squareform, pdist
from pathlib import Path

from utils.plotParams import (
    DATE_RANGE_DOMINANT_WILDTYPE,
    DATE_RANGE_DOMINANT_ALPHA,
    DATE_RANGE_DOMINANT_DELTA,
    DATE_RANGE_DOMINANT_OMICRON,
    SYMPTOMS_DISPLAY_DICT,
)
from utils.dataParams import (
    ROOT_DIR,
    SYMPTOMS,
    SYMPTOM_DEGREES,
    CONTACTS,
    CONTACT_PROF_KIND,
    CONTACT_PRIV_KIND,
    CONTACT_PROF_FREQ,
    VACC_COLS,
    VACC_DATE_COLS,
    INFECTION_DATE_COLS,
    DAYS_PER_MONTH,
)

DATA_DIR = Path(ROOT_DIR / "data")


def returnDataDir():
    return DATA_DIR


def IQR(column):
    q25, q75 = IQRQuartiles(column)
    return q75 - q25


def IQRQuartiles(column):
    return tuple(column.quantile([0.25, 0.75]))


def addJitterCol(df, origCol, jitterCol, sd=0.01):
    """
    Add a column "jitterCol" with random jittering applied to the values in "origCol".
    @param df: A C{pandas DataFrame}.
    @param origCol: The name of the column whose values you want to apply jittering to.
    @param jitterCol: The name of the column to be created, containing the values that
    have been jittered.
    """
    stdev = sd * (df[origCol].max() - df[origCol].min())
    df.loc[:, jitterCol] = df[origCol] + np.random.randn(len(df[origCol])) * stdev


def getPCRsOutsideVariantPrevalentRanges(df):
    """
    Return data with PCRs performed outside time periods where one viral lineage
    dominated (>90% of smaples in Berlin).
    @param df: A C{pd.DataFrame} with rapid test, PCR and metadata.
    @return: A C{pd.DataFrame} with data from tests performed outside time periods where
    one viral lineage was predominant.
    """
    return df[
        (
            (df.pcrDate >= DATE_RANGE_DOMINANT_WILDTYPE[1])
            & (df.pcrDate <= DATE_RANGE_DOMINANT_ALPHA[0])
        )
        | (
            (df.pcrDate >= DATE_RANGE_DOMINANT_ALPHA[1])
            & (df.pcrDate <= DATE_RANGE_DOMINANT_DELTA[0])
        )
        | (
            (df.pcrDate >= DATE_RANGE_DOMINANT_DELTA[1])
            & (df.pcrDate <= DATE_RANGE_DOMINANT_OMICRON[0])
        )
    ]


def returnRtData():
    """
    @return: A C{pd.DataFrame} containing German 7-day Rt values from the RKI.
    """
    data = pd.read_csv(DATA_DIR / "Nowcast_R_aktuell.csv", parse_dates=["Datum"])
    data = data.rename(
        columns={
            "Datum": "date",
            "PS_7_Tage_R_Wert": "rt",
            "UG_PI_7_Tage_R_Wert": "uirt",
            "OG_PI_7_Tage_R_Wert": "lirt",
            "PS_COVID_Faelle": "newCases",
        }
    )

    data["dateCopy"] = data.date
    data.set_index("dateCopy", inplace=True)
    center = True
    data["dateDate"] = data.date.apply(lambda x: x.date())

    data["rtRolling"] = data["rt"].rolling(window="7D", center=center).mean()
    data["rtRolling2"] = data["rt"].rolling(window="14D", center=center).mean()
    data["rtRolling3"] = data["rt"].rolling(window="21D", center=center).mean()
    data["rtRolling4"] = data["rt"].rolling(window="28D", center=center).mean()

    data["rtRolling5"] = data["rt"].rolling(window="35D", center=center).mean()
    data["newCasesRolling"] = (
        data["newCases"].rolling(window="7D", center=center).mean()
    )
    data["newCasesRolling2"] = (
        data["newCases"].rolling(window="14D", center=center).mean()
    )

    data["week"] = pd.to_datetime(data["date"]).dt.to_period(freq="W-MON")
    mediansRt = data.groupby("week").rt.median()
    data = data.merge(
        mediansRt, left_on="week", right_index=True, suffixes=("", "_median")
    )

    return data


def returnIncidenceData():
    """
    @return: A C{pd.DataFrame} containing German incidence values from the RKI.
    """
    data = pd.read_csv(
        DATA_DIR / "fallzahlen_und_indikatoren.csv", sep=";", parse_dates=["Datum"]
    )
    data = data.rename(
        columns={
            "Datum": "date",
            "rel Veraenderung der 7-Tage-Inzidenz": "incidenceChange",
            "7-Tage-Inzidenz": "incidence",
        }
    )
    data["incidenceChange"] = data.incidenceChange.str.replace(
        ",", ".", regex=False
    ).astype(float)
    data["incidenceChangeRaw"] = data["incidenceChange"] / 100
    data["rtApprox"] = 1 + data["incidenceChangeRaw"]
    data["incidence"] = data.incidence.str.replace(",", ".", regex=False).astype(float)
    standardize(data, "incidenceChange", "zIncidenceChange")
    standardize(data, "incidence", "zIncidence")

    data = data[
        (data.date >= datetime(2020, 12, 1)) & (data.date < datetime(2022, 2, 12))
    ]

    center = True
    data["dateCopy"] = data.date
    data.set_index("dateCopy", inplace=True)
    data["rtRolling"] = data["rtApprox"].rolling(window="7D").mean()
    data["rciRolling"] = (
        data["incidenceChange"].rolling(window="7D", center=center).mean()
    )
    data["rciRolling2"] = (
        data["incidenceChange"].rolling(window="14D", center=center).mean()
    )
    data["rciRolling3"] = (
        data["incidenceChange"].rolling(window="21D", center=center).mean()
    )
    data["rciRolling4"] = (
        data["incidenceChange"].rolling(window="28D", center=center).mean()
    )
    data["rciRolling5"] = (
        data["incidenceChange"].rolling(window="35D", center=center).mean()
    )
    data["iRolling"] = data["incidence"].rolling(window="7D", center=center).mean()
    data["iRolling2"] = data["incidence"].rolling(window="14D", center=center).mean()
    data["iRolling3"] = data["incidence"].rolling(window="21D", center=center).mean()
    data["iRolling4"] = data["incidence"].rolling(window="28D", center=center).mean()
    data["iRolling5"] = data["incidence"].rolling(window="35D", center=center).mean()

    return data


def returnDeltaOmicronTransitionPCRs(df):
    """
    Return A data frame with data points from the time of delta/omicron transition.
    """
    return df[
        (df.pcrDate > DATE_RANGE_DOMINANT_DELTA[1])
        & (df.pcrDate < DATE_RANGE_DOMINANT_OMICRON[0])
    ].copy()


def removeDeltaOmicronPCRs(df):
    """
    Remove data points from the time when delta/omicron typing PCRs were done
    (Dec 2021 - Jan 2022) so that they don't give a false impression of viral
    load levels for those variants (typing PCRs in this period were done for
    high viral load level samples only).
    """
    return df[
        (df.pcrDate <= DATE_RANGE_DOMINANT_DELTA[1])
        | (df.pcrDate >= DATE_RANGE_DOMINANT_OMICRON[0])
    ].copy()


def returnWtAlphaTransitionPCRs(df):
    """
    Return data points from the time of wildtype/alpha transition.
    """
    return df[
        (df.pcrDate > DATE_RANGE_DOMINANT_WILDTYPE[1])
        & (df.pcrDate < DATE_RANGE_DOMINANT_ALPHA[0])
    ].copy()


def removeWtAlphaPCRs(df):
    """
    Remove data points from the time when wildtype/alpha typing PCRs were done
    (Feb 2021 - March 2021) so that they don't give a false impression of viral
    load levels for those variants (typing PCRs in this period were done for
    high viral load level samples only).
    """
    return df[
        (df.pcrDate <= DATE_RANGE_DOMINANT_WILDTYPE[1])
        | (df.pcrDate >= DATE_RANGE_DOMINANT_ALPHA[0])
    ].copy()


def removeAlphaDeltaPCRs(df):
    """
    Remove data points from the time when alpha/delta typing PCRs were done
    (May 2021 - July 2021) so that they don't give a false impression of viral
    load levels for those variants (typing PCRs in this period were done for
    high viral load level samples only).
    """
    return df[
        (df.pcrDate <= DATE_RANGE_DOMINANT_ALPHA[1])
        | (df.pcrDate >= DATE_RANGE_DOMINANT_DELTA[0])
    ].copy()


def removeAllUnclearVariantOrTypingPCRs(df):
    """
    Remove tests from C{df} that were performed during viral variant transition times.
    @param df: A C{pd.DataFrame} with Agrdt, PCR and metadata.
    @return: A C{pd.DataFrame} with data from those tests removed that happened during
    a time of variant transition.
    """
    df = removeWtAlphaPCRs(df)
    df = removeAlphaDeltaPCRs(df)
    df = removeDeltaOmicronPCRs(df)
    return df


def removeReleaseTesting(df):
    """
    Remove data points from release testing.
    """
    return df[df.reasonPres != "releaseTesting"].copy()


def standardize(df, origCol, standCol):
    """
    @param df: A C{pandas DataFrame}.
    @param origCol: The name of the column whose values you want to standardize
        (turn into z-scores).
    @param standCol: The name of the column to be created, containing the standardized
    values.
    """
    df[standCol] = (df[origCol] - df[origCol].mean()) / df[origCol].std()


def dataFrameNoRecovered(df):
    """
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    @return: A C{pd.DataFrame} containing no PCRs of people who we know have had a
    previous SARS-CoV-2 infection at the time of testing.
    """
    return df[df.recovered != 1].copy()


def dataFrameFemaleMale(df):
    """
    @param df: A C{pd.DataFrame} with Agrdt, PCR and metadata.
    @return: A C{pd.DataFrame} with datapoints for which the individual's gender is
    either male or female.
    """
    return df[df.gender.isin(("F", "M"))].copy()


def combineToStrVariantSymp(row):
    """
    @param row: A C{pd.Series} with variant and symptom information.
    @return: A string containing the variant and symptom status.
    """
    try:
        return f"{(row.variant)}, {int(row.symptoms)}"
    except ValueError:
        return np.nan


def combineToStrImmunSymp(row):
    """
    @param row: A C{pd.Series} with immunization and symptom information.
    @return: A string containing immunization and symptom status.
    """
    try:
        return f"{int(row.immun2YN)}, {int(row.symptoms)}"
    except ValueError:
        return np.nan


def createDataFramesFigures(df):
    """
    @param df: A C{pd.DataFrame} with Agrdt, PCR and metadata.
    @return: A C{tuple} with C{pd.DataFrames} used for different figures.
    """
    df = removeInvalidAgrdts(df)
    maskRecovered = df.recovered.isna()
    df.loc[maskRecovered, "recovered"] = df.nPrevInfections >= 1
    # Date frame containing only PCRs with corresponding Ag-RDT results.
    dfAgrdtAll = df.dropna(subset=["agrdt"]).copy()
    dfAgrdt = dfAgrdtAll[(df.valid == 1) & (df.dry == 1)].copy()

    # Data frame containing only PCRs with corresponding Ag-RDT results and containing
    # only the first test within an infection.
    dfAgrdtIndInf = dataFrameIndependentInfection(dfAgrdt)

    dfAllInd = dataFrameIndependent(df.dropna(subset=("pcrPositive", "agrdt")))

    # Data frame containing only positive PCRs.
    dfAllPos = dataFramePCRpos(df)

    # Data frames containing only rapid test and PCR data with positive PCRs.
    dfPos = dataFrameAgrdt(dfAllPos)
    dfSymp = dataFrameSymptoms(dfPos)
    dfAsymp = dataFrameNoSymptoms(dfPos)

    # Data frames containing all data with positive PCRs.
    dfAllPosNoRelease = removeReleaseTesting(dfAllPos)

    # Use only the data points corresponding to the first positive PCR of a person
    # (within an infection). Exclude people who presented for release testing.
    dfAllFirstPosPcrsNoRelease = dataFrameFirstPosPcr(dfAllPosNoRelease)
    dfAllFirstPosPcrsSympNoRelease = dataFrameSymptoms(dfAllFirstPosPcrsNoRelease)
    dfAllFirstPosPcrsAsympNoRelease = dataFrameNoSymptoms(dfAllFirstPosPcrsNoRelease)

    # Figure 1
    # Show only symptomatic people, remove tests from release testing, take the first
    # positive PCR of an infection.
    dfFigure1 = dataFrameFirstPosPcr(removeReleaseTesting(dfSymp))
    dfFigure1Asymp = dataFrameFirstPosPcr(removeReleaseTesting(dfAsymp))
    # Only consider PCR tests after Nov 2020 as there are very few data points before
    # and we don't want to display
    # them when stratifying by sampling month (as we start in December).
    dfFigure1 = dfFigure1[dfFigure1.pcrDate >= date(2020, 12, 1)].copy()
    dfFigure1 = dfFigure1.dropna(subset=["samplingMonth2"]).copy()
    dfFigure1Asymp = dfFigure1Asymp[dfFigure1Asymp.pcrDate >= date(2020, 12, 1)].copy()
    dfFigure1Asymp = dfFigure1Asymp.dropna(subset=["samplingMonth2"]).copy()
    # Remove single datapoint in June-July bin.
    # Identify months with only one datapoint.
    dfFigure1_A = (
        dfFigure1.groupby("samplingMonth2").filter(lambda x: len(x) >= 2).copy()
    )
    dfFigure1_B = dfFigure1_A.copy()
    dfFigure1_A_Asymp = (
        dfFigure1Asymp.groupby("samplingMonth2").filter(lambda x: len(x) >= 2).copy()
    )
    dfFigure1_B_Asymp = dfFigure1_A_Asymp.copy()
    # Remove the five datapoints with 4 immunisations (would be hardly visible in the
    # plot).
    dfFigure1_C = dfFigure1_A[dfFigure1_A.immunN < 4].copy()
    dfFigure1_C_Asymp = dfFigure1_A_Asymp[dfFigure1_A_Asymp.immunN < 4].copy()

    # Figure 2
    # PCRs typed during the delta-omicron transition period have higher viral loads
    # (because that's how they are selected for typing), so viral loads from these
    # samples are not representative for VOCs delta/omicron.
    # Also, release testing (done at the end of an infection) is not representative for
    # viral loads --> don't use.
    dfFigure2 = removeAllUnclearVariantOrTypingPCRs(
        dfAllFirstPosPcrsSympNoRelease
    ).copy()
    dfFigure2 = dfFigure2[dfFigure2.pcrDate >= date(2020, 12, 1)].copy()

    dfFigure2_A = dfFigure2.dropna(subset=["daysPostOnset", "variant"])
    dfFigure2_B = dfFigure2.dropna(subset=["variant", "agrdtYN"])

    dfFigure2Asymp = removeAllUnclearVariantOrTypingPCRs(
        dfAllFirstPosPcrsAsympNoRelease
    ).copy()
    dfFigure2Asymp = dfFigure2Asymp[dfFigure2Asymp.pcrDate >= date(2020, 12, 1)].copy()
    dfFigure2_BAsymp = dfFigure2Asymp.dropna(subset=["variant", "agrdtYN"])

    # Figure 3
    # Data frame for Figure 3
    dfFigure3 = dataFrameFirstPosPcr(removeReleaseTesting(dfPos))
    # Create a new column indicating whether people had symptoms and whether they were immunized.
    dfFigure3["immun2YN:symptoms"] = dfFigure3.apply(combineToStrImmunSymp, axis=1)
    dfFigureA2 = dataFrameFirstPosPcr(
        removeReleaseTesting(removeAllUnclearVariantOrTypingPCRs(dfPos))
    )
    # Create a new column indicating whether people had symptoms and which variant they had.
    dfFigureA2["variant:symptoms"] = dfFigureA2.apply(combineToStrVariantSymp, axis=1)

    return (
        dfAgrdt,
        dfAgrdtAll,
        dfAgrdtIndInf,
        dfAllFirstPosPcrsNoRelease,
        dfPos,
        dfAllInd,
        dfFigure1,
        dfFigure1Asymp,
        dfFigure1_A,
        dfFigure1_A_Asymp,
        dfFigure1_B,
        dfFigure1_B_Asymp,
        dfFigure1_C,
        dfFigure1_C_Asymp,
        dfFigure2,
        dfFigure2_A,
        dfFigure2_B,
        dfFigure2_BAsymp,
        dfFigure3,
        dfFigureA2,
    )


def symptomColNames():
    """
    @return: The column names describing the different symptoms (yes/no).
    """
    return SYMPTOMS


def symptomDegrees():
    """
    @return: The degree of "illness", i.e. asymptomatic, mildy symptomatic
    and symptomatic.
    """
    return SYMPTOM_DEGREES


# The date the Roche rapid test started being used (as opposed to the Abbott
# test).
DATE_ROCHE = datetime(2021, 12, 18)

IMMUNIZATION_COLS = (
    VACC_DATE_COLS
    + INFECTION_DATE_COLS
    + (
        "nPrevInfections",
        "datePrevInfection",
        "nInfectionPreOct18",
        "vaccYN",  # "vaccN",
        "vaccNatLeast",
        "immunYN",
        "immun2YN",
        "immunN",
        "immunNatLeast",
        "vaccTime",
        "vaccTime2",
        "vaccTimeAtLeast",
        "recovered",
        "recoveredTime",
        "recoveredRecently",
        "vaccName1",
        "vaccName2",
        "vaccName3",
        "vaccWhat1",
        "vaccWhat2",
        "vaccWhat3",
    )
    + tuple([f"{col}Orig" for col in VACC_DATE_COLS])
    + tuple([f"{col}Orig" for col in INFECTION_DATE_COLS])
)


def removeInvalidAgrdts(df):
    """
    @param df: A C{pd.DataFrame} with Agrdt, PCR and metadata.
    @return: A C{pd.DataFrame} with data associated with valid Ag-RDT testing (test
    device hadn't expired; contained drying material).
    """
    return df[(df.agrdtYN == 0) | ((df.valid == 1) & (df.dry == 1))].copy()


def dataFrameAgrdt(df):
    """
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    @return: A C{pd.DataFrame} containing only tests with an Ag-RDT result available
    that are valid.
    """
    return df.dropna(subset=("agrdt"))


def dataFrameFirstPosPcr(df):
    """
    Return a data frame containing only data points that correspond to the first
    positive PCR of a person within an infection.
    """
    return df[df.isFirstPosPcr == 1].copy()


def dataFrameIndependent(df):
    """
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    @return: A dataframe that only has one test per person (the first) so that test data
    points are independent of each other.
    """

    # Make sure that values are sorted by pcr date so that only the first test will be
    # kept for each person to keep things consistent.
    df = df.sort_values("pcrDate", ignore_index=True)
    return df.drop_duplicates(subset=["personHash"], keep="first")


def dataFrameIndependentInfection(df):
    """
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    @return: A dataframe that only has one test per infection (the first) so that test
    data points are independent of each other.
    """

    # Make sure that values are sorted by pcr date so that only the first test will be
    # kept for each person to keep things consistent.
    df = df.sort_values("pcrDate", ignore_index=True)
    df = df.dropna(subset="infectionKey")
    return df.drop_duplicates(subset=["infectionKey"], keep="first")


def dataFrameSymptoms(df):
    """
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    @return: A dataframe that only has tests in them where the respective person was
    symptomatic.
    """
    return df[df.symptoms == 1].copy()


def dataFrameEmployees(df):
    """
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    @return: A dataframe that only has tests in it where the respective
    person was an employee (not a student).
    """
    return df[df.student != 1].copy()


def dataFrameReasonPresAsymp(df):
    """
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    @return: A dataframe that only has tests in them where the respective
    person most likely showed up because they were asymptomatic.
    """
    return df[df.reasonPres.isin(["outbreak", "no category"])].copy()


def dataFrameReasonPresSymp(df):
    """
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    @return: A dataframe that only has tests in them where the respective
    person showed up because they were symptomatic.
    """
    return df[df.reasonPres == "symptoms"].copy()


def dataFrameSymptomsNoNan(df):
    """
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    @return: A dataframe that only has tests in them where the respective
    person was symptomatic or asymptomatic (does not contain NaN values).
    """

    newDf = df[df["symptoms"].notnull()].copy()
    newDf.symptoms = newDf.symptoms.astype(int)
    return newDf


def dataFrameIndependentSymptomsNoNan(df):
    """
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    @return: A dataframe that only has independent tests in them where the
    respective person was symptomatic or asymptomatic (does not contain NaN
    values).
    """
    return dataFrameIndependent(dataFrameSymptomsNoNan(df))


def dataFrameSymptomsData(df):
    """
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    @return: A dataframe that only has tests in them with detailed information
    available on the symptoms the respective person was having.
    """
    return df[df.surveyData == 1].copy()


def dataFrameIndependentSymptomsData(df):
    """
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    @return: A dataframe that only has independent tests in them with detailed
    information available on the symptoms the respective person was having.
    """
    return dataFrameIndependent(dataFrameSymptomsData(df))


def dataFrameIndependentSymptomsSymptomsData(df):
    """
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    @return: A dataframe that only has independent tests in them with detailed
    information available on the symptoms the respective person was having.
    Excluding asymptomatic people that (wrongly) filled out the questionnaire.
    """
    return dataFrameIndependent(dataFrameSymptoms(dataFrameSymptomsData(df)))


def dataFrameOutbrSymp(df):
    """
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    @return: A dataframe that only has datapoints with the reason for
    presentation being either 'symptoms', 'outbreak'.
    """
    return df[df.reasonPres.isin(("outbreak", "symptoms"))].copy()


def dataFrameOutbrSymp2(df):
    """
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    @return: A dataframe that only has datapoints with either the reason for
    presentation being 'outbreak' or with a symptomsBool value of 1 (this
    will include people that have a reason for presentation other than
    'symptoms' but still specified symptoms in the survey data file).
    """
    return df[(df.reasonPres == "outbreak") | (df.symptoms == 1)].copy()


def dataFrameReasonPresOther(df):
    """
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    @return: A dataframe that only has datapoints with the reason for
    presentation not being symptoms or outbreak, i.e. either
    freetesting, medical student or no category.
    """
    return df[~df.reasonPres.isin(("outbreak", "symptoms"))].copy()


def dataFrameNoStudents(df):
    """
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    @return: A datframe containing only Charité employee data (and no data
    from the student testing).
    """
    return df[df.student == 0].copy()


def dataFrameNoSymptoms(df):
    """
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    @return: A dataframe that only has datapoints with people who were
        asymptomatic at the time of testing.
    """
    return df[df.symptoms == 0].copy()


def dataFramePCRpos(df):
    """
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    @return: A C{pd.DataFrame} containing only data with positive PCR result.
    """
    return df[df.pcrPositive == 1].copy()


def dataFrameReleaseTesting(df):
    """
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    @return: A C{pd.DataFrame} containing only PCRs of people who showed up
    at the CUSMA for release testing.
    """
    return df[df.reasonPres == "freetesting"].copy()


def symptomsDisplayDict():
    """
    @return: A dictionary mapping the names of the columns containing the
    different symptoms to the symptom names that are supposed to be
    displayed.
    """
    return SYMPTOMS_DISPLAY_DICT


def contactColNames():
    """
    @return: The column names of all contact relevant data.
    """
    return CONTACTS


def contactProfFreq():
    return CONTACT_PROF_FREQ


def contactProfKind():
    return CONTACT_PROF_KIND


def contactPrivKind():
    return CONTACT_PRIV_KIND


def vaccColNames():
    return VACC_COLS


def deStandardize(df, col, x):
    """
    @param df: A C{pandas DataFrame} containing C{col} as column.
    @param col: A C{str} specifying the column to which standardization was
    applied to, resulting in C{x}.
    @param x: The C{np.array} or C{pd.Series} of values to apply the
    destandardization to.
    @return: A destandardized version of C{x}.
    """

    st = df[col].std()
    mean = df[col].mean()
    return x * st + mean


# Functions to add columns to the DataFrame.


def _addLongestSymptomTime(df):
    """
    Add a column with the longest duration of symptoms (irrespective of the symptom
    type) until the day of the test. I.e. this specifies the time that has passed
    since the first symptoms showed up.
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    """

    def maxDays(row):
        return row[symptomDaysCols].max()

    symptomDaysCols = [f"{symptom}Days" for symptom in SYMPTOMS]

    df["daysPostOnset"] = df.apply(maxDays, axis=1)


def _addTimeSinceLongestSymptom(df):
    symptomDaysCols = [f"{symptom}Days" for symptom in SYMPTOMS]

    for symptomDaysCol in symptomDaysCols:
        newCol = f"{symptomDaysCol}SinceLongestSymptom"
        df[newCol] = df.daysPostOnset - df[symptomDaysCol]
        assert all(df[newCol].dropna() >= 0)


def _addLongestSymptom(df):
    def longestSymp(row):
        if row[symptomDaysCols].sum() == 0:
            return np.nan
        return row[symptomDaysCols].astype(float).fillna(-1).idxmax()

    symptomDaysCols = [f"{symptom}Days" for symptom in SYMPTOMS]

    df["longestSymptom"] = df.apply(longestSymp, axis=1)


def addSympDegree(df):
    """
    Add symptom degrees (0 to 2) to each datapoint where 0 is asymptomatic,
    1 is mildly symptomatic (i.e. not "feeling ill") and 2 is symptomatic.
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    """
    df["sympDegree"] = np.nan
    df.loc[df.symptoms == 0, "sympDegree"] = 0
    df.loc[((df.symptoms == 1) & (df.ill == 0)), "sympDegree"] = 1
    df.loc[(df.ill == 1), "sympDegree"] = 2


def addDaysPostFirstPosPcrCategories(df):
    # Create two categories for time since first positive PCR:
    # 0 days, >=7 days.
    df["binDaysPostFirstPosPcr"] = np.nan
    maskFirstPosPcr = df.daysSinceFirstPosPcr == 0
    maskMoreThanSixDaysPostFirstPosPcr = df.daysSinceFirstPosPcr >= 7
    df.loc[maskFirstPosPcr, "binDaysPostFirstPosPcr"] = 0
    df.loc[maskMoreThanSixDaysPostFirstPosPcr, "binDaysPostFirstPosPcr"] = 1


def addAgrdt3(df):
    """
    For a DataFrame with NaN values in the C{'agrdt'} column, replace the NaNs
    with 2 (i.e. add a new category for missing antigen test result).
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    """
    df["agrdt3"] = df.agrdt.fillna(2)


def removeDeltaOmicronTyping(df):
    """
    Remove data points from the time when delta/omicron typing PCRs were done
    (Dec 2021 - Jan 2022) so that they don't give a false impression of viral
    load levels for those variants (typing PCRs in this period were done for
    high viral load level samples only).
    """
    return df[
        (df.pcrDate < datetime(2021, 12, 6))
        | (df.pcrDate >= DATE_RANGE_DOMINANT_OMICRON[0])
    ].copy()


def removeStudents(df):
    """
    Remove data from student testing.
    """
    return df[df.reasonPres != "medical student"].copy()


def _addRecoveredRecently(df):
    """
    Add a column specifying whether someone has recovered recently (< 3
    months ago) from a SARS-CoV-2 infection. People in category 0 did not
    have a previous infection.
    @param df: The pandas dataframe returned by calling dataFrameCusma.
    """
    df["recoveredRecently"] = np.nan
    mask = df.recoveredTime == 0
    mask2 = df.recovered == 0
    df.loc[mask, "recoveredRecently"] = 1
    df.loc[mask2, "recoveredRecently"] = 0


def _samplingDateDistance(pcr1, pcr2):
    """
    @param pcr1: A C{PCR} instance.
    @param pcr2: A C{PCR} instance.
    @return: The number of days between pcr1 and pcr2.
    """
    return (pcr2.samplingDate - pcr1.samplingDate).days


def _addSamplingWeek(dfCusma):
    dfCusma["samplingWeek"] = dfCusma.samplingDate.dt.strftime("%Y-W%-V")
    # Handle the first week of 2021 which gets categorized as
    # week 53 of 2020.
    firstWeekMapping = {"2021-W53": "2021-W1"}
    dfCusma.samplingWeek = dfCusma.samplingWeek.replace(firstWeekMapping)


def addSamplingMonthCusma(df, dateCol):
    df["samplingMonth"] = df[dateCol].dt.strftime("%Y-%m")
    # We want to ignore samples from before December 2020 because we'd
    # end up with too little data for that time.
    maskBeforeDec2020 = df[dateCol] < datetime(2020, 12, 1)
    df.loc[maskBeforeDec2020, "samplingMonth"] = np.nan


def addSamplingMonthCusma2(df):
    "Add a column with two-month bins (Dec. 2020 - Jan. 2021...)."
    errorMssg = (
        "Passed data frame must have a column specifying the "
        "sampling month of each test."
    )
    assert "samplingMonth" in df, errorMssg
    # Create a column specifying the sampling month by integers.
    samplingMonthSorted = df.samplingMonth.sort_values().unique()
    samplingMonthToNumeric = {
        samplingMonth: i for i, samplingMonth in enumerate(samplingMonthSorted)
    }
    df["samplingMonthNumeric"] = df.samplingMonth.dropna().map(samplingMonthToNumeric)
    # Create mapping from numeric sampling month to numeric two-month bins.
    mapping = dict(
        zip(
            samplingMonthToNumeric.values(),
            np.repeat(list(samplingMonthToNumeric.values()), 2),
            strict=True,
        )
    )
    df["samplingMonth2"] = df.samplingMonthNumeric.map(mapping)
    df.drop(columns=["samplingMonthNumeric"], inplace=True)


def addSamplingMonthCusma3(df):
    # Same as samplingMonth2 (see above) but merging December '21, January '22 and
    # February '22.
    df["samplingMonth3"] = df["samplingMonth2"].copy()
    # Last three months.
    lastThreeMonths = sorted(
        [month for month in df["samplingMonth2"].unique() if not np.isnan(month)]
    )[-2:]
    lastNewMonth = lastThreeMonths[0]
    df["samplingMonth3"] = df["samplingMonth3"].replace(
        {month2: lastNewMonth for month2 in lastThreeMonths}
    )


def _addDaysPostSymptomOnsetBins(df):
    # Create three bins for time since symptom onset for logistic regression.
    bins1 = [0, 2, 7]
    df["binDaysPostOnset"] = pd.cut(df.daysPostOnset, bins1, include_lowest=True)

    bins2 = list(range(floor(df.daysPostOnset.min()), ceil(df.daysPostOnset.max()) + 1))
    df["binDaysPostOnset2"] = pd.cut(df.daysPostOnset, bins2, include_lowest=True)

    bins3 = [0, 2, 6, ceil(df.daysPostOnset.max())]
    df["binDaysPostOnset3"] = pd.cut(df.daysPostOnset, bins3, include_lowest=True)


def addSymptomsURT(df):
    # Add a new column 'illURT' indicating whether symptoms in the upper
    # respiratory tract (URT) were present at the time of testing.
    df["illURT"] = (
        (df["runnyNose"] == 1)
        | (df["soreThroat"] == 1)
        | (df["cough"] == 1)
        | (df["noSmell"] == 1)
        | (df["noTaste"] == 1)
    )

    return df


def _addTestDevice(dfCusma):
    """
    Add a column specifying which test device (Abbott or Roche) was used. Roche
    started being used on Dec 18th, 2021.
    """
    maskRoche = dfCusma.samplingDate >= DATE_ROCHE
    dfCusma["testDevice"] = "abbott"
    dfCusma.loc[maskRoche, "testDevice"] = "roche"


def returnRecoveredTime(days, timeBins):
    """
    Return a code for time passed since a positive PCR.
    Codes for time since last infection are:
    For 'recoveredTime':
    0: < 3 months
    1: 3-6 months
    2: 6-9 months
    3: 9-12 months
    4: 12-18 months
    5: 18-24 months
    For 'recoveredTime2':
    0: < 2 months
    1: 2-6 months
    2: 6-9 months
    3: 9-12 months
    4: 12-18 months
    5: 18-24 months
    @param days: A C{int} specifying the number of days that have passed since
        a last positive PCR.
    @return: A C{int} code specifying the time since the last positive PCR.
    """
    for i, month in enumerate(timeBins):
        if days < month * DAYS_PER_MONTH:
            return i


def jaccardSim(row1, row2):
    """

    @param row1: A C{pd.Series} containing 0s and 1s.
    @param row2: A C{pd.Series} containing 0s and 1s.
    @return: The jaccard similarity. Computation following:
    https://www.learndatasci.com/glossary/jaccard-similarity/.
    """
    nBothOne = ((row1 == 1) & (row2 == 1)).sum()
    numerator = nBothOne
    denominator = nBothOne + (row1 != row2).sum()
    if denominator:
        return numerator / denominator
    else:
        return 0


def makeSimilarityDf(df):
    """
    @param df:A binary DataFrame where columns contain information on the presence of
    a specific symptom.
    @return: A symmetric matrix of Jaccard similarities between all pairs of columns,
    with diagonal values set to 1.
    """
    jacSim = squareform(pdist(df.T, lambda row1, row2: jaccardSim(row1, row2)))
    jacSim = pd.DataFrame(jacSim, index=df.columns, columns=df.columns)

    # Identical samples are perfectly correlated
    for idx in jacSim.index:
        for col in jacSim.columns:
            if idx == col:
                jacSim.loc[idx, col] = 1

    return jacSim
