import decimal
import re
from collections import defaultdict

from agrdt.dataParams import ROOT_DIR, SYMPTOMS_DISPLAY_DICT
from agrdt.data import dataFramePCRpos, IQRQuartiles
from agrdt.plotParams import (
    getAbbrvsDict,
    getLabels,
    getLabelsLatex,
    PARAM_NAME_MAPPINGS,
    PARAM_NAME_MAPPINGS_TABLE_LATEX,
)

# Contains functions for creating .tsv and latex tables with information on the
# composition of the study population, symptoms and posterior estimates from statistical
# analyses.

OUTPUT_TABLES = ROOT_DIR / "output" / "tables"

OUTPUT_TABLE_EMPLOYEES = OUTPUT_TABLES / "employeeTable.tsv"
OUTPUT_TABLE_EMPLOYEES2 = OUTPUT_TABLES / "employeeTable2.tsv"
OUTPUT_TABLE_EMPLOYEES_LATEX = OUTPUT_TABLES / "employeeTable.txt"
OUTPUT_TABLE_EMPLOYEES_LATEX2 = OUTPUT_TABLES / "employeeTable2.txt"
OUTPUT_TABLE_SYMPTOMS = OUTPUT_TABLES / "symptomsTable.tsv"
OUTPUT_TABLE_SYMPTOMS2 = OUTPUT_TABLES / "symptomsTable2.tsv"
OUTPUT_TABLE_SYMPTOMS_LATEX = OUTPUT_TABLES / "symptomsTable2.txt"
OUTPUT_TABLE_SYMPTOMS_LATEX2 = OUTPUT_TABLES / "symptomsTable3.txt"
OUTPUT_TABLE_SYMPTOM_SEVERITY_LATEX = OUTPUT_TABLES / "symptomsTableSeverity.txt"

TABLE_DIR = OUTPUT_TABLES
TABLE_DIR_IDATA = TABLE_DIR / "iData"

LABEL = getLabels()
LABEL_LATEX = getLabelsLatex()

TABLE_SUFFIX = r"""
\end{tabularx}
\end{table}
"""

EMPLOYEE_TABLE_LATEX_PREFIX = r"""
\begin{table}[h!]
\centering
\small
\caption{\textbf{Characteristics of the study population and SARS-CoV-2 diagnostic test results, considering samples tested by both Ag-RDT and RT-PCR}. Counts (n) and percentages (\%) are provided for different gender categories (top) and indicate the number and fraction of positive test results using RT-PCR or Ag-RDT (bottom). Median ages and corresponding interquartile ranges (IQRs) of individuals at the time of testing are shown (center). Values  are provided according to symptom status ("Symptoms", "No symptoms", "Unknown symptom status") and for all samples ("Total").}
\label{tab:cusma_pop_characteristics}
\begin{tabularx}{\textwidth}{|>{\centering\arraybackslash}p{1.6cm}|>{\RaggedLeft}p{1.5cm}|R|R|R|R|}
\hline
\rowcolor{lightgray}
"""

EMPLOYEE_TABLE_LATEX_SUFFIX = TABLE_SUFFIX

SYMPTOMS_TABLE_LATEX_PREFIX = r"""
\begin{table}[h!]
\centering
\small
\caption{\textbf{Counts and fractions of reported symptoms at the time of testing.} Column 'All tests' reports values considering all tests whereas column 'PCR positive tests reports values for PCR positive tests. Symptoms are ordered by fraction (descending) of PCR positive tests. Note that 'runny nose' was included in the questionnaire only at the beginning of 2020 and therefore is associated with a lower count despite a higher fraction compared to the other symptoms.}
\label{tab:cusma_symptoms}
\begin{tabularx}{\textwidth}{|R|R|R|R|R|}
\cline{1-5}
\noalign{\vskip 0.3pt}
\rowcolor{lightgray}
"""

SYMPTOMS_TABLE_LATEX_SUFFIX = TABLE_SUFFIX

SYMPTOMS_SEVERITY_TABLE_LATEX_PREFIX = r"""
\begin{table}[h!]
\centering
\small
\caption{\textbf{Number of experienced symptoms and perceived severity of illness.} 
The earliest data point associated with a positive RT-PCR result and symptom specific information from each infection was used. The median number of symptoms present, and the severity of symptoms at the time of testing along with IQRs are presented, stratified by gender, age, SARS-CoV-2 variant and immunization status. The number of infections in each subgroup is indicated in the 'Count' column. Severity of illness was reported on an integer scale from 1 to 5, indicating increasing severity.}
\label{tab:cusma_symptoms_severity}
\begin{tabularx}{\textwidth}{|>{\hsize=1.3\hsize}C|R|>{\hsize=0.7\hsize}R|R|R|R|R|}
\hline
\rowcolor{lightgray}
"""

SYMPTOMS_TABLE_LATEX_PREFIX2 = r"""
\begin{table}[h!]
\centering
\small
\caption{\textbf{Reported symptoms and durations at the time of testing}. Days post earliest onset (DPEO) specifies the duration (in days) between the onset of the earliest symptom and the onset of the considered symptom.}
\label{tab:cusma_symptoms}
\begin{tabularx}{\textwidth}{|R|R|R|R|R|R|R|R|R|}
\hline
\rowcolor{lightgray}
"""

SYMPTOMS_TABLE_LATEX_SUFFIX2 = TABLE_SUFFIX

SYMPTOMS_SEVERITY_TABLE_LATEX_SUFFIX = TABLE_SUFFIX


# Note: for now we are specifying the table environment manually in latex
IDATA_TABLE_PREFIX1 = r"""
"""

IDATA_TABLE_PREFIX2 = r"""
"""

IDATA_TABLE_PREFIX3 = r"""
"""

IDATA_TABLE_PREFIX4 = r"""
"""

IDATA_TABLE_PREFIX5 = r"""
"""

IDATA_TABLE_PREFIX6 = r"""
"""

IDATA_TABLE_PREFIX = r"""

"""

IDATA_TABLE_SUFFIX = ""


def roundHalfUp(value, decimals=2):
    """

    @param value: The float value to be rounded.
    @param decimals: The decimal places to round to.
    @return: The rounded value (using commercial rounding).
    """
    with decimal.localcontext() as ctx:
        d = decimal.Decimal(value)
        ctx.rounding = decimal.ROUND_HALF_UP
        return float(round(d, decimals))


def mapRoundHalfUp(values, decimals=2):
    """

    @param values: The float values to be rounded.
    @param decimals: An C{Iterable} with the decimal places to round to.
    @return: The rounded values (using commercial rounding).
    """
    return [roundHalfUp(value, decimals=decimals) for value in values]


def createNewString(feature, stat1, stat2, stat1Value, stat2Value, latex=False):
    """

    @param feature: The feature of interest (e.g. Age of the tested person).
    @param stat1: A string specifying the summary statistic used for computing
    C{stat1Value} (e.g. "sum").
    @param stat2: A string specifying the summary statistic used for computing
    C{stat2Value} (e.g. "mean" or "IQRQuartiles" (Q1 and Q3).
    @param stat1Value: The value of the first summary statistic.
    @param stat2Value: The value of the second summary statistic.
    @param latex: A bool indicating whether the output string is intended for use in
    latex.
    @return: A string representing a table row (with summary statistics for a
    specific subgroup).
    """
    stat1Value = float(stat1Value)
    try:
        stat2Value = float(stat2Value)
    except TypeError:
        stat2Value = tuple(map(float, stat2Value))

    percentSign = "\\%" if latex else "%"
    indentationString = " \\newline " if latex and feature == "Age" else " "

    if stat2 == "mean":
        stat2Value *= 100
        stat2Value = f"({roundHalfUp(stat2Value)}{percentSign})"
    elif stat2 == "IQRQuartiles":
        stat2Value = tuple(map(roundHalfUp, stat2Value))
        stat2Value = f"({stat2Value[0]-stat2Value[1]})"
    else:
        stat2Value = f"({roundHalfUp(stat2Value)})"
    if stat1 == "sum":
        stat1Value = int(stat1Value)
    else:
        stat1Value = roundHalfUp(stat1Value)

    return str(stat1Value) + indentationString + f"{stat2Value}"


def createNewString2(stat, statValue, percent=True, latex=False):
    """

    @param stat: A string specifying the summary statistic used for computing
    C{statValue}.
    @param statValue: The value of the summary statistic.
    @param percent: A bool indicating whether C{statValue} should be presented as a
    percentage value.
    @param latex: A bool indicating whether the output string is intended for use in
    latex.
    @return: C{statValue} as a string (intended for use in a (latex) table).
    """
    try:
        statValue = float(statValue)
    except TypeError:
        statValue = tuple(map(float, statValue))

    percentSign = "\\%" if latex else "%"

    if stat == "mean":
        if percent:
            statValue *= 100
            statValue = f"{roundHalfUp(statValue, decimals=1)}{percentSign}"
        else:
            statValue = roundHalfUp(statValue, decimals=1)
    elif stat == "sum":
        statValue = int(statValue)
    elif stat == "IQRQuartiles":
        statValue = tuple(mapRoundHalfUp(statValue, decimals=1))
        statValue = f"{statValue[0]}, {statValue[1]}"
    else:
        statValue = roundHalfUp(statValue, decimals=1)

    return str(statValue)


def generateTableStringsEmployees(df, colFuncsDict, latex=False):
    """
    Create strings for a table with summary statistics on the study population, grouped
    by symptom status in a multicolumn format if latex==True.

    @param df: A C{pd.DataFrame} with Ag-RDT and metadata.
    @param colFuncsDict: A C{dictionary} with functions to use for computing summary
    statistics.
    @param latex: A bool indicating whether to store table in a latex format (
    default: tsv file).
    @return: A C{dict} with strings that can be used for generating a table.
    """
    summaryDict = dict(df.groupby("symptoms2").agg(colFuncsDict))
    summaryDictTotal = dict(df.agg(colFuncsDict))

    categoryDict = {
        "Female": "Gender",
        "Age": "Age (years)",
        "Cough": "Symptoms",
        "RT-PCR": "Positive test",
    }

    featureGroups = (("Age",), ("Female", "Male", "Unknown"), ("RT-PCR", "Ag-RDT"))

    tableStrings = {}

    if latex:
        header = (
            "& &"
            "\\multicolumn{3}{c|}{\\textbf{Symptoms}} & "
            "\\textbf{Total} \\\\ \n\\cline{3-5} \n\\rowcolor{lightgray} \n& & "
            "\\textbf{Yes} & \\textbf{No} & \\textbf{Unknown} & "
        )
    else:
        header = (
            "Employee tests",
            "",
            f"Symptoms",
            f"No symptoms",
            f"Unknown symptom status",
            f"Total",
        )
    tableStrings["header"] = header
    tableStrings["content"] = []

    for featureGroup in featureGroups:
        for category in featureGroup:
            stats = colFuncsDict[category]
            assert len(stats) == 2
            try:
                featureString = categoryDict[category]
            except KeyError:
                featureString = ""
            printStrings = [featureString, category]

            for symptom_status in ("Symptomatic", "Asymptomatic", "Unknown"):
                stat1, stat2 = stats
                stat1 = "IQRQuartiles" if stat1 == IQRQuartiles else stat1
                stat2 = "IQRQuartiles" if stat2 == IQRQuartiles else stat2
                newString1 = createNewString2(
                    stat1,
                    summaryDict[(category, stat1)][symptom_status],
                    latex=latex,
                )
                newString2 = createNewString2(
                    stat2,
                    summaryDict[(category, stat2)][symptom_status],
                    latex=latex,
                )
                printStrings.append(f"{newString1} ({newString2})")

            # Total counts.
            stat1, stat2 = stats
            stat1 = "IQRQuartiles" if stat1 == IQRQuartiles else stat1
            stat2 = "IQRQuartiles" if stat2 == IQRQuartiles else stat2
            newString1 = createNewString2(
                stat1, summaryDictTotal[category][stat1], latex=latex
            )
            newString2 = createNewString2(
                stat2, summaryDictTotal[category][stat2], latex=latex
            )
            printStrings.append(f"{newString1} ({newString2})")
            tableStrings["content"].append(printStrings)
    return tableStrings


def generateTableSymptomsLatex(outfile, data, includeDays=True):
    """
    Create a latex table with counts and fractions of reported symptoms for all tests
    and only PCR-positive tests.

    @param outfile: A path to the output file.
    @param data: A C{dict} with keys "header" and "content", with values for the table
    header and rows.
    @param includeDays: A bool indicating whether to include the average number of
    days between onset and testing by symptom.
    @return: A string representing a table with counts and fractions of reported
    symptoms.
    """
    latexCode = (
        SYMPTOMS_TABLE_LATEX_PREFIX2 if includeDays else (SYMPTOMS_TABLE_LATEX_PREFIX)
    )
    preHeader = (
        (
            r"\textbf{Symptom} & \multicolumn{4}{|c|}{\textbf{All tests}} & "
            r"\multicolumn{4}{|c|}{\textbf{PCR positive tests}} \\ \hline"
        )
        if includeDays
        else (
            r"\textbf{Symptom} & \multicolumn{2}{|c|}{\textbf{All tests}} & "
            r"\multicolumn{2}{|c|}{\textbf{PCR positive tests}} \\ \hline"
        )
    )
    latexCode += preHeader + "\n"
    latexCode += (
        " & ".join(f"\\textbf{{{headerCell}}}" for headerCell in data["header"])
        + r"  \\ \hline"
        + "\n"
    )
    for row in data["content"]:
        latexCode += " & ".join(row) + r" \\ \hline" + "\n"

    latexCode += (
        SYMPTOMS_TABLE_LATEX_SUFFIX2 if includeDays else SYMPTOMS_TABLE_LATEX_SUFFIX
    )
    with open(outfile, "w") as oEL:
        oEL.write(latexCode)


def writeTableEmployeesLatex(outfile, data):
    """
    Write a table in latex format with summary statistics on the study population,
    grouped by symptom status.

    @param outfile: A path to the output file.
    @param data: A C{dict} with keys "header" and "content", with values for the table
    header and rows.
    """
    latexCode = EMPLOYEE_TABLE_LATEX_PREFIX
    latexCode += data["header"] + r"  \\ " + "\n"
    feature2NumberOfCategories = {"Gender": 3, "Age (years)": 1, "Positive test": 2}
    for row in data["content"]:
        if row[0] in feature2NumberOfCategories:
            latexCode += r"\cline{1-6}" + "\n"
            if row[0] == "Age (years)":
                latexCode += r"\textbf{" + row[0] + r"} & & "
                cellWidths = ("2.4cm", "2.5cm", "2.3cm", "2.5cm")
                newAgeRow = (
                    "\\parbox[t]{" + cellWidth + "}{"
                    "\\RaggedLeft " + "\\\\".join(ageStat.split(" ")) + "}"
                    for cellWidth, ageStat in zip(cellWidths, row[2:])
                )
                latexCode += " & ".join(newAgeRow) + r" \\ \cline{2-6}" + "\n"
            else:
                nCats = feature2NumberOfCategories[row[0]]
                latexCode += (
                    rf"\multirow{{{nCats}}}{{1.4cm}}{{\RaggedLeft \textbf{{"
                    rf"{row[0]}}}}} & "
                )
                latexCode += " & ".join(row[1:]) + r" \\ \cline{2-6}" + "\n"
        else:
            latexCode += "& " + " & ".join(row[1:]) + r" \\ \cline{2-6}" + "\n"

    latexCode += "\n\\cline{1-6}" + EMPLOYEE_TABLE_LATEX_SUFFIX
    with open(outfile, "w") as oEL:
        oEL.write(latexCode)


def writeTableSymptoms(df, colFuncsDictSymptoms, outfile):
    """
    Save a table with counts and fractions of reported symptoms.

    @param df: A C{pd.DataFrame} with information on reported symptoms.
    @param colFuncsDictSymptoms: A C{dictionary} with functions to use for computing
    summary statistics on reported symptoms.
    @param outfile: The path to the output file.
    """
    summaryDict = dict(df[df.symptoms2 == "Symptomatic"].agg(colFuncsDictSymptoms))

    sortingFunc = lambda x: summaryDict[x]["mean"]
    symptomsSorted = sorted(summaryDict, key=sortingFunc, reverse=True)
    sortedDict = {key: colFuncsDictSymptoms[key] for key in symptomsSorted}

    with open(outfile, "w") as oS:
        oS.write("\t".join(["Symptom", "Count (mean)\n"]))
        for feature, stats in sortedDict.items():
            stat1, stat2 = stats
            if stat2 == IQRQuartiles:
                stat2 = "IQRQuartiles"
            newString = createNewString(
                feature,
                stat1,
                stat2,
                summaryDict[feature][stat1],
                summaryDict[feature][stat2],
            )
            oS.write("\t".join([feature, newString]) + "\n")


def generateTableStringsSymptoms2(
    df, dfPos, colFuncsDict, displayDict, latex=False, includeDays=True
):
    """
    Create strings of values to generate a table with symptom counts, fractions and
    statistics on the day of symptom onset relative to the day of onset of the
    earliest symptom.

    @param df: A C{pd.DataFrame} with information on reported symptoms.
    @param dfPos: A C{pd.DataFrame} with information on reported symptoms from tests
    with a positive PCR result.
    @param colFuncsDict: A C{dictionary} with functions to use for computing
    summary statistics on reported symptoms.
    @param displayDict: A C{dict} mapping the symptom name stored in C{df} to that to
    be displayed in the table.
    @param latex: A bool indicating whether to format strings for use in latex.
    @return: A C{dict} with keys "header" and "content", with values for the table
    header and rows.
    """
    summaryDict = dict(df[df.symptoms2 == "Symptomatic"].agg(colFuncsDict))
    summaryDictPos = dict(dfPos[dfPos.symptoms2 == "Symptomatic"].agg(colFuncsDict))
    colFuncsDictSymptoms = {
        feature: stats
        for feature, stats in colFuncsDict.items()
        if not "Days" in feature
    }

    tableStrings = {"content": []}
    header = (
        [
            "",
            "Count",
            "Percentage",
            "DPEO (mean)",
            "DPEO (sd)",
            "Count",
            "Percentage",
            "DPEO (mean)",
            "DPEO (sd)\n",
        ]
        if includeDays
        else ["", "Count", "Percentage", "Count", "Percentage"]
    )
    tableStrings["header"] = header
    sortingFunc = lambda x: summaryDictPos[x]["mean"]
    symptomsSorted = [
        key
        for key in sorted(summaryDictPos, key=sortingFunc, reverse=True)
        if not "Days" in key
    ]
    sortedDict = {key: colFuncsDictSymptoms[key] for key in symptomsSorted}
    for feature, stats in sortedDict.items():
        printStrings = [feature]
        if includeDays:
            featureDays = displayDict[feature]
            statsDays = colFuncsDict[featureDays]
            for summaryDictCurr in (summaryDict, summaryDictPos):
                printStringsDays = []
                printStringsOther = []
                for stat, statDays in zip(stats, statsDays, strict=True):
                    if stat == IQRQuartiles:
                        stat = "IQRQuartiles"
                    if statDays == IQRQuartiles:
                        statDays = "IQRQuartiles"
                    newString1 = createNewString2(
                        stat, summaryDictCurr[feature][stat], latex=latex
                    )
                    newString2 = createNewString2(
                        statDays,
                        summaryDictCurr[featureDays][statDays],
                        percent=False,
                        latex=latex,
                    )
                    printStringsOther.append(newString1)
                    printStringsDays.append(newString2)
                printStrings.extend(printStringsOther + printStringsDays)
        else:
            for summaryDictCurr in (summaryDict, summaryDictPos):
                for stat in stats:
                    if stat == IQRQuartiles:
                        stat = "IQRQuartiles"
                    newString1 = createNewString2(
                        stat, summaryDictCurr[feature][stat], latex=latex
                    )
                    printStrings.append(newString1)

        tableStrings["content"].append(printStrings)

    return tableStrings


def writeTableSymptoms2(
    df, dfPos, colFuncsDict, displayDict, latex, outfile, includeDays=True
):
    """
    Save a table as a tsv file with counts and fractions of reported symptoms,
    and statistics on the day of symptom onset relative to the day of onset of the
    earliest symptom.

    @param df: A C{pd.DataFrame} with information on reported symptoms.
    @param colFuncsDict: A C{dictionary} with functions to use for computing
    summary statistics on reported symptoms.
    @param displayDict: A C{dict} mapping the symptom name stored in C{df} to that to
    be displayed in the table.
    @param latex: A bool indicating whether to format strings for use in latex.
    @param outfile: The path to the output file.
    """

    tableStrings = generateTableStringsSymptoms2(
        df=df,
        dfPos=dfPos,
        colFuncsDict=colFuncsDict,
        displayDict=displayDict,
        latex=latex,
        includeDays=includeDays,
    )

    with open(outfile, "w") as oS:
        oS.write("\t".join(tableStrings["header"]))
        for printStrings in tableStrings["content"]:
            oS.write("\t".join(printStrings) + "\n")


def generateTableEmployeesLatex(
    dfTable, colFuncsDict, outfile=OUTPUT_TABLE_EMPLOYEES_LATEX2
):
    tableStrings = generateTableStringsEmployees(dfTable, colFuncsDict, latex=True)
    writeTableEmployeesLatex(outfile, tableStrings)


def generateTableEmployeesTSV(dfTable, colFuncsDict, outfile=OUTPUT_TABLE_EMPLOYEES2):
    tableStrings = generateTableStringsEmployees(dfTable, colFuncsDict, latex=False)
    with open(outfile, "w") as oE:
        oE.write("\t".join(tableStrings["header"]) + "\n")
        for printString in tableStrings["content"]:
            oE.write("\t".join(printString) + "\n")


def writeTableSymptomsLatex(
    dfTable, dfTablePos, colFuncsDict, displayDict, outfile=OUTPUT_TABLE_SYMPTOMS_LATEX2
):
    tableStrings = generateTableStringsSymptoms2(
        dfTable,
        dfTablePos,
        colFuncsDict,
        displayDict,
        latex=True,
        includeDays=False,
    )
    generateTableSymptomsLatex(outfile, data=tableStrings, includeDays=False)


def writeTableSymptomsLatex2(
    dfTable, dfTablePos, colFuncsDict, displayDict, outfile=OUTPUT_TABLE_SYMPTOMS_LATEX
):
    tableStrings = generateTableStringsSymptoms2(
        dfTable, dfTablePos, colFuncsDict, displayDict, latex=True
    )
    generateTableSymptomsLatex(outfile, data=tableStrings)


def generateTableStringsSymptomSeverityLatex(df, groupbyDict, colFuncsDict):
    """
    Create strings for a table with summary statistics on symptom severity, according
    to subgroups of the study population.

    @param df: A C{pd.DataFrame} with symptom severity data.
    @param groupbyDict: A C{dict} of strings specifying the categories for each
    feature (e.g. SARS-CoV-2 variant: pre-VOC, Alpha, Delta, Omicron).
    @param colFuncsDict: A C{dictionary} with functions to use for computing summary
    statistics.
    @return: A C{dict} with keys "header" and "content", with values for the table
    header and rows.
    """
    summaryDicts = defaultdict(dict)
    for feature in groupbyDict:
        for varname, stats in colFuncsDict.items():
            summaryDict = dict(df.groupby(feature, observed=True)[varname].agg(stats))
            summaryDicts[feature][varname] = summaryDict

    header = ["", "", "Count", "Median", "IQR", "Median", "IQR"]
    tableStrings = {}
    tableStrings["header"] = header
    tableStrings["content"] = []

    for feature in groupbyDict:
        nCats = len(set(groupbyDict[feature]))
        featureString = (
            rf"\specialrule{{0.75pt}}{{0pt}}{{0pt}}"
            + "\n"
            + rf"\multirow{{{nCats}}}{{*}}{{\textbf{{{feature}}}}} "
        )

        for i, category in enumerate(groupbyDict[feature]):
            featureString = "" if i else featureString
            printStrings = [
                featureString,
                category,
                str(int((df[feature] == category).sum())),
            ]
            for varname, stats in colFuncsDict.items():
                for stat in stats:
                    stat = "IQRQuartiles" if stat == IQRQuartiles else stat
                    newString = createNewString2(
                        stat,
                        summaryDicts[feature][varname][str(stat)][category],
                        latex=True,
                    )

                    printStrings.append(newString)
            tableStrings["content"].append(printStrings)

    return tableStrings


def generateTableSymptomSeverityLatex(outfile, data):
    """
    Save a table in latex format with summary statistics on symptom severity by
    subgroup in the study population.

    @param outfile: The path to the output file.
    @param data: A C{dict} with keys "header" and "content", with values for the table
    header and rows.
    """
    linerule = "\cline{2-7}"
    latexCode = SYMPTOMS_SEVERITY_TABLE_LATEX_PREFIX
    preHeader = (
        r" & & & \multicolumn{2}{c|}{\textbf{Number of symptoms}} & \multicolumn{2}{c|}{"
        r"\textbf{Severity of illness}} \\ " + linerule + "\n" + r"\rowcolor{lightgray}"
    )
    latexCode += preHeader + "\n"
    latexCode += (
        " & ".join(f"\\textbf{{{headerCell}}}" for headerCell in data["header"])
        + r"  \\ "
        + linerule
        + "\n"
    )
    for row in data["content"]:
        latexCode += " & ".join(row) + r" \\ " + linerule + "\n"

    latexCode += (
        r"\specialrule{0.75pt}{0pt}{0pt}" + "\n" + SYMPTOMS_SEVERITY_TABLE_LATEX_SUFFIX
    )
    with open(outfile, "w") as oss:
        oss.write(latexCode)


def writeSummaryTablesEmployees(df):
    # Write tables with summary statistics of the study population (as tsv file and in
    # latex format).

    # Prepare a dataframe
    df["symptoms2"] = df.symptoms.replace({0: "Asymptomatic", 1: "Symptomatic"})
    df["symptoms2"] = df["symptoms2"].fillna("Unknown")

    df["Male"] = df.gender == "M"
    df["Female"] = df.gender == "F"
    df["Unknown"] = df.gender == "U"

    colsOfInterest = [
        "Female",
        "Male",
        "Unknown",
        "age",
        "pcrPositive",
        "agrdt",
        "symptoms2",
    ]
    dfTable = df[colsOfInterest].copy()
    dfTable = dfTable.rename(
        columns={"age": "Age", "pcrPositive": "RT-PCR", "agrdt": "Ag-RDT"}
    )

    # Specify which summary values to compute.
    colFuncsFeatures = {
        "Female": ["sum", "mean"],
        "Male": ["sum", "mean"],
        "Unknown": ["sum", "mean"],
        "Age": ["median", IQRQuartiles],
    }
    colFuncsTests = {"RT-PCR": ["sum", "mean"], "Ag-RDT": ["sum", "mean"]}

    colFuncsDict = {**colFuncsFeatures, **colFuncsTests}
    generateTableEmployeesLatex(
        dfTable, colFuncsDict=colFuncsDict, outfile=OUTPUT_TABLE_EMPLOYEES_LATEX2
    )
    generateTableEmployeesTSV(
        dfTable=dfTable, colFuncsDict=colFuncsDict, outfile=OUTPUT_TABLE_EMPLOYEES2
    )


def writeSummaryTablesSymptoms(df):
    """
    Save tables with summary statistics of the symptoms reported at the time of
    testing (as a tsv file and in latex format).

    @param df: A C{pd.DataFrame} with rapid test and metadata, including on reported
    symptoms at the time of testing.
    """

    symptoms = sorted(
        [
            "fatigue",
            "headache",
            "vertigo",
            "melalgia",
            "fever",
            "cough",
            "runnyNose",
            "soreThroat",
            "dyspnea",
            "noSmell",
            "noTaste",
            "nausea",
            "noAppetite",
            "vomiting",
            "diarrhea",
        ]
    )
    symptomDays = [f"{symptom}DaysSinceLongestSymptom" for symptom in symptoms]
    symptomDaysDisplayDict = {
        SYMPTOMS_DISPLAY_DICT[symptom]: (f"{symptom}DaysSinceLongestSymptom")
        for symptom in symptoms
    }

    # Prepare a dataframe
    df["symptoms2"] = df.symptoms.replace({0: "Asymptomatic", 1: "Symptomatic"})
    df["symptoms2"] = df["symptoms2"].fillna("Unknown")
    dfTable = df.rename(columns=SYMPTOMS_DISPLAY_DICT)
    dfTablePos = dataFramePCRpos(dfTable)

    # Specify which summary values to compute.
    colFuncsSymptoms = {
        SYMPTOMS_DISPLAY_DICT[symptom]: ["sum", "mean"] for symptom in sorted(symptoms)
    }
    colFuncsSymptomDays = {symptomDay: ["mean", "std"] for symptomDay in symptomDays}
    colFuncsSymptomsWDays = {**colFuncsSymptoms, **colFuncsSymptomDays}

    writeTableSymptoms(dfTable, colFuncsSymptoms, outfile=OUTPUT_TABLE_SYMPTOMS)

    writeTableSymptoms2(
        dfTable,
        dfTablePos,
        colFuncsSymptomsWDays,
        displayDict=symptomDaysDisplayDict,
        latex=False,
        outfile=OUTPUT_TABLE_SYMPTOMS2,
    )

    writeTableSymptomsLatex(
        dfTable,
        dfTablePos,
        colFuncsSymptoms,
        SYMPTOMS_DISPLAY_DICT,
        outfile=OUTPUT_TABLE_SYMPTOMS_LATEX2,
    )

    writeTableSymptomsLatex2(
        dfTable,
        dfTablePos,
        colFuncsSymptomsWDays,
        symptomDaysDisplayDict,
        outfile=OUTPUT_TABLE_SYMPTOMS_LATEX,
    )


def writeSummaryTableSymptomSeverity(df):
    """
    Write a table (in latex format) of summary statistics on symptom severity.

    @param df: A C{pd.DataFrame} with rapid test and metadata, including on symptom
    severity.
    @return:
    """
    abbrvs = getAbbrvsDict()
    dfTable = df.rename(
        columns={
            "gender": "Gender",
            "ageBin": "Age (years)",
            "variant": "Variant",
            "immun2YN": "Immunized",
            "nSymptoms": "Number of symptoms",
            "illSeverity": "Severity of illness",
        }
    )
    dfTable["Variant"] = dfTable["Variant"].replace(abbrvs["variant"])
    dfTable["Age (years)"] = dfTable["Age (years)"].replace(abbrvs["ageBin"])
    dfTable["Gender"] = dfTable["Gender"].replace(
        {"F": "Female", "M": "Male", "U": "Unknown"}
    )
    dfTable["Immunized"] = dfTable["Immunized"].replace(abbrvs["immun2YN"])
    groupbyDict = {
        "Gender": ["Female", "Male", "Unknown"],
        "Age (years)": abbrvs["ageBin"].values(),
        "Variant": abbrvs["variant"].values(),
        "Immunized": abbrvs["immun2YN"].values(),
    }
    # Specify which summary values to compute.
    colFuncsDict = {
        "Number of symptoms": ["median", IQRQuartiles],
        "Severity of illness": ["median", IQRQuartiles],
    }

    tableStrings = generateTableStringsSymptomSeverityLatex(
        dfTable, groupbyDict, colFuncsDict
    )
    generateTableSymptomSeverityLatex(
        outfile=OUTPUT_TABLE_SYMPTOM_SEVERITY_LATEX, data=tableStrings
    )


def generateIDataSummaryDf(df, varnames, latex=False):
    """
    Generate a table with MCMC summary statistics.

    @param df: A C{pd.DataFrame} with MCMC summary statistics created by az.summary.
    @param varnames: An C{iterable} with variable names of interest (present in index
    column of C{df}).
    @param latex: A C{bool} indicating whether to format cells in the resulting table
    for use in latex.
    @return: A C{pd.DataFrame} that can be used to create the final table.
    """

    label = LABEL_LATEX if latex else LABEL
    paramNameMappings = (
        PARAM_NAME_MAPPINGS_TABLE_LATEX if latex else PARAM_NAME_MAPPINGS
    )
    paramNameMappings["immun2YN[1]"] = label.immun2YN
    paramNameMappings["zAge"] = "Age"

    params = []
    for varname in varnames:
        for param in df.index:
            if re.match(varname, param):
                params.append(param)
    params = [param for param in params if param in paramNameMappings]

    # Extract relevant rows.
    df = df.loc[params]
    df = df[["mean", "sd", "hdi_3%", "hdi_97%"]].round(decimals=2)
    df["hdi_94%"] = df["hdi_3%"].astype(str) + ", " + df["hdi_97%"].astype(str)
    # Extract relevant columns.
    df = df[["mean", "sd", "hdi_3%", "hdi_97%"]]

    if latex:
        columns = {
            "hdi_3%": "3\%",
            "hdi_97%": "97\%",
            "mean": "Mean",
            "sd": "Sd",
        }
        df = df.applymap(lambda x: f"${x}$")  # This is to display minus signs
        # correctly in the latex table
        index = {
            param: paramNameMappings.get(param).replace("-", "--")
            or param.replace("-", "--")
            for param in df.index
        }
    else:
        columns = {
            "hdi_3%": "3%",
            "hdi_97%": "97%",
            "mean": "Mean",
            "sd": "Sd",
        }
        index = {param: paramNameMappings.get(param) or param for param in df.index}

    df.rename(
        index=index,
        columns=columns,
        inplace=True,
    )

    return df


def writeIDataSummaryTable(summaryDf, varnames, outfile):
    """
    Save MCMC summary statistics to a table in tsv format.
    @param summaryDf: A C{pd.DataFrame} with MCMC summary statistics created by
    az.summary.
    @param varnames: An C{iterable} with variable names of interest (present in index
    column of C{df}).
    @param outfile: The path to the output file.
    """
    # Save a table of MCMC summary statistics
    summaryDf = generateIDataSummaryDf(summaryDf, varnames)
    summaryDf.to_csv(outfile, sep="\t")


def writeIDataSummaryTableLatex(summaryDf, varnames, outfile, model_no=1):
    """
    Save MCMC summary statistics to a table in latex format.
    @param summaryDf: A C{pd.DataFrame} with MCMC summary statistics created by
    az.summary.
    @param varnames: An C{iterable} with variable names of interest (present in index
    column of C{df}).
    @param outfile: The path to the output file.
    @param model_no: The model number for which C{summaryDf} contains the MCMC results.
    """
    summaryDf = generateIDataSummaryDf(summaryDf, varnames, latex=True)
    assert model_no in (1, 2, 3, 4, 5, 6, None)
    if model_no == 1:
        latexCode = IDATA_TABLE_PREFIX1
    elif model_no == 2:
        latexCode = IDATA_TABLE_PREFIX2
    elif model_no == 3:
        latexCode = IDATA_TABLE_PREFIX3
    elif model_no == 4:
        latexCode = IDATA_TABLE_PREFIX4
    elif model_no == 5:
        latexCode = IDATA_TABLE_PREFIX5
    elif model_no == 6:
        latexCode = IDATA_TABLE_PREFIX6
    else:
        latexCode = IDATA_TABLE_PREFIX
    latexCode += "Parameter & " + " & ".join(summaryDf.columns) + "\\\\\hline" + "\n"
    for index, row in summaryDf.iterrows():
        index = index.replace("|", "\\textbar")
        latexCode += (
            " & ".join([index] + [str(value) for value in row.values])
            + "\\\\\hline"
            + "\n"
        )
    latexCode += IDATA_TABLE_SUFFIX
    with open(outfile, "w") as out:
        out.write(latexCode)


def returnTableDir():
    return TABLE_DIR


def returnTableDirIData():
    return TABLE_DIR_IDATA
