import marimo

__generated_with = "0.16.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from datetime import date
    from collections import defaultdict

    # Bayesian MCMC analysis
    import arviz as az

    from agrdt.data import createDataFramesFigures
    from agrdt.tables import roundHalfUp, mapRoundHalfUp
    from agrdt.plotting import (annotateWithLetter, annotateWithLetters,
                                returnPlotDirSensitivity, plotFig1_A, plotFig3_B,
                                plotFigA2_B, setCustomTheme, saveFigure)
    from agrdt.plotParams import (getPalettes, getLabels, getLegends, getOrders, getAbbrvsDict, COL_WIDTH, CM,
                                  ANNOTATION_LETTER_SIZE, ANNOTATION_COORDS)
    from agrdt.regression import returnIDataDirSensitivity
    from agrdt.dataParams import ROOT_DIR
    setCustomTheme()
    return (
        ANNOTATION_COORDS,
        ANNOTATION_LETTER_SIZE,
        CM,
        COL_WIDTH,
        ROOT_DIR,
        annotateWithLetter,
        annotateWithLetters,
        az,
        createDataFramesFigures,
        date,
        defaultdict,
        getAbbrvsDict,
        getLabels,
        getLegends,
        getOrders,
        getPalettes,
        mapRoundHalfUp,
        mo,
        np,
        pd,
        plotFig1_A,
        plotFig3_B,
        plotFigA2_B,
        plt,
        returnIDataDirSensitivity,
        returnPlotDirSensitivity,
        roundHalfUp,
        saveFigure,
    )


@app.cell
def _(mo):
    mo.md(r"""# Load""")
    return


@app.cell
def _(mo):
    mo.md(r"""## Plotting parameters""")
    return


@app.cell
def _(
    getLabels,
    getOrders,
    getPalettes,
    returnPlotDirSensitivity,
):
    pal = getPalettes()
    order = getOrders()
    label = getLabels()
    plotDir = returnPlotDirSensitivity()
    return label, order, pal, plotDir


@app.cell
def _(getAbbrvsDict):
    abbrvDictPaper = getAbbrvsDict()
    variantAbbrvs = abbrvDictPaper["variant"]

    return abbrvDictPaper, variantAbbrvs


@app.cell
def _(COL_WIDTH):
    colWidth = COL_WIDTH
    fontSizePlot = 12
    return colWidth, fontSizePlot


@app.cell
def _(mo):
    mo.md(r"""## Data""")
    return


@app.cell
def _(ROOT_DIR, pd):
    # Load data
    target = "T1"
    dataPath = ROOT_DIR / "data" / f"agrdtDataThesis{target}.tsv"
    df = pd.read_csv(dataPath, sep="\t", low_memory=False, parse_dates=['pcrDate'])
    # Turn into datetime.date objects.
    df["pcrDate"] = df["pcrDate"].apply(lambda pcrDate: pcrDate.date())
    return (df,)


@app.cell
def _(mo):
    mo.md(r"""## Dataframes for figures""")
    return


@app.cell
def _(createDataFramesFigures, df):
    dfTuple = createDataFramesFigures(df)
    (dfAgrdt,
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
     dfFigureA2) = dfTuple
    return dfFigure1_A, dfFigure1_A_Asymp, dfFigure3, dfFigureA2, dfPos


@app.cell
def _(mo):
    mo.md(r"""## Labels""")
    return


@app.cell
def _(abbrvDictPaper, date, dfPos):
    month2Labels = list(abbrvDictPaper['samplingMonth2'].values())
    months2 = dfPos[dfPos.pcrDate > date(2020, 11, 30)].samplingMonth2.sort_values().unique()
    month2ToMonths = {monthBin: label for monthBin, label in zip(months2, month2Labels)}
    return (month2Labels,)


@app.cell
def _(mo):
    mo.md(r"""# Plots and stats""")
    return


@app.cell
def _(mo):
    mo.md(r"""## Sensitivity over time""")
    return


@app.cell
def _(defaultdict):
    statsDictTime = defaultdict(lambda: defaultdict(dict))
    statsDictTimeAsymp = defaultdict(lambda: defaultdict(dict))
    return statsDictTime, statsDictTimeAsymp


@app.cell
def _(dfFigure1_A):
    len(dfFigure1_A)
    return


@app.cell
def _(
    ANNOTATION_LETTER_SIZE,
    abbrvDictPaper,
    annotateWithLetters,
    colWidth,
    dfFigure1_A,
    dfFigure1_A_Asymp,
    dfPos,
    label,
    month2Labels,
    order,
    pal,
    plotDir,
    plotFig1_A,
    plt,
    saveFigure,
    statsDictTime,
    statsDictTimeAsymp,
):
    _figA_full, _ax = plt.subplots(3, 1, figsize=(colWidth * 2, colWidth * 2), 
                                   gridspec_kw={'height_ratios': [0.02, 1, 1]}, 
                                   constrained_layout=True)
    modelSensSymp, iDataSensSymp = plotFig1_A(dfFigure1_A, dfPos, timeVar='samplingMonth2', statsDictTime=statsDictTime,
                                              palette=pal, label=label, xLabels=None, abbrvDictPaper=abbrvDictPaper, 
                                              order=order,
                                              ax=_ax[:2])
    modelSensAsymp, iDataSensAsymp = plotFig1_A(dfFigure1_A_Asymp, dfPos, timeVar='samplingMonth2',
                                                statsDictTime=statsDictTimeAsymp, palette=pal, label=label,
                                                xLabels=month2Labels, abbrvDictPaper=abbrvDictPaper, order=order, ax=_ax[2],
                                                showVariantTimes=False)
    plt.tight_layout()
    annotateWithLetters(_ax[1:], size=ANNOTATION_LETTER_SIZE, coords=(-0.2, 1))
    saveFigure(_figA_full, plotDir / 'Figure1-sens-thesis.png')
    saveFigure(_figA_full, plotDir / 'Figure1-sens-thesis.pdf')
    _ax
    return iDataSensAsymp, iDataSensSymp, modelSensAsymp, modelSensSymp


@app.cell
def _(modelSensSymp):
    modelSensSymp
    return


@app.cell
def _(modelSensAsymp):
    modelSensAsymp
    return


@app.cell
def _(dfFigure1_A, dfFigure1_A_Asymp, iDataSensAsymp, iDataSensSymp):
    referenceLevel = iDataSensSymp.posterior.sel({"samplingMonth2__factor_dim": "0"})["1|samplingMonth2"]
    for _month2Cat in dfFigure1_A.samplingMonth2.unique():
        if _month2Cat:
            currentLevel = iDataSensSymp.posterior.sel({"samplingMonth2__factor_dim":str(int(_month2Cat))})["1|samplingMonth2"]
            iDataSensSymp.posterior[f"1|samplingMonth2[0-{int(_month2Cat)}]"] = referenceLevel - currentLevel


    referenceLevel = iDataSensAsymp.posterior.sel({"samplingMonth2__factor_dim": "0"})["1|samplingMonth2"]
    for _month2Cat in dfFigure1_A_Asymp.samplingMonth2.unique():
        if _month2Cat:
            currentLevel = iDataSensAsymp.posterior.sel({"samplingMonth2__factor_dim":str(int(_month2Cat))})["1|samplingMonth2"]
            iDataSensAsymp.posterior[f"1|samplingMonth2[0-{int(_month2Cat)}]"] = referenceLevel - currentLevel
    return


@app.cell
def _(az, iDataSensSymp):
    az.summary(iDataSensSymp, var_names='~p')
    return


@app.cell
def _(az, iDataSensAsymp):
    az.summary(iDataSensAsymp, var_names='~p')
    return


@app.cell
def _(mo):
    mo.md(r"""#### Counts""")
    return


@app.cell
def _(date, dfFigure1_A, dfFigure1_A_Asymp, dfPos):
    _months2 = dfPos[dfPos.pcrDate > date(2020, 11, 30)].samplingMonth2.sort_values().unique()
    countsFigure1A = [f'n={(dfFigure1_A.samplingMonth2 == month2).sum()}' for month2 in _months2]
    countsFigure1A_Asymp = [f'n={(dfFigure1_A_Asymp.samplingMonth2 == month2).sum()}' for month2 in _months2]
    return countsFigure1A, countsFigure1A_Asymp


@app.cell
def _(countsFigure1A):
    countsFigure1A
    return


@app.cell
def _(countsFigure1A_Asymp):
    countsFigure1A_Asymp
    return


@app.cell
def _(mo):
    mo.md(r"""## Ag-RDT sensitivity difference between immunized and immune naive people""")
    return


@app.cell
def _(defaultdict):
    statsDictImmun = defaultdict(lambda: defaultdict(dict))
    statsDictImmunAsymp = defaultdict(lambda: defaultdict(dict))
    return (statsDictImmun,)


@app.cell
def _(
    ANNOTATION_COORDS,
    ANNOTATION_LETTER_SIZE,
    annotateWithLetter,
    colWidth,
    dfFigure3,
    fontSizePlot,
    label,
    plotDir,
    plotFig3_B,
    plt,
    saveFigure,
    statsDictImmun,
):
    fig_sense_immun, _ax = plt.subplots(1, 1, figsize=(colWidth * 2, colWidth * 0.8), constrained_layout=True)
    model_immun, iData_immun = plotFig3_B(dfFigure3, statsDictImmun, label=label, immunSympPalette="black", ax=_ax)

    annotateWithLetter(_ax, "A", size=ANNOTATION_LETTER_SIZE, coords=ANNOTATION_COORDS)
    fig_sense_immun.text(0.34, 1, 'Symptomatic', ha='center', fontdict={'fontsize': fontSizePlot})
    fig_sense_immun.text(0.79, 1, 'Asymptomatic', ha='center', fontdict={'fontsize': fontSizePlot})
    saveFigure(fig_sense_immun, plotDir / 'Figure2A-sens-thesis.pdf')
    saveFigure(fig_sense_immun, plotDir / 'Figure2A-sens-thesis.png')
    _ax
    return (model_immun,)


@app.cell
def _(model_immun):
    model_immun
    return


@app.cell
def _(roundHalfUp, statsDictImmun):
    sympImmuneNaiveMean = statsDictImmun["immun2YN:symptoms"]["0, 1"]["mean"]
    sympImmuneNaiveHdi = statsDictImmun["immun2YN:symptoms"]["0, 1"]["hdi"]

    sympImmunizedMean = statsDictImmun["immun2YN:symptoms"]["1, 1"]["mean"]
    sympImmunizedHdi = statsDictImmun["immun2YN:symptoms"]["1, 1"]["hdi"]

    print("Ag-RDT sensitivity for symptomatic people:")
    print(f"Immune naive: {roundHalfUp(sympImmuneNaiveMean)}, {list(map(roundHalfUp, sympImmuneNaiveHdi))}")
    print(f"Immunized: {roundHalfUp(sympImmunizedMean)}, {list(map(roundHalfUp, sympImmunizedHdi))}")
    return


@app.cell
def _(roundHalfUp, statsDictImmun):
    asympImmuneNaiveMean = statsDictImmun["immun2YN:symptoms"]["0, 0"]["mean"]
    asympImmuneNaiveHdi = statsDictImmun["immun2YN:symptoms"]["0, 0"]["hdi"]

    asympImmunizedMean = statsDictImmun["immun2YN:symptoms"]["1, 0"]["mean"]
    asympImmunizedHdi = statsDictImmun["immun2YN:symptoms"]["1, 0"]["hdi"]

    print("Ag-RDT sensitivity for asymptomatic people:")
    print(f"Immune naive: {roundHalfUp(asympImmuneNaiveMean)}, {list(map(roundHalfUp, asympImmuneNaiveHdi))}")
    print(f"Immunized: {roundHalfUp(asympImmunizedMean)}, {list(map(roundHalfUp, asympImmunizedHdi))}")
    return


@app.cell
def _(az, np, statsDictImmun):
    diffImmunAsymp = statsDictImmun["immun2YN:symptoms"]["0, 0"]["samples"] - statsDictImmun["immun2YN:symptoms"]["1, 0"]["samples"]
    diffImmunSymp = statsDictImmun["immun2YN:symptoms"]["0, 1"]["samples"] - statsDictImmun["immun2YN:symptoms"]["1, 1"]["samples"]

    immunAsympMean = np.mean(diffImmunAsymp)
    immunAsympHdi = az.hdi(diffImmunAsymp, 0.94)
    immunSympMean = np.mean(diffImmunSymp)
    immunSympHdi = az.hdi(diffImmunSymp, 0.94)
    return immunAsympHdi, immunAsympMean, immunSympHdi, immunSympMean


@app.cell
def _(
    immunAsympHdi,
    immunAsympMean,
    immunSympHdi,
    immunSympMean,
    mapRoundHalfUp,
):
    print("Mean sensitivity differences (HPDIs) between immune naive and multiply-immunized people")
    print(f"Asymptomatics: {immunAsympMean:.4f} ({mapRoundHalfUp(immunAsympHdi, 4)})")
    print(f"Symptomatics: {immunSympMean:.4f} ({mapRoundHalfUp(immunSympHdi, 4)})")
    return


@app.cell
def _(mo):
    mo.md(r"""## Ag-RDT sensitivity difference between SARS-CoV-2 variants""")
    return


@app.cell
def _(defaultdict):
    statsDictVariant = defaultdict(lambda: defaultdict(dict))
    statsDictVariantAsymp = defaultdict(lambda: defaultdict(dict))
    return (statsDictVariant,)


@app.cell
def _(
    ANNOTATION_COORDS,
    ANNOTATION_LETTER_SIZE,
    annotateWithLetter,
    colWidth,
    dfFigureA2,
    label,
    plotDir,
    plotFigA2_B,
    plt,
    saveFigure,
    statsDictVariant,
):
    fig_sense_variant, _ax = plt.subplots(1, 1, figsize=(colWidth * 2, colWidth * 0.95), constrained_layout=True)
    model_variant, iData_variant = plotFigA2_B(dfFigureA2, statsDictVariant, variantSympPalette="black", label=label, ax=_ax)

    annotateWithLetter(_ax, "B", size=ANNOTATION_LETTER_SIZE, coords=ANNOTATION_COORDS)
    saveFigure(fig_sense_variant, plotDir / 'Figure2B-sens-thesis.pdf')
    saveFigure(fig_sense_variant, plotDir / 'Figure2B-sens-thesis.png')
    _ax
    return (model_variant,)


@app.cell
def _(model_variant):
    model_variant
    return


@app.cell
def _(roundHalfUp, statsDictVariant):
    sympWtMean = statsDictVariant["variant:symptoms"]["wildtype, 1"]["mean"]
    sympWtHdi = statsDictVariant["variant:symptoms"]["wildtype, 1"]["hdi"]

    sympAlphaMean = statsDictVariant["variant:symptoms"]["alpha, 1"]["mean"]
    sympAlphaHdi = statsDictVariant["variant:symptoms"]["alpha, 1"]["hdi"]

    sympDeltaMean = statsDictVariant["variant:symptoms"]["delta, 1"]["mean"]
    sympDeltaHdi = statsDictVariant["variant:symptoms"]["delta, 1"]["hdi"]

    sympOmicronMean = statsDictVariant["variant:symptoms"]["omicron, 1"]["mean"]
    sympOmicronHdi = statsDictVariant["variant:symptoms"]["omicron, 1"]["hdi"]

    print("Ag-RDT sensitivity in symptomatic people:")
    print(f"Pre-VOC: {roundHalfUp(sympWtMean)}, {list(map(roundHalfUp, sympWtHdi))}")
    print(f"Alpha: {roundHalfUp(sympAlphaMean)}, {list(map(roundHalfUp, sympAlphaHdi))}")
    print(f"Delta: {roundHalfUp(sympDeltaMean)}, {list(map(roundHalfUp, sympDeltaHdi))}")
    print(f"Omicron: {roundHalfUp(sympOmicronMean)}, {list(map(roundHalfUp, sympOmicronHdi))}")
    return


@app.cell
def _(az, mapRoundHalfUp, np, statsDictVariant, variantAbbrvs):
    for _variant, _variantLabel in variantAbbrvs.items():
        _diffVariantAsymp = (statsDictVariant["variant:symptoms"]["wildtype, 0"]["samples"] -
                             statsDictVariant["variant:symptoms"][f"{_variant}, 0"]["samples"])
        _diffVariantSymp = (statsDictVariant["variant:symptoms"]["wildtype, 1"]["samples"] -
                            statsDictVariant["variant:symptoms"][f"{_variant}, 1"]["samples"])
        _variantAsympMean = np.mean(_diffVariantAsymp)
        _variantAsympHdi = az.hdi(_diffVariantAsymp, 0.94)
        _variantSympMean = np.mean(_diffVariantSymp)
        _variantSympHdi = az.hdi(_diffVariantSymp, 0.94)

        print(f"Mean sensitivity differences (HPDIs) between infections with pre-VOC and {_variantLabel}")
        print(f"Asymptomatics: {_variantAsympMean:.4f} ({mapRoundHalfUp(_variantAsympHdi, 4)})")
        print(f"Symptomatics: {_variantSympMean:.4f} ({mapRoundHalfUp(_variantSympHdi, 4)})")
        print()
    return


@app.cell
def _(roundHalfUp, statsDictVariant):
    asympWtMean = statsDictVariant["variant:symptoms"]["wildtype, 0"]["mean"]
    asympWtHdi = statsDictVariant["variant:symptoms"]["wildtype, 0"]["hdi"]

    asympAlphaMean = statsDictVariant["variant:symptoms"]["alpha, 0"]["mean"]
    asympAlphaHdi = statsDictVariant["variant:symptoms"]["alpha, 0"]["hdi"]

    asympDeltaMean = statsDictVariant["variant:symptoms"]["delta, 0"]["mean"]
    asympDeltaHdi = statsDictVariant["variant:symptoms"]["delta, 0"]["hdi"]

    asympOmicronMean = statsDictVariant["variant:symptoms"]["omicron, 0"]["mean"]
    asympOmicronHdi = statsDictVariant["variant:symptoms"]["omicron, 0"]["hdi"]

    print("Ag-RDT sensitivity in symptomatic people:")
    print(f"Pre-VOC: {roundHalfUp(asympWtMean)}, {list(map(roundHalfUp, asympWtHdi))}")
    print(f"Alpha: {roundHalfUp(asympAlphaMean)}, {list(map(roundHalfUp, asympAlphaHdi))}")
    print(f"Delta: {roundHalfUp(asympDeltaMean)}, {list(map(roundHalfUp, asympDeltaHdi))}")
    print(f"Omicron: {roundHalfUp(asympOmicronMean)}, {list(map(roundHalfUp, asympOmicronHdi))}")
    return


if __name__ == "__main__":
    app.run()
