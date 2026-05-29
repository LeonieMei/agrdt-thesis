import marimo

__generated_with = "0.16.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    import matplotlib as mpl
    import seaborn as sns
    from datetime import date, datetime
    from pathlib import Path
    from scipy.stats import skew

    # Bayesian MCMC analysis
    import arviz as az
    import bambi as bmb

    from agrdt.data import createDataFramesFigures, returnRtData, IQR, IQRQuartiles
    from agrdt.tables import roundHalfUp
    from agrdt.plotting import (annotateWithLetter, annotateWithLetters, replaceLegend,
                                setFontSize, returnPlotDirViralLoad, plotFig1_B,
                                plotFig2_B, plotFig3_A, plotFigA2_A, setCustomTheme,
                                saveFigure, plotFigVlResult, plotFigVlTestline)
    from agrdt.plotParams import (getPalettes, getLabels, getLegends, getOrders, getAbbrvsDict, COL_WIDTH,
                                  ANNOTATION_LETTER_SIZE, GGPLOT_PALETTE, ANNOTATION_COORDS)
    from agrdt.regression import SEED, sampleVl, returnIDataDirViralLoad
    from agrdt.dataParams import ROOT_DIR

    setCustomTheme()
    return (
        ANNOTATION_COORDS,
        ANNOTATION_LETTER_SIZE,
        COL_WIDTH,
        GGPLOT_PALETTE,
        IQR,
        IQRQuartiles,
        Path,
        ROOT_DIR,
        SEED,
        annotateWithLetter,
        annotateWithLetters,
        az,
        bmb,
        createDataFramesFigures,
        date,
        datetime,
        getAbbrvsDict,
        getLabels,
        getLegends,
        getOrders,
        getPalettes,
        mo,
        mpl,
        np,
        pd,
        plotFig1_B,
        plotFig2_B,
        plotFig3_A,
        plotFigA2_A,
        plotFigVlResult,
        plotFigVlTestline,
        plt,
        replaceLegend,
        returnIDataDirViralLoad,
        returnPlotDirViralLoad,
        returnRtData,
        roundHalfUp,
        sampleVl,
        saveFigure,
        setFontSize,
        skew,
        sns,
        ticker,
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
    getLegends,
    getOrders,
    getPalettes,
    returnIDataDirViralLoad,
    returnPlotDirViralLoad,
):
    pal = getPalettes()
    legend = getLegends()
    order = getOrders()
    label = getLabels()
    plotDir = returnPlotDirViralLoad()
    iDataDir = returnIDataDirViralLoad()
    return iDataDir, label, legend, order, pal, plotDir


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
def _(ROOT_DIR, datetime, pd, returnRtData):
    # Load data
    target = "T1"
    cusmaMaxDate = datetime(2022, 2, 11)
    cusmaMinDate = datetime(2020, 12, 1)
    dataPath = ROOT_DIR / "data" / f"agrdtDataThesis{target}.tsv"
    df = pd.read_csv(dataPath, sep="\t", low_memory=False, parse_dates=["pcrDate"])

    # Turn into datetime.date objects.
    df["pcrDate"] = df["pcrDate"].apply(lambda pcrDate: pcrDate.date())
    dfRt = returnRtData()
    dfRt = dfRt[(dfRt.date <= cusmaMaxDate) & (dfRt.date >= cusmaMinDate)]
    df = df.merge(dfRt[["dateDate", "rtRolling4", "rt"]], left_on="pcrDate", right_on="dateDate", how="left")
    df["epidemicGrowth"] = df.rt >= 1
    return cusmaMaxDate, cusmaMinDate, df, dfRt


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
    return (
        dfAgrdtIndInf,
        dfAllFirstPosPcrsNoRelease,
        dfFigure1_B,
        dfFigure1_B_Asymp,
        dfFigure2_B,
        dfFigure2_BAsymp,
        dfFigure3,
        dfFigureA2,
        dfPos,
    )


@app.cell
def _(mo):
    mo.md(r"""## Labels""")
    return


@app.cell
def _(abbrvDictPaper, date, dfPos):
    month2Labels = list(abbrvDictPaper["samplingMonth2"].values())
    months2 = dfPos[dfPos.pcrDate > date(2020, 11, 30)].samplingMonth2.sort_values().unique()
    month2ToMonths = {monthBin: label for monthBin, label in zip(months2, month2Labels)}
    return month2Labels, month2ToMonths, months2


@app.cell
def _(dfFigure1_B, dfFigure1_B_Asymp, months2):
    countsFigure1 = [f"n={(dfFigure1_B.samplingMonth2 == month2).sum()}" for month2 in months2]
    countsFigure1Asymp = [f"n={(dfFigure1_B_Asymp.samplingMonth2 == month2).sum()}" for month2 in months2]
    return countsFigure1, countsFigure1Asymp


@app.cell
def _(mo):
    mo.md(r"""# Plots""")
    return


@app.cell
def _(mo):
    mo.md(r"""## Viral load over time""")
    return


@app.cell
def _(
    abbrvDictPaper,
    colWidth,
    dfFigure1_B,
    dfPos,
    label,
    month2Labels,
    month2ToMonths,
    pal,
    plotDir,
    plotFig1_B,
    plt,
    saveFigure,
):
    _figB, _ax = plt.subplots(2, 1, figsize=(colWidth * 2, colWidth * 1.6), gridspec_kw={"height_ratios": [0.02, 1]})
    plotFig1_B(dfFigure1_B, dfPos, timeVar="samplingMonth2", xOrder=month2ToMonths, palette=pal, label=label,
               xLabels=month2Labels, abbrvDictPaper=abbrvDictPaper, ax=_ax)
    saveFigure(_figB, plotDir / "Figure1-B.png")
    saveFigure(_figB, plotDir / "Figure1-B.pdf")
    _ax
    return


@app.cell
def _(mo):
    mo.md(r"""#### Plotting both symptomatic and asymptomatic infections""")
    return


@app.cell
def _(
    ANNOTATION_COORDS,
    ANNOTATION_LETTER_SIZE,
    abbrvDictPaper,
    annotateWithLetters,
    colWidth,
    dfFigure1_B,
    dfFigure1_B_Asymp,
    dfPos,
    label,
    month2Labels,
    month2ToMonths,
    pal,
    plotDir,
    plotFig1_B,
    plt,
    saveFigure,
):
    _figB_full, _ax = plt.subplots(3, 1, figsize=(colWidth * 2, colWidth * 2.3), 
                              gridspec_kw={"height_ratios": [0.02, 1, 1]}, constrained_layout=True)
    plotFig1_B(dfFigure1_B, dfPos, timeVar="samplingMonth2", xOrder=month2ToMonths, palette=pal, label=label, xLabels=None,
               abbrvDictPaper=abbrvDictPaper, ax=_ax[:2])
    plotFig1_B(dfFigure1_B_Asymp, dfPos, timeVar="samplingMonth2", xOrder=month2ToMonths, palette=pal, label=label,
               xLabels=month2Labels, abbrvDictPaper=abbrvDictPaper, ax=_ax[2], showVariantTimes=False)
    annotateWithLetters(_ax[1:], size=ANNOTATION_LETTER_SIZE, coords=ANNOTATION_COORDS)
    saveFigure(_figB_full, plotDir / "Figure1-vl-thesis.png")
    saveFigure(_figB_full, plotDir / "Figure1-vl-thesis.pdf")
    _ax
    return


@app.cell
def _(countsFigure1, month2ToMonths):
    print("Number of tests (symptomatic)")
    for _timePeriod, _nTests in zip(month2ToMonths.values(), countsFigure1):
        print(f"{_timePeriod}: {_nTests}")
    return


@app.cell
def _(dfFigure1_B, month2ToMonths, roundHalfUp):
    print("Viral load medians (symptomatic)\n")
    for _idxMonth2, vlMedian in dfFigure1_B.groupby(by=["samplingMonth2"]).vl.agg("median").items():
        print(f"{month2ToMonths[_idxMonth2]}: {roundHalfUp(vlMedian)}")
    return


@app.cell
def _(IQRQuartiles, dfFigure1_B, month2ToMonths, roundHalfUp):
    print("Viral load IQRs (Q1, Q3) (symptomatic)\n")
    for _idxMonth2, (Q1, Q3) in dfFigure1_B.groupby(by=["samplingMonth2"]).vl.agg(IQRQuartiles).items():
        print(f"{month2ToMonths[_idxMonth2]}: {(roundHalfUp(Q1), roundHalfUp(Q3))}")
    return


@app.cell
def _(countsFigure1Asymp, month2ToMonths):
    print("Number of tests (asymptomatic)")
    for _timePeriod, _nTests in zip(month2ToMonths.values(), countsFigure1Asymp):
        print(f"{_timePeriod}: {_nTests}")
    return


@app.cell
def _(dfFigure1_B_Asymp, month2ToMonths, roundHalfUp):
    print("Viral load medians (asymptomatic)\n")
    for _idxMonth2, vlMedianAsymp in dfFigure1_B_Asymp.groupby(by=["samplingMonth2"]).vl.agg("median").items():
        print(f"{month2ToMonths[_idxMonth2]}: {roundHalfUp(vlMedianAsymp)}")
    return


@app.cell
def _(IQRQuartiles, dfFigure1_B_Asymp, month2ToMonths, roundHalfUp):
    print("Viral load IQRs (Q1, Q3) (asymptomatic)\n")
    for _idxMonth2, (Q1Asymp, Q3Asymp) in dfFigure1_B_Asymp.groupby(by=["samplingMonth2"]).vl.agg(IQRQuartiles).items():
        print(f"{month2ToMonths[_idxMonth2]}: {(roundHalfUp(Q1Asymp), roundHalfUp(Q3Asymp))}")
    return


@app.cell
def _(mo):
    mo.md(r"""## Viral load by subgroup""")
    return


@app.cell
def _(
    ANNOTATION_COORDS,
    ANNOTATION_LETTER_SIZE,
    annotateWithLetter,
    colWidth,
    dfFigure3,
    fontSizePlot,
    label,
    order,
    pal,
    plotDir,
    plotFig3_A,
    plt,
    saveFigure,
):
    fig3A, _ax = plt.subplots(1, 1, figsize=(colWidth * 2, colWidth * 0.8), constrained_layout=True)
    plotFig3_A(dfFigure3, palette=pal, order=order, label=label, xLabels=True, ax=_ax, markersize=2.3)
    annotateWithLetter(_ax, "A", size=ANNOTATION_LETTER_SIZE, coords=ANNOTATION_COORDS)
    fig3A.text(0.34, 1, "Symptomatic", ha="center", fontdict={"fontsize": fontSizePlot})
    fig3A.text(0.79, 1, "Asymptomatic", ha="center", fontdict={"fontsize": fontSizePlot})
    saveFigure(fig3A, plotDir / "Figure2A-vl-thesis.pdf")
    saveFigure(fig3A, plotDir / "Figure2A-vl-thesis.png")
    _ax
    return


@app.cell
def _(IQR, IQRQuartiles, dfFigure3):
    print("Viral load medians (Q1, Q3)")
    medianDictFig3 = dict(dfFigure3[dfFigure3.symptoms==1].groupby("immun2YN").vl.median().items())
    iqrsDictFig3 = dict(dfFigure3[dfFigure3.symptoms==1].groupby("immun2YN").vl.apply(IQR))
    iqrsDict2Fig3 = dict(dfFigure3[dfFigure3.symptoms==1].groupby("immun2YN").vl.apply(IQRQuartiles))

    medianDictFig3Asymp = dict(dfFigure3[dfFigure3.symptoms==0].groupby("immun2YN").vl.median().items())
    iqrsDictFig3Asymp = dict(dfFigure3[dfFigure3.symptoms==0].groupby("immun2YN").vl.apply(IQR))
    iqrsDict2Fig3Asymp = dict(dfFigure3[dfFigure3.symptoms==0].groupby("immun2YN").vl.apply(IQRQuartiles))

    print("Symptomatic people")
    print(f"Not multiply-immunized: {medianDictFig3[0]:.2f} ({iqrsDict2Fig3[0][0]:.2f}-{iqrsDict2Fig3[0][1]:.2f})")
    print(f"Multiply-immunized: {medianDictFig3[1]:.2f} ({iqrsDict2Fig3[1][0]:.2f}-{iqrsDict2Fig3[1][1]:.2f})")
    print()
    print("Asymptomatic people")
    print(f"Not multiply-immunized: {medianDictFig3Asymp[0]:.2f} ({iqrsDict2Fig3Asymp[0][0]:.2f}-{iqrsDict2Fig3Asymp[0][1]:.2f})")
    print(f"Multiply-immunized: {medianDictFig3Asymp[1]:.2f} ({iqrsDict2Fig3Asymp[1][0]:.2f}-{iqrsDict2Fig3Asymp[1][1]:.2f})")
    return


@app.cell
def _(
    ANNOTATION_COORDS,
    ANNOTATION_LETTER_SIZE,
    annotateWithLetter,
    colWidth,
    dfFigureA2,
    label,
    order,
    pal,
    plotDir,
    plotFigA2_A,
    plt,
    saveFigure,
):
    figA2A, _ax = plt.subplots(1, 1, figsize=(colWidth * 2, colWidth * 1), constrained_layout=True)
    plotFigA2_A(dfFigureA2, palette=pal, order=order, label=label, xLabels=True, ax=_ax, markersize=2.3)
    annotateWithLetter(_ax, "B", size=ANNOTATION_LETTER_SIZE, coords=ANNOTATION_COORDS)
    saveFigure(figA2A, plotDir / "Figure2B-vl-thesis.pdf")
    saveFigure(figA2A, plotDir / "Figure2B-vl-thesis.png")
    _ax
    return


@app.cell
def _(IQR, IQRQuartiles, abbrvDictPaper, dfFigureA2, variantAbbrvs):
    print("Viral load medians (Q1, Q3)")
    medianDictFig3B = dict(dfFigureA2[dfFigureA2.symptoms==1].groupby("variant").vl.median().items())
    iqrsDictFig3B = dict(dfFigureA2[dfFigureA2.symptoms==1].groupby("variant").vl.apply(IQR))
    iqrsDict2Fig3B = dict(dfFigureA2[dfFigureA2.symptoms==1].groupby("variant").vl.apply(IQRQuartiles))

    medianDictFig3BAsymp = dict(dfFigureA2[dfFigureA2.symptoms==0].groupby("variant").vl.median().items())
    iqrsDictFig3BAsymp = dict(dfFigureA2[dfFigureA2.symptoms==0].groupby("variant").vl.apply(IQR))
    iqrsDict2Fig3BAsymp = dict(dfFigureA2[dfFigureA2.symptoms==0].groupby("variant").vl.apply(IQRQuartiles))


    i = 0
    _sympStatus = ("Symptomatic", "Asymptomatic")
    for _medianDict, _iqrsDict, in zip((medianDictFig3B, medianDictFig3BAsymp), 
                                       (iqrsDict2Fig3B, iqrsDict2Fig3BAsymp)):
        print()
        print(_sympStatus[i])
        i+=1
        for variant, variantDisplay in variantAbbrvs.items():
            print(f"{abbrvDictPaper["variant"][variant]}: "
                  f"{_medianDict[variant]:.2f} ({_iqrsDict[variant][0]:.2f}-{_iqrsDict[variant][1]:.2f})")
    return


@app.cell
def _(mo):
    mo.md(r"""## Viral load by Ag-RDT administration""")
    return


@app.cell
def _(
    ANNOTATION_COORDS,
    ANNOTATION_LETTER_SIZE,
    abbrvDictPaper,
    annotateWithLetters,
    colWidth,
    dfFigure2_B,
    dfFigure2_BAsymp,
    label,
    order,
    plotDir,
    plotFig2_B,
    plt,
    saveFigure,
):
    _fig2B, _ax = plt.subplots(2, 1, figsize=(colWidth * 2, colWidth * 2.4))
    # _fig2B.subplots_adjust(hspace=-0.1)
    plotFig2_B(dfFigure2_B, palette=None, label=label, order=order, abbrvDictPaper=abbrvDictPaper, ax=_ax[0])
    plotFig2_B(dfFigure2_BAsymp, palette=None, label=label, order=order, abbrvDictPaper=abbrvDictPaper, ax=_ax[1])
    # Remove x-label and xticklabels
    _ax[0].tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)
    _ax[0].set_xlabel("")       # Remove x-axis label
    annotateWithLetters(_ax, coords=ANNOTATION_COORDS, size=ANNOTATION_LETTER_SIZE)
    saveFigure(_fig2B,  plotDir / "Figure3-thesis.pdf")
    saveFigure(_fig2B, plotDir / "Figure3-thesis.png")
    _ax
    return


@app.cell
def _(dfFigure2_B):
    len(dfFigure2_B)
    return


@app.cell
def _(dfFigure2_BAsymp):
    len(dfFigure2_BAsymp)
    return


@app.cell
def _(mo):
    mo.md(r"""## Viral load by test result""")
    return


@app.cell
def _(mo):
    mo.md(r"""### Plots""")
    return


@app.cell
def _(
    ANNOTATION_LETTER_SIZE,
    annotateWithLetters,
    colWidth,
    dfAgrdtIndInf,
    label,
    plotDir,
    plotFigVlResult,
    plotFigVlTestline,
    plt,
    saveFigure,
):
    dfFigureAgrdtResult = dfAgrdtIndInf[dfAgrdtIndInf.vl > 3]

    _figVlResult, _ax = plt.subplots(2, 1, figsize=(colWidth*2, colWidth*1.8), constrained_layout=True)
    plotFigVlResult(dfFigureAgrdtResult, label, ax=_ax[0], markersize=1.75)
    plotFigVlTestline(dfFigureAgrdtResult.dropna(subset="testline"), label, ax=_ax[1], markersize=1.75)
    annotateWithLetters(_ax, size=ANNOTATION_LETTER_SIZE, coords=(-0.2, 1))
    saveFigure(_figVlResult, plotDir / "vl-result-thesis.png")
    _ax
    return


@app.cell
def _(mo):
    mo.md(r"""### Medians and IQRs""")
    return


@app.cell
def _(IQR, IQRQuartiles, dfAgrdtIndInf):
    medianDictTestResults = dict(dfAgrdtIndInf.groupby("testline").vl.median().items())
    iqrsDictTestResults = dict(dfAgrdtIndInf.groupby("testline").vl.apply(IQR))
    iqrsDict2TestResults = dict(dfAgrdtIndInf.groupby("testline").vl.apply(IQRQuartiles))
    return iqrsDict2TestResults, medianDictTestResults


@app.cell
def _(
    IQRQuartiles,
    dfAgrdtIndInf,
    iqrsDict2TestResults,
    medianDictTestResults,
):
    print("Test result")
    for _result in range(2):
        print(f"Test result: {_result}")
        print(f"Median: {dfAgrdtIndInf[dfAgrdtIndInf.agrdt==_result].vl.median():.2f}")
        _iqrsResult = IQRQuartiles(dfAgrdtIndInf[dfAgrdtIndInf.agrdt==_result].vl)
        print(f"IQR: {_iqrsResult[0]:.2f}-{_iqrsResult[1]:.2f}")
        print()

    print("Testline strength")
    for _testlineStrength in range(4):
        print(f"Testline strength: {_testlineStrength}")
        print(f"Median: {medianDictTestResults[_testlineStrength]:.2f}")
        print(f"IQR: {iqrsDict2TestResults[_testlineStrength][0]:.2f}-{iqrsDict2TestResults[_testlineStrength][1]:.2f}")
        print()
    return


@app.cell
def _(mo):
    mo.md(r"""# Regression""")
    return


@app.cell
def _(mo):
    mo.md(r"""## Viral load differences between time periods (Figure 1)""")
    return


@app.cell
def _(mo):
    mo.md(r"""### Symptomatic""")
    return


@app.cell
def _(bmb, dfFigure1_B, sampleVl):
    _priors = {"1|samplingMonth2": bmb.Prior("Normal", mu=0, sigma=bmb.Prior("HalfNormal", sigma=2))}
    modelVl, iDataVl, statsVl = sampleVl(df=dfFigure1_B, catVars=["samplingMonth2"], likelihood="skewnormal", formula="zVl ~ (1|samplingMonth2)", priors=_priors, target_accept=0.98)
    return iDataVl, modelVl, statsVl


@app.cell
def _(modelVl):
    modelVl
    return


@app.cell
def _(iDataDir, iDataVl):
    iDataVl.to_netcdf(iDataDir / "fig1BVl.nc")
    return


@app.cell
def _(az, iDataVl, month2ToMonths, statsVl):
    print("Difference in log10 viral loads between start of study (Dec 20/Jan 21) and remaining time periods.\n")
    _vlStd = statsVl["sd"]
    baseline = iDataVl.posterior.sel({"samplingMonth2__factor_dim": "0"})["1|samplingMonth2"].stack(samples=("chain", "draw")).values.flatten()
    for samplingTime, samplingTimeStr in month2ToMonths.items():
        if samplingTime == 3:
            continue
        print(samplingTimeStr)
        _samples = iDataVl.posterior.sel({"samplingMonth2__factor_dim": f"{int(samplingTime)}"})["1|samplingMonth2"].stack(samples=("chain", "draw")).values.flatten()
        samplesDiff = (_samples - baseline) * _vlStd
        samplingTimeDiffMean = samplesDiff.mean()
        samplingTimeDiffHDI = az.hdi(samplesDiff, 0.94)
        print(f"{samplingTimeStr}:")
        print(f"{samplingTimeDiffMean:.2f} ({samplingTimeDiffHDI[0]:.2f}, {samplingTimeDiffHDI[1]:.2f})")
        print()
    return


@app.cell
def _(mo):
    mo.md(r"""### Asymptomatic""")
    return


@app.cell
def _(bmb, dfFigure1_B_Asymp, sampleVl):
    _priors = {"1|samplingMonth2": bmb.Prior("Normal", mu=0, sigma=bmb.Prior("HalfNormal", sigma=2))}
    modelVlAsymp, iDataVlAsymp, statsVlAsymp = sampleVl(df=dfFigure1_B_Asymp, 
                                                        catVars=["samplingMonth2"], 
                                                        likelihood="skewnormal", 
                                                        formula="zVl ~ (1|samplingMonth2)", 
                                                        priors=_priors, 
                                                        target_accept=0.98)
    return iDataVlAsymp, statsVlAsymp


@app.cell
def _(iDataDir, iDataVlAsymp):
    iDataVlAsymp.to_netcdf(iDataDir / "fig1BVlAsymp.nc")
    return


@app.cell
def _(az, iDataVlAsymp, month2ToMonths, statsVlAsymp):
    print("Difference in log10 viral loads between start of study (Dec 20/Jan 21) and remaining time periods.\n")
    _vlStd = statsVlAsymp["sd"]
    baseline = iDataVlAsymp.posterior.sel({"samplingMonth2__factor_dim": "0"})["1|samplingMonth2"].stack(samples=("chain", "draw")).values.flatten()
    for samplingTime, samplingTimeStr in month2ToMonths.items():
        if samplingTime in (3, 4):
            continue
        print(samplingTimeStr)
        _samples = (iDataVlAsymp.posterior.sel({"samplingMonth2__factor_dim": f"{int(samplingTime)}"})
                    ["1|samplingMonth2"].stack(samples=("chain", "draw")).values.flatten())
        samplesDiff = (_samples - baseline) * _vlStd
        samplingTimeDiffMean = samplesDiff.mean()
        samplingTimeDiffHDI = az.hdi(samplesDiff, 0.94)
        print(f"{samplingTimeStr}:")
        print(f"{samplingTimeDiffMean:.2f} ({samplingTimeDiffHDI[0]:.2f}, {samplingTimeDiffHDI[1]:.2f})")
        print()
    return


@app.cell
def _(mo):
    mo.md(r"""## Viral load difference between negative and positive Ag-RDT and between different testline strengths""")
    return


@app.cell
def _(dfFigure1_B, sampleVl):
    modelVlResult, iDataVlResult, statsVlResult = sampleVl(df=dfFigure1_B, catVars=["agrdt"], likelihood="student-t", 
                                                           formula="zVl ~ agrdt", target_accept=0.98)
    return iDataVlResult, modelVlResult, statsVlResult


@app.cell
def _(modelVlResult):
    modelVlResult
    return


@app.cell
def _(statsVlResult):
    statsVlResult["sd"]
    return


@app.cell
def _(az, iDataVlResult, statsVlResult):
    print("Difference in log10 viral loads comparing samples with positive and negative Ag-RDT result.\n")
    _vlStd = statsVlResult["sd"]
    _samples = (iDataVlResult.posterior.sel({"agrdt_dim": "1"})["agrdt"].stack(samples=("chain", "draw")).values.flatten())

    _samples = _samples * _vlStd
    agrdtResultParamMean = _samples.mean()
    agrdtResultParamHDI = az.hdi(_samples, 0.94)
    print("Difference in viral loads:")
    print(f"{agrdtResultParamMean:.2f} ({agrdtResultParamHDI[0]:.2f}, {agrdtResultParamHDI[1]:.2f})")
    return


@app.cell
def _(dfFigure1_B, sampleVl):
    modelVlTestline, iDataVlTestline, statsVlTestline = sampleVl(df=dfFigure1_B, 
                                                                 catVars=["testline"], 
                                                                 likelihood="student-t", 
                                                                 formula="zVl ~ testline", 
                                                                 target_accept=0.98)
    return iDataVlTestline, modelVlTestline, statsVlTestline


@app.cell
def _(modelVlTestline):
    modelVlTestline
    return


@app.cell
def _(abbrvDictPaper, az, iDataVlTestline, statsVlTestline):
    print("Difference in log10 viral loads comparing samples with difference Ag-RDT testline strengths.\n")
    _vlStd = statsVlTestline["sd"]
    for _testline in range(1, 4):
        _samples = (iDataVlTestline.posterior.sel({"testline_dim": str(_testline)})
                    ["testline"].stack(samples=("chain", "draw")).values.flatten())

        _samples = _samples * _vlStd
        agrdtTestlineParamMean = _samples.mean()
        agrdtTestlineParamHDI = az.hdi(_samples, 0.94)
        print(f"Testline strength: {abbrvDictPaper["testline"][_testline]}")
        print(f"{agrdtTestlineParamMean:.2f} ({agrdtTestlineParamHDI[0]:.2f}, {agrdtTestlineParamHDI[1]:.2f})")
        print()
    return


@app.cell
def _(mo):
    mo.md(r"""## Viral load difference between people tested by Ag-RDT and tested only by RT-PCR""")
    return


@app.cell
def _(dfFigure2_B, dfFigure2_BAsymp):
    variantModels = {}
    variantIData = {}
    variantDfs = {}
    dfAgrdtYNVariant = dfFigure2_B.copy()
    dfAgrdtYNVariantAsymp = dfFigure2_BAsymp.copy()
    return dfAgrdtYNVariant, dfAgrdtYNVariantAsymp


@app.cell
def _(mo):
    mo.md(r"""### Symptomatic""")
    return


@app.cell
def _(abbrvDictPaper, dfAgrdtYNVariant, order):
    for _variant in order.variant:
        means = dfAgrdtYNVariant[dfAgrdtYNVariant.variant == _variant].groupby("agrdtYN").vl.mean()
        print(f"Mean log10 viral load for {abbrvDictPaper["variant"][_variant]} samples.\nNo Ag-RDT performed: {means[0]:.2f}, Ag-RDT performed: {means[1]:.2f}\n")
    return


@app.cell
def _(SEED, bmb, dfAgrdtYNVariant, sampleVl):
    _priors = {"variant": bmb.Prior("Normal", mu=0, sigma=2), "agrdtYN:variant": bmb.Prior("Normal", mu=0, sigma=2)}
    modelAgrdtYNVl, iDataAgrdtYNVl, statsAgrdtYNVl = sampleVl(dfAgrdtYNVariant, catVars=["agrdtYN", "variant"], likelihood="skewnormal", target_accept=0.95, priors=_priors, seed=SEED)
    return iDataAgrdtYNVl, modelAgrdtYNVl, statsAgrdtYNVl


@app.cell
def _(modelAgrdtYNVl):
    modelAgrdtYNVl
    return


@app.cell
def _(Path, iDataAgrdtYNVl, iDataDir):
    iDataAgrdtYNVl.to_netcdf(Path(iDataDir, "fig2BVl.nc"))
    return


@app.cell
def _(abbrvDictPaper, az, iDataAgrdtYNVl, order, statsAgrdtYNVl):
    print("Difference in log10 viral loads when comparing people who were not Ag-RDT tested vs who were.\n")
    _vlStd = statsAgrdtYNVl["sd"]
    for _variant in order.variant:
        _samples = (iDataAgrdtYNVl.posterior.sel({"agrdtYN:variant_dim": f"1, {_variant}"})
                    ["agrdtYN:variant"].stack(samples=("chain", "draw")).values.flatten())
        _samples = _samples * _vlStd
        agrdtYNParamMean = _samples.mean()
        agrdtYNParamHDI = az.hdi(_samples, 0.94)
        print(f"{abbrvDictPaper["variant"][_variant]}:")
        print(f"{agrdtYNParamMean:.2f} ({agrdtYNParamHDI[0]:.2f}, {agrdtYNParamHDI[1]:.2f})")
        # print(f"{abbrvDictPaper["variant"][_variant]}: {agrdtYNParamMean:.3f}")
        # print(f"94% HDI: {agrdtYNParamHDI[0]:.3f}, {agrdtYNParamHDI[1]:.3f}")
        print()
    return


@app.cell
def _(mo):
    mo.md(r"""### Asymptomatic""")
    return


@app.cell
def _(SEED, bmb, dfAgrdtYNVariantAsymp, sampleVl):
    _priors = {"variant": bmb.Prior("Normal", mu=0, sigma=2), "agrdtYN:variant": bmb.Prior("Normal", mu=0, sigma=2)}
    modelAgrdtYNVlAsymp, iDataAgrdtYNVlAsymp, statsAgrdtYNVlAsymp = sampleVl(dfAgrdtYNVariantAsymp, catVars=["agrdtYN", "variant"], likelihood="skewnormal", target_accept=0.95, priors=_priors, seed=SEED)
    return iDataAgrdtYNVlAsymp, modelAgrdtYNVlAsymp, statsAgrdtYNVlAsymp


@app.cell
def _(modelAgrdtYNVlAsymp):
    modelAgrdtYNVlAsymp
    return


@app.cell
def _(Path, iDataAgrdtYNVlAsymp, iDataDir):
    iDataAgrdtYNVlAsymp.to_netcdf(Path(iDataDir, "fig2BVlAsymp.nc"))
    return


@app.cell
def _(abbrvDictPaper, az, iDataAgrdtYNVlAsymp, order, statsAgrdtYNVlAsymp):
    print("Difference in log10 viral loads when comparing people who were not Ag-RDT tested vs who were (asymptomatic infections).\n")
    _vlStd = statsAgrdtYNVlAsymp["sd"]
    for _variant in order.variant:
        _samples = (iDataAgrdtYNVlAsymp.posterior.sel({"agrdtYN:variant_dim": f"1, {_variant}"})
                    ["agrdtYN:variant"].stack(samples=("chain", "draw")).values.flatten())
        _samples = _samples * _vlStd
        _agrdtYNParamMean = _samples.mean()
        _agrdtYNParamHDI = az.hdi(_samples, 0.94)
        print(f"{abbrvDictPaper["variant"][_variant]}:")
        print(f"{_agrdtYNParamMean:.2f} ({_agrdtYNParamHDI[0]:.2f}, {_agrdtYNParamHDI[1]:.2f})")
        print()
    return


@app.cell
def _(mo):
    mo.md(r"""## Viral load difference between immunized and immune naive people (Figure 2)""")
    return


@app.cell
def _(SEED, bmb, dfFigure3, sampleVl):
    _priors = {"symptoms": bmb.Prior("Normal", mu=0, sigma=2), "immun2YN:symptoms": bmb.Prior("Normal", mu=0, sigma=2)}
    modelImmunVlSkewed, iDataImmunVlSkewed, statsImmunVlSkewed = sampleVl(df=dfFigure3, catVars=["immun2YN", "symptoms"], likelihood="skewnormal", priors=_priors, target_accept=0.95, seed=SEED, interaction=True)
    return iDataImmunVlSkewed, modelImmunVlSkewed, statsImmunVlSkewed


@app.cell
def _(modelImmunVlSkewed):
    modelImmunVlSkewed
    return


@app.cell
def _(iDataDir, iDataImmunVlSkewed):
    iDataImmunVlSkewed.to_netcdf(iDataDir / "fig3AVl.nc")
    return


@app.cell
def _(abbrvDictPaper, az, iDataImmunVlSkewed, np, order, statsImmunVlSkewed):
    print("Mean log10 viral loads according to immunization status, depending on whether people had symptoms.\n")
    _vlStd = statsImmunVlSkewed["sd"]
    _vlMean = statsImmunVlSkewed["mean"]
    _zVl_sigma = iDataImmunVlSkewed.posterior["sigma"].stack(samples=("chain", "draw")).values.flatten()
    _zVl_alpha = iDataImmunVlSkewed.posterior["alpha"].stack(samples=("chain", "draw")).values.flatten()
    _commonTermMean = _zVl_sigma * np.sqrt(2 / np.pi) * (_zVl_alpha / np.sqrt(1 + np.power(_zVl_alpha, 2)))
    _intercept = iDataImmunVlSkewed.posterior["Intercept"].stack(samples=("chain", "draw")).values.flatten()
    _betaSymptoms = iDataImmunVlSkewed.posterior["symptoms"].stack(samples=("chain", "draw")).values.flatten()
    for _symptomStatus in order.symp:
        for immunized in (0, 1):
            print(f"{abbrvDictPaper["symptoms"][_symptomStatus]}, immunized: {abbrvDictPaper["immun2YN"][immunized]}")
            if immunized:
                _samples = (iDataImmunVlSkewed.posterior.sel({"immun2YN:symptoms_dim": f"1, {_symptomStatus}"})
                            ["immun2YN:symptoms"].stack(samples=("chain", "draw")).values)
            else:
                _samples = 0
            _samples = (_intercept + _betaSymptoms * _symptomStatus + _samples + _commonTermMean) * _vlStd + _vlMean
            _immun2YNParamMean = _samples.mean()
            _immun2YNParamHDI = az.hdi(_samples, 0.94)
            print(f"{_immun2YNParamMean:.2f} ({_immun2YNParamHDI[0]:.3f}, {_immun2YNParamHDI[1]:.3f})")
            print()
    return


@app.cell
def _(az, iDataImmunVlSkewed, statsImmunVlSkewed):
    print("Difference in log10 viral loads, depending on whether people had symptoms.\n")
    _vlStd = statsImmunVlSkewed["sd"]
    for _symptomStatus in (0, 1):
        _samples = (iDataImmunVlSkewed.posterior.sel({"immun2YN:symptoms_dim": f"1, {_symptomStatus}"})
                    ["immun2YN:symptoms"].stack(samples=("chain", "draw")).values)
        _samples = _samples * _vlStd
        _immun2YNParamMean = _samples.mean()
        _immun2YNParamHDI = az.hdi(_samples, 0.94)
        sympString = "Symptomatic" if _symptomStatus else "Asymptomatic"
        print(sympString)
        print(f"{_immun2YNParamMean:.2f} ({_immun2YNParamHDI[0]:.2f}, {_immun2YNParamHDI[1]:.2f})")
        print()
    return


@app.cell
def _(mo):
    mo.md(r"""## Viral load difference between SARS-CoV-2 variants""")
    return


@app.cell
def _(SEED, bmb, dfFigureA2, order, pd, sampleVl):
    dfVariantVl = dfFigureA2.dropna(subset=["variant"]).copy()
    _priors = {"symptoms": bmb.Prior("Normal", mu=0, sigma=2), "variantCode:symptoms": bmb.Prior("Normal", mu=0, sigma=2)}
    dfVariantVl["variantCode"] = pd.Categorical(dfVariantVl.variant, categories=order.variant).codes
    modelVariantVlSkewed, iDataVariantVlSkewed, statsVariantVlSkewed = sampleVl(df=dfVariantVl, catVars=["variantCode", 
                                                                                                         "symptoms"],
                                                                                likelihood="skewnormal", target_accept=0.95,
                                                                                priors=_priors, seed=SEED)
    return iDataVariantVlSkewed, modelVariantVlSkewed, statsVariantVlSkewed


@app.cell
def _(modelVariantVlSkewed):
    modelVariantVlSkewed
    return


@app.cell
def _(iDataDir, iDataVariantVlSkewed):
    iDataVariantVlSkewed.to_netcdf(iDataDir / "figA2AVl.nc")
    return


@app.cell
def _(abbrvDictPaper, az, iDataVariantVlSkewed, order, statsVariantVlSkewed):
    print("Difference in log10 viral loads (with wildtype being the baseline), depending on whether people had symptoms.\n")
    _vlStd = statsVariantVlSkewed["sd"]
    for _variantCode, _variantStr in zip(range(1, 4), order.variant[1:]):
        for _symptomStatus in order.symp:
            print(abbrvDictPaper["symptoms"][_symptomStatus])
            _samples = (iDataVariantVlSkewed.posterior.sel({"variantCode:symptoms_dim": f"{_variantCode}, {_symptomStatus}"})
                        ["variantCode:symptoms"].stack(samples=("chain", "draw")).values.flatten())
            _samples = _samples * _vlStd
            _variantParamMean = _samples.mean()
            _variantParamHDI = az.hdi(_samples, 0.94)
            print(f"{abbrvDictPaper["variant"][_variantStr]}:")
            print(f"{_variantParamMean:.2f} ({_variantParamHDI[0]:.2f}, {_variantParamHDI[1]:.2f})")
            print()
        print()
    return


@app.cell
def _(mo):
    mo.md(r"""# Correlations between viral load metrics (median and skewness) and different variables""")
    return


@app.cell
def _(mo):
    mo.md(r"""## Prepare data""")
    return


@app.cell
def _(cusmaMaxDate, cusmaMinDate, dfAllFirstPosPcrsNoRelease, pd):
    # # Prepare the data
    dfCusmaData = dfAllFirstPosPcrsNoRelease.copy()
    dfCusmaData = dfCusmaData[(dfCusmaData.pcrDate >= cusmaMinDate.date()) & (dfCusmaData.pcrDate <= cusmaMaxDate.date())]
    dfCusmaData["pcrDate"] = pd.to_datetime(dfCusmaData.pcrDate)
    dfCusmaData["dateNumerical"] = (dfCusmaData.pcrDate - dfCusmaData.pcrDate.min()).dt.days
    dfCusmaData["week"] = pd.to_datetime(dfCusmaData["pcrDate"]).dt.to_period(freq="W-MON").dt.start_time
    dfCusmaData["biweek"] = pd.to_datetime(dfCusmaData["pcrDate"]).dt.to_period(freq="2W-MON").dt.start_time
    dfCusmaData["month"] = pd.to_datetime(dfCusmaData["pcrDate"]).dt.to_period(freq="M").dt.start_time

    dfCusmaDataAsymp = dfCusmaData[dfCusmaData.symptoms==0].copy()
    dfCusmaDataSymp = dfCusmaData[dfCusmaData.symptoms==1].copy()
    dfCusmaDataOutbreak = dfCusmaData[dfCusmaData.reasonPres == "outbreak"].copy()
    dfCusmaDataNoCat = dfCusmaData[dfCusmaData.reasonPres == "no category"].copy()
    vlMediansCusmaData = dfCusmaData.groupby("pcrDate").vl.median()
    vlMediansCusmaData = vlMediansCusmaData.sort_index()
    vlMediansCusmaDataRolling = vlMediansCusmaData.rolling(window=14, center=True).median()

    vlMediansCusmaDataSymp = dfCusmaDataSymp.groupby("pcrDate").vl.median()
    vlMediansCusmaDataSymp = vlMediansCusmaDataSymp.sort_index()
    vlMediansCusmaDataSympRolling = vlMediansCusmaDataSymp.rolling(window=14, center=True).median()

    vlMediansCusmaDataAsymp = dfCusmaDataAsymp.groupby("pcrDate").vl.median()
    vlMediansCusmaDataAsymp = vlMediansCusmaDataAsymp.sort_index()
    vlMediansCusmaDataAsympRolling = vlMediansCusmaDataAsymp.rolling(window=14, center=True).median()

    vlMediansCusmaDataOutbreak = dfCusmaDataOutbreak.groupby("pcrDate").vl.median()
    vlMediansCusmaDataOutbreak = vlMediansCusmaDataOutbreak.sort_index()
    vlMediansCusmaDataOutbreakRolling = vlMediansCusmaDataOutbreak.rolling(window=14, center=True).median()

    vlMediansCusmaDataNoCat = dfCusmaDataNoCat.groupby("pcrDate").vl.median()
    vlMediansCusmaDataNoCat = vlMediansCusmaDataNoCat.sort_index()
    vlMediansCusmaDataNoCatRolling = vlMediansCusmaDataNoCat.rolling(window=14, center=True).median()
    return (
        dfCusmaData,
        dfCusmaDataAsymp,
        dfCusmaDataSymp,
        vlMediansCusmaDataAsympRolling,
        vlMediansCusmaDataSympRolling,
    )


@app.cell
def _(mo):
    mo.md(r"""## Plots""")
    return


@app.cell
def _(colWidth, dfCusmaData, pd, plt, variantAbbrvs):
    # Group the data by variant
    # Step 1: Count occurrences
    dfPies = dfCusmaData.copy()
    dfPies["variant"] = pd.Categorical(dfPies["variant"], categories=variantAbbrvs, ordered=True)
    counts = dfPies.groupby(["reasonPres", "variant"]).size().reset_index(name="count")

    # Step 2: Calculate total count per variant
    variant_totals = counts.groupby("variant")["count"].transform("sum")

    # Step 3: Calculate percentage
    counts["percentage"] = (counts["count"] / variant_totals) * 100

    grouped = counts.groupby("variant")

    _fig, _ax = plt.subplots(2, 2, figsize=(colWidth * 2, colWidth * 2))
    _axes = _ax.flatten()
    # Create a pie plot for each variant
    for i, (variant, group) in enumerate(grouped):
        _axes[i].pie(
            group["count"],
            labels=group["reasonPres"],
            autopct="%1.1f%%",
            startangle=90,
            labeldistance=0.7,
            pctdistance=0.3,
            colors=plt.cm.Pastel1.colors,
        )
        _axes[i].set_title(variantAbbrvs[variant])
    plt.tight_layout()
    _fig
    return


@app.cell
def _(mo):
    mo.md(r"""#### Temporal patterns""")
    return


@app.cell
def _(
    ANNOTATION_LETTER_SIZE,
    GGPLOT_PALETTE,
    annotateWithLetters,
    colWidth,
    dfCusmaData,
    dfRt,
    label,
    legend,
    pal,
    plotDir,
    plt,
    replaceLegend,
    saveFigure,
    setFontSize,
    sns,
    ticker,
    vlMediansCusmaDataAsympRolling,
    vlMediansCusmaDataSympRolling,
):
    dfRt["dateNumerical"] = (dfRt.date - dfRt.date.min()).dt.days
    dfRt["displayDate"] = dfRt.date.dt.strftime("%Y-%m")

    # plot the results
    _fig, _ax = plt.subplots(2, 1, figsize=(colWidth*2, colWidth*1.7), sharex=True, 
                             gridspec_kw={"height_ratios": [0.55, 0.45]}, constrained_layout=True)
    _ax0_2 = _ax[0].twinx()
    _ax1_2 = _ax[1].twinx()

    thin_space = "\u2009"  # narrow no-break space
    _ax1_2.yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda x, pos: f"{int(x):,}".replace(",", thin_space))
    )

    # Plot R value data
    sns.scatterplot(dfCusmaData, x="pcrDate", y="vl", ax=_ax[0], alpha=0.2, s=5, hue="symptoms", palette=pal.symp)
    sns.lineplot(x=vlMediansCusmaDataAsympRolling.index, y=vlMediansCusmaDataAsympRolling, ax=_ax[0], color=pal.symp[0])
    sns.lineplot(x=vlMediansCusmaDataSympRolling.index, y=vlMediansCusmaDataSympRolling, ax=_ax[0], color=pal.symp[1])

    _ax0_2.plot(dfRt.date, dfRt.rt, color="blue", linewidth=1, alpha=0.55)
    _ax0_2.plot(dfRt.date, dfRt.rtRolling4, color="blue", linewidth=1.5)

    _ax0_2.set_ylim(0.327, 1.63)
    _ax0_2.grid(False)

    replaceLegend(_ax[0], handles=legend.sympAsympRt, loc="lower right")

    # Plot case count data
    countsSamples = dfCusmaData.value_counts("pcrDate")
    countsSamples = countsSamples.sort_index()
    countsSamplesRolling = countsSamples.rolling(window=14, center=True).mean()
    _ax[1].plot(countsSamples.index, countsSamplesRolling, linewidth=1.5)
    replaceLegend(_ax[1], handles=legend.newCases, loc="upper center")

    _ax1_2.plot(dfRt.date, dfRt.newCasesRolling2, color=GGPLOT_PALETTE[2], linewidth=1.5)
    _ax1_2.grid(False)
    _ax1_2.set_ylim(-12000, 180000)

    # axis labels
    _ax[0].set_ylabel(f"{label.vl}\n(Charité cohort)", labelpad=1)
    _ax[1].set_ylabel("Number of new cases\n(Charité cohort)", labelpad=10)
    _ax0_2.set_ylabel(r"7-day $R_t$ value (Germany)", rotation=270, labelpad=20)
    _ax1_2.set_ylabel("Number of new cases\n(Germany)", rotation=270, labelpad=30)
    _ax[1].set_xlabel("Time")
    _ax[1].patch.set_visible(False)
    _ax[1].set_xticklabels(_ax[1].get_xticklabels(), rotation=45)

    for a in (_ax[0], _ax[1], _ax0_2, _ax1_2):
        setFontSize(a)

    plt.tight_layout()
    annotateWithLetters(_ax, size=ANNOTATION_LETTER_SIZE, coords=(-0.31, 1))

    saveFigure(_fig, plotDir / "vlVsRt.png", dpi=500)
    _fig
    return countsSamples, countsSamplesRolling


@app.cell
def _(mo):
    mo.md(r"""#### Viral load medians and skewness""")
    return


@app.cell
def _(
    ANNOTATION_COORDS,
    ANNOTATION_LETTER_SIZE,
    annotateWithLetters,
    colWidth,
    dfCusmaDataAsymp,
    dfCusmaDataSymp,
    dfRt,
    label,
    legend,
    mpl,
    pal,
    pd,
    plotDir,
    plt,
    replaceLegend,
    saveFigure,
    setFontSize,
    skew,
    sns,
):
    # Plot skewness vs vl
    colorBy = "variant"#"rt", "immun2YN"
    thresh = 20
    timeframe = "week"

    _fig, _ax = plt.subplots(2, 1, figsize=(colWidth*2, colWidth*2), sharey=True, sharex=True, constrained_layout=True)
    dfRt["week"] = (pd.to_datetime(dfRt["date"]).dt.to_period(freq="W-MON")).dt.start_time
    dfRt["biweek"] = (pd.to_datetime(dfRt["date"]).dt.to_period(freq="2W-MON")).dt.start_time
    dfRt["month"] = (pd.to_datetime(dfRt["date"]).dt.to_period(freq="M")).dt.start_time

    countsPeriodAsymp = dfCusmaDataAsymp[timeframe].value_counts()
    countsPeriodSymp = dfCusmaDataSymp[timeframe].value_counts()
    validPeriodAsymp = countsPeriodAsymp[countsPeriodAsymp >= thresh].index
    validPeriodSymp = countsPeriodSymp[countsPeriodSymp >= thresh].index

    mediansVlAsymp = dfCusmaDataAsymp.groupby(timeframe).vl.median()
    mediansVlSymp = dfCusmaDataSymp.groupby(timeframe).vl.median()

    countsMedianAsymp = dfCusmaDataAsymp.groupby(timeframe).vl.count()
    countsMedianSymp = dfCusmaDataSymp.groupby(timeframe).vl.count()

    mediansSkewnessAsymp = dfCusmaDataAsymp.groupby(timeframe).vl.apply(skew)
    mediansSkewnessSymp = dfCusmaDataSymp.groupby(timeframe).vl.apply(skew)

    mediansRt = dfRt.groupby(timeframe).rt.median()

    majorityVariantAsymp = dfCusmaDataAsymp.groupby(timeframe).variant.agg(lambda x: x.mode())
    majorityVariantSymp = dfCusmaDataSymp.groupby(timeframe).variant.agg(lambda x: x.mode())

    immun2YNAsymp = dfCusmaDataAsymp.groupby(timeframe).immun2YN.mean()
    immun2YNSymp = dfCusmaDataSymp.groupby(timeframe).immun2YN.mean()             


    dfPeriodAsymp = pd.DataFrame({"vl": mediansVlAsymp, "skew": mediansSkewnessAsymp, "rt": mediansRt, 
                                 "count": countsMedianAsymp, "variant": majorityVariantAsymp, 
                                 "immun2YN": immun2YNAsymp})
    dfPeriodSymp = pd.DataFrame({"vl": mediansVlSymp, "skew": mediansSkewnessSymp, "rt": mediansRt, 
                               "count": countsMedianSymp, "variant": majorityVariantSymp, 
                               "immun2YN": immun2YNSymp})
    dfPeriodAsymp = dfPeriodAsymp[dfPeriodAsymp.index.isin(validPeriodAsymp)]
    dfPeriodSymp = dfPeriodSymp[dfPeriodSymp.index.isin(validPeriodSymp)]

    if colorBy == "variant":
        sns.scatterplot(dfPeriodAsymp, y="vl", x="skew", hue="variant", palette=pal.variant, legend=False, size="count",
                        ax=_ax[1], alpha=1)
        sns.scatterplot(dfPeriodSymp, y="vl", x="skew", hue="variant", palette=pal.variant, legend=False, size="count",
                        ax=_ax[0], alpha=1)
        replaceLegend(_ax[0], legend.variantPoints2, loc="lower left")
    elif colorBy == "rt":
        cmap= sns.color_palette("Reds", as_cmap= True)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=mpl.colors.Normalize())
        sns.scatterplot(dfPeriodAsymp, y="vl", x="skew", hue="rt", hue_norm=sm.norm, palette=cmap, legend=False, 
                        size="count", ax=_ax[1], alpha=1)
        sns.scatterplot(dfPeriodSymp, y="vl", x="skew", hue="rt", hue_norm=sm.norm, palette=cmap, legend=False, size="count",
                        ax=_ax[0], alpha=1)
    elif colorBy == "immun2YN":
        cmap= sns.color_palette("viridis", as_cmap= True)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=mpl.colors.Normalize())
        sns.scatterplot(dfPeriodAsymp, y="vl", x="skew", hue="immun2YN",  hue_norm=sm.norm, palette=cmap,
                        legend=False, size="count", ax=_ax[1], alpha=1)
        sns.scatterplot(dfPeriodSymp, y="vl", x="skew", hue="immun2YN", hue_norm=sm.norm, palette=cmap, legend=False, 
                        size="count", ax=_ax[0], alpha=1)

    for _a in _ax:
        _a.set_ylabel(f"Median {label.vl.lower()}")
        _a.set_xlabel("Skewness")
        if colorBy in ("rt", "immun2YN"):
            cbar = plt.colorbar(sm, ax=_a)
            if colorBy == "rt":
                cbar.set_label("Median 7-day $R_t$ value", rotation=270, labelpad=15)
            elif colorBy == "immun2YN":
                cbar.set_label("Multiple immunization rate", rotation=270, labelpad=15)
    setFontSize(_ax[0])
    setFontSize(_ax[1])
    annotateWithLetters(_ax, size=ANNOTATION_LETTER_SIZE, coords=ANNOTATION_COORDS)
    saveFigure(_fig, plotDir / f"vlVsSkewBy{colorBy.capitalize()}{timeframe.capitalize()}.png")
    _fig
    return


@app.cell
def _(countsSamples, countsSamplesRolling, plt):
    _fig, _ax = plt.subplots(1, 1, figsize=(10, 5))
    _ax.plot(countsSamples.index, countsSamplesRolling, linewidth=1.5)
    _ax.vlines(x=85, ymin=0, ymax=11, color="black")
    _ax.vlines(x=155, ymin=0, ymax=11, color="black")
    _ax.vlines()
    return


@app.cell
def _(mo):
    mo.md(r"""## Compute skewness values""")
    return


@app.cell
def _(bmb, dfCusmaData, dfRt, pd, skew):
    dfCusmaDataCorr = dfCusmaData.copy()
    dfRtCorr = dfRt.copy()
    iDataCorrDict = {}
    corrIndexDict = {}
    for freq in ("W-MON", "W-TUE", "W-WED", "W-THU", "W-FRI", "W-SAT", "W-SUN"):
        dfCusmaDataCorr["week"] = (pd.to_datetime(dfCusmaDataCorr["pcrDate"]).dt.to_period(freq=freq))
        dfRtCorr["week"] = (pd.to_datetime(dfRtCorr["date"]).dt.to_period(freq=freq))
        countsWeeksCorr = dfCusmaDataCorr["week"].value_counts()
        validWeeksCorr = countsWeeksCorr[countsWeeksCorr >= 15].index
        mediansVlCorr = dfCusmaDataCorr.groupby("week").vl.median()
        mediansSkewnessCorr = dfCusmaDataCorr.groupby("week").vl.apply(skew)
        mediansRtCorr = dfRtCorr.groupby("week").rt.median()

        dfWeekCorr = pd.DataFrame({"vl": mediansVlCorr, "skew": mediansSkewnessCorr, "rt": mediansRtCorr})
        dfWeekCorr = dfWeekCorr[dfWeekCorr.index.isin(validWeeksCorr)]

        iDataCorrDict[freq] = bmb.Model(formula="rt ~ vl", data=dfWeekCorr).fit()
        corrIndexDict[freq] = dfWeekCorr.rt.corr(dfWeekCorr.vl)
    return corrIndexDict, iDataCorrDict


@app.cell
def _(az, corrIndexDict, iDataCorrDict):
    for _freq in corrIndexDict:
        print(_freq)
        print(f"Correlation: {corrIndexDict[_freq]}")
        print("Parameter estimate vl (mean, 3%, 97%)")
        print(az.summary(iDataCorrDict[_freq])[["mean", "hdi_3%", "hdi_97%"]].loc["vl"].values)
        print()
    return


if __name__ == "__main__":
    app.run()
