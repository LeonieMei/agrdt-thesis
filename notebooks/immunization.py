import marimo

__generated_with = "0.16.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import matplotlib.pyplot as plt

    from agrdt.data import createDataFramesFigures
    from agrdt.plotting import (annotateWithLetters, returnPlotDirImmunization,
                                setCustomTheme, saveFigure)
    from agrdt.plotParams import getLabels, getLegends, getOrders, getAbbrvsDict, COL_WIDTH, CM, ANNOTATION_LETTER_SIZE
    from agrdt.regression import returnIDataDirImmunization
    from agrdt.dataParams import ROOT_DIR

    setCustomTheme()
    return (
        ANNOTATION_LETTER_SIZE,
        CM,
        COL_WIDTH,
        ROOT_DIR,
        annotateWithLetters,
        createDataFramesFigures,
        getAbbrvsDict,
        getLabels,
        getLegends,
        getOrders,
        mo,
        pd,
        plt,
        returnIDataDirImmunization,
        returnPlotDirImmunization,
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
    getLegends,
    getOrders,
    getPalettes,
    returnIDataDirImmunization,
    returnPlotDirImmunization,
):
    pal = getPalettes()
    legend = getLegends()
    order = getOrders()
    label = getLabels()
    plotDir = returnPlotDirImmunization()
    iDataDir = returnIDataDirImmunization()
    return label, legend, order, pal, plotDir


@app.cell
def _(getAbbrvsDict):
    abbrvDictPaper = getAbbrvsDict()

    genderAbbrvs = abbrvDictPaper["gender"]
    symptomAbbrvs = abbrvDictPaper["symptoms"]
    variantAbbrvs = abbrvDictPaper["variant"]
    immunAbbrvs = abbrvDictPaper["immun2YN"]
    daysAbbrvs = abbrvDictPaper["binDaysPostOnset"]
    variants = list(variantAbbrvs.keys())
    return (abbrvDictPaper,)


@app.cell
def _(CM, COL_WIDTH):
    cm = CM
    colWidth = COL_WIDTH
    fontSizePlot = 12
    fontSizePlotSuppl = 10
    return (colWidth,)


@app.cell
def _(abbrvDictPaper):
    month2Labels = list(abbrvDictPaper['samplingMonth2'].values())
    return (month2Labels,)


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
    return dfFigure1_C, dfFigure1_C_Asymp, dfPos


@app.cell
def _(mo):
    mo.md(r"""# Plots""")
    return


@app.cell
def _(
    ANNOTATION_LETTER_SIZE,
    annotateWithLetters,
    colWidth,
    dfFigure1_C,
    dfFigure1_C_Asymp,
    dfPos,
    label,
    legend,
    month2Labels,
    order,
    pal,
    plotDir,
    plotFig1_C,
    plt,
    saveFigure,
):
    _fig1C, _ax = plt.subplots(3, 1, figsize=(colWidth * 2, colWidth * 2.1), gridspec_kw={'height_ratios': [0.02, 1, 1]},)
                               # constrained_layout=True)
    # Ensure original figure size is kept
    _fig1C.subplots_adjust(
        top=0.9,
        bottom=0.1,
        left=0.15,
        right=0.88,
        hspace=0.3
    )
    plotFig1_C(dfFigure1_C, dfPos, timeVar='samplingMonth2', xOrder=order.samplingMonth2, palette=pal, label=label,
               legend=legend, xLabels=None, ax=_ax[:2], legendLoc=(0.07, 0.6))
    plotFig1_C(dfFigure1_C_Asymp, dfPos, timeVar='samplingMonth2', xOrder=order.samplingMonth2, palette=pal, label=label,
               legend=None, xLabels=month2Labels, showVariantTimes=False, ax=_ax[2], legendLoc=(0.07, 0.6))
    annotateWithLetters(_ax[1:], size=ANNOTATION_LETTER_SIZE, coords=(-0.2, 1))
    saveFigure(_fig1C, plotDir / 'Figure1-C-thesis.png')
    saveFigure(_fig1C, plotDir / 'Figure1-C-thesis.png')
    _ax
    return


if __name__ == "__main__":
    app.run()
