import marimo

__generated_with = "0.16.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import numpy as np
    import bambi as bmb
    import arviz as az
    import warnings
    import matplotlib.pyplot as plt
    from collections import defaultdict

    from utils.plotUtils import setCustomTheme, saveFigure, ridgeForestPlot, spaghettiPlotCategorical, returnPlotDirRegression, DPI
    from utils.plotParams import CM, COL_WIDTH
    from utils.regression import predictionsNewData, SEED
    from utils.dataParams import ROOT_DIR
    from utils.tableUtils import writeIDataSummaryTable, writeIDataSummaryTableLatex, returnTableDirIData
    return (
        CM,
        COL_WIDTH,
        DPI,
        ROOT_DIR,
        SEED,
        az,
        bmb,
        defaultdict,
        mo,
        np,
        pd,
        plt,
        predictionsNewData,
        returnPlotDirRegression,
        returnTableDirIData,
        ridgeForestPlot,
        saveFigure,
        setCustomTheme,
        spaghettiPlotCategorical,
        warnings,
        writeIDataSummaryTable,
        writeIDataSummaryTableLatex,
    )


@app.cell
def _(CM, COL_WIDTH, setCustomTheme, warnings):
    setCustomTheme()
    warnings.filterwarnings("ignore")

    cm = CM
    colWidth = COL_WIDTH
    fontSizePlot = 12
    fontSizePlotSuppl = 10
    return (colWidth,)


@app.cell
def _(returnPlotDirRegression, returnTableDirIData):
    plotDir = returnPlotDirRegression()
    tableDirIData = returnTableDirIData()
    return plotDir, tableDirIData


@app.cell
def _(mo):
    mo.md(r"""# Utils""")
    return


@app.cell
def _(pd):
    def fillRNACol(grp):
        nonNAGrp = grp[grp["RNA"].notna()]
        nDilutionsTotal = grp["Dilution"].max() - grp["Dilution"].min() + 1
        nDilutions = max(nonNAGrp["Dilution"]) - min(nonNAGrp["Dilution"])
        RNAmax, RNAmin = max(nonNAGrp["RNA"]), min(nonNAGrp["RNA"])
        dilutionFactor = (RNAmax / RNAmin)**(1/nDilutions) # formula: RNAmin * dilutionFactor^dilutionStep = RNAmax
        RNAimputed = []
        for dilution in range(nDilutionsTotal)[::-1]:
            RNAimputed.append(RNAmin * dilutionFactor**dilution)

        return pd.Series(RNAimputed, index=grp.index)
    return (fillRNACol,)


@app.cell
def _(mo):
    mo.md(r"""# Load data""")
    return


@app.cell
def _(ROOT_DIR, fillRNACol, pd):
    # Note: this file contains only data from experiments where RNA measurements are available
    dataPath = ROOT_DIR / "data" / "rapidTestVariantsSummary.tsv"
    df = pd.read_csv(dataPath, sep="\t")
    df["RNA"] = df.RNA.str.replace(",", ".").astype(float)
    df["testline"] = df.Result.replace({"(1)": 0}).astype(int) # handling very weak line
    df["agrdt"] = df.testline > 0
    df["RNAFull"] = df.groupby(["Variant", "Experiment_no"], group_keys=False).apply(fillRNACol, include_groups=False)
    return (df,)


@app.cell
def _(mo):
    mo.md(r"""# Regression""")
    return


@app.cell
def _(SEED, bmb, df, np, pd):
    # Drop duplicates from repeating application on rapid test devices with the same diluted stock solution
    dfBmb = df[~(df.Variant.isin(("BA.1", "BA.2")) & (df.Experiment_no.isin((2, 3))))].copy()
    dfBmb["VariantCode"] = pd.Categorical(dfBmb.Variant, ordered=True, categories=["WT", "Delta", "BA.1", "BA.2"]).codes
    dfBmb["log10Load"] = np.log10(dfBmb.RNAFull)
    priors = {"VariantCode": bmb.Prior("Normal", mu=0, sigma=1)}
    formula = "agrdt ~ VariantCode + log10Load"
    model = bmb.Model(data=dfBmb, formula=formula, categorical=["VariantCode"], 
                      priors=priors, family="bernoulli")
    iData = model.fit(1000, 1000, target_accept=0.8, seed=SEED)
    return dfBmb, iData, model


@app.cell
def _(model):
    model
    return


@app.cell
def _(az, iData):
    summaryDf = az.summary(iData, var_names="~p")
    summaryDf
    return (summaryDf,)


@app.cell
def _(colWidth, iData, plotDir, plt, ridgeForestPlot, saveFigure):
    _fig, ax = plt.subplots(1, 1, figsize=(colWidth * 1.5, colWidth * 1))
    ridgeForestPlot(iData, varNames=["VariantCode", "log10Load"], addGrid=True, ridge=False, ax=ax)
    saveFigure(_fig, plotDir / "analyticSensVariantsForestPlot.png")
    _fig
    return


@app.cell
def _(defaultdict, dfBmb, iData, model, predictionsNewData):
    predDict = defaultdict(dict)
    predictionsNewData(dfBmb, model, iData, predDict, "VariantCode", (0, 1, 2, 3), vlVarName="log10Load", vlRange=(0, 11))
    return (predDict,)


@app.cell
def _(
    summaryDf,
    tableDirIData,
    writeIDataSummaryTable,
    writeIDataSummaryTableLatex,
):
    writeIDataSummaryTable(summaryDf=summaryDf, varnames=("VariantCode", "log10Load"), outfile=tableDirIData / "iDataSummaryAnalyticSens.tsv")
    writeIDataSummaryTableLatex(summaryDf=summaryDf, varnames=("VariantCode", "log10Load"), outfile=tableDirIData / "iDataSummaryAnalyticSensLatex.txt")
    return


@app.cell
def _(
    DPI,
    colWidth,
    dfBmb,
    plotDir,
    plt,
    predDict,
    saveFigure,
    spaghettiPlotCategorical,
):
    _fig, _ax = plt.subplots(1, 1, figsize=(colWidth * 2, colWidth * 1.15), constrained_layout=True)
    spaghettiPlotCategorical(_ax, "VariantCode", (0, 1, 2, 3), predDict, dfBmb, showXLabel=True, 
                             legendLoc=(0.6, 0.15), vlVarName="log10Load", vlRange=(0, 11))
    _ax
    saveFigure(_fig, plotDir / "analyticSensVariantsSpaghPlot.png", dpi=DPI)
    _fig
    return


@app.cell
def _(colWidth):
    colWidth * 2
    return


if __name__ == "__main__":
    app.run()
