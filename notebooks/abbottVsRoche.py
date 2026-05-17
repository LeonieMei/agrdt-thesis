import marimo

__generated_with = "0.16.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import xarray
    import marimo as mo

    # Bayesian analysis
    import arviz as az
    import bambi as bmb

    from os import path
    from collections import defaultdict

    from utils.dataUtils import standardize
    from utils.plotUtils import (ridgeForestPlot, setFontSize, boldStr, spaghettiPlotCategorical, returnPlotDirRegression, 
                                 setCustomTheme, saveFigure)
    from utils.plotParams import getLabels, getPalettes, getLegends, COL_WIDTH, CM
    from utils.regression import SEED, returnIDataDirRegression, postProcessAbbottRoche, postProcessAbbottRocheHierarchical
    from utils.dataParams import ROOT_DIR
    from utils.tableUtils import writeIDataSummaryTable, writeIDataSummaryTableLatex, returnTableDirIData

    setCustomTheme()
    return (
        CM,
        COL_WIDTH,
        ROOT_DIR,
        SEED,
        az,
        bmb,
        boldStr,
        defaultdict,
        getLabels,
        getLegends,
        getPalettes,
        mo,
        np,
        path,
        pd,
        plt,
        postProcessAbbottRoche,
        postProcessAbbottRocheHierarchical,
        returnIDataDirRegression,
        returnPlotDirRegression,
        returnTableDirIData,
        ridgeForestPlot,
        saveFigure,
        setFontSize,
        sns,
        spaghettiPlotCategorical,
        standardize,
        writeIDataSummaryTable,
        writeIDataSummaryTableLatex,
        xarray,
    )


@app.cell
def _(ROOT_DIR, path, pd, returnTableDirIData):
    tableDirIData = returnTableDirIData()
    DATA_DIR = ROOT_DIR / "data"
    DATA_FILE1 = path.join(DATA_DIR, "abbottVsRocheWildtype.tsv")
    DATA_FILE2 = path.join(DATA_DIR, "abbottVsRocheOmicron.tsv")

    df1 = pd.read_csv(DATA_FILE1, sep="\t")
    df2 = pd.read_csv(DATA_FILE2, sep="\t")
    return df1, df2, tableDirIData


@app.cell
def _(
    getLabels,
    getLegends,
    getPalettes,
    returnIDataDirRegression,
    returnPlotDirRegression,
):
    pal = getPalettes()
    label = getLabels()
    legend = getLegends()
    plotDir = returnPlotDirRegression()
    iDataDir= returnIDataDirRegression()
    return iDataDir, label, legend, pal, plotDir


@app.cell
def _(label, sns):
    pal1 = sns.color_palette("colorblind")
    labelVl = label.vl
    labelRes = "AgPOCT result"
    return labelRes, labelVl, pal1


@app.cell
def _(CM, COL_WIDTH):
    cm = CM
    colWidth = COL_WIDTH
    fontSizePlot = 12
    fontSizePlotSuppl = 10
    return colWidth, fontSizePlot, fontSizePlotSuppl


@app.cell
def _(df1, df2, pd, standardize):
    df = pd.concat([df1, df2])
    standardize(df, "vl", "zVl")
    return (df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Plots""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Pre-VOC samples""")
    return


@app.cell
def _(df1, plt, sns):
    _fig, _ax = plt.subplots(1, 2, figsize=(15, 7), sharey=True)
    df1Abbott = df1[df1.test == 0].copy()
    sns.countplot(data=df1Abbott, x='agrdt', ax=_ax[0])
    _ax[0].set_title('Abbott')
    df1Roche = df1[df1.test == 1].copy()
    sns.countplot(data=df1Roche, x='agrdt', ax=_ax[1])
    _ax[1].set_title('Roche')
    return df1Abbott, df1Roche


@app.cell
def _(df1Abbott, df1Roche, plt, sns):
    _fig, _ax = plt.subplots(1, 2, figsize=(15, 7), sharey=True)
    _min_ = 3
    _max_ = 10
    sns.stripplot(x='agrdt', y='vl', data=df1Abbott, ax=_ax[0])
    _ax[0].set_title('Abbott')
    _ax[0].set_ylim(_min_, _max_)
    sns.stripplot(x='agrdt', y='vl', data=df1Roche, ax=_ax[1])
    _ax[1].set_title('Roche')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Omicron samples""")
    return


@app.cell
def _(np):
    def jitter(pdSeries):
        return np.random.normal(pdSeries, 0.015)
    return (jitter,)


@app.cell
def _(df2):
    # Abbott Omicron data
    df2Abbott = df2[df2.test==0].copy()
    # Roche Omicron data
    df2Roche = df2[df2.test==1].copy()
    return df2Abbott, df2Roche


@app.cell
def _(df2Abbott, df2Roche, plt, sns):
    _fig, _ax = plt.subplots(1, 2, figsize=(15, 7), sharey=True)
    sns.countplot(data=df2Abbott, x='agrdt', ax=_ax[0])
    _ax[0].set_title('Abbott')
    sns.countplot(data=df2Roche, x='agrdt', ax=_ax[1])
    _ax[1].set_title('Roche')
    return


@app.cell
def _(df2Abbott, df2Roche, plt, sns):
    _fig, _ax = plt.subplots(1, 2, figsize=(15, 7), sharey=True)
    _min_ = 3
    _max_ = 10
    sns.stripplot(x='agrdt', y='vl', data=df2Abbott, ax=_ax[0])
    _ax[0].set_title('Abbott')
    _ax[0].set_ylim(_min_, _max_)
    sns.stripplot(x='agrdt', y='vl', data=df2Roche, ax=_ax[1])
    _ax[1].set_title('Roche')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""#### Check for batch effect""")
    return


@app.cell
def _(
    boldStr,
    df2,
    df2Abbott,
    df2Roche,
    jitter,
    labelRes,
    labelVl,
    pal1,
    plt,
    sns,
):
    _, _ax = plt.subplots(2, 1, figsize=(8, 10), sharex=True)
    ax1, ax2 = _ax.flatten()
    sns.scatterplot(x=df2Abbott.vl, y=jitter(df2Abbott.agrdt), hue=df2Abbott.batch, palette={1: pal1[0], 2: pal1[1]}, s=40, alpha=0.5, ax=ax1)
    ax1.set_ylabel(labelRes)
    ax1.set_title(f"{boldStr('Abbott')}\n{(df2Abbott.agrdt & (df2Abbott.batch == 1)).sum()} positive in batch 1, {((df2Abbott.agrdt == 0) & (df2.batch == 1)).sum()} negative in batch 1\n{(df2Abbott.agrdt & (df2Abbott.batch == 2)).sum()} positive in batch 2, {((df2Abbott.agrdt == 0) & (df2Abbott.batch == 2)).sum()} negative in batch 2\n")
    sns.scatterplot(x=df2Roche.vl, y=jitter(df2Roche.agrdt), hue=df2Roche.batch, palette={1: pal1[0], 2: pal1[1]}, s=40, alpha=0.5, ax=ax2)
    ax2.set_ylabel(labelRes)
    ax2.set_title(f"{boldStr('Roche')}\n{(df2Roche.agrdt & (df2Roche.batch == 1)).sum()} positive in batch 1, {((df2Roche.agrdt == 0) & (df2Roche.batch == 1)).sum()} negative in batch 1\n{(df2Roche.agrdt & (df2Roche.batch == 2)).sum()} positive in batch 2, {((df2Roche.agrdt == 0) & (df2Roche.batch == 2)).sum()} negative in batch 2\n")
    ax2.set_xlabel(labelVl)
    plt.tight_layout()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Regression""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## All data (Pre-VOC and Omicron)""")
    return


@app.cell
def _(bmb, df):
    dfBmb = df.copy()
    dfBmb = dfBmb.replace({'variant': {'wt': 0, 'omicron': 1}})
    _formula = 'agrdt ~ variant + test + test:variant + zVl + zVl:variant + zVl:test'
    _priors = {'test:variant': bmb.Prior('Normal', mu=0, sigma=1), 'zVl': bmb.Prior('Lognormal', mu=1, sigma=1), 
               'zVl:test': bmb.Prior('Normal', mu=0, sigma=1), 'zVl:variant': bmb.Prior('Normal', mu=0, sigma=1),}
    bmbModel = bmb.Model(_formula, dfBmb, categorical=['variant', 'test'], family='bernoulli', priors=_priors, dropna=True)
    bmbModel.build()
    return bmbModel, dfBmb


@app.cell
def _(bmbModel):
    bmbModel
    return


@app.cell
def _(bmbModel):
    priorPred = bmbModel.prior_predictive()
    return (priorPred,)


@app.cell
def _(az, priorPred):
    az.plot_dist(priorPred.prior_predictive["agrdt"])
    return


@app.cell
def _(az, bmbModel, iDataDir, postProcessAbbottRoche):
    iDataPath = iDataDir / "abbottRocheWtOmicron.nc"
    if iDataPath.exists():
        iData = az.from_netcdf(iDataPath)
    else:
        iData = bmbModel.fit(target_accept=0.99, draws=5000, tune=4000, idata_kwargs={'log_likelihood': True})
        postProcessAbbottRoche(iData)
        iData.to_netcdf(iDataPath)
    return (iData,)


@app.cell
def _(az, iData, plt):
    az.plot_trace(iData, var_names=["~abbott-roche_variant"], compact=True);
    plt.tight_layout()
    return


@app.cell
def _(az, iData):
    az.summary(iData, "~p")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""##### Plots""")
    return


@app.cell
def _(az, iData):
    az.plot_forest(iData, combined=True)
    return


@app.cell
def _(colWidth, iData, label, plotDir, plt, ridgeForestPlot, saveFigure):
    var_names = ["variant", "roche_abbott_variant", "zVl_variant_test"]
    yTickLabels = ("Omicron - pre-VOC", "Roche - Abbott | pre-VOC", "Roche - Abbott | Omicron",
                   f"{label.vl} | pre-VOC, Abbott", f"{label.vl} | pre-VOC, Roche", f"{label.vl} | Omicron, Abbott", 
                   f"{label.vl} | Omicron, Roche")
    _fig, _ax = plt.subplots(1, 1, figsize=(colWidth * 2, colWidth * 1))
    _ax = ridgeForestPlot(iData, varNames=tuple(var_names), yTickLabels=yTickLabels, addGrid=True, ridge=False, ax=_ax)
    plt.show()
    saveFigure(_fig, plotDir / "abbottVsRocheForestPlot.png")
    return (var_names,)


@app.cell
def _(az, iData):
    summaryDf = az.summary(iData, var_names="~p")
    return (summaryDf,)


@app.cell
def _(summaryDf):
    summaryDf
    return


@app.cell
def _(
    summaryDf,
    tableDirIData,
    var_names,
    writeIDataSummaryTable,
    writeIDataSummaryTableLatex,
):
    writeIDataSummaryTable(summaryDf=summaryDf, varnames=var_names, outfile=tableDirIData / "iDataSummaryAbbottRoche.tsv")
    writeIDataSummaryTableLatex(summaryDf=summaryDf, varnames=var_names, outfile=tableDirIData / "iDataSummaryAbbottRocheLatex.txt")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""#### Posterior predictions""")
    return


@app.cell
def _(df, np):
    vlMean = df.vl.mean()
    vlSd = df.vl.std()
    nCats = 2

    vl = np.arange(3, 11, step=0.1)
    zVl = (vl - vlMean) / vlSd
    nVl = len(zVl)
    n = nCats * nVl

    dtype = np.int32
    return dtype, n, nCats, nVl, vl, zVl


@app.cell
def _(dtype, n, nCats, nVl, np, pd, zVl):
    # Displaying difference for variant:
    newData = pd.DataFrame({
        "variant": np.repeat([0, 1], nVl),
        "test": np.zeros(n, dtype=dtype),
        "zVl": np.tile(zVl, nCats),
    })

    # Displaying difference for rapid test:
    newData2 = pd.DataFrame({
        "variant": np.ones(n, dtype=dtype),
        "test": np.repeat([0, 1], nVl),
        "zVl": np.tile(zVl, nCats),
    })
    return (newData,)


@app.cell
def _(bmbModel, iData, newData):
    bmbModel.predict(iData, data=newData)
    return


@app.cell
def _(defaultdict, iData, newData, vl):
    predDict = defaultdict(dict)
    predDict["variant"]["vl"] = vl
    predDict["variant"]["data"] = newData
    predDict["variant"]["ppSamples"] = iData.posterior["p"].stack(samples=("chain", "draw")).values
    return (predDict,)


@app.cell
def _(
    colWidth,
    dfBmb,
    fontSizePlot,
    legend,
    pal,
    plotDir,
    plt,
    predDict,
    saveFigure,
    setFontSize,
    spaghettiPlotCategorical,
):
    _fig, _ax = plt.subplots(1, 1, figsize=(colWidth * 2, colWidth * 1.1), constrained_layout=True)
    palette = {i: pal.variant[variant] for i, variant in enumerate(('wildtype', 'omicron'))}
    spaghettiPlotCategorical(_ax, 'variant', (0, 1), predDict, dfBmb, palette=palette, legend=legend.variantLinesWtOmicron)
    setFontSize(_ax, size=fontSizePlot)
    # plt.tight_layout()
    saveFigure(_fig, plotDir / "abbottVsRocheSpaghettiPlot.png")
    _ax
    return


@app.cell
def _(mo):
    mo.md(r"""### Hierarchical model""")
    return


@app.cell
def _(bmb, df):
    dfBmbHr = df.copy()
    dfBmbHr = dfBmbHr.replace({'variant': {'wt': 0, 'omicron': 1}})
    _formula = 'agrdt ~ variant + test + (0 + test|variant) + zVl + (0 + zVl|variant + test)'
    _priors = {'test|variant': bmb.Prior('Normal', mu=0, sigma=bmb.Prior('HalfNormal', sigma=1)), 'zVl': bmb.Prior('Lognormal', mu=1, sigma=1), 'zVl|test': bmb.Prior('Normal', mu=0, sigma=bmb.Prior('HalfNormal', sigma=1)), 
              'zVl|variant': bmb.Prior('Normal', mu=0, sigma=bmb.Prior('HalfNormal', sigma=1))}
    bmbModelHr = bmb.Model(_formula, dfBmbHr, categorical=['variant', 'test'], family='bernoulli', priors=_priors, 
                           dropna=True)
    bmbModelHr.build()
    return (bmbModelHr,)


@app.cell
def _(az, bmbModelHr, iDataDir, postProcessAbbottRocheHierarchical):
    iDataHrPath = iDataDir / "abbottRocheWtOmicronHierarchical.nc"
    if iDataHrPath.exists():
        iDataHr = az.from_netcdf(iDataHrPath)
    else:
        iDataHr = bmbModelHr.fit(target_accept=0.99, draws=5000, tune=4000, idata_kwargs={'log_likelihood': True})
        postProcessAbbottRocheHierarchical(iDataHr)
        iDataHr.to_netcdf(iDataHrPath)
    return (iDataHr,)


@app.cell
def _(mo):
    mo.md(r"""### Compare fixed effects and mixed effects model""")
    return


@app.cell
def _(az, iData, iDataHr):
    compareDict = {"fixed_effects": iData,
                   "mixed_effects": iDataHr}
    az.compare(compareDict, ic="loo")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Only Omicron samples""")
    return


@app.cell
def _(bmb, df2):
    _formula = 'agrdt ~ test + zVl + (0 + zVl|test)'
    _priors = {'zVl': bmb.Prior('Lognormal', mu=1, sigma=1), 'zVl|test': bmb.Prior('Normal', mu=0, sigma=bmb.Prior('HalfNormal', sigma=1))}
    bmbModel2 = bmb.Model(_formula, df2, categorical=['test'], family='bernoulli', priors=_priors, dropna=True)
    bmbModel2.build()
    return (bmbModel2,)


@app.cell
def _(bmbModel2):
    bmbModel2
    return


@app.cell
def _(SEED, bmbModel2):
    iData2_1 = bmbModel2.fit(target_accept=0.99, draws=10000, tune=4000, random_seed=SEED)
    return (iData2_1,)


@app.cell
def _(iData2_1, iDataDir):
    iData2_1.to_netcdf(iDataDir + '/figA1.nc')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""##### Process zVl|test - compute differences between tests""")
    return


@app.cell
def _(iData2_1):
    postZVl_1 = iData2_1.posterior['zVl']
    return (postZVl_1,)


@app.cell
def _(iData2_1, postZVl_1):
    abbott_zVl_1 = postZVl_1 + iData2_1.posterior.sel(test__factor_dim='0')['zVl|test']
    roche_zVl_1 = postZVl_1 + iData2_1.posterior.sel(test__factor_dim='1')['zVl|test']
    return abbott_zVl_1, roche_zVl_1


@app.cell
def _(abbott_zVl_1):
    abbott_zVl_1['test__factor_dim'] = 'abbott'
    abbott_zVl_1['test__factor_dim'] = 'roche'
    return


@app.cell
def _(abbott_zVl_1, roche_zVl_1, xarray):
    zVlAbbottRoche_1 = xarray.concat([abbott_zVl_1, roche_zVl_1], dim='new_dim')
    zVlAbbottRoche_1.name = 'zVl_roche_abbott'
    return (zVlAbbottRoche_1,)


@app.cell
def _(iData2_1, zVlAbbottRoche_1):
    iData2_1.posterior['zVl_roche_abbott'] = zVlAbbottRoche_1
    return


@app.cell
def _(az, iData2_1):
    az.summary(iData2_1)
    return


@app.cell
def _(
    colWidth,
    colWidthAppendix,
    fontSizePlotSuppl,
    iData2_1,
    label,
    plotDir,
    plt,
    ridgeForestPlot,
):
    _ax = ridgeForestPlot(iData2_1, varNames=('test', 'zVl_roche_abbott'), yTickLabels=('Roche - Abbott', f'{label.vl} | Abbott', f'{label.vl} | Roche'), figsize=(colWidthAppendix, colWidth), fontsize=fontSizePlotSuppl)
    plt.tight_layout()

    plt.savefig(f'{plotDir}/FigureA1.png', dpi=600, bbox_inches='tight')
    plt.savefig(f'{plotDir}/FigureA1.pdf', dpi=600, bbox_inches='tight')
    return


if __name__ == "__main__":
    app.run()
