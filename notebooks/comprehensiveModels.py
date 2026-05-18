import marimo

__generated_with = "0.16.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import matplotlib.pyplot as plt
    from collections import defaultdict
    from pathlib import Path

    # Bayesian MCMC analysis
    import arviz as az
    import bambi as bmb

    from agrdt.data import createDataFramesFigures, returnRtData, removeReleaseTesting, addSymptomsURT
    from agrdt.tables import (returnTableDirIData, writeIDataSummaryTable,
                              writeIDataSummaryTableLatex)
    from agrdt.plotting import (annotateWithLetters, annotateWithLetter, ridgeForestPlot,
                                setFontSize, spaghettiPlotCategorical, returnPlotDirRegression,
                                plotDataPointsRegression, setCustomTheme, saveFigure)
    from agrdt.plotParams import (getPalettes, getLabels, getLegends, getOrders, getAbbrvsDict, COL_WIDTH, DINA4_HEIGHT, CM,
                                  ANNOTATION_LETTER_SIZE, ANNOTATION_COORDS)
    from agrdt.regression import (SEED, logisticRegressionVars, logisticRegressionDf, predictionsNewData, returnIDataDirRegression,
                                  postProcessModel1, postProcessModel2, postProcessModel3, postProcessModel4,
                                  postProcessModel5, postProcessModel6)
    from agrdt.dataParams import ROOT_DIR

    setCustomTheme()
    pd.set_option('display.max_columns', 30)

    # I had problems with C code compilation when trying to run bambi models (after updating to a newer macOS version (Tahoe))
    # See https://discourse.pymc.io/t/environment-not-working-anymore-on-macos/14210/16?page=2
    # None of the the above worked.

    # This hopefully works:
    # from https://github.com/pymc-devs/pytensor/issues/1342
    # Problematic path in conda environment
    # (pymc_env) $ which clang++
    # /Users/USER/miniconda3/envs/pymc_env/bin/clang++  ❌

    # # Remove conda's compiler toolchain
    # conda remove --force clang clangxx -n pymc_env

    # # Verification
    # (pymc_env) $ which clang++
    # /usr/bin/clang++  ✅
    return (
        ANNOTATION_COORDS,
        ANNOTATION_LETTER_SIZE,
        CM,
        COL_WIDTH,
        DINA4_HEIGHT,
        Path,
        ROOT_DIR,
        SEED,
        addSymptomsURT,
        annotateWithLetter,
        annotateWithLetters,
        az,
        bmb,
        createDataFramesFigures,
        defaultdict,
        getAbbrvsDict,
        getLabels,
        getLegends,
        getOrders,
        getPalettes,
        logisticRegressionDf,
        logisticRegressionVars,
        mo,
        pd,
        plotDataPointsRegression,
        plt,
        postProcessModel1,
        postProcessModel2,
        postProcessModel3,
        postProcessModel4,
        postProcessModel5,
        postProcessModel6,
        predictionsNewData,
        removeReleaseTesting,
        returnIDataDirRegression,
        returnPlotDirRegression,
        returnRtData,
        returnTableDirIData,
        ridgeForestPlot,
        saveFigure,
        setFontSize,
        spaghettiPlotCategorical,
        writeIDataSummaryTable,
        writeIDataSummaryTableLatex,
    )


@app.cell(hide_code=True)
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
    returnIDataDirRegression,
    returnPlotDirRegression,
    returnTableDirIData,
):
    pal = getPalettes()
    legend = getLegends()
    order = getOrders()
    label = getLabels()
    plotDir = returnPlotDirRegression()
    iDataDir = returnIDataDirRegression()
    tableDirIData = returnTableDirIData()
    return iDataDir, label, legend, pal, plotDir, tableDirIData


@app.cell
def _(getAbbrvsDict):
    abbrvDictPaper = getAbbrvsDict()

    genderAbbrvs = abbrvDictPaper["gender"]
    symptomAbbrvs = abbrvDictPaper["symptoms"]
    variantAbbrvs = abbrvDictPaper["variant"]
    immunAbbrvs = abbrvDictPaper["immun2YN"]
    daysAbbrvs = abbrvDictPaper["binDaysPostOnset"]
    return


@app.cell
def _(CM, COL_WIDTH):
    cm = CM
    colWidth = COL_WIDTH
    fontSizePlot = 12
    fontSizePlotSuppl = 10
    return colWidth, fontSizePlot


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Data""")
    return


@app.cell
def _(Path, ROOT_DIR, pd, returnRtData):
    # Load data
    target="T1"
    dataPath = Path(ROOT_DIR, "data", f"agrdtDataThesis{target}.tsv")
    dataPath2 = Path(ROOT_DIR, "data", f"agrdtDataThesisT2.tsv")
    df = pd.read_csv(dataPath, sep="\t", low_memory=False, parse_dates=['pcrDate'])
    df2 = pd.read_csv(dataPath2, sep="\t", low_memory=False, parse_dates=['pcrDate'])
    # Turn into datetime.date objects.
    df["pcrDate"] = df["pcrDate"].apply(lambda pcrDate: pcrDate.date())
    df["binDaysPostOnset4"] = pd.cut(df.daysPostOnset, bins=(0, 1, 7), include_lowest=True)
    dfRt = returnRtData()
    df = df.merge(dfRt[["dateDate", "rtRolling4", "rt"]], 
                  left_on="pcrDate", right_on="dateDate", how="left")
    df["epidemicGrowth"] = df.rt >= 1
    return df, target


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Data frames for figures""")
    return


@app.cell
def _(addSymptomsURT, createDataFramesFigures, df):
    dfTuple = createDataFramesFigures(df)
    (dfAgrdt,
     dfAgrdtAll,
     dfAgrdtIndInf,
     dfAllFirstPosPcrsSympNoRelease,
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
     dfFigureA2,) = dfTuple

    dfPos = addSymptomsURT(dfPos)
    return (dfPos,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Regression""")
    return


@app.cell
def _(mo):
    mo.md(r"""## Model parameter specifications""")
    return


@app.cell
def _():
    varNames1 = ('zAge', 'gender', 'testDevice', 'variant', 'binDaysPostOnset', 'zVl')
    varNames1Days4 = ('zAge', 'gender', 'testDevice', 'variant', 'binDaysPostOnset4', 'zVl')
    varNames1InteractionVariantDays = ('zAge', 'gender', 'testDevice', 'variant', 'daysVariant', 'zVl')
    varNames1InteractionVariantDays4 = ('zAge', 'gender', 'testDevice', 'variant', 'days4Variant', 'zVl')
    varNames1InteractionVariantDays4SymptomsURT = ('zAge', 'gender', 'testDevice', 'variant', 'illURT', 'days4Variant',
                                                   'zVl')
    varNames1SymptomsURTDays4Interaction = ('zAge', 'gender', 'testDevice', 'variant', 'illURT', 'days4SymptomsURT', 'zVl')

    varNames2 = ('zAge', 'gender', 'testDevice', 'variant', 'immun2YN', 'recovered', 'binDaysPostOnset', 'zVl')
    varNames2InteractionVariantDays = ('zAge', 'gender', 'testDevice', 'variant', 'immun2YN', 'recovered', 'daysVariant',
                                       'zVl')
    varNames2InteractionVariantDaysReverseDays = ('zAge', 'gender', 'testDevice', 'variant', 'immun2YN', 'recovered',
                                                  'daysReverseVariant', 'zVl')
    varNames2InteractionVariantDays4 = ('zAge', 'gender', 'testDevice', 'variant', 'immun2YN', 'recovered', 'days4Variant', 
                                        'zVl')
    varNames2InteractionVariantDays4ReverseDays = ('zAge', 'gender', 'testDevice', 'variant', 'immun2YN', 'recovered', 
                                                   'days4ReverseVariant', 'zVl')

    varNames2SymptomsURT = ('zAge', 'gender', 'testDevice', 'variant', 'immun2YN', 'recovered','illURT', 
                            'binDaysPostOnset', 'zVl')
    varNames2InteractionVariantDaysSymptomsURT = ('zAge', 'gender', 'testDevice', 'variant', 'immun2YN', 'recovered',
                                                'illURT', 'daysVariant', 'zVl')
    varNames2InteractionVariantDaysReverseDaysSymptomsURT = ('zAge', 'gender', 'testDevice', 'variant', 'immun2YN', 
                                                           'recovered', 'illURT', 'daysReverseVariant', 'zVl')
    varNames2InteractionVariantDays4SymptomsURT = ('zAge', 'gender', 'testDevice', 'variant', 'immun2YN', 'recovered', 
                                                   'illURT', 'days4Variant', 'zVl')


    varNames3 = ('zAge', 'gender', 'testDevice', 'immun2YN', 'binDaysPostOnset', 'zVl')
    varNames3Immun2YNDaysInteraction = ('zAge', 'gender', 'testDevice', 'immun2YN', 'recovered', 'daysImmun2YN', 'zVl')
    varNames3Immun2YNDays4InteractionSymptomsURT = ('zAge', 'gender', 'testDevice', 'immun2YN', 'recovered', 'illURT',
                                                    'days4Immun2YN', 'zVl')
    varNames3Immun2YNDays4Interaction = ('zAge', 'gender', 'testDevice', 'immun2YN', 'recovered',
                                         'days4Immun2YN', 'zVl')

    varNames4 = ('zAge', 'gender', 'testDevice', 'variantSymptoms', 'symptoms', 'zVl')

    varNames5 = ('zAge',  'gender', 'testDevice', 'variant', 'symptoms', 'immun2YNSymptoms', 'recoveredSymptoms', 'zVl')
    varNames5VariantSymptomsInteraction = ('zAge',  'gender', 'testDevice', 'variantSymptoms', 'symptoms', 
                                            'immun2YNSymptoms', 'recoveredSymptoms', 'zVl')

    varNames6 = ('zAge',  'gender', 'testDevice', 'symptoms', 'immun2YNSymptoms', 'recoveredSymptoms', 'zVl')
    varNames6SymptomsURT = ('zAge',  'gender', 'testDevice', 'illURT', 'immun2YNSymptomsURT', 'recoveredSymptomsURT', 'zVl')

    modelNoMapping = {1: 12, 2: 13, 3: 14, 4: 15, 5: 16, 6: 17}
    return (
        modelNoMapping,
        varNames1,
        varNames1Days4,
        varNames1InteractionVariantDays,
        varNames1InteractionVariantDays4,
        varNames1InteractionVariantDays4SymptomsURT,
        varNames1SymptomsURTDays4Interaction,
        varNames2,
        varNames2InteractionVariantDays,
        varNames2InteractionVariantDays4,
        varNames2InteractionVariantDays4ReverseDays,
        varNames2InteractionVariantDays4SymptomsURT,
        varNames2InteractionVariantDaysReverseDays,
        varNames2InteractionVariantDaysReverseDaysSymptomsURT,
        varNames2InteractionVariantDaysSymptomsURT,
        varNames2SymptomsURT,
        varNames3,
        varNames3Immun2YNDays4Interaction,
        varNames3Immun2YNDays4InteractionSymptomsURT,
        varNames3Immun2YNDaysInteraction,
        varNames4,
        varNames5,
        varNames5VariantSymptomsInteraction,
        varNames6,
        varNames6SymptomsURT,
    )


@app.cell
def _(mo):
    mo.md(r"""## Only symptomatic infections""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### Include variable for SARS-CoV-2 variant (Model 1)""")
    return


@app.cell
def _(bmb, dfPos, logisticRegressionDf):
    familyModel1 = "cumulative" # bernoulli
    binDaysPostOnsetVar1 = "binDaysPostOnset4"
    _outcome = "testline" if familyModel1 == "cumulative" else "agrdt"

    interactionVariantDays1 = True
    includeSymptomsURT1 = False
    interactionSymptomsDays1 = False
    assert (interactionSymptomsDays1 + interactionVariantDays1) < 2

    _indVars = ("vl", "age", "testline") if familyModel1=="cumulative" else ("vl", "age")
    _indCatVars = (("testDevice", "variant", "gender", binDaysPostOnsetVar1, "illURT") if 
                   includeSymptomsURT1 else 
                   ("testDevice", "variant", "gender", binDaysPostOnsetVar1))

    formulaSymptomsSuffix1 = "+ illURT" if includeSymptomsURT1 else ""

    _daysVariantVar = f"{binDaysPostOnsetVar1}:variant" if interactionVariantDays1 else binDaysPostOnsetVar1
    _daysSymptomsVar = f"{binDaysPostOnsetVar1}:illURT" if interactionSymptomsDays1 else binDaysPostOnsetVar1
    _daysVar = _daysVariantVar if interactionVariantDays1 else _daysSymptomsVar

    _formulaAllInd = f'{_outcome} ~ testDevice + variant + gender + {_daysVar} + zVl + zAge {formulaSymptomsSuffix1}'

    dfLogisticBmb1 = logisticRegressionDf(dfPos, indCatVars=_indCatVars, indVars=_indVars)

    _priors = {'zVl': bmb.Prior('Lognormal', mu=1, sigma=1),
               'variant': bmb.Prior('Normal', mu=0, sigma=1),
               f'{binDaysPostOnsetVar1}:variant': bmb.Prior('Normal', mu=0, sigma=1), 
               f'{binDaysPostOnsetVar1}:illURT': bmb.Prior('Normal', mu=0, sigma=1), 
               f'{binDaysPostOnsetVar1}': bmb.Prior('Normal', mu=0, sigma=1),
               'testDevice': bmb.Prior('Normal', mu=-0.5, sigma=0.5), 
               'zAge': bmb.Prior('Normal', mu=0, sigma=1), 
               'gender': bmb.Prior('Normal', mu=0, sigma=1),
               'illURT': bmb.Prior('Normal', mu=0, sigma=1),
               'recovered': bmb.Prior('Normal', mu=0, sigma=1)}
    bmbModel1 = bmb.Model(_formulaAllInd, dfLogisticBmb1, priors=_priors, family=familyModel1, categorical=list(_indCatVars), 
                          noncentered=True)
    bmbModel1.build()
    return (
        binDaysPostOnsetVar1,
        bmbModel1,
        dfLogisticBmb1,
        familyModel1,
        includeSymptomsURT1,
        interactionSymptomsDays1,
        interactionVariantDays1,
    )


@app.cell
def _(bmbModel1):
    bmbModel1
    return


@app.cell
def _(bmbModel1):
    bmbModel1.graph()
    return


@app.cell
def _(
    Path,
    SEED,
    az,
    binDaysPostOnsetVar1,
    bmbModel1,
    familyModel1,
    iDataDir,
    includeSymptomsURT1,
    interactionSymptomsDays1,
    interactionVariantDays1,
    postProcessModel1,
    target,
):
    binDaysPostOnsetVarSuffix1 = binDaysPostOnsetVar1[0].capitalize() + binDaysPostOnsetVar1[1:]
    interactionVariantDaysSuffix1 = "InteractionVariantDays" if interactionVariantDays1 else ""
    interactionSymptomsDaysSuffix1 = "InteractionSymptomsURTDays" if interactionSymptomsDays1 else ""
    symptomsVarSuffix1 = "SymptomsURT" if includeSymptomsURT1 else ""
    iDataModel1File = Path(iDataDir, f"model1{target}{binDaysPostOnsetVarSuffix1}{familyModel1.capitalize()}"
                           f"{interactionVariantDaysSuffix1}{interactionSymptomsDaysSuffix1}{symptomsVarSuffix1}.nc")
    if iDataModel1File.exists():
        iDataBmb1 = az.from_netcdf(iDataModel1File)
    else:
        iDataBmb1 = bmbModel1.fit(target_accept=0.99, draws=10000, tune=4000, random_seed=SEED, 
                                  idata_kwargs={'log_likelihood': True})
        postProcessModel1(iDataBmb1, interactionVariantDays=interactionVariantDays1, 
                          interactionSymptomsDays=interactionSymptomsDays1,
                          binDaysPostOnsetVar=binDaysPostOnsetVar1)
        iDataBmb1.to_netcdf(iDataModel1File)
    return (
        binDaysPostOnsetVarSuffix1,
        iDataBmb1,
        interactionSymptomsDaysSuffix1,
        interactionVariantDaysSuffix1,
        symptomsVarSuffix1,
    )


@app.cell
def _(az, iDataBmb1):
    az.summary(iDataBmb1)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""#### Posterior analysis""")
    return


@app.cell
def _(az, iDataBmb1, plt):
    az.plot_trace(iDataBmb1)
    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(
    binDaysPostOnsetVar,
    binDaysPostOnsetVar1,
    includeSymptomsURT1,
    interactionSymptomsDays1,
    interactionVariantDays1,
    varNames1,
    varNames1Days4,
    varNames1InteractionSymptomsURTDays,
    varNames1InteractionVariantDays,
    varNames1InteractionVariantDays4,
    varNames1InteractionVariantDays4SymptomsURT,
    varNames1SymptomsURTDays4Interaction,
):
    if interactionVariantDays1:
        if binDaysPostOnsetVar1 == "binDaysPostOnset":
            varNames1Tmp = varNames1InteractionVariantDays
        elif binDaysPostOnsetVar1 == "binDaysPostOnset4":
            varNames1Tmp = varNames1InteractionVariantDays4
            if includeSymptomsURT1:
                varNames1Tmp = varNames1InteractionVariantDays4SymptomsURT
    elif interactionSymptomsDays1:
        if binDaysPostOnsetVar1 == "binDaysPostOnset":
            varNames1Tmp = varNames1InteractionSymptomsURTDays
        elif binDaysPostOnsetVar1 == "binDaysPostOnset4":
            varNames1Tmp = varNames1SymptomsURTDays4Interaction
    else:
        if binDaysPostOnsetVar == "binDaysPostOnset4":
            varNames1Tmp = varNames1Days4
        else:
            varNames1Tmp = varNames1
    return (varNames1Tmp,)


@app.cell
def _(
    binDaysPostOnsetVarSuffix1,
    colWidth,
    familyModel1,
    iDataBmb1,
    interactionSymptomsDaysSuffix1,
    interactionVariantDaysSuffix1,
    modelNoMapping,
    plotDir,
    plt,
    ridgeForestPlot,
    saveFigure,
    symptomsVarSuffix1,
    target,
    varNames1Tmp,
):
    _fig, _ax = plt.subplots(1, 1, figsize=(colWidth * 2, colWidth * 1.3))
    ridgeForestPlot(iDataBmb1, varNames=varNames1Tmp, addGrid=True, ax=_ax)
    saveFigure(_fig, f'{plotDir}/forestPlotModel{modelNoMapping[1]}{target}{familyModel1.capitalize()}'
               f'{interactionVariantDaysSuffix1}{binDaysPostOnsetVarSuffix1}{interactionSymptomsDaysSuffix1}'
               f'{symptomsVarSuffix1}.pdf')
    saveFigure(_fig, f'{plotDir}/forestPlotModel{modelNoMapping[1]}{target}{familyModel1.capitalize()}'
               f'{interactionVariantDaysSuffix1}{binDaysPostOnsetVarSuffix1}{interactionSymptomsDaysSuffix1}'
               f'{symptomsVarSuffix1}.png')
    _ax
    return


@app.cell
def _(
    az,
    iDataBmb1,
    modelNoMapping,
    tableDirIData,
    varNames1Tmp,
    writeIDataSummaryTable,
    writeIDataSummaryTableLatex,
):
    summaryDf1 = az.summary(iDataBmb1, var_names=['~agrdt_mean', '~p'])
    writeIDataSummaryTable(summaryDf1, varnames=varNames1Tmp, 
                           outfile=tableDirIData / f'iDataSummaryModel{modelNoMapping[1]}.tsv')
    writeIDataSummaryTableLatex(summaryDf1, varnames=varNames1Tmp, 
                                outfile=tableDirIData / f'iDataSummaryModel{modelNoMapping[1]}Latex.txt', 
                                model_no=1)
    return (summaryDf1,)


@app.cell
def _(summaryDf1):
    summaryDf1
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""##### Spaghetti plots""")
    return


@app.cell
def _(defaultdict):
    predDictBmb1 = defaultdict(dict)
    return (predDictBmb1,)


@app.cell
def _(bmbModel1, dfLogisticBmb1, iDataBmb1, predDictBmb1, predictionsNewData):
    predictionsNewData(dfLogisticBmb1, bmbModel1, iDataBmb1, predDictBmb1, "variant", (0, 1, 2, 3))
    predictionsNewData(dfLogisticBmb1, bmbModel1, iDataBmb1, predDictBmb1, "binDaysPostOnset4", (0, 1))
    return


@app.cell
def _(
    ANNOTATION_LETTER_SIZE,
    annotateWithLetters,
    colWidth,
    dfLogisticBmb1,
    modelNoMapping,
    plotDir,
    plt,
    predDictBmb1,
    saveFigure,
    spaghettiPlotCategorical,
):
    figModel1Spagh, _ax = plt.subplots(2, 1, figsize=(colWidth * 2, colWidth * 2.5), sharex=True)
    _ax1, _ax2 = _ax.flatten()
    spaghettiPlotCategorical(_ax1, feature='variant', featureLevels=(0, 1, 2, 3), predDict=predDictBmb1, 
                             df=dfLogisticBmb1, showXLabel=False, 
                             legendLoc=(0.565, 0.15))
    spaghettiPlotCategorical(_ax2, feature='binDaysPostOnset4', featureLevels=(0, 1), predDict=predDictBmb1, 
                             df=dfLogisticBmb1)
    plt.tight_layout()

    annotateWithLetters((_ax1, _ax2), size=ANNOTATION_LETTER_SIZE, coords=(-0.15, 1))
    saveFigure(figModel1Spagh, plotDir / f'model{modelNoMapping[1]}Spagh.png')
    saveFigure(figModel1Spagh, plotDir / f'model{modelNoMapping[1]}Spagh.pdf')
    return


@app.cell
def _(mo):
    mo.md(r"""### Include both variant and immunization status (Model 2)""")
    return


@app.cell
def _(logisticRegressionVars):
    _indVars, _zIndVars, _ = logisticRegressionVars(immunVar='immun2YN')
    return


@app.cell
def _(bmb, dfPos, logisticRegressionDf):
    includeSymptomsURT2 = False

    binDaysPostOnsetVar2 = "binDaysPostOnset4"
    familyModel2 = "cumulative" # bernoulli
    _indVars = ("vl", "age", "testline") if familyModel2=="cumulative" else ("vl", "age")
    interactionVariantDays2 = True
    _indCatVars = (("testDevice", "immun2YN", "variant", "gender", binDaysPostOnsetVar2, "illURT", "recovered") 
                   if includeSymptomsURT2 else 
                   ("testDevice", "immun2YN", "variant", "gender", binDaysPostOnsetVar2, "recovered"))

    formulaSymptomsSuffix2 = "+ illURT" if includeSymptomsURT2 else ""
    _outcome = "testline" if familyModel2 == "cumulative" else "agrdt"
    _formulaAllInd = (f'{_outcome} ~ testDevice + immun2YN + variant + gender + {binDaysPostOnsetVar2}:variant '
                      f'+ zVl + zAge {formulaSymptomsSuffix2} + recovered' if interactionVariantDays2 else 
                      f'{_outcome} ~ testDevice + immun2YN + variant + gender + {binDaysPostOnsetVar2} + zVl + '
                      f'zAge {formulaSymptomsSuffix2} + recovered')

    dfLogisticBmb2 = logisticRegressionDf(dfPos, immunVar='immun2YN', indCatVars=_indCatVars, indVars=_indVars)

    sigma2 = 1
    _priors = {'zVl': bmb.Prior('Lognormal', mu=1, sigma=1), 
               'binDaysPostOnset': bmb.Prior('Normal', mu=0, sigma=sigma2),
               'binDaysPostOnset4': bmb.Prior('Normal', mu=0, sigma=sigma2),
               'binDaysPostOnsetReverse': bmb.Prior('Normal', mu=0, sigma=sigma2),
               'binDaysPostOnset4Reverse': bmb.Prior('Normal', mu=0, sigma=sigma2),
               'binDaysPostOnset:variant': bmb.Prior('Normal', mu=0, sigma=sigma2),
               'binDaysPostOnset4:variant': bmb.Prior('Normal', mu=0, sigma=sigma2),
               'binDaysPostOnsetReverse:variant': bmb.Prior('Normal', mu=0, sigma=sigma2),
               'binDaysPostOnset4Reverse:variant': bmb.Prior('Normal', mu=0, sigma=sigma2),
               'testDevice': bmb.Prior('Normal', mu=-0.5, sigma=0.5), 
               'zAge': bmb.Prior('Normal', mu=0, sigma=1), 
               'gender': bmb.Prior('Normal', mu=0, sigma=1),
               'variant': bmb.Prior('Normal', mu=0, sigma=sigma2),
               'immun2YN': bmb.Prior('Normal', mu=0, sigma=sigma2),
               'illURT': bmb.Prior('Normal', mu=0, sigma=sigma2),
               'recovered': bmb.Prior('Normal', mu=0, sigma=sigma2),}
    bmbModel2 = bmb.Model(_formulaAllInd, dfLogisticBmb2, priors=_priors, family=familyModel2, 
                          categorical=list(_indCatVars), noncentered=True)
    bmbModel2.build()
    return (
        binDaysPostOnsetVar2,
        bmbModel2,
        dfLogisticBmb2,
        familyModel2,
        includeSymptomsURT2,
        interactionVariantDays2,
        sigma2,
    )


@app.cell
def _(bmbModel2):
    bmbModel2
    return


@app.cell
def _(dfLogisticBmb2):
    maskNegPoints = (dfLogisticBmb2.immun2YN==0) & dfLogisticBmb2.variant.isin((2, 3)) & (dfLogisticBmb2.agrdt==0)
    negPointsYCoords = dfLogisticBmb2[maskNegPoints].vl
    return (negPointsYCoords,)


@app.cell
def _(dfLogisticBmb2):
    dfPrevocAlpha = dfLogisticBmb2[dfLogisticBmb2.variant.isin((0, 1))].copy()
    dfDeltaOmicron_1 = dfLogisticBmb2[dfLogisticBmb2.variant.isin((2, 3))]
    return dfDeltaOmicron_1, dfPrevocAlpha


@app.cell
def _(dfDeltaOmicron_1, dfPrevocAlpha):
    print(f'Pre-Voc and Alpha (pos, neg): {dfPrevocAlpha.agrdt.value_counts().loc[1]}, {dfPrevocAlpha.agrdt.value_counts().loc[0]}')
    print(f'Delta and Omicron, immune naive (pos, neg): {dfDeltaOmicron_1[dfDeltaOmicron_1.immun2YN == 0].agrdt.value_counts().loc[1]}, {dfDeltaOmicron_1[dfDeltaOmicron_1.immun2YN == 0].agrdt.value_counts().loc[0]}')
    print(f'Delta and Omicron, immunized (pos, neg): {dfDeltaOmicron_1[dfDeltaOmicron_1.immun2YN == 1].agrdt.value_counts().loc[1]}, {dfDeltaOmicron_1[dfDeltaOmicron_1.immun2YN == 1].agrdt.value_counts().loc[0]}')
    return


@app.cell
def _(
    colWidth,
    dfLogisticBmb2,
    fontSizePlot,
    label,
    pal,
    plotDir,
    plotFigA4_left,
    plt,
    saveFigure,
):
    fig_A4_left, _ax = plt.subplots(1, 1, figsize=(colWidth, colWidth))
    plotFigA4_left(df=dfLogisticBmb2, palette=pal, label=label, ax=_ax, fontSizePlot=fontSizePlot)
    saveFigure(fig_A4_left, plotDir / 'FigureA4_left.png')
    saveFigure(fig_A4_left, plotDir / 'FigureA4_left.pdf')
    return


@app.cell
def _(
    colWidth,
    dfLogisticBmb2,
    fontSizePlot,
    legend,
    negPointsYCoords,
    pal,
    plotDir,
    plotFigA4_right,
    plt,
    saveFigure,
):
    fig_A4_right, _ax = plt.subplots(1, 1, figsize=(colWidth, colWidth))
    plotFigA4_right(dfLogisticBmb2, negPointsYCoords, palette=pal, legend=legend, ax=_ax, fontSizePlot=fontSizePlot)
    _ax.set_yticklabels('')
    saveFigure(fig_A4_right, plotDir / 'FigureA4_right.png')
    saveFigure(fig_A4_right, plotDir / 'FigureA4_right.pdf')
    return


@app.cell
def _(
    colWidth,
    dfLogisticBmb2,
    fontSizePlot,
    label,
    legend,
    negPointsYCoords,
    pal,
    plotDataPointsRegression,
    plotDir,
    plt,
    saveFigure,
):
    figA4, _ax = plt.subplots(1, 2, figsize=(2 * colWidth, colWidth), sharey=True)
    _ax1, _ax2 = _ax.flatten()
    plotDataPointsRegression(dfLogisticBmb2, negPointsYCoords, _ax1, _ax2, pal=pal, label=label, legend=legend, 
                             fontSizePlot=fontSizePlot, fig=figA4)
    # Remove x-tick marks
    _ax1.tick_params(axis="x", which="both", bottom=False, top=False)
    _ax2.tick_params(axis="x", which="both", bottom=False, top=False)
    plt.tight_layout()
    saveFigure(figA4, plotDir / 'FigureA4.png')
    saveFigure(figA4, plotDir / 'FigureA4.pdf')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""#### Using all data points""")
    return


@app.cell
def _(
    Path,
    SEED,
    az,
    binDaysPostOnsetVar2,
    bmbModel2,
    familyModel2,
    iDataDir,
    includeSymptomsURT2,
    interactionVariantDays2,
    postProcessModel2,
    sigma2,
    target,
):
    interactionVariantDaysSuffix2 = "InteractionVariantDays" if interactionVariantDays2 else ""
    binDaysPostOnsetVarSuffix2 = binDaysPostOnsetVar2[0].capitalize() + binDaysPostOnsetVar2[1:]
    symptomsVarSuffix2 = "SymptomsURT" if includeSymptomsURT2 else ""
    iDataModel2File = Path(iDataDir, 
                           f"model2{target}{familyModel2.capitalize()}{interactionVariantDaysSuffix2}"
                           f"{binDaysPostOnsetVarSuffix2}{symptomsVarSuffix2}Sigma{sigma2}.nc")
    if iDataModel2File.exists():
        iDataBmb2 = az.from_netcdf(iDataModel2File)
    else:
        iDataBmb2 = bmbModel2.fit(target_accept=0.99, draws=10000, tune=4000, random_seed=SEED, 
                                  idata_kwargs={'log_likelihood': True})
        postProcessModel2(iDataBmb2, interactionVariantDays=interactionVariantDays2, 
                          binDaysPostOnsetVar=binDaysPostOnsetVar2)
        iDataBmb2.to_netcdf(iDataModel2File)
    return (
        binDaysPostOnsetVarSuffix2,
        iDataBmb2,
        interactionVariantDaysSuffix2,
        symptomsVarSuffix2,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""##### Posterior analysis""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""##### Trace plot""")
    return


@app.cell
def _(az, iDataBmb2, plt):
    az.plot_trace(iDataBmb2, combined=True, var_names=["~zVlImmun2YN", "~daysImmun2YN", "~agrdt_mean"])
    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(pd):
    pd.set_option("display.max_rows", 50)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""##### Forest plot""")
    return


@app.cell
def _(
    binDaysPostOnsetVar2,
    includeSymptomsURT2,
    interactionVariantDays2,
    varNames2,
    varNames2InteractionVariantDays,
    varNames2InteractionVariantDays4,
    varNames2InteractionVariantDays4ReverseDays,
    varNames2InteractionVariantDays4ReverseDaysSymptomsURT,
    varNames2InteractionVariantDays4SymptomsURT,
    varNames2InteractionVariantDaysReverseDays,
    varNames2InteractionVariantDaysReverseDaysSymptomsURT,
    varNames2InteractionVariantDaysSymptomsURT,
    varNames2SymptomsURT,
):
    if interactionVariantDays2:
        if binDaysPostOnsetVar2 == "binDaysPostOnsetReverse":
            varNames2Tmp = (varNames2InteractionVariantDaysReverseDaysSymptomsURT if includeSymptomsURT2 else
                            varNames2InteractionVariantDaysReverseDays)
        elif binDaysPostOnsetVar2 == "binDaysPosOnset":
            varNames2Tmp = (varNames2InteractionVariantDaysSymptomsURT if includeSymptomsURT2 else 
                            varNames2InteractionVariantDays)
        if binDaysPostOnsetVar2 == "binDaysPostOnset4":
            varNames2Tmp = (varNames2InteractionVariantDays4SymptomsURT if includeSymptomsURT2 else 
                            varNames2InteractionVariantDays4)
        elif binDaysPostOnsetVar2 == "binDaysPostOnset4Reverse":
            varNames2Tmp = (varNames2InteractionVariantDays4ReverseDaysSymptomsURT if includeSymptomsURT2 else 
                            varNames2InteractionVariantDays4ReverseDays)
    else:
        varNames2Tmp = varNames2SymptomsURT if includeSymptomsURT2 else varNames2
    return (varNames2Tmp,)


@app.cell
def _(
    binDaysPostOnsetVarSuffix2,
    colWidth,
    familyModel2,
    iDataBmb2,
    interactionVariantDaysSuffix2,
    modelNoMapping,
    plotDir,
    plt,
    ridgeForestPlot,
    sigma2,
    symptomsVarSuffix2,
    target,
    varNames2Tmp,
):
    _ax = ridgeForestPlot(iDataBmb2, varNames=varNames2Tmp, figsize=(colWidth * 2, colWidth * 1.3), xlim=(-6, 6), 
                          addGrid=True)
    plt.savefig(f'{plotDir}/forestPlotModel{modelNoMapping[2]}{target}{familyModel2.capitalize()}'
                f'{interactionVariantDaysSuffix2}{binDaysPostOnsetVarSuffix2}{symptomsVarSuffix2}'
                f'Sigma{sigma2}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{plotDir}/forestPlotModel{modelNoMapping[2]}{target}{familyModel2.capitalize()}'
                f'{interactionVariantDaysSuffix2}{binDaysPostOnsetVarSuffix2}{symptomsVarSuffix2}Sigma{sigma2}.pdf', 
                dpi=300, bbox_inches='tight')
    _ax
    return


@app.cell
def _(bmbModel2, defaultdict, dfLogisticBmb2, iDataBmb2, predictionsNewData):
    predDictBmb2 = defaultdict(dict)
    predictionsNewData(dfLogisticBmb2, bmbModel2, iDataBmb2, predDictBmb2, "immun2YN", (0, 1))
    predictionsNewData(dfLogisticBmb2, bmbModel2, iDataBmb2, predDictBmb2, "binDaysPostOnset", (0, 1))
    predictionsNewData(dfLogisticBmb2, bmbModel2, iDataBmb2, predDictBmb2, "variant", (0, 1, 2, 3))
    return (predDictBmb2,)


@app.cell
def _(
    ANNOTATION_LETTER_SIZE,
    DINA4_HEIGHT,
    annotateWithLetters,
    colWidth,
    dfLogisticBmb2,
    modelNoMapping,
    plotDir,
    plt,
    predDictBmb2,
    saveFigure,
    spaghettiPlotCategorical,
):
    model2Spagh, _ax = plt.subplots(3, 1, figsize=(colWidth * 2, DINA4_HEIGHT), sharex=True)
    _ax1, _ax2, _ax3 = _ax.flatten()
    spaghettiPlotCategorical(_ax1, 'variant', (0, 1, 2, 3), predDictBmb2, dfLogisticBmb2, showXLabel=False, 
                             legendLoc=(0.565, 0.15))
    spaghettiPlotCategorical(_ax2, 'binDaysPostOnset', (0, 1), predDictBmb2, dfLogisticBmb2, showXLabel=False)
    plt.tight_layout()
    spaghettiPlotCategorical(_ax3, 'immun2YN', (0, 1), predDictBmb2, dfLogisticBmb2)

    annotateWithLetters((_ax1, _ax2, _ax3), size=ANNOTATION_LETTER_SIZE, coords=(-0.15, 1))
    saveFigure(model2Spagh, plotDir / f'model{modelNoMapping[2]}Spagh.png')
    saveFigure(model2Spagh, plotDir / f'model{modelNoMapping[2]}Spagh.pdf')
    return


@app.cell
def _(
    az,
    iDataBmb2,
    modelNoMapping,
    tableDirIData,
    varNames2Tmp,
    writeIDataSummaryTable,
    writeIDataSummaryTableLatex,
):
    summaryDf2 = az.summary(iDataBmb2, var_names=['~agrdt_mean', '~p'])
    writeIDataSummaryTable(summaryDf2, varnames=varNames2Tmp, 
                           outfile=tableDirIData / f'iDataSummaryModel{modelNoMapping[2]}.tsv')
    writeIDataSummaryTableLatex(summaryDf2, varnames=varNames2Tmp, 
                                outfile=tableDirIData / f'iDataSummaryModel{modelNoMapping[2]}Latex.txt', model_no=2)
    return


@app.cell
def _(mo):
    mo.md(r"""### Include immunization but no variant variable (Model 3)""")
    return


@app.cell
def _(bmb, dfPos, logisticRegressionDf):
    includeSymptomsURT3 = False
    binDaysPostOnsetVar3 = "binDaysPostOnset4"
    familyModel3 = "cumulative" # bernoulli
    _outcome = "testline" if familyModel3 == "cumulative" else "agrdt"

    _indVars = ("vl", "age", "testline") if familyModel3=="cumulative" else ("vl", "age")
    _indCatVars = (("testDevice", "immun2YN", "gender", binDaysPostOnsetVar3, "recovered", "illURT") if 
                   includeSymptomsURT3 else ("testDevice", "immun2YN", "gender", binDaysPostOnsetVar3, "recovered"))

    immun2YNDaysInteraction3=True
    _daysVar = f"{binDaysPostOnsetVar3}:immun2YN" if immun2YNDaysInteraction3 else ""

    formulaSymptomsSuffix3 = "+ illURT" if includeSymptomsURT3 else ""
    _formulaAllInd = (f'{_outcome} ~ testDevice + immun2YN + gender + {binDaysPostOnsetVar3} + {_daysVar} + zVl + '
                      f'zAge + recovered {formulaSymptomsSuffix3}')

    dfLogisticBmb3 = logisticRegressionDf(dfPos, immunVar='immun2YN', indCatVars=_indCatVars, indVars=_indVars)

    _priors = {'zVl': bmb.Prior('Lognormal', mu=1, sigma=1), 
               'binDaysPostOnset:immun2YN': bmb.Prior('Normal', mu=0, sigma=1), 
               'binDaysPostOnset4:immun2YN': bmb.Prior('Normal', mu=0, sigma=1),
               'binDaysPostOnset4': bmb.Prior('Normal', mu=0, sigma=1),
               'binDaysPostOnset': bmb.Prior('Normal', mu=0, sigma=1),
               'testDevice': bmb.Prior('Normal', mu=-0.5, sigma=0.5), 
               'immun2YN': bmb.Prior('Normal', mu=0, sigma=1),
               'illURT': bmb.Prior('Normal', mu=0, sigma=1),
               'zAge': bmb.Prior('Normal', mu=0, sigma=1), 
               'gender': bmb.Prior('Normal', mu=0, sigma=1),
               'recovered': bmb.Prior('Normal', mu=0, sigma=1),}
    bmbModel3 = bmb.Model(_formulaAllInd, dfLogisticBmb3, priors=_priors, family=familyModel3, categorical=list(_indCatVars), 
                          noncentered=True)
    bmbModel3.build()
    return (
        binDaysPostOnsetVar3,
        bmbModel3,
        dfLogisticBmb3,
        familyModel3,
        immun2YNDaysInteraction3,
        includeSymptomsURT3,
    )


@app.cell
def _(bmbModel3):
    bmbModel3
    return


@app.cell
def _(bmbModel3):
    bmbModel3.graph()
    return


@app.cell
def _(
    Path,
    SEED,
    az,
    binDaysPostOnsetVar3,
    bmbModel3,
    familyModel3,
    iDataDir,
    immun2YNDaysInteraction3,
    includeSymptomsURT3,
    postProcessModel3,
    target,
):
    immun2YNDaysInteractionSuffix3 = "Immun2YNDaysInteraction" if immun2YNDaysInteraction3 else ""
    binDaysPostOnsetVarSuffix3 = binDaysPostOnsetVar3[0].capitalize() + binDaysPostOnsetVar3[1:]
    symptomsVarSuffix3 = "SymptomsURT" if includeSymptomsURT3 else ""
    iDataModel3File = Path(iDataDir, 
                           f"model3{target}{familyModel3.capitalize()}{immun2YNDaysInteractionSuffix3}"
                           f"{binDaysPostOnsetVarSuffix3}{symptomsVarSuffix3}.nc")

    if iDataModel3File.exists():
        iDataBmb3 = az.from_netcdf(iDataModel3File)
    else:
        iDataBmb3 = bmbModel3.fit(target_accept=0.99, draws=10000, tune=4000, random_seed=SEED, 
                                  idata_kwargs={'log_likelihood': True})
        postProcessModel3(iDataBmb3, binDaysPostOnsetVar=binDaysPostOnsetVar3)
        iDataBmb3.to_netcdf(iDataModel3File)
    return (
        binDaysPostOnsetVarSuffix3,
        iDataBmb3,
        immun2YNDaysInteractionSuffix3,
        symptomsVarSuffix3,
    )


@app.cell
def _(az, iDataBmb3):
    az.summary(iDataBmb3)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""##### Trace plot""")
    return


@app.cell
def _(az, iDataBmb3, plt):
    az.plot_trace(iDataBmb3, combined=True);
    plt.tight_layout()
    plt.show()
    return


@app.cell
def _(pd):
    pd.set_option("display.max_rows", 160)
    return


@app.cell
def _(
    az,
    iDataBmb3,
    modelNoMapping,
    tableDirIData,
    varNames3Tmp,
    writeIDataSummaryTable,
    writeIDataSummaryTableLatex,
):
    summaryDf3 = az.summary(iDataBmb3, var_names=['~agrdt_mean', '~p'])
    writeIDataSummaryTable(summaryDf3, varnames=varNames3Tmp, 
                           outfile=tableDirIData / f'iDataSummaryModel{modelNoMapping[3]}.tsv')
    writeIDataSummaryTableLatex(summaryDf3, varnames=varNames3Tmp, 
                                outfile=tableDirIData / f'iDataSummaryModel{modelNoMapping[3]}Latex.txt',
                                model_no=3)
    return (summaryDf3,)


@app.cell
def _(summaryDf3):
    summaryDf3
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""##### Forest plot""")
    return


@app.cell
def _(
    binDaysPostOnsetVar3,
    immun2YNDaysInteraction3,
    includeSymptomsURT2,
    includeSymptomsURT3,
    varNames3,
    varNames3Immun2YNDays4Interaction,
    varNames3Immun2YNDays4InteractionReverseDays,
    varNames3Immun2YNDays4InteractionReverseDaysSymptomsURT,
    varNames3Immun2YNDays4InteractionSymptomsURT,
    varNames3Immun2YNDaysInteraction,
    varNames3Immun2YNDaysInteractionReverseDays,
    varNames3Immun2YNDaysInteractionReverseDaysSymptomsURT,
    varNames3Immun2YNDaysInteractionSymptomsURT,
    varNames3SymptomsURT,
):
    if immun2YNDaysInteraction3:
        if binDaysPostOnsetVar3 == "binDaysPostOnsetReverse":
            varNames3Tmp = (varNames3Immun2YNDaysInteractionReverseDaysSymptomsURT if includeSymptomsURT3 else
                            varNames3Immun2YNDaysInteractionReverseDays)
        elif binDaysPostOnsetVar3 == "binDaysPostOnset":
            varNames3Tmp = (varNames3Immun2YNDaysInteractionSymptomsURT if includeSymptomsURT3 else 
                            varNames3Immun2YNDaysInteraction)
        if binDaysPostOnsetVar3 == "binDaysPostOnset4":
            varNames3Tmp = (varNames3Immun2YNDays4InteractionSymptomsURT if includeSymptomsURT2 else 
                            varNames3Immun2YNDays4Interaction)
        elif binDaysPostOnsetVar3 == "binDaysPostOnset4Reverse":
            varNames3Tmp = (varNames3Immun2YNDays4InteractionReverseDaysSymptomsURT if includeSymptomsURT3 else 
                            varNames3Immun2YNDays4InteractionReverseDays)
    else:
        varNames3Tmp = varNames3SymptomsURT if includeSymptomsURT3 else varNames3
    return (varNames3Tmp,)


@app.cell
def _(
    binDaysPostOnsetVarSuffix3,
    colWidth,
    familyModel3,
    iDataBmb3,
    immun2YNDaysInteractionSuffix3,
    modelNoMapping,
    plotDir,
    plt,
    ridgeForestPlot,
    symptomsVarSuffix3,
    target,
    varNames3Tmp,
):
    _ax = ridgeForestPlot(iDataBmb3, varNames=varNames3Tmp, figsize=(colWidth * 2, colWidth * 1), addGrid=True)
    _ax.tick_params(axis='y', which='both', labelleft=False, labelright=True)

    plt.savefig(f'{plotDir}/forestPlotModel{modelNoMapping[3]}{target}{familyModel3.capitalize()}'
                f'{immun2YNDaysInteractionSuffix3}{binDaysPostOnsetVarSuffix3}{symptomsVarSuffix3}.png', 
                dpi=300, bbox_inches='tight')
    plt.savefig(f'{plotDir}/forestPlotModel{modelNoMapping[3]}{target}{familyModel3.capitalize()}'
                f'{immun2YNDaysInteractionSuffix3}{binDaysPostOnsetVarSuffix3}{symptomsVarSuffix3}.pdf', 
                dpi=300, bbox_inches='tight')
    _ax
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""##### Spaghetti plots (Figure 5)""")
    return


@app.cell
def _(defaultdict):
    predDictBmb3 = defaultdict(dict)
    return (predDictBmb3,)


@app.cell
def _(bmbModel3, dfLogisticBmb3, iDataBmb3, predDictBmb3, predictionsNewData):
    predictionsNewData(dfLogisticBmb3, bmbModel3, iDataBmb3, predDictBmb3, "immun2YN", (0, 1))
    predictionsNewData(dfLogisticBmb3, bmbModel3, iDataBmb3, predDictBmb3, "binDaysPostOnset", (0, 1))
    return


@app.cell
def _(
    annotateWithLetter,
    colWidth,
    dfLogisticBmb3,
    fontSizePlot,
    modelNoMapping,
    plotDir,
    plt,
    predDictBmb3,
    saveFigure,
    spaghettiPlotCategorical,
):
    _figModel3A, _ax = plt.subplots(1, 1, figsize=(colWidth * 2, colWidth * 1.25))
    spaghettiPlotCategorical(_ax, 'binDaysPostOnset', (0, 1), predDictBmb3, dfLogisticBmb3)
    annotateWithLetter(_ax, "B", size=fontSizePlot, coords=(-0.15, 1))
    saveFigure(_figModel3A, plotDir / f'model{modelNoMapping[3]}Spagh-B.png')
    saveFigure(_figModel3A, plotDir / f'model{modelNoMapping[3]}Spagh-B.pdf')
    return


@app.cell
def _(
    annotateWithLetter,
    colWidth,
    dfLogisticBmb3,
    fontSizePlot,
    modelNoMapping,
    plotDir,
    plt,
    predDictBmb3,
    saveFigure,
    spaghettiPlotCategorical,
):
    _figModel3B, _ax = plt.subplots(1, 1, figsize=(colWidth * 2, colWidth * 1.25))
    spaghettiPlotCategorical(_ax, 'immun2YN', (0, 1), predDictBmb3, dfLogisticBmb3, legendLoc=(0.55, 0.15))
    annotateWithLetter(_ax, "A", size=fontSizePlot, coords=(-0.15, 1))
    saveFigure(_figModel3B, plotDir / f'model{modelNoMapping[3]}Spagh-A.png')
    saveFigure(_figModel3B, plotDir / f'model{modelNoMapping[3]}Spagh-A.pdf')
    return


@app.cell
def _(
    ANNOTATION_COORDS,
    annotateWithLetters,
    colWidth,
    dfLogisticBmb3,
    fontSizePlot,
    modelNoMapping,
    plotDir,
    plt,
    predDictBmb3,
    saveFigure,
    setFontSize,
    spaghettiPlotCategorical,
):
    fig5, _ax = plt.subplots(2, 1, figsize=(colWidth * 2, colWidth * 2.5), constrained_layout=True)
    _ax1, _ax2 = _ax.flatten()
    spaghettiPlotCategorical(_ax2, 'binDaysPostOnset', (0, 1), predDictBmb3, dfLogisticBmb3)
    spaghettiPlotCategorical(_ax1, 'immun2YN', (0, 1), predDictBmb3, dfLogisticBmb3, legendLoc=(0.55, 0.15))
    for _a in (_ax1, _ax2):
        setFontSize(_a, size=fontSizePlot)
    annotateWithLetters(_ax, size=fontSizePlot, coords=ANNOTATION_COORDS)
    saveFigure(fig5, plotDir / f'model{modelNoMapping[3]}Spagh.png')
    saveFigure(fig5, plotDir / f'model{modelNoMapping[3]}Spagh.pdf')
    return


@app.cell
def _(mo):
    mo.md(r"""## Include asymptomatic infections""")
    return


@app.cell
def _(mo):
    mo.md(r"""### Include variable for SARS-CoV-2 variant (Model 4)""")
    return


@app.cell
def _(bmb, dfPos, logisticRegressionDf, removeReleaseTesting):
    _indCatVars = ("testDevice", "gender", "variant", "symptoms")

    familyModel4 = "cumulative" # bernoulli
    _indVars = ("vl", "age", "testline") if familyModel4=="cumulative" else ("vl", "age")


    _formulaAllIndBernoulli = f'agrdt ~ testDevice + variant:symptoms + gender + symptoms + zVl + zAge'
    _formulaAllIndCumulative = f'testline ~ testDevice + variant:symptoms + gender + symptoms + zVl + zAge'
    _formulaAllInd = _formulaAllIndCumulative if familyModel4 == "cumulative" else _formulaAllIndBernoulli

    dfLogisticBmb4 = logisticRegressionDf(removeReleaseTesting(dfPos), indCatVars=_indCatVars,
                                           includeAsymptomatic=True, indVars=_indVars)
    sigma4 = 1
    _priors = {'zVl': bmb.Prior('Lognormal', mu=1, sigma=1),
               'symptoms': bmb.Prior('Normal', mu=0, sigma=sigma4),
               'variant:symptoms': bmb.Prior('Normal', mu=0, sigma=sigma4),
               'testDevice': bmb.Prior('Normal', mu=-0.5, sigma=0.5),
               'zAge': bmb.Prior('Normal', mu=0, sigma=1),
               'gender': bmb.Prior('Normal', mu=0, sigma=1),
                }
    bmbModel4 = bmb.Model(_formulaAllInd, dfLogisticBmb4, priors=_priors, family=familyModel4,
                          categorical=list(_indCatVars), noncentered=True)
    bmbModel4.build()
    return bmbModel4, dfLogisticBmb4, familyModel4


@app.cell
def _(bmbModel4):
    bmbModel4
    return


@app.cell
def _(
    Path,
    SEED,
    az,
    bmbModel4,
    familyModel4,
    iDataDir,
    postProcessModel4,
    target,
):
    iDataModel4File = Path(iDataDir, f"model4{target}{familyModel4.capitalize()}.nc")
    if iDataModel4File.exists():
        iDataBmb4 = az.from_netcdf(iDataModel4File)
    else:
        iDataBmb4 = bmbModel4.fit(target_accept=0.99, draws=10000, tune=4000, random_seed=SEED,
                                    idata_kwargs={'log_likelihood': True})
        postProcessModel4(iDataBmb4)
        iDataBmb4.to_netcdf(iDataModel4File)
    return (iDataBmb4,)


@app.cell
def _(az, iDataBmb4):
    az.summary(iDataBmb4)
    return


@app.cell
def _(varNames4):
    varNames4Tmp = varNames4
    return (varNames4Tmp,)


@app.cell
def _(
    colWidth,
    familyModel4,
    iDataBmb4,
    modelNoMapping,
    plotDir,
    plt,
    ridgeForestPlot,
    target,
    varNames4Tmp,
):
    _ax = ridgeForestPlot(iDataBmb4, varNames=varNames4Tmp, figsize=(colWidth * 2, colWidth * 1.2), addGrid=True)
    _ax.tick_params(axis='y', which='both', labelleft=False, labelright=True)

    plt.savefig(f'{plotDir}/forestPlotModel{modelNoMapping[4]}{target}{familyModel4.capitalize()}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{plotDir}/forestPlotModel{modelNoMapping[4]}{target}{familyModel4.capitalize()}.pdf', dpi=300, bbox_inches='tight')

    _ax
    return


@app.cell
def _(
    az,
    iDataBmb4,
    modelNoMapping,
    tableDirIData,
    varNames4Tmp,
    writeIDataSummaryTable,
    writeIDataSummaryTableLatex,
):
    summaryDf4 = az.summary(iDataBmb4, var_names=['~agrdt_mean', '~p'])
    writeIDataSummaryTable(summaryDf4, varnames=varNames4Tmp,
                           outfile=tableDirIData / f'iDataSummaryModel{modelNoMapping[4]}.tsv')
    writeIDataSummaryTableLatex(summaryDf4, varnames=varNames4Tmp,
                                outfile=tableDirIData / f'iDataSummaryModel{modelNoMapping[4]}Latex.txt')
    return


@app.cell
def _(bmbModel4, defaultdict, dfLogisticBmb4, iDataBmb4, predictionsNewData):
    predDictBmb4Symp = defaultdict(dict)
    predDictBmb4Asymp = defaultdict(dict)
    predictionsNewData(dfLogisticBmb4, bmbModel4, iDataBmb4, predDictBmb4Asymp, "variant", (0, 1, 2, 3),
                       otherFeaturesDict={"symptoms":0})
    predictionsNewData(dfLogisticBmb4, bmbModel4, iDataBmb4, predDictBmb4Symp, "variant", (0, 1, 2, 3),
                       otherFeaturesDict={"symptoms":1})
    return predDictBmb4Asymp, predDictBmb4Symp


@app.cell
def _(
    ANNOTATION_COORDS,
    ANNOTATION_LETTER_SIZE,
    annotateWithLetters,
    colWidth,
    dfLogisticBmb4,
    modelNoMapping,
    pal,
    plotDir,
    plt,
    predDictBmb4Asymp,
    predDictBmb4Symp,
    saveFigure,
    spaghettiPlotCategorical,
):
    figModel4Spagh, _ax = plt.subplots(2, 1, figsize=(colWidth * 2, colWidth * 2), sharex=True, constrained_layout=True)
    _ax1, _ax2 = _ax.flatten()
    spaghettiPlotCategorical(_ax2, 'variant', (0, 1, 2, 3), predDictBmb4Asymp,
                             dfLogisticBmb4[dfLogisticBmb4.symptoms==0], showXLabel=True,
                             legendLoc=(0.565, 0.15), showLegend=False, palette=list(pal.variant.values()))
    spaghettiPlotCategorical(_ax1, 'variant', (0, 1, 2, 3), predDictBmb4Symp, dfLogisticBmb4[dfLogisticBmb4.symptoms==1],
                             showXLabel=False, legendLoc=(0.565, 0.15), palette=list(pal.variant.values()))

    annotateWithLetters((_ax1, _ax2), size=ANNOTATION_LETTER_SIZE, coords=ANNOTATION_COORDS)
    saveFigure(figModel4Spagh, plotDir / f'model{modelNoMapping[4]}Spagh.png')
    saveFigure(figModel4Spagh, plotDir / f'model{modelNoMapping[4]}Spagh.pdf')
    figModel4Spagh
    return


@app.cell
def _(mo):
    mo.md(r"""### Include both variant and immunization status (Model 5)""")
    return


@app.cell
def _(bmb, dfPos, logisticRegressionDf, removeReleaseTesting):
    familyModel5 = "cumulative" # bernoulli
    interactionVariantSymptoms5 = True
    interactionImmun2YNSymptoms5 = True
    symptomsVar5 = "symptoms"#"illURT"#"illURTReverse"

    _variantVar = f"variant:{symptomsVar5}" if interactionVariantSymptoms5 else "variant"
    _immun2YNVar = f"immun2YN:{symptomsVar5}" if interactionImmun2YNSymptoms5 else "immun2YN"

    _indCatVars = ("testDevice", "gender", "variant", "immun2YN", "symptoms", "recovered")
    _indVars = ("vl", "age", "testline") if familyModel5=="cumulative" else ("vl", "age")

    _formulaAllIndBernoulli = (f'agrdt ~ testDevice + {_variantVar} + {_immun2YNVar } + gender + {symptomsVar5} + zAge + '
                               f'zVl + recovered:symptoms')
    _formulaAllIndCumulative = (f'testline ~ testDevice + {_variantVar} + {_immun2YNVar } + gender + {symptomsVar5} + zAge '
                                f'+ recovered:symptoms + zVl')
    _formulaAllInd = _formulaAllIndCumulative if familyModel5 == "cumulative" else _formulaAllIndBernoulli
    dfLogisticBmb5 = logisticRegressionDf(removeReleaseTesting(dfPos), immunVar='immun2YN', indCatVars=_indCatVars,
                                           indVars=_indVars, includeAsymptomatic=True)

    sigma5 = 1
    _priors = {'zVl': bmb.Prior('Lognormal', mu=1, sigma=1),
               'immun2YN': bmb.Prior('Normal', mu=0, sigma=sigma5),
               'immun2YN:symptoms': bmb.Prior('Normal', mu=0, sigma=sigma5),
               'variant:symptoms': bmb.Prior('Normal', mu=0, sigma=sigma5),
               'testDevice': bmb.Prior('Normal', mu=-0.5, sigma=0.5),
               'zAge': bmb.Prior('Normal', mu=0, sigma=1),
               'gender': bmb.Prior('Normal', mu=0, sigma=1),
               'symptoms': bmb.Prior('Normal', mu=0, sigma=sigma5),
               'recovered': bmb.Prior('Normal', mu=0, sigma=sigma5),
               'recovered:symptoms': bmb.Prior('Normal', mu=0, sigma=sigma5),}
    bmbModel5 = bmb.Model(_formulaAllInd, dfLogisticBmb5, priors=_priors, family=familyModel5,
                           categorical=list(_indCatVars), noncentered=True, dropna=True)
    bmbModel5.build()
    return (
        bmbModel5,
        dfLogisticBmb5,
        familyModel5,
        interactionImmun2YNSymptoms5,
        interactionVariantSymptoms5,
        sigma5,
    )


@app.cell
def _(bmbModel5):
    bmbModel5
    return


@app.cell
def _(
    Path,
    SEED,
    az,
    bmbModel5,
    familyModel5,
    iDataDir,
    interactionImmun2YNSymptoms5,
    interactionVariantSymptoms5,
    postProcessModel5,
    sigma5,
    target,
):
    interactionVariantSymptomsSuffix5 = "InteractionVariantSymptoms" if interactionVariantSymptoms5 else ""
    interactionImmun2YNSymptomsSuffix5 = "InteractionImmun2YNSymptoms" if interactionImmun2YNSymptoms5 else ""

    iDataModel5File = Path(iDataDir, f"model5{target}{familyModel5.capitalize()}{interactionVariantSymptomsSuffix5}"
                            f"{interactionImmun2YNSymptomsSuffix5}Sigma{sigma5}InteractionRecoveredSymptoms.nc")
    if iDataModel5File.exists():
        iDataBmb5 = az.from_netcdf(iDataModel5File)
    else:
        iDataBmb5 = bmbModel5.fit(target_accept=0.99, draws=10000, tune=4000, random_seed=SEED,
                                    idata_kwargs={'log_likelihood': True})
        postProcessModel5(iDataBmb5, interactionVariantSymptoms=interactionVariantSymptoms5,
                          interactionImmun2YNSymptoms=interactionImmun2YNSymptoms5)
        iDataBmb5.to_netcdf(iDataModel5File)
    return (
        iDataBmb5,
        interactionImmun2YNSymptomsSuffix5,
        interactionVariantSymptomsSuffix5,
    )


@app.cell
def _(az, iDataBmb5):
    az.summary(iDataBmb5)
    return


@app.cell
def _(
    interactionVariantSymptoms5,
    varNames5,
    varNames5VariantSymptomsInteraction,
):
    varNames5Tmp = varNames5VariantSymptomsInteraction if interactionVariantSymptoms5 else varNames5
    return (varNames5Tmp,)


@app.cell
def _(
    colWidth,
    familyModel5,
    iDataBmb5,
    interactionImmun2YNSymptomsSuffix5,
    interactionVariantSymptomsSuffix5,
    modelNoMapping,
    plotDir,
    plt,
    ridgeForestPlot,
    sigma5,
    target,
    varNames5Tmp,
):
    _ax = ridgeForestPlot(iDataBmb5, varNames=varNames5Tmp, figsize=(colWidth * 2, colWidth * 1.7), addGrid=True)
    _ax.tick_params(axis='y', which='both', labelleft=False, labelright=True)

    plt.savefig(f'{plotDir}/forestPlotModel{modelNoMapping[5]}{target}{familyModel5.capitalize()}'
                f'{interactionVariantSymptomsSuffix5}{interactionImmun2YNSymptomsSuffix5}Sigma{sigma5}'
                f'InteractionRecoveredSymptoms.png', dpi=300,
                bbox_inches='tight')
    plt.savefig(f'{plotDir}/forestPlotModel{modelNoMapping[5]}{target}{familyModel5.capitalize()}'
                f'{interactionVariantSymptomsSuffix5}{interactionImmun2YNSymptomsSuffix5}Sigma{sigma5}'
                f'InteractionRecoveredSymptoms.pdf', dpi=300,
                bbox_inches='tight')
    _ax
    return


@app.cell
def _(
    colWidth,
    dfLogisticBmb5,
    fontSizePlot,
    legend,
    negPointsYCoords,
    pal,
    plotFigA4_right,
    plt,
):
    fig_A4_right_asymp_symp, _ax = plt.subplots(1, 1, figsize=(colWidth, colWidth))
    plotFigA4_right(dfLogisticBmb5, negPointsYCoords, palette=pal, legend=legend, ax=_ax, fontSizePlot=fontSizePlot,
                    markersize=3)
    _ax
    return


@app.cell
def _(dfLogisticBmb5):
    ((dfLogisticBmb5.symptoms == 0) & (dfLogisticBmb5.immun2YN == 1)).sum()
    return


app._unparsable_cell(
    r"""
    ((dfLogisticBmb5.symptoms == 0) & (dfLogisticBmb5.variant == )).sum()
    """,
    name="_"
)


@app.cell
def _(dfLogisticBmb5):
    ((dfLogisticBmb5.symptoms == 0) & (dfLogisticBmb5.immun2YN == 0)).sum()
    return


@app.cell
def _(dfLogisticBmb5):
    (dfLogisticBmb5.symptoms == 0).sum()
    return


@app.cell
def _(
    colWidth,
    dfLogisticBmb5,
    fontSizePlot,
    label,
    legend,
    negPointsYCoords,
    pal,
    plotDataPointsRegression,
    plt,
):
    figA4AsympSymp, _ax = plt.subplots(1, 2, figsize=(2 * colWidth, colWidth), sharey=True)
    _ax1, _ax2 = _ax.flatten()
    plotDataPointsRegression(dfLogisticBmb5, negPointsYCoords, _ax1, _ax2, pal=pal, label=label, legend=legend,
                             fig=figA4AsympSymp, fontSizePlot=fontSizePlot, markersize=3)
    # Remove x-tick marks
    _ax1.tick_params(axis="x", which="both", bottom=False, top=False)
    _ax2.tick_params(axis="x", which="both", bottom=False, top=False)
    plt.tight_layout()
    _ax
    # saveFigure(figA4, plotDir / 'FigureA4.png')
    # saveFigure(figA4, plotDir / 'FigureA4.pdf')
    return


@app.cell
def _(bmbModel5, defaultdict, dfLogisticBmb5, iDataBmb5, predictionsNewData):
    predDictBmb5Symp  = defaultdict(dict)
    predDictBmb5Asymp = defaultdict(dict)
    predictionsNewData(dfLogisticBmb5, bmbModel5, iDataBmb5, predDictBmb5Asymp, "immun2YN", (0, 1), 
                       otherFeaturesDict={"symptoms":0, "recovered": 0})
    predictionsNewData(dfLogisticBmb5, bmbModel5, iDataBmb5, predDictBmb5Symp , "immun2YN", (0, 1), 
                       otherFeaturesDict={"symptoms":1, "recovered": 0})
    predictionsNewData(dfLogisticBmb5, bmbModel5, iDataBmb5, predDictBmb5Symp, "variant", (0, 1, 2, 3),
                       otherFeaturesDict={"symptoms":1, "recovered": 0})
    return predDictBmb5Asymp, predDictBmb5Symp


@app.cell
def _(
    ANNOTATION_LETTER_SIZE,
    DINA4_HEIGHT,
    annotateWithLetters,
    colWidth,
    dfLogisticBmb5,
    modelNoMapping,
    plotDir,
    plt,
    predDictBmb5Asymp,
    predDictBmb5Symp,
    saveFigure,
    spaghettiPlotCategorical,
):
    figModel5Spagh, _ax = plt.subplots(3, 1, figsize=(colWidth * 2, DINA4_HEIGHT), sharex=True)
    _ax1, _ax2, _ax3 = _ax.flatten()
    spaghettiPlotCategorical(_ax1, 'variant', (0, 1, 2, 3), predDictBmb5Symp , dfLogisticBmb5, showXLabel=False,
                             legendLoc=(0.565, 0.15))
    spaghettiPlotCategorical(_ax2, 'immun2YN', (0, 1), predDictBmb5Asymp, dfLogisticBmb5, showXLabel=False)
    spaghettiPlotCategorical(_ax3, 'immun2YN', (0, 1), predDictBmb5Symp , dfLogisticBmb5)
    plt.tight_layout()

    annotateWithLetters((_ax1, _ax2, _ax3), size=ANNOTATION_LETTER_SIZE, coords=(-0.15, 1))
    saveFigure(figModel5Spagh, plotDir / f'model{modelNoMapping[5]}Spagh.png')
    saveFigure(figModel5Spagh, plotDir / f'model{modelNoMapping[5]}Spagh.pdf')
    return


@app.cell
def _(
    az,
    iDataBmb5,
    modelNoMapping,
    tableDirIData,
    varNames5Tmp,
    writeIDataSummaryTable,
    writeIDataSummaryTableLatex,
):
    summaryDf5 = az.summary(iDataBmb5, var_names=['~agrdt_mean', '~p'])
    writeIDataSummaryTable(summaryDf5, varnames=varNames5Tmp,
                           outfile=tableDirIData / f'iDataSummaryModel{modelNoMapping[5]}.tsv')
    writeIDataSummaryTableLatex(summaryDf5, varnames=varNames5Tmp,
                                outfile=tableDirIData / f'iDataSummaryModel{modelNoMapping[5]}Latex.txt',
                                model_no=10)
    return


@app.cell
def _(mo):
    mo.md(r"""### Include immunization but no variant variable (Model 6)""")
    return


@app.cell
def _(bmb, dfPos, logisticRegressionDf, removeReleaseTesting):
    # Note: we remove people from release testing to not capture asymptomatic people who are at the end of their infection (would be unproportionally many)
    # We don't do that when regressing on only symptomatic people because there we can control for the time of testing (day 1-7)
    interactionImmun2YNSymptoms6 = True
    familyModel6 = "cumulative" # bernoulli
    symptomsVar6 = "symptoms"#"illURT"#"illURTReverse"
    _indVars = ("vl", "age", "testline") if familyModel6=="cumulative" else ("vl", "age")
    _indCatVars = ("testDevice", "gender", "immun2YN", symptomsVar6, "recovered")
    _immun2YNSymptomsVar = f"immun2YN:{symptomsVar6}" if interactionImmun2YNSymptoms6 else "immun2YN"

    _formulaAllIndBernoulli = (f'agrdt ~ testDevice + {_immun2YNSymptomsVar} + gender + {symptomsVar6} + zAge + zVl + '
                               f'recovered:{symptomsVar6}')
    _formulaAllIndCumulative = (f'testline ~ testDevice + {_immun2YNSymptomsVar} + gender + {symptomsVar6} + zAge + zVl + '
                                f'recovered:{symptomsVar6}')
    _formulaAllInd = _formulaAllIndCumulative if familyModel6 == "cumulative" else _formulaAllIndBernoulli

    dfLogisticBmb6 = logisticRegressionDf(removeReleaseTesting(dfPos), immunVar='immun2YN', indCatVars=_indCatVars, 
                                          indVars=_indVars, includeAsymptomatic=True)
    sigma6 = 1
    _priors = {'zVl': bmb.Prior('Lognormal', mu=1, sigma=1), 
               'immun2YN:symptoms': bmb.Prior('Normal', mu=0, sigma=sigma6), 
               'immun2YN:illURT': bmb.Prior('Normal', mu=0, sigma=sigma6), 
               'testDevice': bmb.Prior('Normal', mu=-0.5, sigma=0.5), 
               'zAge': bmb.Prior('Normal', mu=0, sigma=1), 
               'gender': bmb.Prior('Normal', mu=0, sigma=1),
               'symptoms': bmb.Prior('Normal', mu=0, sigma=sigma6),
               'illURT': bmb.Prior('Normal', mu=0, sigma=sigma6),
               'recovered': bmb.Prior('Normal', mu=0, sigma=sigma6),
               'recovered:symptoms': bmb.Prior('Normal', mu=0, sigma=sigma6)}
    bmbModel6 = bmb.Model(_formulaAllInd, dfLogisticBmb6, priors=_priors, family=familyModel6, 
                          categorical=list(_indCatVars), noncentered=True)
    bmbModel6.build()
    return (
        bmbModel6,
        dfLogisticBmb6,
        familyModel6,
        interactionImmun2YNSymptoms6,
        symptomsVar6,
    )


@app.cell
def _(bmbModel6):
    bmbModel6
    return


@app.cell
def _(
    Path,
    SEED,
    az,
    bmbModel6,
    familyModel6,
    iDataDir,
    interactionImmun2YNSymptoms6,
    postProcessModel6,
    symptomsVar6,
    target,
):
    interactionImmun2YNSymptomsSuffix6 = "InteractionImmun2YNSymptoms" if interactionImmun2YNSymptoms6 else ""
    iDataModel6Path = Path(iDataDir, f"model6{target}{familyModel6.capitalize()}{interactionImmun2YNSymptomsSuffix6}"
                           f"{symptomsVar6.capitalize()}InteractionRecoveredSymptoms.nc")

    if iDataModel6Path.exists():
        iDataBmb6 = az.from_netcdf(iDataModel6Path)
    else:
        iDataBmb6 = bmbModel6.fit(target_accept=0.99, draws=10000, tune=4000, random_seed=SEED, 
                                  idata_kwargs={'log_likelihood': True})
        postProcessModel6(iDataBmb6, symptomsVar=symptomsVar6, interactionRecoveredSymptoms=True)
        iDataBmb6.to_netcdf(iDataModel6Path)
    return iDataBmb6, interactionImmun2YNSymptomsSuffix6


@app.cell
def _(az, iDataBmb6):
    az.summary(iDataBmb6)
    return


@app.cell
def _(symptomsVar6, varNames6, varNames6SymptomsURT):
    varNames6Tmp = varNames6SymptomsURT if symptomsVar6 == "illURT" else varNames6
    return (varNames6Tmp,)


@app.cell
def _(
    colWidth,
    familyModel6,
    iDataBmb6,
    interactionImmun2YNSymptomsSuffix6,
    modelNoMapping,
    plotDir,
    plt,
    ridgeForestPlot,
    symptomsVar6,
    target,
    varNames6Tmp,
):
    _ax = ridgeForestPlot(iDataBmb6, varNames=varNames6Tmp, figsize=(colWidth * 2, colWidth * 1.2), addGrid=True)
    _ax.tick_params(axis='y', which='both', labelleft=False, labelright=True)
    #annotateWithLetter(_ax, "E", coords=(-0.15, 1.05), size=fontSizePlot)

    plt.savefig(f'{plotDir}/forestPlotModel{modelNoMapping[6]}{target}{familyModel6.capitalize()}'
                f'{interactionImmun2YNSymptomsSuffix6}{symptomsVar6.capitalize()}InteractionRecoveredSymptoms.png', 
                dpi=300, bbox_inches='tight')
    plt.savefig(f'{plotDir}/forestPlotModel{modelNoMapping[6]}{target}{familyModel6.capitalize()}'
                f'{interactionImmun2YNSymptomsSuffix6}{symptomsVar6.capitalize()}InteractionRecoveredSymptoms.pdf', 
                dpi=300, bbox_inches='tight')

    _ax
    return


@app.cell
def _(
    az,
    iDataBmb6,
    modelNoMapping,
    tableDirIData,
    varNames6Tmp,
    writeIDataSummaryTable,
    writeIDataSummaryTableLatex,
):
    summaryDf6 = az.summary(iDataBmb6, var_names='~agrdt_mean')
    writeIDataSummaryTable(summaryDf6, varnames=varNames6Tmp, 
                           outfile=tableDirIData / f'iDataSummaryModel{modelNoMapping[6]}.tsv')
    writeIDataSummaryTableLatex(summaryDf6, varnames=varNames6Tmp, 
                                outfile=tableDirIData / f'iDataSummaryModel{modelNoMapping[6]}Latex.txt', 
                                model_no=6)
    return


@app.cell
def _(bmbModel6, defaultdict, dfLogisticBmb6, iDataBmb6, predictionsNewData):
    predDictBmb6Symp = defaultdict(dict)
    predDictBmb6Asymp = defaultdict(dict)
    predictionsNewData(dfLogisticBmb6, bmbModel6, iDataBmb6, predDictBmb6Asymp, "immun2YN", (0, 1), 
                       otherFeaturesDict={"symptoms":0, "recovered": 0})
    predictionsNewData(dfLogisticBmb6, bmbModel6, iDataBmb6, predDictBmb6Symp, "immun2YN", (0, 1), 
                       otherFeaturesDict={"symptoms":1, "recovered": 0})
    return predDictBmb6Asymp, predDictBmb6Symp


@app.cell
def _(
    ANNOTATION_LETTER_SIZE,
    annotateWithLetters,
    colWidth,
    dfLogisticBmb6,
    modelNoMapping,
    plotDir,
    plt,
    predDictBmb6Asymp,
    predDictBmb6Symp,
    saveFigure,
    spaghettiPlotCategorical,
):
    figModel6Spagh, _ax = plt.subplots(2, 1, figsize=(colWidth * 2, colWidth * 2), sharex=True)
    _ax1, _ax2= _ax.flatten()
    spaghettiPlotCategorical(_ax2, 'immun2YN', (0, 1), predDictBmb6Asymp, dfLogisticBmb6[dfLogisticBmb6.symptoms==0], 
                             showXLabel=True, showLegend=False)
    spaghettiPlotCategorical(_ax1, 'immun2YN', (0, 1), predDictBmb6Symp, dfLogisticBmb6[dfLogisticBmb6.symptoms==1], 
                            showXLabel=False)
    plt.tight_layout()

    annotateWithLetters((_ax1, _ax2), size=ANNOTATION_LETTER_SIZE, coords=(-0.15, 1))
    saveFigure(figModel6Spagh, plotDir / f'model{modelNoMapping[6]}Spagh.png')
    saveFigure(figModel6Spagh, plotDir / f'model{modelNoMapping[6]}Spagh.pdf')
    figModel6Spagh
    return


@app.cell
def _(mo):
    mo.md(r"""## Summary forest plots""")
    return


@app.cell
def _(
    ANNOTATION_LETTER_SIZE,
    DINA4_HEIGHT,
    annotateWithLetters,
    colWidth,
    iDataBmb1,
    iDataBmb2,
    iDataBmb3,
    plotDir,
    plt,
    ridgeForestPlot,
    saveFigure,
    varNames1Tmp,
    varNames2Tmp,
    varNames3Tmp,
):
    _fig, _ax = plt.subplots(3, 1, figsize=(colWidth * 2, DINA4_HEIGHT), sharex=True) #height = colWidth * 0.9 * 3
    _axes = _ax.flatten()

    _varNamesList = (varNames1Tmp, varNames2Tmp, varNames3Tmp)
    _iDataList = (iDataBmb1, iDataBmb2, iDataBmb3)
    _counter = 0
    _highlightColor = "darkred"
    for _varNames, _iData, _ax in zip(_varNamesList, _iDataList, _axes):
        _ax1 = ridgeForestPlot(_iData, varNames=_varNames, xlim=(-3, 3), ax=_ax, fontsize=10, addGrid=True, ridge=False)
        _ax1.set_title(None)
        _counter += 1
        if _counter < len(_varNamesList):
            _ax1.set_xlabel("")         # removes x-axis label text
            _ax1.set_xticks([])         # removes tick marks
            _ax1.set_xticklabels([])    # removes tick labels
    annotateWithLetters(_axes, coords=(-0.07, 1), size=ANNOTATION_LETTER_SIZE)
    saveFigure(_fig, plotDir / "forestPlotModelSummaryA.png")
    saveFigure(_fig, plotDir / "forestPlotModelSummaryA.pdf")
    _fig
    return


@app.cell
def _(
    ANNOTATION_LETTER_SIZE,
    DINA4_HEIGHT,
    annotateWithLetters,
    colWidth,
    iDataBmb4,
    iDataBmb5,
    iDataBmb6,
    plotDir,
    plt,
    ridgeForestPlot,
    saveFigure,
    varNames4Tmp,
    varNames5Tmp,
    varNames6Tmp,
):
    _fig, _ax = plt.subplots(3, 1, figsize=(colWidth * 2, DINA4_HEIGHT), sharex=True)
    _axes = _ax.flatten()

    _highlightColor = "darkred"
    _varNamesList = (varNames4Tmp, varNames5Tmp, varNames6Tmp)
    _iDataList = (iDataBmb4, iDataBmb5, iDataBmb6)
    _counter = 0
    for _varNames, _iData, _ax in zip(_varNamesList, _iDataList, _axes):
        _ax1 = ridgeForestPlot(_iData, varNames=_varNames, xlim=(-3, 3), ax=_ax, fontsize=9, addGrid=True, ridge=False)
        _ax1.set_title(None)
        _counter += 1
        if _counter < len(_varNamesList):
            _ax1.set_xlabel("")         # removes x-axis label text
            _ax1.set_xticks([])         # removes tick marks
            _ax1.set_xticklabels([])    # removes tick labels
    annotateWithLetters(_axes, coords=(-0.1, 1), size=ANNOTATION_LETTER_SIZE)
    saveFigure(_fig, plotDir / "forestPlotModelSummaryB.png")
    saveFigure(_fig, plotDir / "forestPlotModelSummaryB.pdf")
    _fig
    return


if __name__ == "__main__":
    app.run()
