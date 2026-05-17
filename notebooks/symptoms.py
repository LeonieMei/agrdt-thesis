import marimo

__generated_with = "0.16.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    import marimo as mo
    import bambi as bmb
    import arviz as az

    from scipy.stats import skew
    from utils.dataUtils import (dataFrameIndependentSymptomsSymptomsData, dataFramePCRpos, addJitterCol,
                                 IQRQuartiles, dataFrameSymptomsData)
    from utils.tableUtils import writeSummaryTablesSymptoms, writeSummaryTableSymptomSeverity
    from utils.plotUtils import (setCustomTheme, plotSymptomHeatmap, annotateWithLetter,
                                 returnPlotDirSymptoms, replaceLegend, saveFigure, boxNSwarmplot, setFontSize)
    from utils.plotParams import (getAbbrvsDict, getPalettes, getLegends, getLabels,
                                  getOrders, COL_WIDTH, ANNOTATION_LETTER_SIZE)
    from utils.dataParams import SYMPTOMS
    from utils.regression import SEED
    from utils.dataParams import ROOT_DIR

    setCustomTheme()
    return (
        ANNOTATION_LETTER_SIZE,
        COL_WIDTH,
        IQRQuartiles,
        ROOT_DIR,
        SEED,
        SYMPTOMS,
        addJitterCol,
        annotateWithLetter,
        az,
        bmb,
        boxNSwarmplot,
        dataFrameIndependentSymptomsSymptomsData,
        dataFramePCRpos,
        dataFrameSymptomsData,
        getAbbrvsDict,
        getLabels,
        getLegends,
        getOrders,
        getPalettes,
        mo,
        np,
        pd,
        plotSymptomHeatmap,
        plt,
        replaceLegend,
        returnPlotDirSymptoms,
        saveFigure,
        setFontSize,
        skew,
        sns,
        writeSummaryTableSymptomSeverity,
        writeSummaryTablesSymptoms,
    )


@app.cell
def _(mo):
    mo.md("""# Load""")
    return


@app.cell
def _(mo):
    mo.md(r"""## Plotting parameters""")
    return


@app.cell
def _(
    getAbbrvsDict,
    getLabels,
    getLegends,
    getOrders,
    getPalettes,
    returnPlotDirSymptoms,
):
    plotDir = returnPlotDirSymptoms()
    abbrvs = getAbbrvsDict()
    pal = getPalettes()
    legend = getLegends()
    label = getLabels()
    variants = list(abbrvs["variant"].keys())
    order = getOrders()
    return abbrvs, label, legend, order, pal, plotDir, variants


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
def _(
    SYMPTOMS,
    dataFrameIndependentSymptomsSymptomsData,
    dataFramePCRpos,
    dataFrameSymptomsData,
    df,
    pd,
):
    dfSympDatSymp = dataFrameSymptomsData(df)
    dfSympDatSympPos = dataFramePCRpos(dfSympDatSymp)
    dfSympDatSympInd = dataFrameIndependentSymptomsSymptomsData(df)
    dfSympDatSympIndPos = dataFramePCRpos(dfSympDatSympInd)
    dfSympDatSympIndPos["ageBin"] = pd.cut(dfSympDatSympIndPos.age, bins=[17, 30, 50, 70], include_lowest=True)
    dfSympDatSympIndPos["nSymptoms"] = dfSympDatSympIndPos[[symptom for symptom in SYMPTOMS if symptom != "ill"]].sum(axis=1)
    return dfSympDatSympInd, dfSympDatSympIndPos


@app.cell
def _(dataFrameSymptomsData, df):
    print(f"Total number of data points with symptom data available: {len(dataFrameSymptomsData(df))}")
    return


@app.cell
def _(mo):
    mo.md(r"""# Tables""")
    return


@app.cell
def _(
    df,
    dfSympDatSympIndPos,
    writeSummaryTableSymptomSeverity,
    writeSummaryTablesSymptoms,
):
    writeSummaryTablesSymptoms(df)
    writeSummaryTableSymptomSeverity(dfSympDatSympIndPos)
    return


@app.cell
def _(mo):
    mo.md(r"""# Plots and stats""")
    return


@app.cell
def _(mo):
    mo.md(r"""## By variant""")
    return


@app.cell
def _(mo):
    mo.md(r"""## Correlations between symptoms""")
    return


@app.cell
def _(
    ANNOTATION_LETTER_SIZE,
    abbrvs,
    annotateWithLetter,
    dfSympDatSympIndPos,
    plotDir,
    plotSymptomHeatmap,
    plt,
    saveFigure,
    sns,
    variants,
):
    # Correlation between different symptoms (heatmap)
    # distance measure: jaccard similarity
    with sns.axes_style("white"):
        fig1, axesOrig1 = plt.subplots(2, 1, figsize=(7, 12))
        fig2, axesOrig2 = plt.subplots(2, 1, figsize=(7, 12))


    axes = list(axesOrig1.flatten()) + list(axesOrig2.flatten())
    lastIdx = len(variants) - 1
    cbar = True
    letters = ("A", "B", "C", "D")

    jaccSimVariants = {}

    for i, (variant, ax) in enumerate(zip(variants, axes)):
        dfCurr = dfSympDatSympIndPos[dfSympDatSympIndPos.variant==variant].copy()
        jaccSim = plotSymptomHeatmap(dfCurr, ax, cbar=cbar, title=abbrvs["variant"][variant])
        jaccSimVariants[variant] = jaccSim
        annotateWithLetter(ax, letters[i], coords=(-0.4, 1.07), size=ANNOTATION_LETTER_SIZE)
        if i in (1, 3):
            suffix = "WtAlpha" if i == 1 else "DeltaOmicron"
            fig = fig1 if i == 1 else fig2
            fig.tight_layout();
            saveFigure(fig, plotDir / f"symptomHeatmaps{suffix}.png")
    return (jaccSimVariants,)


@app.cell
def _(mo):
    mo.md(r"""### Stats""")
    return


@app.cell
def _(abbrvs, dfSympDatSympIndPos, variants):
    for _variant in variants:
        _df = dfSympDatSympIndPos[dfSympDatSympIndPos.variant==_variant].copy()
        print(abbrvs["variant"][_variant])
        print(f"Frequency of loss of smell: {_df.noSmell.mean():.2f}")

    print()
    for _variant in variants:
        _df = dfSympDatSympIndPos[dfSympDatSympIndPos.variant==_variant].copy()
        print(abbrvs["variant"][_variant])
        print(f"Frequency of fever: {_df.fever.mean():.2f}")

    print()
    for _variant in variants:
        _df = dfSympDatSympIndPos[dfSympDatSympIndPos.variant==_variant].copy()
        print(abbrvs["variant"][_variant])
        print(f"Frequency of dyspnea: {_df.dyspnea.mean():.2f}")
    return


@app.cell
def _(SYMPTOMS, dfSympDatSympIndPos, jaccSimVariants, variants):
    symptomsCurr = [symptom for symptom in SYMPTOMS if symptom != "ill"]
    for variantCurr in variants:
        corrVariant = jaccSimVariants[variantCurr]
        dfCurrVariant = dfSympDatSympIndPos[dfSympDatSympIndPos.variant==variantCurr].copy()
        mostCommonSymptoms = []
        mostCommonSymptom = ""
        symptomPercentageMax = 0

        print(variantCurr)
        print(f"Number of data points: {len(dfCurrVariant)}")
        for symptom in symptomsCurr:
            symptomPercentage = dfCurrVariant[symptom].mean()
            if symptomPercentage >= symptomPercentageMax:
                mostCommonSymptom = symptom
                symptomPercentageMax = symptomPercentage
        for symptomNew in symptomsCurr:
            symptomPercentageNew = dfCurrVariant[symptomNew].mean()
            if abs(symptomPercentageNew - symptomPercentageMax) < 0.01:
                mostCommonSymptoms.append(symptomNew)
        print(f"Most common symptoms:")
        for mostCommonSymptom in mostCommonSymptoms:
            print(f"{mostCommonSymptom} ({dfCurrVariant[mostCommonSymptom].mean()})")
        print()
        seen = set()
        for symptom1 in corrVariant.index:
            for symptom2 in corrVariant.columns:
                if symptom1 == symptom2: continue
                symptomTupleSorted = tuple(sorted((symptom1, symptom2)))
                if (corrVariant.loc[symptom1, symptom2] > 0.5) and not symptomTupleSorted in seen:
                    print(f"{symptomTupleSorted[0]}, {symptomTupleSorted[1]}: {corrVariant.loc[symptomTupleSorted[0], 
                    symptomTupleSorted[1]]}")
                    seen.add(symptomTupleSorted)
        print()
    return


@app.cell
def _(mo):
    mo.md(r"""#### Numbers of reported symptoms""")
    return


@app.cell
def _(IQRQuartiles, dfSympDatSympIndPos):
    # Number of reported symptoms
    print("Overall median")
    print(dfSympDatSympIndPos.nSymptoms.median())
    print("Overall IQR")
    print(IQRQuartiles(dfSympDatSympIndPos.nSymptoms))
    print()
    print("median")
    print(dfSympDatSympIndPos.groupby("variant", observed=True).nSymptoms.agg("median"))
    print("IQR")
    print(dfSympDatSympIndPos.groupby("variant", observed=True).nSymptoms.agg(IQRQuartiles))
    return


@app.cell
def _(dfSympDatSympIndPos, sns):
    sns.boxplot(dfSympDatSympIndPos, x="variant", y="nSymptoms")
    return


@app.cell
def _(mo):
    mo.md(r"""### Symptom severity""")
    return


@app.cell
def _(dfSympDatSympIndPos, sns):
    sns.boxplot(data=dfSympDatSympIndPos, x="variant", y="illSeverity")
    return


@app.cell
def _(IQRQuartiles, dfSympDatSympIndPos, rint):
    rint("Overall median")
    print(dfSympDatSympIndPos.illSeverity.median())
    print("Overall IQR")
    print(IQRQuartiles(dfSympDatSympIndPos.illSeverity))
    print()
    print("median")
    print(dfSympDatSympIndPos.groupby("variant").illSeverity.agg("median"))
    print("IQR")
    print(dfSympDatSympIndPos.groupby("variant").illSeverity.agg(IQRQuartiles))
    return


@app.cell
def _(mo):
    mo.md(r"""## By age""")
    return


@app.cell
def _(mo):
    mo.md(r"""### Correlations between symptoms""")
    return


@app.cell
def _(
    ANNOTATION_LETTER_SIZE,
    annotateWithLetter,
    dfSympDatSympIndPos,
    plotDir,
    plotSymptomHeatmap,
    plt,
    saveFigure,
    sns,
):
    # Correlation between different symptoms (heatmap) for different ages
    # distance measure: jaccard similarity
    with sns.axes_style("white"):
        fig1Age, axesOrig1Age = plt.subplots(2, 1, figsize=(7, 12))
        fig2Age, axesOrig2Age = plt.subplots(1, 1, figsize=(7, 6))


    axesAge = list(axesOrig1Age.flatten())
    axesAge.append(axesOrig2Age)
    ageBins = dfSympDatSympIndPos.ageBin.cat.categories
    cbarAge = True
    lettersAge = ("A", "B", "C")

    for iAge, (ageBin, axAge) in enumerate(zip(ageBins, axesAge)):
        dfCurrAge = dfSympDatSympIndPos[dfSympDatSympIndPos.ageBin == ageBin].copy()
        plotSymptomHeatmap(dfCurrAge, axAge, cbar=cbarAge, title=ageBin)
        annotateWithLetter(axAge, lettersAge[iAge], coords=(-0.4, 1.07), size=ANNOTATION_LETTER_SIZE)
        if iAge in (1, 2):
            if iAge == 1:
                suffixAge = ageBins[:2]
                figAge = fig1Age
            else:
                suffixAge = ageBin
                figAge = fig2Age
            figAge.tight_layout();
            saveFigure(figAge, plotDir / f"symptomHeatmaps{suffixAge}.png")
    return


@app.cell
def _(mo):
    mo.md("""### Stats""")
    return


@app.cell
def _(dfSympDatSympIndPos, sns):
    sns.boxplot(dfSympDatSympIndPos, x="ageBin", y="nSymptoms")
    return


@app.cell
def _(IQRQuartiles, dfSympDatSympIndPos):
    print("median")
    print(dfSympDatSympIndPos.groupby("ageBin", observed=True).nSymptoms.agg("median"))
    print("IQR")
    print(dfSympDatSympIndPos.groupby("ageBin", observed=True).nSymptoms.agg(IQRQuartiles))
    return


@app.cell
def _(mo):
    mo.md("""### Symptom severity""")
    return


@app.cell
def _(IQRQuartiles, dfSympDatSympIndPos):
    print("median")
    print(dfSympDatSympIndPos.groupby("ageBin", observed=True).illSeverity.agg("median"))
    print("IQR")
    print(dfSympDatSympIndPos.groupby("ageBin", observed=True).illSeverity.agg(IQRQuartiles))
    return


@app.cell
def _(mo):
    mo.md("""## By gender""")
    return


@app.cell
def _(mo):
    mo.md(r"""### Stats""")
    return


@app.cell
def _(dfSympDatSympIndPos, sns):
    sns.boxplot(dfSympDatSympIndPos, x="gender", y="nSymptoms")
    return


@app.cell
def _(IQRQuartiles, dfSympDatSympIndPos):
    print("median")
    print(dfSympDatSympIndPos.groupby("gender", observed=True).nSymptoms.agg("median"))
    print("IQR")
    print(dfSympDatSympIndPos.groupby("gender", observed=True).nSymptoms.agg(IQRQuartiles))
    return


@app.cell
def _(mo):
    mo.md(r"""### Symptom severity""")
    return


@app.cell
def _(dfSympDatSympIndPos, sns):
    sns.boxplot(dfSympDatSympIndPos, x="gender", y="illSeverity")
    return


@app.cell
def _(IQRQuartiles, dfSympDatSympIndPos):
    print("median")
    print(dfSympDatSympIndPos.groupby("gender", observed=True).illSeverity.agg("median"))
    print("IQR")
    print(dfSympDatSympIndPos.groupby("gender", observed=True).illSeverity.agg(IQRQuartiles))
    return


@app.cell
def _(mo):
    mo.md(r"""## Loss of taste and smell""")
    return


@app.cell
def _(SEED, bmb, dfSympDatSympIndPos, order, pd):
    # Logistic regression to compare percentages between different variants

    dfNoTasteOrSmell = dfSympDatSympIndPos.copy()
    dfNoTasteOrSmell = dfNoTasteOrSmell.dropna(subset=("noTaste", "noSmell", "variant"))
    dfNoTasteOrSmell["noTasteOrSmell"] = (dfNoTasteOrSmell.noTaste == 1) | (dfNoTasteOrSmell.noSmell == 1)
    dfNoTasteOrSmell['variantCode'] = pd.Categorical(dfNoTasteOrSmell.variant, categories=order.variant).codes
    modelNoTasteOrSmell = bmb.Model(formula = "noTasteOrSmell ~ variantCode", family="bernoulli", data=dfNoTasteOrSmell, 
                                    categorical="variantCode", dropna=True)
    iDataNoTasteOrSmell = modelNoTasteOrSmell.fit(draws=1000, tune=1000, seed=SEED)
    return iDataNoTasteOrSmell, modelNoTasteOrSmell


@app.cell
def _(modelNoTasteOrSmell):
    modelNoTasteOrSmell
    return


@app.cell
def _(iDataNoTasteOrSmell, modelNoTasteOrSmell, np, pd):
    newDataNoTasteOrSmell = pd.DataFrame({"variantCode": np.array([0, 1, 2, 3])})
    # We have to make posterior predictions to get probabilities.
    modelNoTasteOrSmell.predict(iDataNoTasteOrSmell, data=newDataNoTasteOrSmell)
    return (newDataNoTasteOrSmell,)


@app.cell
def _(az, iDataNoTasteOrSmell):
    az.summary(iDataNoTasteOrSmell)
    return


@app.cell
def _(az, iDataNoTasteOrSmell, newDataNoTasteOrSmell, np):
    samplesNoTasteOrSmellVariant = {}
    meanDiffNoTasteOrSmellVariant = {}
    hdiDiffNoTasteOrSmellVariant = {}
    for variantCode in range(4):
        _idx = newDataNoTasteOrSmell.index[newDataNoTasteOrSmell.variantCode == variantCode].tolist()
        samplesNoTasteOrSmellVariant[variantCode] = iDataNoTasteOrSmell.posterior.p.stack(sample=("chain", "draw"))[_idx].values.flatten()

    for variantCode in range(1, 4):
        _diff = samplesNoTasteOrSmellVariant[variantCode] - samplesNoTasteOrSmellVariant[0]
        meanDiffNoTasteOrSmellVariant[variantCode] = np.mean(_diff)
        hdiDiffNoTasteOrSmellVariant[variantCode] = az.hdi(_diff, 0.94)
    return (hdiDiffNoTasteOrSmellVariant,)


@app.cell
def _(hdiDiffNoTasteOrSmellVariant):
    hdiDiffNoTasteOrSmellVariant
    return


@app.cell
def _(mo):
    mo.md("""## Viral load""")
    return


@app.cell
def _(mo):
    mo.md(r"""### Splines model""")
    return


@app.cell
def _(plotVlDaysSpDf):
    plotVlDaysSpDf.value_counts("variant")
    return


@app.cell
def _(bmb, dfSympDatSympInd, np, pd):
    plotVlDaysSpDf = dfSympDatSympInd.dropna(subset="daysPostOnset").copy()
    plotVlDaysSpDf = plotVlDaysSpDf[plotVlDaysSpDf.infectionKey.notna()]
    plotVlDaysSpDf.vl = plotVlDaysSpDf.vl.fillna(3)
    plotVlDaysSpDf["vlCens"] = np.where(plotVlDaysSpDf.vl <= 3, "left", "none")
    plotVlDaysSpDf = plotVlDaysSpDf.dropna(subset=["variant", "daysPostOnset", "vl"])
    plotVlDaysSpDf["variant"] = (pd.Categorical(plotVlDaysSpDf.variant, categories=["wildtype", "alpha", "delta", 
                                                                                    "omicron"]).codes)
    degreesFreedom = 5
    # Default degree (of the piecewise polynomial - not the same as degrees of freedom) is 3 for cubic splies
    # from https://github.com/bambinos/formulae/blob/5c28351e5c429e367008a43a1ad7042509e6c5e6/formulae/transforms.py#L228-L361
    # order = degree + 1
    # if df is not None:
    #     n_inner_knots = df - order
    # --> we define the inner knot instead (setting it to 1.5); the placement of the knot is therefore fixed by us and does not say anything about the timing of symptom onset relative to peak viral load between the different variants; but if symptom onset occurred after peak viral load in one of the variants, we would not see a bump.

    iknots=np.array([1.5])

    priors = {
        f"bs(daysPostOnset, df={degreesFreedom}, intercept=True, knots=iknots)": bmb.Prior("Normal", mu=0, sigma=2),
        "Intercept": bmb.Prior("Normal", mu=0, sigma=2)
    }
    formula = bmb.Formula(
        f"censored(vl, vlCens) ~ 0 + bs(daysPostOnset, df={degreesFreedom}, intercept=True, knots=iknots)",
    )

    modelsVsDaysSp = {}
    iDataVlDaysSp = {}
    for _variant in range(4):
        _data = plotVlDaysSpDf[plotVlDaysSpDf.variant == _variant]
        modelsVsDaysSp[_variant] = modelVlDaysSp = bmb.Model(formula=formula, data=_data, dropna=True, family="t", 
                                                             priors=priors)
    return iDataVlDaysSp, modelsVsDaysSp, plotVlDaysSpDf


@app.cell
def _(modelsVsDaysSp):
    modelsVsDaysSp[0]
    return


@app.cell
def _(
    COL_WIDTH,
    abbrvs,
    iDataVlDaysSp,
    label,
    legend,
    modelsVsDaysSp,
    np,
    pal,
    pd,
    plotDir,
    plt,
    replaceLegend,
    saveFigure,
    setFontSize,
):
    q = [0.03, 0.97]
    dims = ("chain", "draw")

    means = {}
    meanIntervals = {}
    yIntervals = {}

    days = np.linspace(0, 10, num=100)

    for _variant in (0, 1, 2, 3):
        newData = pd.DataFrame({"daysPostOnset": days, "variant": np.repeat(_variant, 100)})
        modelsVsDaysSp[_variant].predict(iDataVlDaysSp[_variant], data=newData, kind="response")

        means[_variant] = iDataVlDaysSp[_variant].posterior["mu"].mean(dims).to_numpy()
        meanIntervals[_variant] = iDataVlDaysSp[_variant].posterior["mu"].quantile(q, dims).to_numpy()


    _fig, _ax = plt.subplots(figsize=(COL_WIDTH * 2, COL_WIDTH * 1.3), constrained_layout=True)
    for _variant in (0, 2, 3):
        _color = pal.variant[list(abbrvs["variant"])[_variant]]
        _ax.plot(days, means[_variant], color=_color)
        _ax.fill_between(days, meanIntervals[_variant][0], meanIntervals[_variant][1], alpha=0.5, color=_color)


    _ax.set_ylim(2, 11)
    _ax.set_xlim(0, 10)
    _ax.set_xlabel(label.days)
    _ax.set_ylabel(label.vl)
    replaceLegend(_ax, legend.variantLines[:2] + legend.variantLines[3:], loc=(0.55, 0.7))
    setFontSize(_ax)

    saveFigure(_fig, plotDir / "vlDaysPostOnset.png")
    _fig
    return


@app.function
def mapDays(day):
    if day <= 1:
        return 1
    elif day <= 2:
        return 2
    else:
        return 3


@app.cell
def _(plotVlDaysSpDf):
    plotVlDaysSpDf["daysPostOnsetCoarse"] = plotVlDaysSpDf.daysPostOnset.map(mapDays)
    return


@app.cell
def _(plotVlDaysSpDf, plt, sns):
    _fig, _ax = plt.subplots(1, 1, figsize=(20, 5))
    sns.swarmplot(plotVlDaysSpDf[plotVlDaysSpDf.vl > 3], x="daysPostOnsetCoarse", y="vl", hue="variant", dodge=True, ax=_ax, size=4)
    return


@app.cell
def _(pd, plotVlDaysSpDf, skew):
    variantSkewness = (pd.DataFrame(plotVlDaysSpDf[plotVlDaysSpDf.vl > 3].groupby(["variant", "daysPostOnsetCoarse"]).vl.apply(skew)).reset_index())
    return (variantSkewness,)


@app.cell
def _(variantSkewness):
    variantSkewness
    return


@app.cell
def _(mo):
    mo.md(r"""### Linear regression""")
    return


@app.cell
def _(SEED, bmb, dfSympDatSympInd, np, pd):
    plotVlDaysDf = dfSympDatSympInd.dropna(subset="daysPostOnset").copy()
    plotVlDaysDf = plotVlDaysDf[plotVlDaysDf.infectionKey.notna()]
    plotVlDaysDf.vl = plotVlDaysDf.vl.fillna(3)
    plotVlDaysDf["vlCens"] = np.where(plotVlDaysDf.vl <= 3, "left", "none")

    _priors = {"Intercept": bmb.Prior("Normal", mu=7, sigma=1),
               "daysPostOnset": bmb.Prior("Normal", mu=0, sigma=1.5),
               "sigma": bmb.Prior("HalfStudentT", nu=4, sigma=3)}
    modelVlDays = bmb.Model(data=plotVlDaysDf, formula="censored(vl, vlCens) ~ daysPostOnset", family="t", priors=_priors)
    iDataVlDays = modelVlDays.fit(tune=1000, draws=2000, random_seed=SEED, target_accept=0.9)

    daysRange = np.linspace(min(plotVlDaysDf.daysPostOnset), max(plotVlDaysDf.daysPostOnset), 200)
    predictDf = pd.DataFrame({"daysPostOnset": daysRange})

    modelVlDays.predict(iDataVlDays, data=predictDf, kind="response")
    return daysRange, iDataVlDays, modelVlDays, plotVlDaysDf


@app.cell
def _(modelVlDays):
    modelVlDays
    return


@app.cell
def _(az, iDataVlDays):
    az.summary(iDataVlDays, var_names="~mu")
    return


@app.cell
def _(COL_WIDTH, addJitterCol, plotVlDaysDf, plt, sns):
    _, _axVl = plt.subplots(1, 1, figsize=(COL_WIDTH * 2, COL_WIDTH * 1.2))

    addJitterCol(plotVlDaysDf, "daysPostOnset", "daysPostOnsetJitter", sd=0.0075)
    sns.swarmplot(data=plotVlDaysDf[plotVlDaysDf.variant.isin(("wildtype", "alpha"))], x="daysPostOnset", y="vl", 
                    alpha=1, edgecolor="none", size=2.5, ax=_axVl)
    return


@app.cell
def _(COL_WIDTH, addJitterCol, pal, plotVlDaysDf, plt, sns):
    _, _axVl = plt.subplots(1, 1, figsize=(COL_WIDTH * 2, COL_WIDTH * 1.2))

    addJitterCol(plotVlDaysDf, "daysPostOnset", "daysPostOnsetJitter", sd=0.0075)
    sns.swarmplot(data=plotVlDaysDf[plotVlDaysDf.variant.isin(("delta", "omicron"))], x="daysPostOnset", y="vl", alpha=1,
                  edgecolor="none", hue="variant", size=2.5, palette=pal.variant, ax=_axVl)
    return


@app.cell
def _(
    COL_WIDTH,
    addJitterCol,
    az,
    daysRange,
    iDataVlDays,
    label,
    plotDir,
    plotVlDaysDf,
    plt,
    saveFigure,
    setFontSize,
    sns,
):
    figVl, axVl = plt.subplots(1, 1, figsize=(COL_WIDTH * 2, COL_WIDTH * 1.2), constrained_layout=True)

    addJitterCol(plotVlDaysDf, "daysPostOnset", "daysPostOnsetJitter", sd=0.0075)
    sns.scatterplot(data=plotVlDaysDf, x="daysPostOnsetJitter", y="vl", alpha=0.5,
                    edgecolor="none", ax=axVl)

    iDataVlDaysStacked = az.extract(iDataVlDays)
    # Plot recovered robust linear regression
    axVl.plot(daysRange, iDataVlDaysStacked.mu.mean(axis=1), linestyle="-")
    # Plot HDIs
    for interval in [0.94]:
        az.plot_hdi(daysRange, iDataVlDaysStacked.mu.T, 
                    hdi_prob=interval, color="red", fill_kwargs={"alpha": 0.1},
                    ax=axVl)
    # Customize the plot to set alpha (transparency)
    for artist in axVl.get_children():
        if isinstance(artist, plt.Polygon):
            artist.set_alpha(0.2)  # Set the desired alpha value


    axVl.set_xlabel(label.daysPostOnset, labelpad=10)
    axVl.set_ylabel(label.vl, labelpad=2)
    setFontSize(axVl)
    saveFigure(figVl, plotDir / "vlDaysPostOnset.png")
    axVl
    return


@app.cell
def _(mo):
    mo.md(r"""### Swarmplots""")
    return


@app.cell
def _(COL_WIDTH, pal, plotDf, plt, sns):
    _, axVlImmun = plt.subplots(2, 1, figsize=(COL_WIDTH * 2, 6), sharex=True, sharey=True)
    for axImmun, immunVl in zip(axVlImmun.flatten(), (0, 1), strict=True):
        sns.scatterplot(data=plotDf[(plotDf.immun2YN==immunVl)], x="daysPostOnset", 
                        y="vl", alpha=0.5, ax=axImmun, color=pal.immun2YN[immunVl])
    axVlImmun
    return


@app.cell
def _(COL_WIDTH, boxNSwarmplot, dfSympDatSympIndPos, label, plt):
    figVlSeverity, axVlSeverity = plt.subplots(1, 1, figsize=(COL_WIDTH * 2, 4))

    boxNSwarmplot(dfSympDatSympIndPos, x="illSeverity", y="vl", xLabel="Severity of illness", yLabel=label.vl,
                  xCats=sorted(dfSympDatSympIndPos.illSeverity.dropna().unique()), palBoxplot="none", markercolor="#F8766D",
                  markersize=5, ax=axVlSeverity)
    axVlSeverity
    return


@app.cell
def _(COL_WIDTH, abbrvs, boxNSwarmplot, dfSympDatSympIndPos, label, pal, plt):
    figVlEarlyVariant, axVlEarlyVariant = plt.subplots(1, 1, figsize=(COL_WIDTH * 2, 4))

    boxNSwarmplot(dfSympDatSympIndPos[(dfSympDatSympIndPos.daysPostOnset < 2)], x="variant", y="vl", xLabel=label.variant,
                  yLabel=label.vl, xCats=abbrvs["variant"], palBoxplot="none", palSwarmplot=pal.variant,
                  markersize=5, ax=axVlEarlyVariant)
    axVlEarlyVariant
    return


@app.cell
def _(COL_WIDTH, abbrvs, boxNSwarmplot, dfSympDatSympIndPos, label, pal, plt):
    figVlEarlyImmun, axVlEarlyImmun = plt.subplots(1, 1, figsize=(COL_WIDTH * 2, 4))

    boxNSwarmplot(dfSympDatSympIndPos[dfSympDatSympIndPos.daysPostOnset < 1], x="immun2YN", y="vl", xLabel=label.immun2YN,
                  yLabel=label.vl, xCats=abbrvs["immun2YN"], palBoxplot="none", palSwarmplot=pal.immun2YN,
                  markersize=5, ax=axVlEarlyImmun)
    axVlEarlyImmun
    return


if __name__ == "__main__":
    app.run()
