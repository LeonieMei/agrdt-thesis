import marimo

__generated_with = "0.16.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    from datetime import date
    from collections import defaultdict

    # Bayesian MCMC analysis
    import arviz as az

    from agrdt.data import (createDataFramesFigures, returnWtAlphaTransitionPCRs,
                            returnDeltaOmicronTransitionPCRs)
    from agrdt.tables import writeSummaryTablesEmployees, returnTableDirIData
    from agrdt.plotting import replaceLegend, returnPlotDirRegression, setCustomTheme
    from agrdt.plotParams import (getPalettes, getLabels, getLegends, getOrders,
                                  getAbbrvsDict, COL_WIDTH, CM)
    from agrdt.regression import sampleAgrdtSpec, returnIDataDirRegression
    from agrdt.dataParams import ROOT_DIR

    setCustomTheme()
    pd.set_option('display.max_columns', 30)
    pd.set_option('display.float_format', '{:.2f}'.format)
    return (
        ROOT_DIR,
        az,
        createDataFramesFigures,
        date,
        defaultdict,
        getAbbrvsDict,
        getLabels,
        getLegends,
        getOrders,
        getPalettes,
        mo,
        np,
        pd,
        plt,
        replaceLegend,
        returnDeltaOmicronTransitionPCRs,
        returnWtAlphaTransitionPCRs,
        sampleAgrdtSpec,
        sns,
        writeSummaryTablesEmployees,
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
def _(getAbbrvsDict, getLabels, getLegends, getOrders, getPalettes):
    pal = getPalettes()
    legend = getLegends()
    order = getOrders()
    label = getLabels()
    abbrvDictPaper = getAbbrvsDict()
    return abbrvDictPaper, label, legend, order, pal


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Data""")
    return


@app.cell
def _(ROOT_DIR, pd):
    # Load data
    target="T1"
    dataPath = ROOT_DIR / "data" / f"agrdtDataThesis{target}.tsv"
    df = pd.read_csv(dataPath, sep="\t", low_memory=False, parse_dates=['pcrDate'])
    # Turn into datetime.date objects.
    df["pcrDate"] = df["pcrDate"].apply(lambda pcrDate: pcrDate.date())
    return (df,)


@app.cell
def _(createDataFramesFigures, df):
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
    return (
        dfAgrdt,
        dfAllInd,
        dfFigure1,
        dfFigure1_A,
        dfFigure1_C,
        dfFigure2,
        dfFigure2_A,
        dfFigure2_B,
        dfFigure3,
        dfFigureA2,
        dfPos,
    )


@app.cell
def _(mo):
    mo.md(r"""# Counts""")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Figures""")
    return


@app.cell
def _(
    date,
    dfFigure1,
    dfFigure2,
    dfFigure2_A,
    dfFigure2_B,
    dfFigure3,
    dfFigureA2,
    dfPos,
    order,
):
    _months2 = dfPos[dfPos.pcrDate > date(2020, 11, 30)].samplingMonth2.sort_values().unique()
    countsPcr = [f'{((dfFigure2_B.samplingMonth2 == month2) & (dfFigure2_B.agrdtYN == 1)).sum()}{(dfFigure2_B.samplingMonth2 == month2).sum()}' for month2 in _months2]
    countsFigure1 = [f'n={(dfFigure1.samplingMonth2 == month2).sum()}' for month2 in _months2]
    countsFigure1VaccInfoAvailable = [f'n={((dfFigure1.samplingMonth2 == month2) & dfFigure1.immunN.notna()).sum()}' for month2 in _months2]
    countsFigure2_A = [f'n={(dfFigure2_A.variant == _variant).sum()}' for _variant in order.variant]
    countsFigure2_B = [f'No Agrdt: n={((dfFigure2.variant == _variant) & (dfFigure2.agrdtYN == 0)).sum()}, Agrdt: n={((dfFigure2.variant == _variant) & (dfFigure2.agrdtYN == 1)).sum()}' for _variant in order.variant]
    countsFigure3_Asymp = [f'n={((dfFigure3.symptoms == 0) & (dfFigure3.immun2YN == immun)).sum()}' for immun in order.immun2YN]
    countsFigure3_Symp = [f'n={((dfFigure3.symptoms == 1) & (dfFigure3.immun2YN == immun)).sum()}' for immun in order.immun2YN]
    countsFigureA2_Asymp = [f'n={((dfFigureA2.symptoms == 0) & (dfFigureA2.variant == _variant)).sum()}' for _variant in order.variant]
    countsFigureA2_Symp = [f'n={((dfFigureA2.symptoms == 1) & (dfFigureA2.variant == _variant)).sum()}' for _variant in order.variant]
    return


@app.cell
def _(dfFigure3, dfFigureA2, order):
    countsFigure3Verbose_Asymp = [f'({immun}, {res}): {((dfFigure3.symptoms == 0) & (dfFigure3.immun2YN == immun) & (dfFigure3.agrdt == res)).sum()}' for immun in order.immun2YN for res in (0, 1)]
    countsFigure3Verbose_Symp = [f'({immun}, {res}): {((dfFigure3.symptoms == 1) & (dfFigure3.immun2YN == immun) & (dfFigure3.agrdt == res)).sum()}' for immun in order.immun2YN for res in (0, 1)]
    countsFigureA2Verbose_Asymp = [f'({_variant}, {res}): {((dfFigureA2.symptoms == 0) & (dfFigureA2.variant == _variant) & (dfFigureA2.agrdt == res)).sum()}' for _variant in order.variant for res in (0, 1)]
    countsFigureA2Verbose_Symp = [f'({_variant}, {res}): {((dfFigureA2.symptoms == 1) & (dfFigureA2.variant == _variant) & (dfFigureA2.agrdt == res)).sum()}' for _variant in order.variant for res in (0, 1)]
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Overall""")
    return


@app.cell
def _(df, dfAgrdt, dfFigure1_A, dfFigure1_C, dfPos):
    print(f"Total number of PCR test performed at CUSMA test centres: {len(df)}")
    print(f"Total number of employees tested by PCR at CUSMA test centres: {df.personHash.nunique()}")
    print(f"Total number of Ag-RDTs performed: {len(dfAgrdt)}")
    print(f"Number of Abbott tests performed: {(dfAgrdt.testDevice == 'abbott').sum()}")
    print(f"Number of Roche tests performed: {(dfAgrdt.testDevice == 'roche').sum()}")
    print(f"Fraction of Roche tests performed: {(dfAgrdt.testDevice == 'roche').sum() / dfAgrdt.testDevice.notna().sum()}")
    print(f"Fraction of Abbott tests performed: {(dfAgrdt.testDevice == 'abbott').sum() / dfAgrdt.testDevice.notna().sum()}")
    print(f"Start date of Roche test usage: {dfAgrdt[dfAgrdt.testDevice == 'roche'].pcrDate.min()}")
    print(f"Total number of people tested by Ag-RDT: {dfAgrdt.personHash.nunique()}")
    print()

    print(f"Number of tests with questionnaire info: {(dfAgrdt.surveyData==1).sum()}")
    print(f"Number of tests for which time post symptom onset was available: {(dfAgrdt.daysPostOnset.notna().sum())}")
    print(f"Number of tests for which vaccination data is available: {dfAgrdt.vaccN.notna().sum()}")
    print(f"Number of PCR positive tests with questionnaire info: {(dfPos.surveyData==1).sum()}")
    print(f"Number of PCR positive tests for which vaccination data is available: {dfPos.vaccN.notna().sum()}")
    print()

    print(f"Number of tests where person was symptomatic: {(dfAgrdt.symptoms == 1).sum()}")
    print(f"Number of tests where person was asymptomatic: {(dfAgrdt.symptoms == 0).sum()}")
    print(f"Percentage of symptomatic people: {dfAgrdt.symptoms.mean():.2f}")
    print(f"Number of tested asymptomatic people: {dfAgrdt[dfAgrdt.symptoms == 0].personHash.nunique()}")
    print(f"Number of tested symptomatic people: {dfAgrdt[dfAgrdt.symptoms == 1].personHash.nunique()}")
    print()

    print(f"Number of positive PCRs: {(dfAgrdt.pcrPositive).sum()}")
    print(f"Number of positive PCRs (where symptom status is known): {(dfAgrdt.symptoms.notna() & dfAgrdt.pcrPositive).sum()}")
    print(f"Number of tests with known symptom status and PCR result: {(dfAgrdt.symptoms.notna() & dfAgrdt.pcrPositive.notna()).sum()}")
    print()

    print(f"Number of positive Ag-RDTs: {(dfAgrdt.agrdt).sum()}")
    print(f"Number of positive Ag-RDTs (where symptom status is known): {(dfAgrdt.symptoms.notna() & dfAgrdt.agrdt).sum()}")
    print(f"Number of tests with known symptom status and Ag-RDT result: {(dfAgrdt.symptoms.notna() & dfAgrdt.agrdt.notna()).sum()}")
    print()

    # Note: we are considering all tests, irrespective of known or unknown symptom status
    print(f"Fraction of typed PCRs (only tests with corresponding rapid test): {(dfAgrdt.pcrPositive & (dfAgrdt.hasTyping == 1) & dfAgrdt.variant.notna()).sum()/(dfAgrdt.pcrPositive & dfAgrdt.agrdt.notna()).sum():.3f}")
    print(f"Number of typed PCRs (only tests with corresponding rapid test): {(dfAgrdt.pcrPositive & (dfAgrdt.hasTyping == 1) & dfAgrdt.variant.notna()).sum()}")
    print(f"Fraction of typed PCRs (all tests): {df.hasTyping.mean():.2f}")
    infectionsWithTypingPCR = df[df.hasTyping == 1].infectionKey
    maskPcrNoTyping = dfAgrdt.pcrPositive & dfAgrdt.variant.notna() & ((dfAgrdt.hasTyping == 0) | dfAgrdt.hasTyping.isna())
    print(f"Number of PCRs where variant was assigned because typed PCR is in the same infection: {(maskPcrNoTyping & dfAgrdt.infectionKey.isin(infectionsWithTypingPCR)).sum()}")
    print(f"Number of PCRs were variant was assigned by date: {(maskPcrNoTyping & ~dfAgrdt.infectionKey.isin(infectionsWithTypingPCR)).sum():.2f}")
    print()

    print(f"Number of positive PCRs for which we know the person had a prior infection: {dfPos.recovered.sum():.2f}")
    print(f"Fraction of positive PCRs for which we know the person had a prior infection: {dfPos.recovered.mean():.2f}")
    print()

    dfRecovered = dfPos[dfPos.recovered==1]
    print(f"Number of positive PCRs for which the person had at least two infections (and no vaccinations): {((dfRecovered.vaccN==0) & (dfRecovered.immunN >= 2)).sum():.2f}")
    print(f"Fraction of positive PCRs for which we know the person had a prior infection and at least one vaccination: {(dfRecovered.vaccNatLeast > 0).mean():.2f}")
    print(f"Number of positive PCRs for which we know the person had a prior infection and at least one vaccination: {(dfRecovered.vaccNatLeast > 0).sum():.2f}")
    print()

    print(f"Ratio of PCR positive people for whom we know vaccination/recovery status: "
          f"{dfPos.immunN.notna().sum()/len(dfPos):.2f}")
    print(f"Ratio of PCR positive people for whom we know vaccination/recovery status (of those referred to in Figure 1): "
          f"{dfFigure1_A.immunN.notna().sum()/len(dfFigure1_A):.2f}")
    print(f"Ratio of people for whom we know vaccination/recovery status (immun2YN): {dfPos.immun2YN.notna().sum()/len(dfPos):.2f}")
    print(f"Fraction of immunized (at least once) people that had a prior vaccination (as opposed to a prior infection): {1 - dfPos[dfPos.immunYN == 1].recovered.mean():.2f}")
    print(f"Count of immunized (at least once) people that had a prior vaccination (as opposed to a prior infection): {(dfPos.immunYN == 1).sum() - ((dfPos.immunYN == 1) & (dfPos.recovered == 1)).sum()}")
    print()

    # Note: we're using Figure1_A here to keep the people that had 4 immunizations.
    print(f"Number of tests/tested employees with known number of immunizations in Dec 20/Jan 21: {(dfFigure1_A.immunN.notna() & (dfFigure1_A.samplingMonth2 == 0)).sum()}")
    print(f"Number of tests/tested employees with at least two immunizations in Dec 20/Jan 21: {((dfFigure1_A.samplingMonth2 == 0) & dfFigure1_A.immun2YN).sum()}")
    print()


    print(f"Number of tests/tested employees with known number of immunizations in Oct 21/Nov 21: {(dfFigure1_A.immunN.notna() & (dfFigure1_A.samplingMonth2 == 5)).sum()}")
    print(f"Number of tests/tested employees with at least two immunizations in Oct 21/Nov 21: {((dfFigure1_A.samplingMonth2 == 5) & dfFigure1_A.immun2YN).sum()}")
    print()

    # Note: we don't know the exact number of vaccinations for all employees, but for some we know the minimum number (these are considered in 
    # columns "immunNatLeast" and "immun2YN").Therefore, the given fractions are lower bounds.
    dfDec21Jan22 = dfFigure1_A[dfFigure1_A.samplingMonth2 == 6]
    print(f"Fraction (lower bound) of employees with third vaccination in Dec 21/Jan 22: {(dfDec21Jan22.immunNatLeast >= 3).sum()/len(dfDec21Jan22):.2f}")
    dfWtAlpha = dfFigure1_C[dfFigure1_C.samplingMonth2 < 4]
    dfDeltaOmicron = dfFigure1_C[dfFigure1_C.samplingMonth2 >= 4]

    print(f"Number of employees with at least two immunizations during wildtype/alpha waves: {(dfWtAlpha.immun2YN).sum()}")
    print(f"Fraction (lower bound) of employees with at least two immunizations during wildtype/alpha waves: {(dfWtAlpha.immun2YN).sum()/len(dfWtAlpha):.2f}")
    print(f"Fraction (lower bound) of employees with at least two immunizations during delta/omicron waves: {(dfDeltaOmicron.immun2YN).sum()/len(dfDeltaOmicron):.2f}")
    print(f"Fraction (upper bound) of employees with at least two immunizations during delta/omicron waves: {(dfDeltaOmicron.immun2YN).sum()/(dfDeltaOmicron.immun2YN.notna()).sum():.2f}")

    print(f"Number of employees with prior infection: {dfFigure1_A.recovered.sum()}")
    print(f"Fraction of employees with prior infection: {dfFigure1_A.recovered.mean()}")
    return


@app.cell
def _():
    # posCount = 0
    # sampleCount = 0
    # for samplingMonthCurrent in dfDeltaOmicron.samplingMonth2.unique():
    #     print(samplingMonthCurrent)
    #     dfDeltaOmicronCurrent = dfDeltaOmicron[dfDeltaOmicron.samplingMonth2 == samplingMonthCurrent]
    #     print(dfDeltaOmicronCurrent.immun2YN.sum()/ len(dfDeltaOmicronCurrent))
    #     posCount += dfDeltaOmicronCurrent.immun2YN.sum()
    #     sampleCount += len(dfDeltaOmicronCurrent)
    #     print(dfDeltaOmicronCurrent.immun2YN.sum())
    #     print(len(dfDeltaOmicronCurrent))
    #     print()

    # print(dfDeltaOmicron.immun2YN.sum() / len(dfDeltaOmicron))
    return


@app.cell
def _(dfPos):
    # Fraction of tests with unknown exact number of immunizations but information on vaccinations doses at least received
    (dfPos.immunN.isna() & (dfPos.immunYN.notna() | dfPos.immun2YN.notna())).sum() / len(dfPos)
    return


@app.cell
def _(mo):
    mo.md(r"""## Test centers""")
    return


@app.cell
def _(df, dfAgrdt):
    print(f"PCR counts by test center: {df.testCentre.value_counts()}")
    print(f"Ag-RDT counts by test center: {dfAgrdt.testCentre.value_counts()}")
    print(f"Survey completed by test center: {dfAgrdt.groupby("testCentre").surveyData.agg("sum")}")
    return


@app.cell
def _(mo):
    mo.md(r"""## Summary tables""")
    return


@app.cell
def _(dfAgrdt, writeSummaryTablesEmployees):
    writeSummaryTablesEmployees(dfAgrdt)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Exploring PCRs performed during variant transition times""")
    return


@app.cell
def _(dfPos, returnDeltaOmicronTransitionPCRs, returnWtAlphaTransitionPCRs):
    dfWtAlphaTrans = returnWtAlphaTransitionPCRs(dfPos[dfPos.agrdt.notna()])
    dfWtAlphaTransVariant = dfWtAlphaTrans[dfWtAlphaTrans.variant.notna()]

    dfDeltaOmicronTrans = returnDeltaOmicronTransitionPCRs(dfPos[dfPos.agrdt.notna()])
    dfDeltaOmicronTransVariant = dfDeltaOmicronTrans[dfDeltaOmicronTrans.variant.notna()]
    return (
        dfDeltaOmicronTrans,
        dfDeltaOmicronTransVariant,
        dfWtAlphaTrans,
        dfWtAlphaTransVariant,
    )


@app.cell
def _(
    dfDeltaOmicronTrans,
    dfDeltaOmicronTransVariant,
    dfWtAlphaTrans,
    dfWtAlphaTransVariant,
):
    print(dfWtAlphaTrans.vl.median() - dfWtAlphaTransVariant.vl.median())
    print(dfWtAlphaTrans.vl.mean() - dfWtAlphaTransVariant.vl.mean())
    print()
    print(dfDeltaOmicronTrans.vl.median() - dfDeltaOmicronTransVariant.vl.median())
    print(dfDeltaOmicronTrans.vl.mean() - dfDeltaOmicronTransVariant.vl.mean())
    return


@app.cell
def _(
    dfDeltaOmicronTrans,
    dfDeltaOmicronTransVariant,
    dfWtAlphaTrans,
    dfWtAlphaTransVariant,
    label,
    plt,
    sns,
):
    _, _ax = plt.subplots(2, 2, figsize=(10, 10), sharey=True)
    _ax = _ax.flatten()
    sns.boxplot(data=dfWtAlphaTrans, y='vl', ax=_ax[0])
    sns.boxplot(data=dfWtAlphaTransVariant, y='vl', ax=_ax[1])
    sns.boxplot(data=dfDeltaOmicronTrans, y='vl', ax=_ax[2])
    sns.boxplot(data=dfDeltaOmicronTransVariant, y='vl', ax=_ax[3])
    for _a in _ax:
        _a.set_ylabel(label.vl)
    return


@app.cell
def _(dfWtAlphaTrans, dfWtAlphaTransVariant, label, plt, sns):
    _, _ax = plt.subplots(2, 1, figsize=(5, 5), sharex=True)
    sns.histplot(dfWtAlphaTrans.vl, bins=20, ax=_ax[0])
    sns.histplot(dfWtAlphaTransVariant.vl, bins=20, ax=_ax[1])
    _ax[1].set_xlabel(label.vl)
    return


@app.cell
def _(dfDeltaOmicronTrans, dfDeltaOmicronTransVariant, label, plt, sns):
    _, _ax = plt.subplots(2, 1, figsize=(5, 5), sharex=True)
    sns.histplot(dfDeltaOmicronTrans.vl, bins=20, ax=_ax[0])
    sns.histplot(dfDeltaOmicronTransVariant.vl, bins=20, ax=_ax[1])
    _ax[1].set_xlabel(label.vl)
    return


@app.cell
def _(df, label, legend, pal, plt, replaceLegend, sns):
    _, _ax = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
    sns.histplot(data=df[df.hasTyping == 0], x='vl', hue='variant', multiple='stack', palette=pal.variant, ax=_ax[0])
    _ax[0].set_xlabel(label.vl)
    replaceLegend(_ax[0], legend.variantPatch)
    sns.histplot(data=df[df.hasTyping == 1], x='vl', hue='variant', multiple='stack', palette=pal.variant, ax=_ax[1])
    _ax[1].set_xlabel(label.vl)
    _ax[1].legend(handles=legend.variantPatch)
    replaceLegend(_ax[1], legend.variantPatch)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Specificity""")
    return


@app.cell
def _(dfAllInd, np, order, pd, sampleAgrdtSpec):
    n = 1000
    newData = pd.DataFrame({'samplingMonth2': np.concatenate(list((np.repeat(months, n) for months in order.samplingMonth2)))})
    model, iData = sampleAgrdtSpec(dfAllInd, var='samplingMonth2', cats=sorted(dfAllInd.samplingMonth2.unique()), 
                                   newData=newData, target_accept=0.99)
    return iData, newData


@app.cell
def _(az, defaultdict, iData, newData, np, order):
    statsDictSpec = defaultdict(lambda: defaultdict(dict))
    for _samplingMonth2 in order.samplingMonth2:
        idx = newData.index[newData.samplingMonth2 == _samplingMonth2].tolist()
        _samples = iData.posterior.p.stack(sample=('chain', 'draw'))[idx].values.flatten()
        mean = np.mean(_samples)
        hdi = az.hdi(_samples, 0.94)
        statsDictSpec['samplingMonth2'][_samplingMonth2]['mean'] = mean
        statsDictSpec['samplingMonth2'][_samplingMonth2]['hdi'] = hdi
    return (statsDictSpec,)


@app.cell
def _(abbrvDictPaper, order, statsDictSpec):
    for _samplingMonth2 in order.samplingMonth2:
        print(abbrvDictPaper['samplingMonth2'][_samplingMonth2])
        print(f"Specificity (mean, [94% hdi]): {statsDictSpec['samplingMonth2'][_samplingMonth2]['mean']:.3f}{[statsDictSpec['samplingMonth2'][_samplingMonth2]['hdi']]}")
        print()
    return


if __name__ == "__main__":
    app.run()
