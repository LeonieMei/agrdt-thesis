import pandas as pd
import seaborn as sns
from datetime import date

from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

DATE_RANGE_DOMINANT_WILDTYPE = date(2019, 11, 1), date(2021, 2, 7)
DATE_RANGE_DOMINANT_ALPHA = date(2021, 3, 29), date(2021, 5, 23)
DATE_RANGE_DOMINANT_DELTA = date(2021, 7, 12), date(2021, 12, 19)
DATE_RANGE_DOMINANT_OMICRON = date(2022, 1, 17), date.today()

DPI = 300
CM = 1 / 2.54
COL_WIDTH = 5.91 / 2  # 3.54
DINA4_HEIGHT = 8.35

# Define size constants
TITLE_SIZE = 14
TEXT_SIZE = 12
LEGEND_TITLE_SIZE = 12
LEGEND_TEXT_SIZE = 10
LEGEND_KEY_SIZE = 20
PLOT_TITLE_SIZE = 16
ANNOTATION_LETTER_SIZE = 14
ANNOTATION_COORDS = (-0.175, 1)

MARKERSIZE_LEGEND_VARIANT_POINTS = 6

_GRAY = "0.6"

_DAY_LEVELS2 = list(range(-2, 17))
_DAY_LEVELS2.append(19)
_DAY_LEVELS = [0, 2, 4, 6, 8]

_AGE_BRACKETS = ((16.999, 30), (30, 50), (50, 70))

_AGE_BINS = [pd.Interval(low, up, closed="right") for low, up in _AGE_BRACKETS]

_AGE_BINS_STR = (
    f"17 - {_AGE_BRACKETS[0][1]}",
    f"{_AGE_BRACKETS[1][0]} - " f"{_AGE_BRACKETS[1][1]}",
    f"{_AGE_BRACKETS[2][0]} - " f"{_AGE_BRACKETS[2][1]}",
)

# Define brackets for days post symptom onset and corresponding labels.
_SYMP_DAYS_BRACKETS = ((-0.001, 2), (2, 7))
_SYMP_DAYS4_BRACKETS = ((-0.001, 1), (1, 7))

_SYMP_DAYS_BINS = [
    pd.Interval(low, up, closed="right") for low, up in _SYMP_DAYS_BRACKETS
]

_SYMP_DAYS4_BINS = [
    pd.Interval(low, up, closed="right") for low, up in _SYMP_DAYS4_BRACKETS
]

_SYMP_DAYS_BINS_STR = (
    f"0 - {_SYMP_DAYS_BRACKETS[0][1]} days",
    f"{_SYMP_DAYS_BRACKETS[0][1]} - {_SYMP_DAYS_BRACKETS[1][1]} days",
)

_SYMP_DAYS4_BINS_STR = (
    f"0 - {_SYMP_DAYS4_BRACKETS[0][1]} day",
    f"{_SYMP_DAYS4_BRACKETS[0][1]} - {_SYMP_DAYS4_BRACKETS[1][1]} days",
)

GGPLOT_PALETTE = [
    "#F8766D",
    "#7CAE00",
    "#00BFC4",
    "#C77CFF",
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#D55E00",
    "#CC79A7",
    "#999999",
]

_COLORS_VARIANTS = GGPLOT_PALETTE[:4]
_COLORS_VARIANTS2 = [GGPLOT_PALETTE[0], GGPLOT_PALETTE[2]] + GGPLOT_PALETTE[4:6]
_COLORS_GENDER = GGPLOT_PALETTE[8:10]
_COLORS_IMMUN = sns.color_palette("rocket_r") + ["black"]
_COLORS_IMMUN_2YN = GGPLOT_PALETTE[4:6]
_COLORS_RES = [GGPLOT_PALETTE[0], GGPLOT_PALETTE[2], _GRAY]
_COLORS_DAY = dict(
    zip(_DAY_LEVELS, sns.color_palette("magma", n_colors=len(_DAY_LEVELS)))
)


FEMALE = 0
MALE = 1

WT = 0
ALPHA = 1
DELTA = 2
OMICRON = 3
DELTA2 = 1
BA1 = 2
BA2 = 3

WT_LABEL = "Pre-VOC"
ALPHA_LABEL = "Alpha"
DELTA_LABEL = "Delta"
OMICRON_LABEL = "Omicron"
OMICRON_BA1_LABEL = "Omicron BA.1"
OMICRON_BA2_LABEL = "Omicron BA.2"

NEG = 0
POS = 1
MISSING = 2

CCM = 0
CVK = 1
CBF = 2

ROCHE = 0
ABBOTT = 1

# _COLOR_DAYS_POST_FIRST_POS_PCR_BIN = sns.color_palette()


def _combine(list1, list2):
    return list1 + list2


SYMPTOMS_DISPLAY_DICT = {
    "ill": "Feeling ill",
    "fatigue": "Fatigue",
    "headache": "Headache",
    "vertigo": "Vertigo",
    "melalgia": "Melalgia",
    "fever": "Fever",
    "cough": "Cough",
    "runnyNose": "Runny nose",
    "soreThroat": "Sore throat",
    "dyspnea": "Dyspnea",
    "noSmell": "Loss of smell",
    "noTaste": "Loss of taste",
    "nausea": "Nausea",
    "noAppetite": "No appetite",
    "vomiting": "Vomiting",
    "diarrhea": "Diarrhea",
}

# Palettes for plotting.
_PALETTES = {
    "symp": {0: "green", 1: "orangered"},
    "res": {NEG: _COLORS_RES[NEG], POS: _COLORS_RES[POS]},
    "resPCR": {NEG: _COLORS_RES[NEG], POS: _COLORS_RES[POS]},
    "agrdt": {NEG: "black", POS: "white"},
    "agrdtYN": {0: "gray", 1: "darkred"},
    "gender": {"F": _COLORS_GENDER[FEMALE], "M": _COLORS_GENDER[MALE]},
    "variant": {
        "wildtype": _COLORS_VARIANTS[WT],
        "alpha": _COLORS_VARIANTS[ALPHA],
        "delta": _COLORS_VARIANTS[DELTA],
        "omicron": _COLORS_VARIANTS[OMICRON],
    },
    "variantCodes": {i: _COLORS_VARIANTS[i] for i in range(4)},
    "variantCode": {i: _COLORS_VARIANTS[i] for i in range(4)},
    "VariantCode": {i: _COLORS_VARIANTS2[i] for i in range(4)},
    "immunN": dict(zip(range(5), _COLORS_IMMUN)),
    # "immunN": dict(zip((0.0, 1.0, 2.0, 3.0, 4.0), _COLORS_IMMUN)),
    "immun2YN": {
        0.0: _COLORS_IMMUN_2YN[0],
        1.0: _COLORS_IMMUN_2YN[1],
        "0.0": _COLORS_IMMUN_2YN[0],
        "1.0": _COLORS_IMMUN_2YN[1],
    },
    "binDaysPostOnset": {
        sympDaysBin: col
        for col, sympDaysBin in zip(GGPLOT_PALETTE[7:9], _SYMP_DAYS_BINS, strict=True)
    },
    "binDaysPostOnsetCode": {i: GGPLOT_PALETTE[7:9][i] for i in range(2)},
    "binDaysPostOnsetReverseCode": {i: GGPLOT_PALETTE[7:9][i][::-1] for i in (1, 0)},
    "binDaysPostOnset4": {
        sympDaysBin: col
        for col, sympDaysBin in zip(GGPLOT_PALETTE[7:9], _SYMP_DAYS4_BINS, strict=True)
    },
    "binDaysPostOnsetReverse": {
        sympDaysBin: col
        for col, sympDaysBin in zip(
            GGPLOT_PALETTE[7:9][::-1], _SYMP_DAYS_BINS[::-1], strict=True
        )
    },
    "binDaysPostOnset4Reverse": {
        sympDaysBin: col
        for col, sympDaysBin in zip(
            GGPLOT_PALETTE[7:9][::-1], _SYMP_DAYS4_BINS[::-1], strict=True
        )
    },
    # "binDaysPostFirstPosPcr": {
    #     0: _COLOR_DAYS_POST_FIRST_POS_PCR_BIN[0],
    #     1: _COLOR_DAYS_POST_FIRST_POS_PCR_BIN[1],
    # },
    "cellculture": {
        "BA.1_CaCo2": "firebrick",
        "BA.1_VeroE6": "maroon",
        "BA.2_CaCo2": "turquoise",
        "BA.2_VeroE6": "lightseagreen",
        "WT_CaCo2": "hotpink",
        "WT_VeroE6": "palevioletred",
        "Mock_CaCo2": "skyblue",
        "Mock_VeroE6": "deepskyblue",
    },
    "day": _COLORS_DAY,
}


_ABBRVS = {
    "agrdt": {0.0: "negative", 1.0: "positive"},
    "agrdtYN": {False: "no", True: "yes"},
    "testline": {0.0: 0, 1.0: 1, 2.0: 2, 3.0: 3},
    "pcrPositive": {False: "negative", True: "positive"},
    "symptoms": {0: "asymptomatic", 1: "symptomatic"},
    "gender": {"F": "women", "M": "men"},
    "variant": {
        "wildtype": WT_LABEL,
        "alpha": ALPHA_LABEL,
        "delta": DELTA_LABEL,
        "omicron": OMICRON_LABEL,
    },
    "ageBin": dict(zip(_AGE_BINS, _AGE_BINS_STR)),
    "binDaysPostOnset": dict(zip(_SYMP_DAYS_BINS, _SYMP_DAYS_BINS_STR)),
    "binDaysPostOnset4": dict(zip(_SYMP_DAYS4_BINS, _SYMP_DAYS4_BINS_STR)),
    "binDaysPostFirstPosPcr": {0: "0", 1: ">= 7"},
    "vaccN": dict(zip(range(4), range(4))),
    "immunN": dict(zip(range(4), range(4))),
    "immun2YN": {0.0: "no", 1.0: "yes"},
    "symptom": {0: "no", 1: "yes"},
    "recoveredRecently": {0: "no", 1: "yes"},
    # First month is November 2020 but we actually only have antigen test
    # data for this month and no PCR data (missing PCR ids).
    "samplingMonth": dict(
        zip(
            [float(i) for i in range(15)],
            (
                "Dec '20",
                "Jan '21",
                "Feb '21",
                "March '21",
                "April '21",
                "May '21",
                "June '21",
                "July '21",
                "Aug '21",
                "Sep '21",
                "Oct '21",
                "Nov '21",
                "Dec '21",
                "Jan '22",
                "Feb '22",
            ),
            strict=True,
        )
    ),
    "samplingMonth2": dict(
        zip(
            [float(i) for i in range(8)],
            (
                "Dec '20 - Jan '21",
                "Feb '21 - March '21",
                "April '21 - May '21",
                "June '21 - July '21",
                "Aug '21 - Sep '21",
                "Oct '21 - Nov '21",
                "Dec '21 - Jan '22",
                "Feb '22",
            ),
            strict=True,
        )
    ),
    "cellculture": {
        "BA.1_CaCo2": "BA.1, Caco-2",
        "BA.1_VeroE6": "BA.1, Vero E6",
        "BA.2_CaCo2": "BA.2, Caco-2",
        "BA.2_VeroE6": "BA.2, Vero E6",
        "WT_CaCo2": "Wildtype, Caco-2",
        "WT_VeroE6": "Wildtype, Vero E6",
        "Mock_CaCo2": "Mock, Caco-2",
        "Mock_VeroE6": "Mock, Vero E6",
    },
    "day": dict(zip(_DAY_LEVELS, _DAY_LEVELS)),
}


_LABELS = {
    "hdi": "94% HPDI",
    "sens": "Ag-RDT sensitivity",
    "sensPercent": "Ag-RDT sensitivity (%)",
    "res": "Ag-RDT result",
    "resPCR": "PCR result",
    "samplingMonth2": "Sampling period",
    "testline": "Testline strength",
    "vl": r"$\mathregular{Log_{10}}$ viral load",
    "zVl": r"$\mathregular{Log_{10}}$ viral load",
    "coi": r"$\mathregular{Log_{10}}$ coi",
    "log10Coi": r"$\mathregular{Log_{10}}$ coi",
    "N_RNA": r"$\mathregular{Log_{10}}$ coi/viral load",
    # This is the format that gets properly displayed when saving a C{
    # pd.DataFrame} (it doesn't work in C{matplotlib} figures though.
    "vlRNA": r"$\mathregular{Log_{10}}$ SARS-CoV-2 mRNA/mL",
    "vlDf": "log\u2081\u2080 viral load",
    "zVlDf": "log\u2081\u2080 viral load",
    "zAge": "Age",
    "days": "Day post symptom onset",
    "daysPostOnset": "Day post symptom onset",
    "binDaysPostOnset": "Day post symptom onset",
    "binDaysPostOnsetReverse": "Day post symptom onset",
    "binDaysPostOnset4": "Day post symptom onset",
    "binDaysPostOnset4Reverse": "Day post symptom onset",
    "day": "Day post symptom onset",
    "dayLegend": "Days post\nsymptom\nonset",
    "daysSinceFirstPosPcr": "Day since first positive PCR",
    "pcr": "PCR positive rate",
    "symptoms": "Symptomatic status",
    "gender": "Sex",
    "variant": "SARS-CoV-2 variant",
    "immunN": "Sum of vaccinations and\nprior infections",
    "immunNLegend": "Sum of vaccinations\nand prior infections",
    "agrdtYN": "Ag-RDT performed",
    "immun2YN": "Multiply immunized",
    "immun2YNLegend": "Multiply immunized",
    "zIncidenceChange": "Relative change in 7-day incidence",
    "zRtRolling4": "7-day R value",
}

_LABELS_LATEX = _LABELS.copy()
_LABELS_LATEX["vl"] = "Log$_{10}$ viral load"
_LABELS_LATEX["zVl"] = "Log$_{10}$ viral load"

_ABBRV_TO_COLOR = {
    "binDaysPostOnset": tuple(
        zip(
            _ABBRVS["binDaysPostOnset"].values(), _PALETTES["binDaysPostOnset"].values()
        )
    ),
    "binDaysPostOnset4": tuple(
        zip(
            _ABBRVS["binDaysPostOnset4"].values(),
            _PALETTES["binDaysPostOnset4"].values(),
        )
    ),
    # "binDaysPostFirstPosPcr": tuple(
    #     zip(
    #         _ABBRVS["binDaysPostFirstPosPcr"].values(),
    #         _PALETTES["binDaysPostFirstPosPcr"].values(),
    #     )
    # ),
    "immunN": tuple(zip(_ABBRVS["immunN"].values(), _PALETTES["immunN"].values())),
    "immun2YN": tuple(
        zip(_ABBRVS["immun2YN"].values(), _PALETTES["immun2YN"].values())
    ),
    "agrdtYN": tuple(zip(_ABBRVS["agrdtYN"].values(), _PALETTES["agrdtYN"].values())),
    "cellculture": tuple(
        zip(_ABBRVS["cellculture"].values(), _PALETTES["cellculture"].values())
    ),
}


# Legends for plotting.
_LEGENDS = {
    "sympAsymp": [
        Line2D(
            [0],
            [0],
            marker="o",
            color="orangered",
            label="symptomatic",
            markerfacecolor="orangered",
            markersize=4,
            linestyle=None,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="green",
            label="asymptomatic",
            markerfacecolor="green",
            markersize=4,
            linestyle=None,
        ),
    ],
    "sympAsympRt": [
        Line2D(
            [0],
            [0],
            marker="o",
            color="orangered",
            label="symptomatic",
            markerfacecolor="orangered",
            markersize=4,
            linestyle=None,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="green",
            label="asymptomatic",
            markerfacecolor="green",
            markersize=4,
            linestyle=None,
        ),
        Line2D([0], [0], color="blue", linestyle="-", label=r"$R_t$"),
    ],
    "sympAsympPatch": [
        Patch(color="tab:blue", label="asymptomatic"),
        Patch(color="tab:orange", label="symptomatic"),
    ],
    "newCases": [
        Line2D([0], [0], color=_COLORS_RES[NEG], linestyle="-", label="Charité"),
        Line2D([0], [0], color=_COLORS_RES[POS], linestyle="-", label="Germany"),
    ],
    "resPatch": [
        Rectangle(
            (0, 0), 0, 0, label=_LABELS["res"], facecolor="none", edgecolor="none"
        ),
        Patch(color=_COLORS_RES[POS], label="positive", alpha=1),
        Patch(color=_COLORS_RES[NEG], label="negative", alpha=1),
    ],
    "resPCRPatch": [
        Rectangle(
            (0, 0), 0, 0, label=_LABELS["resPCR"], facecolor="none", edgecolor="none"
        ),
        Patch(color=_COLORS_RES[POS], label="positive", alpha=1),
        Patch(color=_COLORS_RES[NEG], label="negative", alpha=1),
    ],
    "resPoints": [
        Rectangle(
            (0, 0), 0, 0, label=_LABELS["res"], facecolor="none", edgecolor="none"
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color=_PALETTES["res"][POS],
            label="positive",
            markerfacecolor=_PALETTES["res"][POS],
            markeredgecolor="none",
            markersize=8,
            linestyle="None",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color=_PALETTES["res"][NEG],
            label="negative",
            markerfacecolor=_PALETTES["res"][NEG],
            markeredgecolor="none",
            markersize=8,
            linestyle="None",
        ),
    ],
    "agrdtYNPoints": _combine(
        [
            Rectangle(
                (0, 0),
                0,
                0,
                label=_LABELS["agrdtYN"],
                facecolor="none",
                edgecolor="none",
            )
        ],
        [
            Line2D(
                [0],
                [0],
                marker="o",
                color=col,
                label=cat,
                markerfacecolor=col,
                markersize=8,
                linestyle="None",
            )
            for cat, col in _ABBRV_TO_COLOR["agrdtYN"]
        ],
    ),
    "variantLines": [
        Rectangle(
            (0, 0), 0, 0, label="SARS-CoV-2 variant", facecolor="none", edgecolor="none"
        ),
        Line2D(
            [0],
            [0],
            marker=None,
            color=_COLORS_VARIANTS[WT],
            label=WT_LABEL,
            markerfacecolor=_COLORS_VARIANTS[WT],
            markersize=8,
            linestyle="-",
        ),
        Line2D(
            [0],
            [0],
            marker=None,
            color=_COLORS_VARIANTS[ALPHA],
            label=ALPHA_LABEL,
            markerfacecolor=_COLORS_VARIANTS[ALPHA],
            markersize=8,
            linestyle="-",
        ),
        Line2D(
            [0],
            [0],
            marker=None,
            color=_COLORS_VARIANTS[DELTA],
            label=DELTA_LABEL,
            markerfacecolor=_COLORS_VARIANTS[DELTA],
            markersize=8,
            linestyle="-",
        ),
        Line2D(
            [0],
            [0],
            marker=None,
            color=_COLORS_VARIANTS[OMICRON],
            label=OMICRON_LABEL,
            markerfacecolor=_COLORS_VARIANTS[OMICRON],
            markersize=8,
            linestyle="-",
        ),
    ],
    "variantCodeLines": [
        Rectangle(
            (0, 0), 0, 0, label="SARS-CoV-2 variant", facecolor="none", edgecolor="none"
        ),
        Line2D(
            [0],
            [0],
            marker=None,
            color=_COLORS_VARIANTS[WT],
            label=WT_LABEL,
            markerfacecolor=_COLORS_VARIANTS[WT],
            markersize=8,
            linestyle="-",
        ),
        Line2D(
            [0],
            [0],
            marker=None,
            color=_COLORS_VARIANTS[ALPHA],
            label=ALPHA_LABEL,
            markerfacecolor=_COLORS_VARIANTS[ALPHA],
            markersize=8,
            linestyle="-",
        ),
        Line2D(
            [0],
            [0],
            marker=None,
            color=_COLORS_VARIANTS[DELTA],
            label=DELTA_LABEL,
            markerfacecolor=_COLORS_VARIANTS[DELTA],
            markersize=8,
            linestyle="-",
        ),
        Line2D(
            [0],
            [0],
            marker=None,
            color=_COLORS_VARIANTS[OMICRON],
            label=OMICRON_LABEL,
            markerfacecolor=_COLORS_VARIANTS[OMICRON],
            markersize=8,
            linestyle="-",
        ),
    ],
    "VariantCodeLines": [
        Rectangle(
            (0, 0), 0, 0, label="SARS-CoV-2 variant", facecolor="none", edgecolor="none"
        ),
        Line2D(
            [0],
            [0],
            marker=None,
            color=_COLORS_VARIANTS2[WT],
            label=WT_LABEL,
            markerfacecolor=_COLORS_VARIANTS2[WT],
            markersize=8,
            linestyle="-",
        ),
        Line2D(
            [0],
            [0],
            marker=None,
            color=_COLORS_VARIANTS2[DELTA2],
            label=DELTA_LABEL,
            markerfacecolor=_COLORS_VARIANTS2[DELTA2],
            markersize=8,
            linestyle="-",
        ),
        Line2D(
            [0],
            [0],
            marker=None,
            color=_COLORS_VARIANTS2[BA1],
            label=OMICRON_BA1_LABEL,
            markerfacecolor=_COLORS_VARIANTS2[BA1],
            markersize=8,
            linestyle="-",
        ),
        Line2D(
            [0],
            [0],
            marker=None,
            color=_COLORS_VARIANTS2[BA2],
            label=OMICRON_BA2_LABEL,
            markerfacecolor=_COLORS_VARIANTS2[BA2],
            markersize=8,
            linestyle="-",
        ),
    ],
    "variantPoints": [
        Rectangle(
            (0, 0), 0, 0, label="SARS-CoV-2 variant", facecolor="none", edgecolor="none"
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color=_COLORS_VARIANTS[WT],
            label=WT_LABEL,
            markerfacecolor=_COLORS_VARIANTS[WT],
            markersize=MARKERSIZE_LEGEND_VARIANT_POINTS,
            linestyle="None",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color=_COLORS_VARIANTS[ALPHA],
            label=ALPHA_LABEL,
            markerfacecolor=_COLORS_VARIANTS[ALPHA],
            markersize=MARKERSIZE_LEGEND_VARIANT_POINTS,
            linestyle="None",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color=_COLORS_VARIANTS[DELTA],
            label=DELTA_LABEL,
            markerfacecolor=_COLORS_VARIANTS[DELTA],
            markersize=MARKERSIZE_LEGEND_VARIANT_POINTS,
            linestyle="None",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color=_COLORS_VARIANTS[OMICRON],
            label=OMICRON_LABEL,
            markerfacecolor=_COLORS_VARIANTS[OMICRON],
            markersize=MARKERSIZE_LEGEND_VARIANT_POINTS,
            linestyle="None",
        ),
    ],
    "variantPoints2": [
        Rectangle(
            (0, 0), 0, 0, label="SARS-CoV-2 variant", facecolor="none", edgecolor="none"
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color=_COLORS_VARIANTS[WT],
            label=WT_LABEL,
            markerfacecolor=_COLORS_VARIANTS[WT],
            markersize=MARKERSIZE_LEGEND_VARIANT_POINTS,
            linestyle="None",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color=_COLORS_VARIANTS[DELTA],
            label=DELTA_LABEL,
            markerfacecolor=_COLORS_VARIANTS[DELTA],
            markersize=MARKERSIZE_LEGEND_VARIANT_POINTS,
            linestyle="None",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color=_COLORS_VARIANTS[OMICRON],
            label=OMICRON_LABEL,
            markerfacecolor=_COLORS_VARIANTS[OMICRON],
            markersize=MARKERSIZE_LEGEND_VARIANT_POINTS,
            linestyle="None",
        ),
    ],
    "VariantCodePoints": [
        Rectangle(
            (0, 0), 0, 0, label="SARS-CoV-2 variant", facecolor="none", edgecolor="none"
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color=_COLORS_VARIANTS2[WT],
            label=WT_LABEL,
            markerfacecolor=_COLORS_VARIANTS2[WT],
            markersize=MARKERSIZE_LEGEND_VARIANT_POINTS,
            linestyle="None",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color=_COLORS_VARIANTS2[DELTA],
            label=DELTA_LABEL,
            markerfacecolor=_COLORS_VARIANTS2[DELTA2],
            markersize=MARKERSIZE_LEGEND_VARIANT_POINTS,
            linestyle="None",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color=_COLORS_VARIANTS2[BA1],
            label=OMICRON_BA1_LABEL,
            markerfacecolor=_COLORS_VARIANTS2[BA1],
            markersize=MARKERSIZE_LEGEND_VARIANT_POINTS,
            linestyle="None",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color=_COLORS_VARIANTS2[BA2],
            label=OMICRON_BA2_LABEL,
            markerfacecolor=_COLORS_VARIANTS2[BA2],
            markersize=MARKERSIZE_LEGEND_VARIANT_POINTS,
            linestyle="None",
        ),
    ],
    "variantPointsLines": [
        Rectangle(
            (0, 0), 0, 0, label="SARS-CoV-2 variant", facecolor="none", edgecolor="none"
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color=_COLORS_VARIANTS[WT],
            label=WT_LABEL,
            markerfacecolor=_COLORS_VARIANTS[WT],
            markersize=8,
            linestyle="-",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color=_COLORS_VARIANTS[ALPHA],
            label=ALPHA_LABEL,
            markerfacecolor=_COLORS_VARIANTS[ALPHA],
            markersize=8,
            linestyle="-",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color=_COLORS_VARIANTS[DELTA],
            label=DELTA_LABEL,
            markerfacecolor=_COLORS_VARIANTS[DELTA],
            markersize=8,
            linestyle="-",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color=_COLORS_VARIANTS[OMICRON],
            label=OMICRON_LABEL,
            markerfacecolor=_COLORS_VARIANTS[OMICRON],
            markersize=8,
            linestyle="-",
        ),
    ],
    "variantPatch": [
        Rectangle(
            (0, 0), 0, 0, label="SARS-CoV-2 variant", facecolor="none", edgecolor="none"
        ),
        Patch(color=_COLORS_VARIANTS[WT], label=WT_LABEL),
        Patch(color=_COLORS_VARIANTS[ALPHA], label=ALPHA_LABEL),
        Patch(color=_COLORS_VARIANTS[DELTA], label=DELTA_LABEL),
        Patch(color=_COLORS_VARIANTS[OMICRON], label=OMICRON_LABEL),
    ],
    "binDaysPostOnset": _combine(
        [
            Rectangle(
                (0, 0),
                0,
                0,
                label=_LABELS["binDaysPostOnset"],
                facecolor="none",
                edgecolor="none",
            )
        ],
        [
            Line2D(
                [0],
                [0],
                marker="o",
                color=col,
                label=cat,
                markerfacecolor=col,
                markersize=7,
            )
            for cat, col in _ABBRV_TO_COLOR["binDaysPostOnset"]
        ],
    ),
    "binDaysPostOnset4": _combine(
        [
            Rectangle(
                (0, 0),
                0,
                0,
                label=_LABELS["binDaysPostOnset4"],
                facecolor="none",
                edgecolor="none",
            )
        ],
        [
            Line2D(
                [0],
                [0],
                marker="o",
                color=col,
                label=cat,
                markerfacecolor=col,
                markersize=7,
            )
            for cat, col in _ABBRV_TO_COLOR["binDaysPostOnset4"]
        ],
    ),
    "binDaysPostOnsetLines": _combine(
        [
            Rectangle(
                (0, 0),
                0,
                0,
                label=_LABELS["binDaysPostOnset"],
                facecolor="none",
                edgecolor="none",
            )
        ],
        [
            Line2D(
                [0],
                [0],
                marker=None,
                color=col,
                label=cat,
                markerfacecolor=col,
                markersize=7,
            )
            for cat, col in _ABBRV_TO_COLOR["binDaysPostOnset"]
        ],
    ),
    "binDaysPostOnset4Lines": _combine(
        [
            Rectangle(
                (0, 0),
                0,
                0,
                label=_LABELS["binDaysPostOnset4"],
                facecolor="none",
                edgecolor="none",
            )
        ],
        [
            Line2D(
                [0],
                [0],
                marker=None,
                color=col,
                label=cat,
                markerfacecolor=col,
                markersize=7,
            )
            for cat, col in _ABBRV_TO_COLOR["binDaysPostOnset4"]
        ],
    ),
    "immunNPatch": _combine(
        [
            Rectangle(
                (0, 0),
                0,
                0,
                label=_LABELS["immunNLegend"],
                facecolor="none",
                edgecolor="none",
            )
        ],
        [Patch(color=col, label=cat) for cat, col in _ABBRV_TO_COLOR["immunN"]],
    ),
    "immun2YNPatch": _combine(
        [
            Rectangle(
                (0, 0),
                0,
                0,
                label=_LABELS["immun2YN"],
                facecolor="none",
                edgecolor="none",
            )
        ],
        [
            Patch(color=_COLORS_IMMUN_2YN[int(i)], label=label)
            for i, label in _ABBRVS["immun2YN"].items()
        ],
    ),
    "immun2YNLines": _combine(
        [
            Rectangle(
                (0, 0),
                0,
                0,
                label=_LABELS["immun2YNLegend"],
                facecolor="none",
                edgecolor="none",
            )
        ],
        [
            Line2D(
                [0],
                [0],
                marker=None,
                color=col,
                label=cat,
                markerfacecolor=col,
                markersize=7,
            )
            for cat, col in _ABBRV_TO_COLOR["immun2YN"]
        ],
    ),
    "immun2YNPoints": _combine(
        [
            Rectangle(
                (0, 0),
                0,
                0,
                label=_LABELS["immun2YNLegend"],
                facecolor="none",
                edgecolor="none",
            )
        ],
        [
            Line2D(
                [0],
                [0],
                marker="o",
                color=col,
                label=cat,
                markerfacecolor=col,
                markersize=7,
                linestyle="none",
            )
            for cat, col in _ABBRV_TO_COLOR["immun2YN"]
        ],
    ),
    "binDaysPostOnsetPointsLines": [
        Rectangle(
            (0, 0),
            0,
            0,
            label=_LABELS["binDaysPostOnset"],
            facecolor="none",
            edgecolor="none",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="gray",
            label=_SYMP_DAYS_BINS_STR[0],
            markerfacecolor="gray",
            markersize=6,
            linestyle="-",
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            color="gray",
            label=_SYMP_DAYS_BINS_STR[1],
            markerfacecolor="gray",
            markersize=6,
            linestyle="-",
        ),
    ],
    "binDaysPostOnset4PointsLines": [
        Rectangle(
            (0, 0),
            0,
            0,
            label=_LABELS["binDaysPostOnset4"],
            facecolor="none",
            edgecolor="none",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="gray",
            label=_SYMP_DAYS_BINS_STR[0],
            markerfacecolor="gray",
            markersize=6,
            linestyle="-",
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            color="gray",
            label=_SYMP_DAYS_BINS_STR[1],
            markerfacecolor="gray",
            markersize=6,
            linestyle="-",
        ),
    ],
    "cellculture": [
        Line2D(
            [0],
            [0],
            marker=None,
            color=col,
            label=cat,
            markerfacecolor=col,
            markersize=8,
            linestyle="-",
        )
        for cat, col in _ABBRV_TO_COLOR["cellculture"]
    ],
    "variantLinesWtOmicron": [
        Rectangle(
            (0, 0), 0, 0, label="SARS-CoV-2 variant", facecolor="none", edgecolor="none"
        ),
        Line2D(
            [0],
            [0],
            marker=None,
            color=_COLORS_VARIANTS[WT],
            label=WT_LABEL,
            markerfacecolor=_COLORS_VARIANTS[WT],
            markersize=8,
            linestyle="-",
        ),
        Line2D(
            [0],
            [0],
            marker=None,
            color=_COLORS_VARIANTS[OMICRON],
            label=OMICRON_LABEL,
            markerfacecolor=_COLORS_VARIANTS[OMICRON],
            markersize=8,
            linestyle="-",
        ),
    ],
    "NandRNA": [
        Line2D(
            [0],
            [0],
            marker=None,
            color="blue",
            label=_LABELS["vl"],
            markerfacecolor="blue",
            markersize=8,
            linestyle="-",
        ),
        Line2D(
            [0],
            [0],
            marker=None,
            color="red",
            label=_LABELS["coi"],
            markerfacecolor="red",
            markersize=8,
            linestyle="-",
        ),
    ],
    "firstSecondThirdInfection": [
        Line2D(
            [0],
            [0],
            marker=None,
            color="blue",
            label="First infection",
            markerfacecolor="blue",
            markersize=8,
            linestyle="-",
        ),
        Line2D(
            [0],
            [0],
            marker=None,
            color="green",
            label="Second or third infection",
            markerfacecolor="green",
            markersize=8,
            linestyle="-",
        ),
    ],
    "dayLines": [
        Rectangle(
            (0, 0),
            0,
            0,
            label=_LABELS["dayLegend"],
            facecolor="none",
            edgecolor="none",
        )
    ]
    + [
        Line2D(
            [0],
            [0],
            marker=None,
            color=_COLORS_DAY[i],
            label=_ABBRVS["day"][i],
            markerfacecolor=_COLORS_DAY[i],
            markersize=8,
            linestyle="-",
        )
        for i in _DAY_LEVELS
    ],
}

_ORDERS = {
    "agrdt": tuple(_PALETTES["agrdt"]),
    "agrdtYN": tuple(_PALETTES["agrdtYN"]),
    "pcrPositive": tuple(_PALETTES["resPCR"]),
    "variant": tuple(_PALETTES["variant"]),
    "gender": tuple(_PALETTES["gender"]),
    "symp": tuple(_PALETTES["symp"]),
    "binDaysPostOnset": tuple(_PALETTES["binDaysPostOnset"]),
    "binDaysPostOnsetReverse": tuple(_PALETTES["binDaysPostOnset"])[::-1],
    "binDaysPostOnset4": tuple(_PALETTES["binDaysPostOnset4"]),
    "binDaysPostOnset4Reverse": tuple(_PALETTES["binDaysPostOnset4"])[::-1],
    "immunN": tuple(_PALETTES["immunN"]),
    "immun2YN": tuple(_PALETTES["immun2YN"]),
    "samplingMonth": (
        "2020-11",
        "2020-12",
        "2021-01",
        "2021-02",
        "2021-03",
        "2021-04",
        "2021-05",
        "2021-06",
        "2021-07",
        "2021-08",
        "2021-09",
        "2021-10",
        "2021-11",
        "2021-12",
        "2022-01",
        "2022-02",
    ),
    "samplingMonth2": tuple(_ABBRVS["samplingMonth2"]),
    "cellculture": tuple(_ABBRVS["cellculture"]),
}


IMMUN_SYMP_DICT = {
    "0, 1": "no",
    "1, 1": "yes",
    "dummy": "",
    "0, 0": "no",
    "1, 0": "yes",
}

VARIANT_SYMP_DICT = {
    "wildtype, 1": _ABBRVS["variant"]["wildtype"],
    "alpha, 1": _ABBRVS["variant"]["alpha"],
    "delta, 1": _ABBRVS["variant"]["delta"],
    "omicron, 1": _ABBRVS["variant"]["omicron"],
    "dummy": "",
    "wildtype, 0": _ABBRVS["variant"]["wildtype"],
    "alpha, 0": _ABBRVS["variant"]["alpha"],
    "delta, 0": _ABBRVS["variant"]["delta"],
    "omicron, 0": _ABBRVS["variant"]["omicron"],
}

VARIANT_SYMP_PAL = dict(
    zip(
        VARIANT_SYMP_DICT,
        list(_PALETTES["variant"].values())
        + ["k"]
        + list(_PALETTES["variant"].values()),
    )
)
IMMUN_SYMP_PAL = dict(
    zip(
        IMMUN_SYMP_DICT,
        list(_PALETTES["immun2YN"].values())[:2]
        + ["k"]
        + list(_PALETTES["immun2YN"].values())[:2],
    )
)

_PARAM_NAME_MAPPINGS = {
    "zAge": "Age",
    "gender[1]": "Male",
    "testDevice[1]": "Roche device",
    "immun2YN[1]": _LABELS["immun2YN"],
    "variant[1]": f"{ALPHA_LABEL}",
    "variant[2]": f"{DELTA_LABEL}",
    "variant[3]": f"{OMICRON_LABEL}",
    "VariantCode[1]": f"{DELTA_LABEL}",
    "VariantCode[2]": f"{OMICRON_BA1_LABEL}",
    "VariantCode[3]": f"{OMICRON_BA2_LABEL}",
    "binDaysPostOnsetReverse[1]": "0-2 d.p.o.",
    "binDaysPostOnset[1]": "2-7 d.p.o.",
    "daysVariant[0]": ("2-7 d.p.o. | Pre-VOC"),
    "daysVariant[1]": ("2-7 d.p.o. | Alpha"),
    "daysVariant[2]": ("2-7 d.p.o. | Delta"),
    "daysVariant[3]": ("2-7 d.p.o. | Omicron"),
    "daysReverseVariant[0]": "0-2 d.p.o. | Pre-VOC",
    "daysReverseVariant[1]": ("0-2 d.p.o. | Alpha"),
    "daysReverseVariant[2]": ("0-2 d.p.o. | Delta"),
    "daysReverseVariant[3]": ("0-2 d.p.o. | Omicron"),
    "days4Variant[0]": (f"1-7 d.p.o. | Pre-VOC"),
    "days4Variant[1]": ("1-7 d.p.o. | Alpha"),
    "days4Variant[2]": ("1-7 d.p.o. | Delta"),
    "days4Variant[3]": ("1-7 d.p.o. | Omicron"),
    "days4ReverseVariant[0]": "0-1 d.p.o. | Pre-VOC",
    "days4ReverseVariant[1]": ("0-1 d.p.o. | Alpha"),
    "days4ReverseVariant[2]": ("0-1 d.p.o. | Delta"),
    "days4ReverseVariant[3]": ("0-1 d.p.o. | Omicron"),
    "days4Symptoms3[1]": (f"1-7 d.p.o. | symptoms"),
    "days4Symptoms3[2]": ("1-7 d.p.o. | URT symptoms"),
    "days4ReverseSymptoms3[1]": "0-1 d.p.o. | symptoms",
    "days4ReverseSymptoms3[2]": ("0-1 d.p.o. | URT symptoms"),
    "daysImmun2YN[0]": "2-7 d.p.o. | immune naive",
    "daysImmun2YN[1]": ("2-7 d.p.o. | multiply immunized"),
    "days4Immun2YN[0]": "1-7 d.p.o. | immune naive",
    "days4Immun2YN[1]": ("1-7 d.p.o. | multiply immunized"),
    "zMaxSympDays[0]": _LABELS["days"],
    "zVl": _LABELS["vl"],
    "log10Load": _LABELS["vl"],
    "zVlVariant[0]": _LABELS["vl"] + " | Pre-VOC",
    "zVlVariant[1]": _LABELS["vl"] + " | Alpha",
    "zVlVariant[2]": _LABELS["vl"] + " | Delta",
    "zVlVariant[3]": _LABELS["vl"] + " | Omicron",
    "zVlSymptoms3[0]": _LABELS["vl"] + " | no symptoms",
    "zVlSymptoms3[1]": _LABELS["vl"] + " | symptoms",
    "zVlSymptoms3[2]": _LABELS["vl"] + " | URT symptoms",
    "zVlSymptoms3Reverse[0]": _LABELS["vl"] + " | URT symptoms",
    "zVlSymptoms3Reverse[1]": _LABELS["vl"] + " | symptoms",
    "zVlSymptoms3Reverse[2]": _LABELS["vl"] + " | no symptoms",
    "zVl_variant_test[wt_abbott, 1, 1]": (_LABELS["vl"] + " | Pre-VOC, Abbott"),
    "zVl_variant_test[wt_roche, 1, 1]": (_LABELS["vl"] + " | Pre-VOC, Roche"),
    "zVl_variant_test[omi_abbott, 1, 1]": (_LABELS["vl"] + " | Omicron, Abbott"),
    "zVl_variant_test[omi_roche, 1, 1]": (_LABELS["vl"] + " | Omicron, Roche"),
    "roche_abbott_variant[0]": "Roche device | Pre-VOC",
    "roche_abbott_variant[1]": "Roche device | Omicron",
    "zVlImmun2YN[0]": (_LABELS["vl"] + " | immune naive"),
    "zVlImmun2YN[1]": (_LABELS["vl"] + " | multiply immunized"),
    "zVlSymptoms[0]": (_LABELS["vl"] + " | asymptomatic"),
    "zVlSymptoms[1]": (_LABELS["vl"] + " | symptomatic"),
    "zVlBinDaysPostOnset[0]": (_LABELS["vl"] + " | 0-2 d.p.o."),
    "zVlBinDaysPostOnset[1]": (_LABELS["vl"] + " | 2-7 d.p.o."),
    "zIncidenceChange": _LABELS["zIncidenceChange"],
    "zRtRolling4": _LABELS["zRtRolling4"],
    "symptoms[1]": "Symptomatic",
    # "symptomaticImmun2YN[0]": _LABELS["immun2YN"] + " | asymptomatic",
    # "symptomaticImmun2YN[1]": _LABELS["immun2YN"] + " | symptomatic",
    "immun2YNSymptoms[0]": _LABELS["immun2YN"] + " | asymptomatic",
    "immun2YNSymptoms[1]": _LABELS["immun2YN"] + " | symptomatic",
    "immun2YNSymptoms3[0]": _LABELS["immun2YN"] + " | no symptoms",
    "immun2YNSymptoms3[1]": _LABELS["immun2YN"] + " | symptoms",
    "immun2YNSymptoms3[2]": _LABELS["immun2YN"] + " | URT symptoms",
    "immun2YNSymptoms3Reverse[0]": _LABELS["immun2YN"] + " | URT symptoms",
    "immun2YNSymptoms3Reverse[1]": _LABELS["immun2YN"] + " | symptoms",
    "immun2YNSymptoms3Reverse[2]": _LABELS["immun2YN"] + " | no symptoms",
    "symptoms3[1]": "Symptoms",
    "symptoms3[2]": "URT symptoms",
    "symptoms3Reverse[1]": "No URT symptoms",
    "symptoms3Reverse[2]": "No symptoms",
    "illURT": "URT symptoms",
    "illSeverityURT": "Severity of URT symptoms",
    "variantSymptoms[0]": f"{ALPHA_LABEL} | asymptomatic",
    "variantSymptoms[1]": f"{DELTA_LABEL} | asymptomatic",
    "variantSymptoms[2]": f"{OMICRON_LABEL} | asymptomatic",
    "variantSymptoms[3]": f"{ALPHA_LABEL} | symptomatic",
    "variantSymptoms[4]": f"{DELTA_LABEL} | symptomatic",
    "variantSymptoms[5]": f"{OMICRON_LABEL} | symptomatic",
    "recovered[1]": "Prior infection",
    "recoveredSymptoms[0]": "Prior infection | asymptomatic",
    "recoveredSymptoms[1]": "Prior infection | symptomatic",
    # "symptomaticVariant[0]": ""
}

_DAY_DICT = {f"day[{i}]": f"Day {i} p.o." for i in _DAY_LEVELS}
_DAY_DICT2 = {f"1|day[{i}]": f"Day {i} p.o." for i in _DAY_LEVELS2}


PARAM_NAME_MAPPINGS = {**_PARAM_NAME_MAPPINGS, **_DAY_DICT}

_PARAM_NAME_MAPPINGS_TABLE_LATEX = {
    "zAge": "Age",
    "gender[1]": "Male",
    "testDevice[1]": "Roche device",
    "immun2YN[1]": _LABELS_LATEX["immun2YN"],
    "variant[1]": f"{ALPHA_LABEL}",
    "variant[2]": f"{DELTA_LABEL}",
    "variant[3]": f"{OMICRON_LABEL}",
    "VariantCode[1]": f"{DELTA_LABEL}",
    "VariantCode[2]": f"{OMICRON_BA1_LABEL}",
    "VariantCode[3]": f"{OMICRON_BA2_LABEL}",
    "zVl_variant_test[wt_abbott, 1, 1]": (_LABELS["vl"] + " | Pre-VOC, Abbott"),
    "zVl_variant_test[wt_roche, 1, 1]": (_LABELS["vl"] + " | Pre-VOC, Roche"),
    "zVl_variant_test[omi_abbott, 1, 1]": (_LABELS["vl"] + " | Omicron, Abbott"),
    "zVl_variant_test[omi_roche, 1, 1]": (_LABELS["vl"] + " | Omicron, Roche"),
    "roche_abbott_variant[0]": "Roche device | pre-VOC",
    "roche_abbott_variant[1]": "Roche device | Omicron",
    "binDaysPostOnset[1]": (f"2-7 d.p.o."),
    "binDaysPostOnsetReverse[1]": (f"0-2 d.p.o."),
    "binDaysPostOnset4[1]": (f"1-7 d.p.o."),
    "binDaysPostOnsetReverse4[1]": (f"0-1 d.p.o."),
    "daysVariant[0]": (f"2-7 d.p.o. | Pre-VOC"),
    "daysVariant[1]": ("2-7 d.p.o. | Alpha"),
    "daysVariant[2]": ("2-7 d.p.o. | Delta"),
    "daysVariant[3]": ("2-7 d.p.o. | Omicron"),
    "daysReverseVariant[0]": f"0-2 d.p.o. | Pre-VOC",
    "daysReverseVariant[1]": ("0-2 d.p.o. | Alpha"),
    "daysReverseVariant[2]": ("0-2 d.p.o. | Delta"),
    "daysReverseVariant[3]": ("0-2 d.p.o. | Omicron"),
    "days4Variant[0]": (f"1-7 d.p.o. | Pre-VOC"),
    "days4Variant[1]": ("1-7 d.p.o. | Alpha"),
    "days4Variant[2]": ("1-7 d.p.o. | Delta"),
    "days4Variant[3]": ("1-7 d.p.o. | Omicron"),
    "days4ReverseVariant[0]": "0-1 d.p.o. | Pre-VOC",
    "days4ReverseVariant[1]": ("0-1 d.p.o. | Alpha"),
    "days4ReverseVariant[2]": ("0-1 d.p.o. | Delta"),
    "days4ReverseVariant[3]": ("0-1 d.p.o. | Omicron"),
    "days4Symptoms3[1]": (f"1-7 d.p.o. | symptoms"),
    "days4Symptoms3[2]": ("1-7 d.p.o. | URT symptoms"),
    "days4ReverseSymptoms3[1]": "0-1 d.p.o. | symptoms",
    "days4ReverseSymptoms3[2]": ("0-1 d.p.o. | URT symptoms"),
    "daysImmun2YN[0]": ("2-7 d.p.o. | immune naive"),
    "daysImmun2YN[1]": ("2-7 d.p.o. | multiply immunized"),
    "days4Immun2YN[0]": "1-7 d.p.o. | immune naive",
    "days4Immun2YN[1]": ("1-7 d.p.o. | multiply immunized"),
    "zVlSymptoms[0]": (_LABELS["vl"] + " | asymptomatic"),
    "zVlSymptoms[1]": (_LABELS["vl"] + " | symptomatic"),
    "zVlBinDaysPostOnset[0]": (_LABELS["vl"] + " | 0-2 d.p.o."),
    "zVlBinDaysPostOnset[1]": (_LABELS["vl"] + " | 2-7 d.p.o."),
    "zMaxSympDays[0]": _LABELS_LATEX["days"],
    "vlAdj": _LABELS_LATEX["vl"],
    "zVlAdj": _LABELS_LATEX["vl"],
    "day": "Day post symptom onset",
    "ID": "Infection ID",
    "zVl": _LABELS_LATEX["vl"],
    "log10Load": _LABELS_LATEX["vl"],
    "zVlVariant[0]": _LABELS_LATEX["vl"] + " | Pre-VOC",
    "zVlVariant[1]": _LABELS_LATEX["vl"] + " | Alpha",
    "zVlVariant[2]": _LABELS_LATEX["vl"] + " | Delta",
    "zVlVariant[3]": _LABELS_LATEX["vl"] + " | Omicron",
    "zVlImmun2YN[0]": (_LABELS_LATEX["vl"] + " | immune naive"),
    "zVlImmun2YN[1]": (_LABELS_LATEX["vl"] + " | multiply immunized"),
    "zVlSymptoms3[0]": _LABELS["vl"] + " | no symptoms",
    "zVlSymptoms3[1]": _LABELS["vl"] + " | symptoms",
    "zVlSymptoms3[2]": _LABELS["vl"] + " | URT symptoms",
    "zVlSymptoms3Reverse[0]": _LABELS["vl"] + " | URT symptoms",
    "zVlSymptoms3Reverse[1]": _LABELS["vl"] + " | symptoms",
    "zVlSymptoms3Reverse[2]": _LABELS["vl"] + " | no symptoms",
    "zIncidenceChange": _LABELS["zIncidenceChange"],
    "zRtRolling4": _LABELS["zRtRolling4"],
    "symptoms[1]": "Symptomatic",
    # "symptomaticImmun2YN[0]": _LABELS["immun2YN"] + " | asymptomatic",
    # "symptomaticImmun2YN[1]": _LABELS["immun2YN"] + " | symptomatic",
    "immun2YNSymptoms[0]": _LABELS["immun2YN"] + " | asymptomatic",
    "immun2YNSymptoms[1]": _LABELS["immun2YN"] + " | symptomatic",
    "immun2YNSymptoms3[0]": _LABELS["immun2YN"] + " | no symptoms",
    "immun2YNSymptoms3[1]": _LABELS["immun2YN"] + " | symptoms",
    "immun2YNSymptoms3[2]": _LABELS["immun2YN"] + " | URT symptoms",
    "immun2YNSymptoms3Reverse[0]": _LABELS["immun2YN"] + " | URT symptoms",
    "immun2YNSymptoms3Reverse[1]": _LABELS["immun2YN"] + " | symptoms",
    "immun2YNSymptoms3Reverse[2]": _LABELS["immun2YN"] + " | no symptoms",
    "symptoms3[1]": "Symptoms",
    "symptoms3[2]": "URT symptoms",
    "illURT": "URT symptoms",
    "symptoms3Reverse[1]": "symptoms",
    "symptoms3Reverse[2]": "no symptoms",
    "illSeverityURT": "severity of URT symptoms",
    "variantSymptoms[0]": f"{ALPHA_LABEL} | asymptomatic",
    "variantSymptoms[1]": f"{DELTA_LABEL} | asymptomatic",
    "variantSymptoms[2]": f"{OMICRON_LABEL} | asymptomatic",
    "variantSymptoms[3]": f"{ALPHA_LABEL} | symptomatic",
    "variantSymptoms[4]": f"{DELTA_LABEL} | symptomatic",
    "variantSymptoms[5]": f"{OMICRON_LABEL} | symptomatic",
    "recovered[1]": "Prior infection",
    "recovere": "Prior infection",
    "recoveredSymptoms[0]": "Prior infection | asymptomatic",
    "recoveredSymptoms[1]": "Prior infection | symptomatic",
}

PARAM_NAME_MAPPINGS_TABLE_LATEX = {
    **_PARAM_NAME_MAPPINGS_TABLE_LATEX,
    **_DAY_DICT,
    **_DAY_DICT2,
}

# PARAM_ORDER = list(PARAM_NAME_MAPPINGS.keys())


class PlotParams:
    def __init__(self, paramDict):
        self.__dict__.update(paramDict)


def getAbbrvsDict():
    return _ABBRVS


def getLabels():
    return PlotParams(_LABELS)


def getLabelsLatex():
    return PlotParams(_LABELS_LATEX)


def getLegends():
    return PlotParams(_LEGENDS)


def getOrders():
    return PlotParams(_ORDERS)


def getPalettes():
    return PlotParams(_PALETTES)
