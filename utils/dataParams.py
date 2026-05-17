from pathlib import Path
from datetime import datetime

import utils
ROOT_DIR = Path(utils.__file__).resolve().parent.parent

# The date the Roche rapid test started being used (as opposed to the Abbott
# test).
DATE_ROCHE = datetime(2021, 12, 18)
DAYS_PER_MONTH = 30.5

SYMPTOMS = [
    "ill",
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

SYMPTOM_DAYS = [f"{symptom}Days" for symptom in SYMPTOMS]


SYMPTOM_DEGREES = ["asymptomatic", "mildly symptomatic", "symptomatic"]

CONTACTS = [
    "contactConf",
    "contactConfProfFreq",
    "contactConfProf",
    "contactConfPriv",
    "contactConfProf1",
    "contactConfProf2",
    "contactConfProf3",
    "contactConfProf4",
    "contactConfPriv1",
    "contactConfPriv2",
    "contactConfPriv3",
    "contactConfPriv4",
    "contactConfPrivHouse",
    "contactSusp",
    "contactSuspProf",
    "contactSuspPriv",
    "contactSuspProfFreq",
    "contactSuspProf1",
    "contactSuspProf2",
    "contactSuspProf3",
    "contactSuspProf4",
    "contactSuspPriv1",
    "contactSuspPriv2",
    "contactSuspPriv3",
    "contactSuspPriv4",
    "contactSuspPrivHouse",
    "contactsValid",
]

CONTACT_PROF_FREQ = {
    0: "Never",
    1: "Rarely",
    2: "Several times a week",
    3: "Every day",
}

CONTACT_PROF_KIND = {
    1: "Close distance (< 1.5 m) multiple times or > 15 min",
    2: "No adequate face mask",
    3: "Being in the same room (>1h)",
    4: "Exposition to aerosols",
}

CONTACT_PRIV_KIND = {
    1: "Close distance (< 1.5 m) multiple times or > 15 min",
    2: "No adequate face mask",
    3: "Being in the same room (>1h)",
    4: "Red message in corona warn app",
}

SYMPTOMS = [
    "ill",
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
    "noSmell": "No smell",
    "noTaste": "No taste",
    "nausea": "Nausea",
    "noAppetite": "No appetite",
    "vomiting": "Vomiting",
    "diarrhea": "Diarrhea",
}


VACC_DATE_COLS = ("vaccDate1", "vaccDate2", "vaccDate3", "vaccDate4")
INFECTION_DATE_COLS = (
    "infectionDate1",
    "infectionDate2",
    "infectionDate3",
    # "infectionDate4",
)

IMMUNIZATION_COLS = (
    VACC_DATE_COLS
    + INFECTION_DATE_COLS
    + (
        "nPrevInfections",
        "datePrevInfection",
        "nInfectionPreOct18",
        "vaccYN",
        # "vaccN",
        "vaccNatLeast",
        "immunYN",
        "immun2YN",
        "immunN",
        "immunNatLeast",
        "vaccTime",
        "vaccTime2",
        "vaccTimeAtLeast",
        # "immunTime",
        # "immunTimeAtLeast",
        "recovered",
        "recoveredTime",
        "recoveredRecently",
        "vaccName1",
        "vaccName2",
        "vaccName3",
        "vaccWhat1",
        "vaccWhat2",
        "vaccWhat3",
    )
    + tuple([f"{col}Orig" for col in VACC_DATE_COLS])
    + tuple([f"{col}Orig" for col in INFECTION_DATE_COLS])
)


# columns containing numeric datatypes
VACC_COLS_NUM = [
    "vaccYN",
    "vaccYN_2",
    "vaccYN2",
    "vaccYN4",
    "vaccYN2_2",
    "vaccYN4_2",
    "vaccN",
    "vaccN_2",
    "vaccTime",
    "vaccTime2",
    "vaccN",
    "vaccStatus1",
    "vaccStatus2",
    "vaccStatus3",
    "vaccStatus4",
    "vaccStatus5",
    "vaccStatus1_2",
    "vaccStatus3_2",
    "vaccStatus5_2",
    "vaccStatus4_2",
]

# columns containing strings
VACC_COLS_STR = [
    "vaccName1",
    "vaccWhat1",
    "vaccName2",
    "vaccWhat2",
    "vaccName3",
    "vaccWhat3",
    "vaccDate1",
    "vaccDate2",
    "vaccDate3",
]

VACC_COLS = VACC_COLS_NUM + VACC_COLS_STR
