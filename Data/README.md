# Restricted inputs

No participant data are included. Supply the two original analysis inputs through `run_analysis.py --data-dir`; do not place them in this code package.

| File | Role | Unit of observation |
|---|---|---|
| `full_RTG_data.csv` | Task behaviour, diagnosis, role, and administration mode | One round of a participant's 10-round game |
| `demographics.csv` | Age and gender linked to the task data | Participant |

The source is the combined dataset analysed in the accompanying manuscript. The code expects the original analysis version, including exclusions handled during preprocessing. `data_dictionary.csv` documents the fields used by the code. Identifiers are required for joining records and recognising repeated rounds; no actual identifier values are supplied.

The manuscript states that de-identified behavioural data may be requested from the corresponding author, Ismail Guennouni (ismail.guennouni@iwr.uni-heidelberg.de), subject to data-custodian approval and applicable ethics and information-governance requirements. A request is not a guarantee of access. Access conditions and timing are determined by the relevant data custodians.

No data redistribution licence is granted by this code package. Identifiable data cannot be shared. Fitted models and saved R objects may embed original observations and are also excluded.
