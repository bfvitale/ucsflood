"""Models states of America."""

# TODO(vitale): consider https://pypi.org/project/us/

import dataclasses
import enum

@dataclasses.dataclass
class AmericaState:
    """A record representing an American State or Territory."""
    abbrev: str  # two-letter postal abbreviation
    name: str
    fips_number: int  # FIPS code
    region: str

state_dict = {}

def _read_states():
    for record in _data:
        abbrev, name, fips_number, region = record
        state = AmericaState(abbrev, name, fips_number, region)
        state_dict[abbrev] = state

# TODO(vitale): consider https://docs.python.org/3/library/enum.html

AK = 'AK'
AL = 'AL'
AR = 'AR'
AS = 'AS'
AZ = 'AZ'
CA = 'CA'
CO = 'CO'
CT = 'CT'
DC = 'DC'
DE = 'DE'
FL = 'FL'
GA = 'GA'
HI = 'HI'
IA = 'IA'
ID = 'ID'
IL = 'IL'
IN = 'IN'
KS = 'KS'
KY = 'KY'
LA = 'LA'
MA = 'MA'
MD = 'MD'
ME = 'ME'
MI = 'MI'
MN = 'MN'
MO = 'MO'
MS = 'MS'
MT = 'MT'
NC = 'NC'
ND = 'ND'
NE = 'NE'
NH = 'NH'
NJ = 'NJ'
NM = 'NM'
NV = 'NV'
NY = 'NY'
OH = 'OH'
OK = 'OK'
OR = 'OR'
PA = 'PA'
PR = 'PR'
RI = 'RI'
SC = 'SC'
SD = 'SD'
TN = 'TN'
TX = 'TX'
UT = 'UT'
VA = 'VA'
VI = 'VI'
VT = 'VT'
VI = 'VI'
WA = 'WA'
WI = 'WI'
WV = 'WV'
WY = 'WY'

SOUTHEAST = 'southeast'
WEST = 'west'
NORTHEAST = 'northeast'
SOUTHWEST = 'southwest'
MIDWEST = 'midwest'
NONCONTIG = 'noncontig'

# abbrev, name, fips_number, region
_data = [
    (AK, 'Alaska', 2, NONCONTIG),
    (AL, 'Alabama', 1, SOUTHEAST), 
    (AR, 'Arkansas', 5, SOUTHEAST), 
    (AS, 'American Samoa', 60, NONCONTIG), 
    (AZ, 'Arizona', 4, SOUTHWEST), 
    (CA, 'California', 6, WEST), 
    (CO, 'Colorado', 8, WEST), 
    (CT, 'Connecticut', 9, NORTHEAST), 
    (DC, 'District of Columbia', 11, NORTHEAST), 
    (DE, 'Delaware', 10, SOUTHEAST), 
    (FL, 'Florida', 12, SOUTHEAST), 
    (GA, 'Georgia', 13, SOUTHEAST), 
    (HI, 'Hawaii', 15, NONCONTIG), 
    (IA, 'Iowa', 19, MIDWEST), 
    (ID, 'Idaho', 16, WEST), 
    (IL, 'Illinois', 17, MIDWEST), 
    (IN, 'Indiana', 18, MIDWEST), 
    (KS, 'Kansas', 20, MIDWEST), 
    (KY, 'Kentucky', 21, SOUTHEAST), 
    (LA, 'Louisiana', 22, SOUTHEAST), 
    (MA, 'Massachusetts', 25, NORTHEAST), 
    (MD, 'Maryland', 24, SOUTHEAST), 
    (ME, 'Maine', 23, NORTHEAST), 
    (MI, 'Michigan', 26, MIDWEST), 
    (MN, 'Minnesota', 27, MIDWEST), 
    (MO, 'Missouri', 29, MIDWEST), 
    (MS, 'Mississippi', 28, SOUTHEAST), 
    (MT, 'Montana', 30, WEST), 
    (NC, 'North Carolina', 37, SOUTHEAST), 
    (ND, 'North Dakota', 38, MIDWEST), 
    (NE, 'Nebraska', 31, MIDWEST), 
    (NH, 'New Hampshire', 33, NORTHEAST), 
    (NJ, 'New Jersey', 34, NORTHEAST), 
    (NM, 'New Mexico', 35, SOUTHWEST), 
    (NV, 'Nevada', 32, WEST), 
    (NY, 'New York', 36, NORTHEAST), 
    (OH, 'Ohio', 39, MIDWEST), 
    (OK, 'Oklahoma', 40, SOUTHWEST), 
    (OR, 'Oregon', 41, WEST), 
    (PA, 'Pennsylvania', 42, NORTHEAST), 
    (PR, 'Puerto Rico', 72, NONCONTIG), 
    (RI, 'Rhode Island', 44, NORTHEAST), 
    (SC, 'South Carolina', 45, SOUTHEAST), 
    (SD, 'South Dakota', 46, MIDWEST), 
    (TN, 'Tennessee', 47, SOUTHEAST), 
    (TX, 'Texas', 48, SOUTHWEST), 
    (UT, 'Utah', 49, WEST), 
    (VA, 'Virginia', 51, SOUTHEAST), 
    (VI, 'US Virgin Islands', 78, NONCONTIG), 
    (VT, 'Vermont', 50, NORTHEAST), 
    (WA, 'Washington', 53, WEST), 
    (WI, 'Wisconsin', 55, MIDWEST), 
    (WV, 'West Virginia', 54, SOUTHEAST), 
    (WY, 'Wyoming', 56, WEST)
]

_read_states()
