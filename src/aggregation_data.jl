"""
185 → 25 sector aggregation mapping for INDOTERM CGE.

The 25 aggregated sectors:

 1  Crops        — All crop agriculture (incl. AgricSvc)
 2  Livestock    — Animal husbandry & products
 3  Forestry     — Forestry & logging
 4  Fisheries    — Fishing & aquaculture
 5  Coal         — Coal & lignite mining
 6  OilGas       — Crude oil, natural gas extraction
 7  OthMining    — Other mining & quarrying
 8  FoodProc     — Food processing
 9  BevTobText   — Beverages, tobacco & textiles
10  Apparel      — Apparel & leather products
11  WoodPaper    — Wood, paper & publishing
12  Chemicals    — Chemicals, pharma, rubber & plastics
13  NonMetalPrd  — Non‑metallic mineral products (glass, cement)
14  BasMetals    — Basic metals & foundry
15  MetalMach    — Metal products, machinery & equipment
16  TranspEquip  — Transport equipment
17  OthManuf     — Other manufacturing incl. repair
18  Utilities    — Electricity, gas, water & waste
19  Construction — Construction
20  Trade        — Trade & margins (CarTrading, OthTrade)
21  Transport    — Transport services (all modes)
22  InfoComm     — Postal, publishing, telecom & IT
23  HotelsRest   — Hotels & restaurants
24  FinBusSvc    — Financial & business services
25  OthServices  — Government & other services

The 9 margin commodities (MAR set) map as:
  CarTrading (154), OthTrade (156)  → Trade  (20)
  RailTransprt (157) … TransportSvc (162) → Transport (21)
  Postal (163) → InfoComm (22)
"""
const AGGCOM = [
    "Crops",
    "Livestock",
    "Forestry",
    "Fisheries",
    "Coal",
    "OilGas",
    "OthMining",
    "FoodProc",
    "BevTobText",
    "Apparel",
    "WoodPaper",
    "Chemicals",
    "NonMetalPrd",
    "BasMetals",
    "MetalMach",
    "TranspEquip",
    "OthManuf",
    "Utilities",
    "Construction",
    "Trade",
    "Transport",
    "InfoComm",
    "HotelsRest",
    "FinBusSvc",
    "OthServices",
]

# SEC_MAP_185_to_25[i] = aggregated sector index (1-25) for COM[i]
const SEC_MAP_185_to_25 = [
    # 1–25  Crops
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1,
    # 26–29 Livestock
    2, 2, 2, 2,
    # 30    AgricSvc → Crops
    1,
    # 31–32 Forestry
    3, 3,
    # 33–36 Fisheries
    4, 4, 4, 4,
    # 37    Coal
    5,
    # 38–39 OilGas
    6, 6,
    # 40–50 OthMining
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    # 51    PetrolNatGas → OilGas
    6,
    # 52    OthMiningQry → OthMining
    7,
    # 53–72 FoodProc
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    # 73–80 BevTobText
    9, 9, 9, 9, 9, 9, 9, 9,
    # 81–85 Apparel
    10, 10, 10, 10, 10,
    # 86–93 WoodPaper
    11, 11, 11, 11, 11, 11, 11, 11,
    # 94–110 Chemicals
    12, 12, 12, 12, 12, 12, 12, 12,
    12, 12, 12, 12, 12, 12, 12, 12, 12,
    # 111–113 NonMetalPrd
    13, 13, 13,
    # 114–117 BasMetals
    14, 14, 14, 14,
    # 118–130 MetalMach
    15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15,
    # 131–136 TranspEquip
    16, 16, 16, 16, 16, 16,
    # 137–144 OthManuf
    17, 17, 17, 17, 17, 17, 17, 17,
    # 145–148 Utilities
    18, 18, 18, 18,
    # 149–153 Construction
    19, 19, 19, 19, 19,
    # 154–156 Trade
    20, 20, 20,
    # 157–162 Transport
    21, 21, 21, 21, 21, 21,
    # 163    Postal → InfoComm
    22,
    # 164–165 HotelsRest
    23, 23,
    # 166–169 InfoComm
    22, 22, 22, 22,
    # 170–176 FinBusSvc
    24, 24, 24, 24, 24, 24, 24,
    # 177–185 OthServices
    25, 25, 25, 25, 25, 25, 25, 25, 25,
]

"""
34 → 6 region aggregation for the Phase-1 validation model.

Collapses the 34 provinces of `REG` (see prepare_sets.jl) into 6 island groups.
This shrinks the levels model from ~2.4M variables to ~80k — inside the regime
where Ipopt+MUMPS already factors cleanly — so the levels equations can be
validated against a benchmark solve before tackling full-scale factorization.

The 6 island groups:
 1  Sumatra      — NAD, SumUt, SumBar, RiauProv, KepRi, Jambi, SumSel, BaBel, Bengkulu, Lampung  (REG 1–10)
 2  Java         — DKI, JaBar, Banten, JaTeng, DIY, JaTim                                          (REG 11–16)
 3  Kalimantan   — KalBar, KalTeng, KalSel, KalTim, KalUt                                          (REG 17–21)
 4  Sulawesi     — SulUt, Gorontalo, SulTeng, SulaSel, SulBar, SulTra                              (REG 22–27)
 5  BaliNusa     — Bali, NTB, NTT                                                                  (REG 28–30)
 6  MalukuPapua  — Maluku, MalUt, PapuaBar, PapuaProv                                              (REG 31–34)
"""
const REG6 = [
    "Sumatra",
    "Java",
    "Kalimantan",
    "Sulawesi",
    "BaliNusa",
    "MalukuPapua",
]

# REG_MAP_34_to_6[i] = island-group index (1-6) for REG[i]
const REG_MAP_34_to_6 = [
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   # 1–10  Sumatra
    2, 2, 2, 2, 2, 2,               # 11–16 Java
    3, 3, 3, 3, 3,                  # 17–21 Kalimantan
    4, 4, 4, 4, 4, 4,               # 22–27 Sulawesi
    5, 5, 5,                        # 28–30 BaliNusa
    6, 6, 6, 6,                     # 31–34 MalukuPapua
]
