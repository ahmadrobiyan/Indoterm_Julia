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
