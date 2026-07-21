"""
Base (pre-aggregation) sets, as declared in `TERM.TAB`'s `Set` block and confirmed against the
185-commodity / 185-industry / 34-region national/regsupp/DISTGONE data (element lists and file
order cross-checked against `1LAB`, `1MAR`, `R001`, `sec.agg`, `reg.agg`).

`COM` and `IND` share one 185-element label set (ORANI-G convention: an industry's characteristic
output commodity has the same name as the industry, e.g. `"Abattoir"` is both a commodity and an
industry code). `sec.agg`/`reg.agg` later aggregate `COM`/`IND` 185 -> 25 and confirm `REG` is
already a 34 -> 34 identity map -- i.e. the base sets below ARE the full base-data sets; regional
aggregation happens post-Step-2 (RAS pipeline), sectoral aggregation in Step 3.
"""

const COM = [
    "Rice", "Corn", "SweetPotato", "Cassava", "OthTubers", "Peanuts", "Soy", "OthNuts",
    "OthGrains", "Vegetables", "DecorPlants", "Cane", "Tobacco", "PlantFiber", "OthPlantaton",
    "Fruits", "PlantBiophrm", "Rubber", "Coconut", "PalmOil", "Coffee", "Tea", "Cocoa", "Clove",
    "Cashew", "Livestock", "FreshMilk", "PoultryEggs", "OthAnimalPrd", "AgricSvc", "Wood",
    "OthForestPrd", "Fish", "Shrimps", "OthAquatic", "Seaweed", "CoalLignite", "CrudeOil",
    "NaturalGas", "IronOre", "TinOre", "BauxiteOre", "CopperOre", "NickelOre", "OthMetalMine",
    "GoldOre", "SilverOre", "GalianPrd", "NonMetlMinrl", "CoarseSalt", "PetrolNatGas",
    "OthMiningQry", "Abattoir", "MeatProcessg", "DriedFish", "FishProcess", "VegFruitProc",
    "EdibleOils", "Copra", "DairyPrd", "OthFlour", "WheatFlour", "RiceMilling", "BreadBiscuit",
    "Sugar", "Confectionry", "PastaNoodle", "CoffeeProc", "TeaProc", "SoyProcessed", "OthFood",
    "AnimalFeed", "AlcoBeverage", "NonAlcoBev", "Cigarettes", "TobaccoPrd", "Yarn", "Textile",
    "RopesCarpets", "OthTextilPrd", "KnittedPrd", "Apparel", "Tanning", "LeatherPrd", "Footwear",
    "SawMill", "Plywood", "WoodBuildMat", "OthWoodPrdn", "PaperPulp", "Paper", "PaperPrd",
    "PrintedPrd", "OthNMmetlPrd", "OilAndGasRef", "BasChemicals", "Fertilizer", "Plastics",
    "Pesticide", "Paints", "Varnish", "Soaps", "Cosmetics", "OthChemPrd", "Pharmaceutcl",
    "TradMedicine", "Tire", "SmokedRubber", "OthRubberPrd", "PlasticPrd", "Glass", "ClayCermcPrd",
    "Cement", "BasIronSteel", "NonFerrMetal", "FoundryPrd", "FabMetalPrd", "WeaponsAmmo",
    "DomMetalPrd", "OthMetalPrd", "ElectroncPrd", "Instruments", "ElecMotorGen", "ElectEqp",
    "Batteries", "OthElecEqp", "DomElectEqp", "PrimaryMover", "OfficeMachin", "OthMachinEqp",
    "MotorVehicle", "Ships", "RailwayEqp", "Aircraft", "OthTrnsEqp", "Motorcycle", "Furniture",
    "Jewelry", "MusicInstrum", "SportsEquip", "GamesAndToys", "MedicalDvces", "OthIndustPrd",
    "ManMtlRepair", "Electricity", "GasDist", "WaterSupply", "WasteManage", "ResidBuildng",
    "EleGasInfstr", "AgricInfstrc", "RoadsBridges", "OthBuildings", "CarTrading", "CarRepair",
    "OthTrade", "RailTransprt", "LandTransprt", "SeaTransport", "RiverTrnsprt", "AirTransport",
    "TransportSvc", "Postal", "Hotels", "Restaurants", "Publishing", "Broadcasting",
    "Telecommunic", "InformatTech", "FinancialSvc", "InsuranceSvc", "PensionSvc", "OthFinancSvc",
    "RealEstatSvc", "ProfSciTech", "RentalSvc", "GenGovernmet", "GovEducation", "GovHealth",
    "OthGovSvc", "PrvEducation", "PrvHealth", "ArtsEntrtain", "HholdRepairs", "OthSvc",
]
const IND = COM

const SRC = ["dom", "imp"]

const OCC = ["Agric", "OtherManual", "Clerical", "Managerial"]

const MAR = [
    "CarTrading", "OthTrade", "RailTransprt", "LandTransprt", "SeaTransport", "RiverTrnsprt",
    "AirTransport", "TransportSvc", "Postal",
]

const REG = [
    "NAD", "SumUt", "SumBar", "RiauProv", "KepRi", "Jambi", "SumSel", "BaBel", "Bengkulu",
    "Lampung", "DKI", "JaBar", "Banten", "JaTeng", "DIY", "JaTim", "KalBar", "KalTeng", "KalSel",
    "KalTim", "KalUt", "SulUt", "Gorontalo", "SulTeng", "SulaSel", "SulBar", "SulTra", "Bali",
    "NTB", "NTT", "Maluku", "MalUt", "PapuaBar", "PapuaProv",
]
# TERM.TAB's regional structure indexes the SAME 34-region list three times under different
# roles (destination / origin / production-of-margins); kept as distinct aliases so later
# equation code can name a role without implying it is a materially different set.
const DST = REG
const ORG = REG
const PRD = REG

const HOU = ["AllHou"]

@assert length(COM) == 185
@assert length(REG) == 34
@assert length(OCC) == 4
@assert length(MAR) == 9
@assert length(HOU) == 1
@assert length(unique(COM)) == length(COM)
@assert length(unique(REG)) == length(REG)
