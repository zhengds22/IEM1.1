import numpy as np

def get_supply(CFw0, CFs0, supply_mix, demand0, country):
    CFw = CFw0[country].values
    CFs = CFs0[country].values
    #CF_hydro = CF_hydro0[:, index]
    nuclear_ratio = float(supply_mix["Nuclear"][supply_mix["Country"] == country])
    hydro_ratio = float(supply_mix["Hydro"][supply_mix["Country"] == country])
    coal_ratio = float(supply_mix["Coal"][supply_mix["Country"] == country])
    oil_ratio = float(supply_mix["Oil"][supply_mix["Country"] == country])
    gas_ratio = float(supply_mix["Gas"][supply_mix["Country"] == country])

    demand = demand0[country].values
    demandm = np.mean(demand)
    demandmax = np.max(demand)
    return CFw, CFs, nuclear_ratio, hydro_ratio,\
           coal_ratio,oil_ratio,gas_ratio,demand,demandm,demandmax





def supply0(CFw0, CFs0, supply_mix, demand0, c, index):
    CFw = CFw0[:, index]
    CFs = CFs0[:, index]
    #CF_hydro = CF_hydro0[:, index]
    nuclear_ratio = supply_mix["Nuclear"].values[c]
    hydro_ratio = supply_mix["Hydro"].values[c]
    coal_ratio = supply_mix["Coal"].values[c]
    oil_ratio = supply_mix["Oil"].values[c]
    gas_ratio = supply_mix["Gas"].values[c]
    demand = demand0[:, index]
    demandm = np.mean(demand0[:, index])
    demandmax = np.max(demand)
    return CFw, CFs, nuclear_ratio, hydro_ratio,\
           coal_ratio,oil_ratio,gas_ratio,demand,demandm,demandmax

