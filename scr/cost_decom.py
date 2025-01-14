import os, time, matplotlib
from pandas import read_csv
from scipy import stats
import numpy as np, pandas as pd, xarray as xr
from joblib import Parallel, delayed

main_path = os.getcwd() + '/'
input_path = '../input/'


hours = 8760

cost179 = pd.read_excel(input_path + "cost_GHG2020.xlsx")
ratiosw = pd.read_excel(input_path + "supply_mix2020.xlsx")

def cost_hourly(supply_hour, cost_entire):
    order = np.array([i + 1 for i in range(8760)]).reshape(8760, 1)

    supply_hour2 = np.concatenate([order, supply_hour], axis=1)
    supply_hour2 = pd.DataFrame(supply_hour2, columns=["order", "supply_hour"])
    supply_hour3 = supply_hour2.sort_values(by='supply_hour', ascending=True)

    cost_hour = np.empty(shape=(8760, 1))
    for d in range(8760):
        if d == 0:
            cost_hour[d, 0] = supply_hour3.values[d, 1] / (8760 - d)
        else:
            cost_hour[d, 0] = cost_hour[d - 1, 0] + (supply_hour3.values[d, 1]
                            - supply_hour3.values[d - 1, 1]) / (8760 - d)
    supply_hour4 = np.concatenate([supply_hour3.values, cost_hour], axis=1)
    supply_hour4 = pd.DataFrame(supply_hour4, columns=["order", "supply_hour", "cost_hour"])

    supply_hour5 = supply_hour4.sort_values(by="order", ascending=True)
    supply_hour6 = supply_hour5.values[:, 2]

    supply_hour7 = (supply_hour6/np.sum(supply_hour6))*cost_entire
    return supply_hour7


def get_invc(cost_cap, discount_rate, life_wind):
    #Capital recovery factor (CRF)
    CRF = discount_rate * ((1 + discount_rate) ** life_wind) / ((1 + discount_rate) ** life_wind - 1)
    cost_invwind = cost_cap * 1000 * CRF
    return cost_invwind

cost_mean = []
cost_90value = []
cost_95value = []
cost_extreme90 = []
cost_normal90 = []
cost_extreme95 = []
cost_normal95 = []
cost_extreme2 = []
cost_normal2 = []
cost_extreme1_5 = []
cost_normal1_5 = []
number_hour2 = []
number_hour1_5 = []


def get_cost(c, results_hourly, results_cap):
    info_hourly = pd.read_excel('../output/' + model + "/hourly/" + results_hourly[c])
    info_cap = pd.read_excel('../output/' + model + "/country/" + results_cap[c])
    start, end = len(model) + 5, -12
    nation = results_hourly[c][start:end]

    discount_rate = 0.07
    life_storage1, life_storage2 = 15, 15
    life_wind, life_solar = 25, 25
    life_hydro, life_nuclear = 40, 40
    life_coal, life_oil, life_gas, life_CCUS = 40, 40, 40, 40


    inv_wind = get_invc(cost179["inv_wind"][cost179["Country"] == nation].values[0], discount_rate, life_wind)
    inv_solar = get_invc(cost179["inv_solar"][cost179["Country"] == nation].values[0], discount_rate, life_solar)
    inv_hydro = get_invc(cost179["inv_hydro"][cost179["Country"] == nation].values[0], discount_rate, life_hydro)
    inv_nuclear = get_invc(cost179["inv_nuclear"][cost179["Country"] == nation].values[0], discount_rate, life_nuclear)
    inv_coal = get_invc(cost179["inv_coal"][cost179["Country"] == nation].values[0], discount_rate, life_coal)
    inv_oil = get_invc(cost179["inv_oil"][cost179["Country"] == nation].values[0], discount_rate, life_oil)
    inv_gas = get_invc(cost179["inv_gas"][cost179["Country"] == nation].values[0], discount_rate, life_gas)
    inv_power1 = get_invc(cost179["inv_power1"][cost179["Country"] == nation].values[0], discount_rate, life_storage1)
    inv_power2 = get_invc(cost179["inv_power2"][cost179["Country"] == nation].values[0], discount_rate, life_storage2)
    inv_stor1 = get_invc(cost179["inv_stor1"][cost179["Country"] == nation].values[0], discount_rate, life_storage1)
    inv_stor2 = get_invc(cost179["inv_stor2"][cost179["Country"] == nation].values[0], discount_rate, life_storage2)


    vari_wind = cost179["vari_wind"][cost179["Country"] == nation].values[0]
    vari_solar = cost179["vari_solar"][cost179["Country"] == nation].values[0]
    vari_hydro = cost179["vari_hydro"][cost179["Country"] == nation].values[0]
    vari_nuclear = cost179["vari_nuclear"][cost179["Country"] == nation].values[0]
    vari_coal = cost179["vari_coal"][cost179["Country"] == nation].values[0]
    vari_oil = cost179["vari_oil"][cost179["Country"] == nation].values[0]
    vari_gas = cost179["vari_gas"][cost179["Country"] == nation].values[0]
    vari_storage1 = cost179["vari_storage1"][cost179["Country"] == nation].values[0]
    vari_storage2 = cost179["vari_storage2"][cost179["Country"] == nation].values[0]

    ratio_coal = ratiosw["Coal2020"][ratiosw["Country"] == nation].values[0]
    ratio_oil = ratiosw["Oil2020"][ratiosw["Country"] == nation].values[0]
    ratio_gas = ratiosw["Gas2020"][ratiosw["Country"] == nation].values[0]
    thermal_max = np.max([ratio_coal, ratio_oil, ratio_gas])

    if ratio_coal == thermal_max:
        inv_thermal = inv_coal
        vari_thermal = vari_coal
    elif ratio_oil == thermal_max:
        inv_thermal = inv_oil
        vari_thermal = vari_oil
    else:
        inv_thermal = inv_gas
        vari_thermal = vari_gas


    cap_wind = info_cap["cap_windx"].values[0]
    cap_solar = info_cap["cap_solarx"].values[0]
    cap_hydro = info_cap["cap_hydrox"].values[0]
    cap_nuclear = info_cap["cap_nuclearx"].values[0]
    cap_thermal = info_cap["cap_thermalx"].values[0]
    cap_storage1 = info_cap["cap_spower1x"].values[0]
    cap_storage2 = info_cap["cap_spower2x"].values[0]
    cap_energy1 = info_cap["cap_senergy1x"].values[0]
    cap_energy2 = info_cap["cap_senergy2x"].values[0]


    supply_wind = np.array(info_hourly["supply_wind"].values).reshape(8760, 1)
    curtail_wind = 0 - np.array(info_hourly["curtail_wind"].values).reshape(8760, 1)
    supply_solar = np.array(info_hourly["supply_solar"].values).reshape(8760, 1)
    curtail_solar = 0 - np.array(info_hourly["curtail_solar"].values).reshape(8760, 1)
    supply_hydro = np.array(info_hourly["supply_hydro"].values).reshape(8760, 1)
    supply_nuclear = np.array(info_hourly["supply_nuclear"].values).reshape(8760, 1)
    supply_thermal = np.array(info_hourly["supply_thermal"].values).reshape(8760, 1)
    demand = np.array(info_hourly["demand"].values).reshape(8760, 1)

    supply_storage10 =  np.array(info_hourly["storage_chadis1x"].values).reshape(hours, 1)
    supply_storage20 =  np.array(info_hourly["storage_chadis2x"].values).reshape(hours, 1)

    supply_storage1 =  np.where(supply_storage10>0, supply_storage10, 0)
    supply_storage2 =  np.where(supply_storage20>0, supply_storage20, 0)

    CFw = (supply_wind + curtail_wind) / cap_wind 
    CFs = (supply_solar + curtail_solar) / cap_solar 


    capw_use = (supply_wind) / CFw
    caps_use = (supply_solar) / CFs

    capw_use[np.isnan(capw_use)] = 0
    caps_use[np.isnan(caps_use)] = 0


    cost_wind1 = cost_hourly(capw_use, cap_wind * inv_wind)/1000
    cost_wind2 = np.array(supply_wind * vari_wind).reshape(8760)/1000
    cost_wind = cost_wind1 + cost_wind2
    cost_wind[np.isnan(cost_wind)] = 0

    cost_solar1 = cost_hourly(caps_use, cap_solar * inv_solar)/1000
    cost_solar2 = np.array(supply_solar * vari_solar).reshape(8760)/1000
    cost_solar = cost_solar1 + cost_solar2
    cost_solar[np.isnan(cost_solar)] = 0

    cost_hydro1 = cost_hourly(supply_hydro, cap_hydro * inv_hydro) / 1000
    cost_hydro2 = np.array(supply_hydro * vari_hydro).reshape(8760) / 1000
    cost_hydro = cost_hydro1 + cost_hydro2
    cost_hydro[np.isnan(cost_hydro)] = 0

    cost_nuclear1 = np.repeat(cap_nuclear * inv_nuclear / 8760, 8760)/1000
    cost_nuclear2 = np.array(supply_nuclear * vari_nuclear).reshape(8760) / 1000
    cost_nuclear = cost_nuclear1 + cost_nuclear2
    cost_nuclear[np.isnan(cost_nuclear)] = 0

    cost_thermal1 = cost_hourly(supply_thermal, cap_thermal * inv_thermal)/1000
    cost_thermal2 = np.array(supply_thermal * vari_thermal).reshape(8760)/1000
    cost_thermal = cost_thermal1 + cost_thermal2
    cost_thermal[np.isnan(cost_thermal)] = 0

    cost_storage_pow1 = cost_hourly(supply_storage1, cap_storage1 * inv_power1)/1000
    cost_storage_ene1 = cost_hourly(supply_storage1, cap_energy1 * inv_stor1)/1000
    cost_storage_pow1a = np.array(supply_storage1 * vari_storage1).reshape(8760)/1000

    cost_storage_pow2 = cost_hourly(supply_storage2, cap_storage2 * inv_power2)/1000
    cost_storage_ene2 = cost_hourly(supply_storage2, cap_energy2 * inv_stor2)/1000
    cost_storage_pow2a = np.array(supply_storage2 * vari_storage2).reshape(8760)/1000

    cost_storage1 = cost_storage_ene1 + cost_storage_pow1a
    cost_storage2 = cost_storage_ene2 + cost_storage_pow2a

    cost_storage1[np.isnan(cost_storage1)] = 0
    cost_storage2[np.isnan(cost_storage2)] = 0

    cost = cost_wind + cost_solar + cost_hydro + cost_nuclear + \
           cost_thermal + cost_storage1 + cost_storage2

    columns_all = ["cost_wind", "cost_solar", "cost_hydro", "cost_nuclear",
                   "cost_thermal", "cost_storage1", "cost_storage2", "cost"]

    cost_all = np.concatenate([cost_wind, cost_solar, cost_hydro, cost_nuclear,
                         cost_thermal, cost_storage1, cost_storage2, cost],
                         axis=0).reshape(len(columns_all), 8760)

    cost_all1 = pd.DataFrame(np.transpose(cost_all), columns=columns_all)


    cost_all2 = pd.concat([info_hourly, cost_all1], axis=1)
    cost_all2.to_excel(folder_path + results_hourly[c])


model_all = ["BCC-CSM2-MR", "CanESM5", "MIROC-ES2H", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0", "NESM3"]


for model in model_all:
    folder_path = '../output/' + model + "/hourly2/"
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)


    results_hourly = sorted(os.listdir('../output/' + model + "/hourly/"))
    results_cap = sorted(os.listdir('../output/' + model + "/country/"))

    start, end = 0, len(results_hourly)
    
    results = Parallel(n_jobs = -1)(delayed(get_cost)(c, results_hourly, results_cap) for c in range(start, end))



