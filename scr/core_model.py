
import os, time, matplotlib
from pandas import read_csv
import numpy as np, pandas as pd, gurobipy as gp
import cvxpy as cvx
from country_run2 import supply_pro0331
from supply1 import get_supply
from get_costGHG import get_costGHG0528
from get_ratio import get_ratio
from scipy.interpolate import interp1d
import time
from joblib import Parallel, delayed



main_path = os.getcwd() + '/'
path_input = '../input/'

#parent_path = os.path.dirname(main_path)
#grandparent_path = os.path.dirname(parent_path)
#path_input = grandparent_path + '/input'

solu = ["base_future"]
#solu = ["base_future", "combination", "expanded_connection"]
#solu = ["base_future", "combination", "demand_response", "expanded_connection"]
#solution =  ["base_future", "combination", "demand_response", "expanded_connection", "enhanced_storage", "improved_thermal"]

year = "2056"
scenario = "126"
technologies = ["wind", "solar", "hydro", "nuclear", "thermal", "storage1", "storage2"]
models = ["CanESM5", "BCC-CSM2-MR", "MIROC-ES2H", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0", "NESM3"]
hours = 8760

country_name = pd.read_excel(path_input + 'country_name.xlsx')
supply_mix2020 = pd.read_excel(path_input + 'supply_mix2020.xlsx')
supply_mix = supply_pro0331(supply_mix2020)
demand0 = read_csv(path_input + 'demand'+scenario+'.csv')
cost_GHG = pd.read_excel(path_input + 'cost_GHG2020.xlsx')
ratio_sw = pd.read_excel(path_input + 'ratio_sw.xlsx')
CFw0 = read_csv(path_input + 'CFw.csv')
CFs0 = read_csv(path_input + 'CFs.csv')

run_number = len(supply_mix2020["Status"][supply_mix2020["Status"]==1])


def model_run(
    # capacity
    cap_wind, cap_solar, cap_nuclear, cap_hydro, cap_coal, cap_oil, cap_gas,
    cap_thermal, cap_imp, cap_senergy1, cap_senergy2,
    # generation
    gen_wind, gen_solar, gen_hydro, gen_nuclear, gen_coal, gen_oil, gen_gas,
    gen_thermal, gen_imp, storage_cha1, storage_cha2, storage_dis1, storage_dis2,
    SOC1, SOC2, supply, res_spin1, res_spin2, \
    # objective function
    sys_GHG, sys_cost1, objects, prob, constraints, \
    # cost_inv information
    cost_invwind, cost_invsolar, cost_invhydro, cost_invnuclear, cost_invcoal, cost_invoil,
    cost_invgas, cost_invthermal, cost_invimp, cost_invstor1, cost_invstor2,
    # cost_OM information
    cost_OMwind, cost_OMsolar, cost_OMhydro, cost_OMnuclear, cost_OMcoal, cost_OMoil,
    cost_OMgas, cost_OMthermal, cost_OMimp, cost_chadis1, cost_chadis2,
    # GHG information
    GHG_wind, GHG_solar, GHG_hydro, GHG_nuclear, GHG_coal, GHG_oil, GHG_gas,
    GHG_thermal, GHG_imp, GHG_stor1, GHG_stor2,
    # supply share
    ratio_windsolar, share_sw, wind_ratio, solar_ratio, hydro_ratio, nuclear_ratio,
    coal_ratio, oil_ratio, gas_ratio, duration1, duration2,
    # system information
    i, CFw, CFs, demand, demand_resp, demandm, demandmax, country, index, technologies):


    share_sw.value = i / 100
    start_time = time.time()
    prob.solve(solver="GUROBI", verbose=True)  #####, enforce_dpp=True

    print(i)
    print(country)
    print("@author_Tong's_Lab")
    

    if "wind" in technologies:
        cap_windx = cap_wind.value[0]
        gen_windx = gen_wind.value
    else:
        cap_windx = cap_wind
        gen_windx = gen_wind

    if "solar" in technologies:
        cap_solarx = cap_solar.value[0]
        gen_solarx = gen_solar.value
    else:
        cap_solarx = cap_solar
        gen_solarx = gen_solar

    if "hydro" in technologies:
        cap_hydrox = cap_hydro.value[0]
        gen_hydrox = gen_hydro.value
    else:
        cap_hydrox = cap_hydro
        gen_hydrox = gen_hydro

    if "nuclear" in technologies:
        cap_nuclearx = cap_nuclear.value[0]
        gen_nuclearx = gen_nuclear.value
    else:
        cap_nuclearx = cap_nuclear
        gen_nuclearx = gen_nuclear

    if "coal" in technologies:
        cap_coalx = cap_coal.value[0]
        gen_coalx = gen_coal.value
    else:
        cap_coalx = cap_coal
        gen_coalx = gen_coal

    if "oil" in technologies:
        cap_oilx = cap_oil.value[0]
        gen_oilx = gen_oil.value
    else:
        cap_oilx = cap_oil
        gen_oilx = gen_oil

    if "gas" in technologies:
        cap_gasx = cap_gas.value[0]
        gen_gasx = gen_gas.value
    else:
        cap_gasx = cap_gas
        gen_gasx = gen_gas

    if "thermal" in technologies:
        cap_thermalx = cap_thermal.value[0]
        gen_thermalx = gen_thermal.value
    else:
        cap_thermalx = cap_thermal
        gen_thermalx = gen_thermal

    if "expanded_connection" in solu:
        cap_impx = cap_imp.value[0]
        gen_impx = gen_imp.value
    else:
        cap_impx = cap_imp
        gen_impx = gen_imp

    if "storage1" in technologies:
        cap_senergy1x = cap_senergy1.value[0]
        cap_spower1x = cap_senergy1x / duration1
        storage_cha1x = 0 - storage_cha1.value
        storage_dis1x = storage_dis1.value
        storage_chadis1x = storage_cha1x + storage_dis1x
        SOC1x = SOC1.value
    else:
        cap_spower1x = 0
        cap_senergy1x = 0
        storage_cha1x = np.zeros(hours)
        storage_dis1x = np.zeros(hours)
        storage_chadis1x = np.zeros(hours)
        SOC1x = np.zeros(hours)

    if "storage2" in technologies:
        cap_senergy2x = cap_senergy2.value[0]
        cap_spower2x = cap_senergy2x / duration2
        storage_cha2x = 0 - storage_cha2.value
        storage_dis2x = storage_dis2.value
        storage_chadis2x = storage_cha2x + storage_dis2x
        SOC2x = SOC2.value
    else:
        cap_spower2x = 0
        cap_senergy2x = 0
        storage_cha2x = np.zeros(hours)
        storage_dis2x = np.zeros(hours)
        storage_chadis2x = np.zeros(hours)
        SOC2x = np.zeros(hours)

    if "demand_response" in solu:
        demand_cha = demand_resp.value
    else:
        demand_cha = np.zeros(hours)

    curtail_windx = 0 - (cap_windx * CFw - gen_windx)
    curtail_solarx = 0 - (cap_solarx * CFs - gen_solarx)

    disp_hourly = [gen_windx, gen_solarx, gen_hydrox, gen_nuclearx, gen_coalx, gen_oilx, gen_gasx,
                   gen_thermalx, gen_impx, storage_chadis1x, storage_chadis2x,
                   curtail_windx, curtail_solarx, demand, demand_cha, SOC1x, SOC2x
                   ]

    disp_hourly1 = np.array(disp_hourly).reshape(len(disp_hourly), 8760)
    disp_hourly2 = np.transpose(disp_hourly1)

    supply_sum2 = pd.DataFrame(disp_hourly2,
                               columns=["supply_wind", "supply_solar", "supply_hydro",
                                        "supply_nuclear", "supply_coal", "supply_oil", "supply_gas",
                                        "supply_thermal", "supply_imp","storage_chadis1x", "storage_chadis2x",
                                        "curtail_wind", "curtail_solar", "demand", "demand_cha",
                                        "SOC1", "SOC2"])

    supply_sum2.to_excel(folder_path1 + model + year + '_' + country + '_CFsw' + str(i) + '.xlsx', index=False)

    supply_wind11 = np.sum(gen_windx)
    supply_solar11 = np.sum(gen_solarx)
    supply_hydro11 = np.sum(gen_hydrox)
    supply_nuclear11 = np.sum(gen_nuclearx)
    supply_coal11 = np.sum(gen_coalx)
    supply_oil11 = np.sum(gen_oilx)
    supply_gas11 = np.sum(gen_gasx)
    supply_thermal11 = np.sum(gen_thermalx)
    supply_imp11 = np.sum(gen_impx)
    curtail_windx1 = np.sum(curtail_windx)
    curtail_solarx1 = np.sum(curtail_solarx)

    battarycc11, battarycc22 = np.sum(storage_cha1x), np.sum(storage_cha2x)
    battarydd11, battarydd22 = np.sum(storage_dis1x), np.sum(storage_dis2x)
    GHG_totalx = sys_GHG.value
    cost_totalx = sys_cost1.value[0]

    capw_accumi, caps_accumi = 0, 0


    cost_wind = (cap_windx * cost_invwind + supply_wind11 * cost_OMwind)/np.sum(demand)/1000
    cost_solar = (cap_solarx * cost_invsolar + supply_solar11 * cost_OMsolar)/np.sum(demand)/1000
    cost_hydro = (cap_hydrox * cost_invhydro + supply_hydro11 * cost_OMhydro)/np.sum(demand)/1000
    cost_nuclear = (cap_nuclearx * cost_invnuclear + supply_nuclear11 * cost_OMnuclear)/np.sum(demand)/1000
    cost_coal = (cap_coalx * cost_invcoal + supply_coal11 * cost_OMcoal)/np.sum(demand)/1000
    cost_oil = (cap_oilx * cost_invoil + supply_oil11 * cost_OMoil)/np.sum(demand)/1000
    cost_gas = (cap_gasx * cost_invgas + supply_gas11 * cost_OMgas)/np.sum(demand)/1000
    cost_thermal = (cap_thermalx * cost_invthermal + supply_thermal11 * cost_OMthermal)/np.sum(demand)/1000
    cost_imp = (cap_impx * cost_invimp + supply_imp11 * cost_OMimp)/np.sum(demand)/1000

    cost_storage1 = (cap_senergy1x * cost_invstor1
                    + cost_chadis1 * (abs(battarycc11) + battarydd11))/np.sum(demand)/1000

    cost_storage2 = (cap_senergy2x * cost_invstor2
                    + cost_chadis2 * (abs(battarycc22) + battarydd22))/np.sum(demand)/1000
 
    battarycc = storage_cha1x + storage_cha2x
    battarydd = storage_dis1x + storage_dis2x
    battarycf = np.sum(np.where(battarycc == 0, 0, 1)) 
    battarydf = np.sum(np.where(battarydd == 0, 0, 1)) 


    therm_max = np.max([coal_ratio, oil_ratio, gas_ratio])
    if therm_max == coal_ratio:
        peakf = cap_thermalx * 0.50 * 1.2
    elif therm_max == oil_ratio:
        peakf = cap_thermalx * 0.45 * 1.2
    else:
        peakf = cap_thermalx * 0.40 * 1.2

    peakh = cap_hydrox * 0.3
    peak_fire = np.sum(np.where(gen_thermalx <= peakf, 1, 0)) 
    peak_hydro = np.sum(np.where(gen_hydrox <= peakh, 1, 0))  
    peak_hg = peak_fire + peak_hydro  


    battarycd = battarycc + battarydd
    battarycd1 = battarycd
    for h in range(8760):
        battarycd1[h] = np.where(battarycd1[h] == 0, battarycd1[h - 1], battarycd1[h])
    battarycd2 = np.where(battarycd1 > 0, 1, np.where(battarycd1 < 0, -1, 0))  # >0的赋值为1，<0赋值为-1
    battary_cyc = np.zeros(8760)
    for h in range(8760):
        battary_cyc[h] = np.where(battarycd2[h - 1] < battarycd2[h], 1, 0)
    battary_freq = np.sum(battary_cyc)


    load_net = demand - gen_windx - gen_solarx
    ramp_rate = np.empty(shape=(8760,))  
    for t in range(8760):
        ramp_rate[t] = (load_net[t] - load_net[t - 1]) / load_net[t]
    slope = abs(ramp_rate)
    slope1 = slope[slope > 0]
    resp = np.percentile(slope1, 50)  


    load_net2 = np.array(load_net).reshape(365, 24)
    peak_d = np.empty(shape=(365, 1)) 
    for d in range(365):
        peak_d[d] = (np.max(load_net2[d]) - np.min(load_net2[d])) / np.mean(demand)
    stab = np.percentile(peak_d, 50)  


    ramp_up = ramp_rate[ramp_rate > 0]  
    ramp_90up = np.percentile(ramp_up, 90)  
    extrem_90up_all = ramp_up[ramp_up > ramp_90up]  
    extrem = np.percentile(extrem_90up_all, 50)


    flu = np.std(load_net)/np.mean(load_net)
    flu_ws = np.std(gen_windx + gen_solarx)/np.mean(gen_windx + gen_solarx)
    gen_windsolar = gen_windx + gen_solarx
    ramp_ratews = np.empty(shape=(8760,))
    for t in range(8760):
        ramp_ratews[t] = (gen_windsolar[t] - gen_windsolar[t - 1]) / np.mean(gen_windsolar)
    slope_ws = abs(ramp_ratews)
    slop90 = np.percentile(slope_ws, 90)
    slop90a = slope_ws[slope_ws >= slop90]
    flu_ws_ext = np.sum(slop90a)


    end_time = time.time()
    gap_time = end_time - start_time

    results = [
        cap_windx, cap_solarx, float(capw_accumi), float(caps_accumi), cap_hydrox,
        cap_nuclearx, cap_coalx, cap_oilx, cap_gasx, cap_thermalx, cap_impx,
        cap_spower1x, cap_spower2x, cap_senergy1x, cap_senergy2x,

        GHG_totalx,
        cost_totalx, cost_wind, cost_solar, cost_hydro, cost_nuclear, cost_coal,
        cost_oil, cost_gas, cost_thermal, cost_imp, cost_storage1, cost_storage2,

        supply_wind11, supply_solar11, supply_hydro11, supply_nuclear11, supply_coal11,
        supply_oil11, supply_gas11, supply_thermal11, supply_imp11, curtail_windx1, curtail_solarx1,
        battarycc11, battarydd11, battarycc22, battarydd22,

        battarycf, battarydf, peak_hg, battary_freq, resp, stab, extrem,
        flu, flu_ws, flu_ws_ext, gap_time
        ]

    print(i)
    print(country)
    print(gap_time)
    print("@author_Tong's_Lab")
    return results


def country_run(c):
    #run selected countries
    country = supply_mix.iloc[c, 1]
    index = supply_mix.iloc[c, 0]

    ratio_windsolar = get_ratio(ratio_sw, index, scenario)
    wind_ratio = float(ratio_sw["wind"+scenario][ratio_sw["Country"] == country])
    solar_ratio = float(ratio_sw["solar"+scenario][ratio_sw["Country"] == country])
    i = round((wind_ratio + solar_ratio) * 100)


    CFw, CFs, nuclear_ratio, hydro_ratio, coal_ratio, oil_ratio, gas_ratio, \
    demand, demandm, demandmax = get_supply(CFw0, CFs0, supply_mix, demand0, country)

    if "base_future" in solu:
        name1 = country_name["country"].values
        name2 = country_name["abbr"].values
        country2 = name2[name1 == country][0]

        CFw = CFw00[country2].values
        CFs = CFs00[country2].values


    cost_invwind, cost_invsolar, cost_invhydro, cost_invnuclear, \
    cost_invstor1, cost_invstor2, cost_invimp, cost_invthermal, \
    cost_invcoal, cost_invoil, cost_invgas, \
    cost_OMwind, cost_OMsolar, cost_OMhydro, cost_OMnuclear, \
    cost_chadis1, cost_chadis2, cost_OMimp, \
    cost_OMthermal, cost_OMcoal, cost_OMoil, cost_OMgas, \
    GHG_wind, GHG_solar, GHG_hydro, GHG_nuclear, GHG_stor1, GHG_stor2, \
    GHG_thermal, GHG_imp, GHG_coal, GHG_oil, GHG_gas, \
    duration1, duration2 = get_costGHG0528(cost_GHG, index, coal_ratio, oil_ratio, gas_ratio)

    if "improved_thermal" in solu:
        # improved_thermal measure
        coal_min = 0.30
        oil_min = 0.25
        gas_min = 0.20
        cha_thermal = 1.20
        cost_invthermal, cost_invcoal, cost_invoil, cost_invgas \
        = cost_invthermal * cha_thermal, \
          cost_invcoal * cha_thermal, cost_invoil * cha_thermal, \
          cost_invgas * cha_thermal
    else:
        coal_min = 0.50
        oil_min = 0.45
        gas_min = 0.40

    if "enhanced_storage" in solu:
        eff_dis1 = 0.99
        eff_ch1 = 0.99
        decay_rate1 = 0.00000114  
        eff_dis2 = 0.95
        eff_ch2 = 0.95
        decay_rate2 = 0.0000000114  
    else:
        eff_dis1 = 0.91
        eff_ch1 = 0.91
        decay_rate1 = 0.00000114  
        eff_dis2 = 0.80
        eff_ch2 = 0.80
        decay_rate2 = 0.0000000114  


    constraints = []
    sys_cost = 0
    sys_GHG = 0
    ##------------------------- wind ---------------------------
    if "wind" in technologies:
        cap_wind = cvx.Variable(1)
        gen_wind = cvx.Variable(hours)
        constraints += [
            cap_wind >= 0
            ]
        constraints += [
            gen_wind >= 0,
            gen_wind <= cap_wind * CFw
            ]
        sys_cost += cap_wind * cost_invwind + cvx.sum(gen_wind) * cost_OMwind
        sys_GHG += cvx.sum(gen_wind) * GHG_wind
    else:
        cap_wind = 0
        gen_wind = np.zeros(hours)

    ##------------------------- solar ---------------------------
    if "solar" in technologies:
        cap_solar = cvx.Variable(1)
        gen_solar = cvx.Variable(hours)
        constraints += [
            cap_solar >= 0
            ]
        constraints += [
            gen_solar >= 0,
            gen_solar <= cap_solar * CFs
            ]
        sys_cost += cap_solar * cost_invsolar + cvx.sum(gen_solar) * cost_OMsolar
        sys_GHG += cvx.sum(gen_solar) * GHG_solar
    else:
        cap_solar = 0
        gen_solar = np.zeros(hours)

    ##------------------------- hydro ---------------------------
    if "hydro" in technologies:
        cap_hydro = cvx.Variable(1)
        gen_hydro = cvx.Variable(hours)
        constraints += [
            cap_hydro >= 0
            ]
        constraints += [
            gen_hydro >= 0,
            gen_hydro <= cap_hydro
            ]
        sys_cost += cap_hydro * cost_invhydro + cvx.sum(gen_hydro) * cost_OMhydro
        sys_GHG += cvx.sum(gen_hydro) * GHG_hydro
    else:
        cap_hydro = 0
        gen_hydro = np.zeros(hours)

    ##------------------------- nuclear ---------------------------
    if "nuclear" in technologies:
        cap_nuclear = cvx.Variable(1)
        gen_nuclear = cvx.Variable(hours)
        constraints += [
            cap_nuclear >= 0
            ]
        constraints += [
            gen_nuclear >= 0,
            gen_nuclear == cap_nuclear
            ]
        sys_cost += cap_nuclear * cost_invnuclear + cvx.sum(gen_nuclear) * cost_OMnuclear
        sys_GHG += cvx.sum(gen_nuclear) * GHG_nuclear
    else:
        cap_nuclear = 0
        gen_nuclear = np.zeros(hours)

    ##------------------------- coal ---------------------------
    if "coal" in technologies:
        cap_coal = cvx.Variable(1)
        gen_coal = cvx.Variable(hours)
        constraints += [
            cap_coal >= 0
            ]
        constraints += [
            gen_coal >= 0, 
            gen_coal >= coal_min * cap_coal,
            gen_coal <= cap_coal
            ]
        sys_cost += cap_coal * cost_invcoal + cvx.sum(gen_coal) * cost_OMcoal
        sys_GHG += cvx.sum(gen_coal) * GHG_coal
    else:
        cap_coal = 0
        gen_coal = np.zeros(hours)

    ##------------------------- oil ---------------------------
    if "oil" in technologies:
        cap_oil = cvx.Variable(1)
        gen_oil = cvx.Variable(hours)
        constraints += [
            cap_oil >= 0
            ]
        constraints += [
            gen_oil >= 0,
            gen_oil >= gas_min * cap_oil,
            gen_oil <= cap_oil
            ]
        sys_cost += cap_oil * cost_invoil+ cvx.sum(gen_oil) * cost_OMoil
        sys_GHG += cvx.sum(gen_oil) * GHG_oil
    else:
        cap_oil = 0
        gen_oil = np.zeros(hours)

    ##------------------------- gas ---------------------------
    if "gas" in technologies:
        cap_gas = cvx.Variable(1)
        gen_gas = cvx.Variable(hours)
        constraints += [
            cap_gas >= 0
            ]
        constraints += [
            gen_gas >= 0,
            gen_gas >= 0.40 * cap_gas,
            gen_gas <= cap_gas
            ]
        sys_cost += cap_gas * cost_invgas + cvx.sum(gen_gas) * cost_OMgas
        sys_GHG += cvx.sum(gen_gas) * GHG_gas
    else:
        cap_gas = 0
        gen_gas = np.zeros(hours)

    ##------------------------- thermal ---------------------------
    if "thermal" in technologies:
        cap_thermal = cvx.Variable(1)
        gen_thermal = cvx.Variable(hours)
        constraints += [
            cap_thermal >= 0
        ]
        constraints += [
            gen_thermal >= 0,
            gen_thermal <= cap_thermal
        ]
        therm_max = np.max([coal_ratio, oil_ratio, gas_ratio])
        if therm_max == coal_ratio:
            constraints += [
                gen_thermal >= coal_min * cap_thermal
            ]
            thermal_min = coal_min
        elif therm_max == oil_ratio:
            constraints += [
                gen_thermal >= oil_min * cap_thermal
            ]
            thermal_min = oil_min
        else:
            constraints += [
                gen_thermal >= gas_min * cap_thermal
            ]
            thermal_min = gas_min
        sys_cost += cap_thermal * cost_invthermal + cvx.sum(gen_thermal) * cost_OMthermal
        sys_GHG += cvx.sum(gen_thermal) * GHG_thermal
    else:
        cap_thermal = 0
        gen_thermal = np.zeros(hours)

    ##------------------------- ccus ---------------------------
    #if "ccus" in technologies:
    #    cap_ccus = cvx.Variable(1)
    #    gen_ccus = cvx.Variable(hours)
    #    constraints += [
    #        cap_ccus >= 0
    #        ]
    #    constraints += [
    #        gen_ccus >= 0,
    #        gen_ccus >= thermal_min * cap_ccus,
    #        gen_ccus <= cap_ccus
    #        ]

    #    sys_cost += cap_ccus * cost_invccus + cvx.sum(gen_ccus) * cost_OMccus
    #    sys_GHG += cvx.sum(gen_ccus) * GHG_ccus
    #else:
    #    cap_ccus = 0
    #    gen_ccus = np.zeros(hours)

    ##------------------------- power import ---------------------------
    if "expanded_connection" in solu:
        imp_min = 0 
        imp_max = 0.2
        cha_cost = 2 

        cap_imp = cvx.Variable(1)
        gen_imp = cvx.Variable(hours)
        constraints += [
            cap_imp >= 0
            ]
        constraints += [
            gen_imp >= 0,
            gen_imp >= imp_min * cap_imp,
            gen_imp <= imp_max * demand,
            gen_imp <= cap_imp
            ]
        sys_cost += cap_imp * cost_invimp + cvx.sum(gen_imp) * cost_OMimp*cha_cost
        sys_GHG += cvx.sum(gen_imp) * GHG_imp
    else:
        cap_imp = 0
        gen_imp = np.zeros(hours)


    ##------------------------- short-duration storage ---------------------------
    if "storage1" in technologies:
        cap_senergy1 = cvx.Variable(1)
        SOC1 = cvx.Variable(hours)
        storage_cha1 = cvx.Variable(hours)
        storage_dis1 = cvx.Variable(hours)


        constraints += [
            cap_senergy1 >= 0,
            SOC1 >= 0,
            storage_cha1 >= 0,
            storage_dis1 >= 0,
            storage_cha1 <= cap_senergy1 / duration1,
            storage_dis1 <= cap_senergy1 / duration1,
            SOC1 <= cap_senergy1
        ]


        for t in range(hours):
            constraints += [
                SOC1[(t + 1) % hours] == SOC1[t]
                + storage_cha1[t] * eff_ch1
                - storage_dis1[t] / eff_dis1
                - SOC1[t] * decay_rate1
            ]

        sys_cost += cap_senergy1 * cost_invstor1 \
                 + (cvx.sum(storage_cha1) + cvx.sum(storage_dis1)) * cost_chadis1

        sys_GHG += (cvx.sum(storage_cha1) + cvx.sum(storage_dis1)) * GHG_stor1
    else:
        cap_senergy1 = 0
        SOC1 = np.zeros(hours)
        storage_cha1 = np.zeros(hours)
        storage_dis1 = np.zeros(hours)

    ##------------------------- long-duration storage  ---------------------------
    if "storage2" in technologies:
        cap_senergy2 = cvx.Variable(1)
        SOC2 = cvx.Variable(hours)
        storage_cha2 = cvx.Variable(hours)
        storage_dis2 = cvx.Variable(hours)

        constraints += [
            cap_senergy2 >= 0,
            SOC2 >= 0,
            storage_cha2 >= 0,
            storage_dis2 >= 0,

            storage_cha2 <= cap_senergy2 / duration2,
            storage_dis2 <= cap_senergy2 / duration2,
            SOC2 <= cap_senergy2
        ]

        for t in range(hours):
            constraints += [
                SOC2[(t + 1) % hours] == SOC2[t]
                + storage_cha2[t] * eff_ch2
                - storage_dis2[t] / eff_dis2
                - SOC2[t] * decay_rate2
            ]

        sys_cost += cap_senergy2 * cost_invstor2 \
                    + (cvx.sum(storage_cha2) + cvx.sum(storage_dis2)) * cost_chadis2

        sys_GHG += (cvx.sum(storage_cha2) + cvx.sum(storage_dis2)) * GHG_stor2
    else:
        cap_senergy2 = 0
        SOC2 = np.zeros(hours)
        storage_cha2 = np.zeros(hours)
        storage_dis2 = np.zeros(hours)

    ##------------------------- spinning reserves ---------------------------
    if "res_spin" in technologies:
        res_spin1 = cvx.Variable(hours)
        res_spin2 = cvx.Variable(hours)
        constraints += [
            storage_dis1 + res_spin1 <= cap_senergy1 / duration1,
            storage_dis2 + res_spin2 <= cap_senergy2 / duration2
            ]
        constraints += [
            SOC1 - (storage_dis1 + res_spin1) / eff_dis1
            - SOC1 * decay_rate1 >= 0,
            SOC2 - (storage_dis2 + res_spin2) / eff_dis2
            - SOC2 * decay_rate2 >= 0
            ]
        constraints += [
            cvx.sum(res_spin1 + res_spin2) >= 0.01 * cvx.sum(gen_solar + gen_wind)
            ]
    else:
        res_spin1 = np.zeros(hours)
        res_spin2 = np.zeros(hours)

    ##------------------------- generation share constraints ---------------------------
    share_sw = cvx.Parameter(nonneg=True, value=0.01)
    supply = cvx.Variable(hours)
    constraints += [
        supply == gen_wind
        + gen_solar
        + gen_hydro
        + gen_nuclear
        + gen_coal
        + gen_oil
        + gen_gas
        + gen_thermal
        + gen_imp
        ]

    # wind-solar share constraint
    if "wind" in technologies:
        constraints += [
            cvx.sum(gen_solar + gen_wind) == share_sw * cvx.sum(supply),
        ]

    # wind/solar ratio constraint
    #if "solar" in technologies:
    #    constraints += [
    #            cvx.sum(gen_wind) == ratio_windsolar * cvx.sum(gen_solar)
    #       ]


    if "hydro" in technologies:
        constraints += [
            cvx.sum(gen_hydro) == (1 - share_sw) * hydro_ratio * cvx.sum(supply)
        ]

    if "nuclear" in technologies:
            constraints += [
                cvx.sum(gen_nuclear) == (1 - share_sw) * nuclear_ratio * cvx.sum(supply)
        ]


    if "coal" in technologies:
        constraints += [
            cvx.sum(gen_coal) == (1 - share_sw) * coal_ratio * cvx.sum(supply)
        ]


    if "oil" in technologies:
        constraints += [
            cvx.sum(gen_oil) == (1 - share_sw) * oil_ratio * cvx.sum(supply)
        ]


    if "gas" in technologies:
        constraints += [
            cvx.sum(gen_gas) == (1 - share_sw) * gas_ratio * cvx.sum(supply)
            ]

    # thermal share constraint
    #thermal_ratio = np.sum(coal_ratio + oil_ratio + gas_ratio)
    #if "thermal" in technologies:
        #     constraints += [
        #          cvx.sum(gen_thermal) == (1 - share_sw) * thermal_ratio * cvx.sum(supply)
    #    ]

    #ccus share constraint
    #if "ccus" in technologies:
        #    if "thermal" in technologies:
        #        constraints += [
        #            cvx.sum(gen_ccus) <= cvx.sum(gen_thermal)
        #        ]
        #    else:
        #        constraints += [
        #        cvx.sum(gen_ccus) <= cvx.sum(gen_coal)+cvx.sum(gen_oil)+cvx.sum(gen_gas)
        #     ]



    ##--------------------------demand response----------------------------
    # Peak response coefficient
    share_max = 0.05
    share_flexi = 0.01 
    share_DR = 0.01
    

    if "demand_response" in solu:
        demand_resp = cvx.Variable(hours)  
        demand2 = cvx.Variable(hours)  
        constraints += [
            demand2 == demand + demand_resp,
            cvx.abs(demand_resp) <= share_max * demandmax,
            cvx.abs(cvx.sum(demand_resp)) <= share_flexi * np.sum(demand),
            cvx.sum(cvx.abs(demand_resp)) <= share_DR * np.sum(demand)
        ]

        sys_cost += cvx.sum(cvx.abs(demand_resp)) * cost_OMthermal * 0.1


        ##------------------------- supply/demand balance constraint with demand response measure---------------------------
        constraints += [
            demand2 ==
            gen_wind
            + gen_solar
            + gen_hydro
            + gen_nuclear
            + gen_coal
            + gen_oil
            + gen_gas
            + gen_thermal
            + gen_imp
            + storage_dis1
            + storage_dis2
            - storage_cha1
            - storage_cha2
        ]

    else:
        demand_resp = np.zeros(hours)
        constraints += [
            demand ==
            gen_wind
            + gen_solar
            + gen_hydro
            + gen_nuclear
            + gen_coal
            + gen_oil
            + gen_gas
            + gen_thermal
            + gen_imp
            + storage_dis1
            + storage_dis2
            - storage_cha1
            - storage_cha2
        ]

    ##------------------------- objective function ---------------------------
    sys_cost1 = sys_cost
    objects = cvx.Minimize(sys_cost1)
    prob = cvx.Problem(objects, constraints)
    assert prob.is_dcp(dpp=True)

    start_time = time.time()
    results = []


    results += model_run(
        # capacity
        cap_wind, cap_solar, cap_nuclear, cap_hydro, cap_coal, cap_oil, cap_gas,
        cap_thermal, cap_imp, cap_senergy1, cap_senergy2,
        # generation
        gen_wind, gen_solar, gen_hydro, gen_nuclear, gen_coal, gen_oil, gen_gas,
        gen_thermal, gen_imp, storage_cha1, storage_cha2, storage_dis1, storage_dis2,
        SOC1, SOC2, supply, res_spin1, res_spin2, \
        # objective function
        sys_GHG, sys_cost1, objects, prob, constraints, \
        # cost_inv information
        cost_invwind, cost_invsolar, cost_invhydro, cost_invnuclear, cost_invcoal, cost_invoil,
        cost_invgas, cost_invthermal, cost_invimp, cost_invstor1, cost_invstor2,
        # cost_OM information
        cost_OMwind, cost_OMsolar, cost_OMhydro, cost_OMnuclear, cost_OMcoal, cost_OMoil,
        cost_OMgas, cost_OMthermal, cost_OMimp, cost_chadis1, cost_chadis2,
        # GHG information
        GHG_wind, GHG_solar, GHG_hydro, GHG_nuclear, GHG_coal, GHG_oil, GHG_gas,
        GHG_thermal, GHG_imp, GHG_stor1, GHG_stor2,
        # supply share
        ratio_windsolar, share_sw, wind_ratio, solar_ratio, hydro_ratio, nuclear_ratio,
        coal_ratio, oil_ratio, gas_ratio, duration1, duration2,
        # system information
        i, CFw, CFs, demand, demand_resp, demandm, demandmax, country, index, technologies
        )

    results1 = np.array(results).reshape(1, len(results))

    results2 = pd.DataFrame(results1, columns=[

        "cap_windx", "cap_solarx", "capw_accumi", "caps_accumi", "cap_hydrox",
        "cap_nuclearx", "cap_coalx", "cap_oilx", "cap_gasx", "cap_thermalx", "cap_impx",
        "cap_spower1x", "cap_spower2x", "cap_senergy1x", "cap_senergy2x",

        "GHG_totalx",
        "cost_totalx", "cost_wind", "cost_solar", "cost_hydro", "cost_nuclear", "cost_coal",
        "cost_oil", "cost_gas", "cost_thermal", "cost_imp", "cost_storage1", "cost_storage2",

        "supply_wind11", "supply_solar11", "supply_hydro11", "supply_nuclear11", "supply_coal11",
        "supply_oil11", "supply_gas11", "supply_thermal11", "supply_imp11", "curtail_windx1", "curtail_solarx1",
        "battarycc11", "battarydd11", "battarycc22", "battarydd22",

        "battarycf", "battarydf", "peak_hg", "battary_freq", "resp", "stab", "extrem",
        "flu", "flu_ws", "flu_ws_ext", "gap_time"
        ])

    results2.insert(0, "order", i)
    
    results2.to_excel(folder_path2 + model + year + "_" + country + ".xlsx", index=False)

    end_time = time.time()
    gap_time = end_time - start_time
    print(gap_time)
    print(country)



if "base_future" in solu:
    for model in models:
        folder_path1 = '../output/' + model + "/hourly/"
        folder_path2 = '../output/' + model + "/country/"
        if not os.path.exists(folder_path1):
            os.makedirs(folder_path1)
        if not os.path.exists(folder_path2):
            os.makedirs(folder_path2)
  
        CFw00 = read_csv(path_input + 'CFsw/' + 'CFw'+scenario+'/' + model + '/csv_CF_wind_' + model + '_ssp'+scenario+'_adjusted12h_' + year + '.csv')
        CFs00 = read_csv(path_input +'CFsw/' +  'CFs'+scenario+'/' + model + '/csv_CF_solar2x_1hr_' + model + '_ssp'+scenario+'_adjusted_' + year + '.csv')
        
        CFw00 = CFw00.fillna(0)
        CFs00 = CFs00.fillna(0)

        demand0 = read_csv(path_input + 'demand/' +  'demand'+scenario+'/' + model + '/demand_' + model + year + '.csv')

        for c in range(run_number):
            country_run(c)
        #Parallel(n_jobs = -1)(delayed(country_run)(c) for c in range(run_number))

