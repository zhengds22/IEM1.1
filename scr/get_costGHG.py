import numpy as np, pandas as pd


def get_costGHG0528(cost_GHG, index, coal_ratio, oil_ratio, gas_ratio):
    #overnight capital cost, discount rate, and life time are required
    #Cfixed_cost = Ccapital_cost * CRF + CO&M
    #catpital recovery factor (CRF), discount rate (DR), life time (LT)
    #CRF = DR * ((1 + DR)^LT) / ((1 + DR) ** LT - 1)

    discount_rate = 0.07
    life_storage1, life_storage2 = 15, 15
    life_wind, life_solar = 25, 25
    life_hydro, life_nuclear = 40, 40
    life_coal, life_oil, life_gas, life_imp = 40, 40, 40, 40

    def get_invc(cost_cap, life_wind):
        #Capital recovery factor (CRF)
        CRF = discount_rate * ((1 + discount_rate) ** life_wind) / ((1 + discount_rate) ** life_wind - 1)
        cost_invwind = cost_cap * 1000 * CRF
        return cost_invwind

    # fixed cost
    cost_invwind = get_invc(cost_GHG["inv_wind"][index], life_wind) + cost_GHG["OMF_wind"][index] * 1000
    cost_invsolar = get_invc(cost_GHG["inv_solar"][index], life_solar) + cost_GHG["OMF_solar"][index] * 1000
    cost_invhydro = get_invc(cost_GHG["inv_hydro"][index], life_hydro) + cost_GHG["OMF_hydro"][index] * 1000
    cost_invnuclear = get_invc(cost_GHG["inv_nuclear"][index], life_nuclear) + cost_GHG["OMF_nuclear"][index] * 1000

    cost_invcoal = get_invc(cost_GHG["inv_coal"][index], life_coal) + cost_GHG["OMF_coal"][index] * 1000
    cost_invgas = get_invc(cost_GHG["inv_gas"][index], life_gas) + cost_GHG["OMF_gas"][index] * 1000
    cost_invoil = get_invc(cost_GHG["inv_oil"][index], life_oil) + cost_GHG["OMF_oil"][index] * 1000
    cost_invimp = get_invc(cost_GHG["inv_imp"][index], life_imp) + cost_GHG["OMF_imp"][index] * 1000

    cost_invstor1 = get_invc(cost_GHG["inv_stor1"][index], life_storage1) + cost_GHG["OMF_stor1"][index] * 1000
    cost_invstor2 = get_invc(cost_GHG["inv_stor2"][index], life_storage2) + cost_GHG["OMF_stor2"][index] * 1000
    # storage system duration
    duration1, duration2 = cost_GHG["duration1"][index], cost_GHG["duration2"][index]

    # variable O&M
    cost_OMwind, cost_OMsolar = cost_GHG["vari_wind"][index], cost_GHG["vari_solar"][index]
    cost_OMhydro, cost_OMnuclear = cost_GHG["vari_hydro"][index], cost_GHG["vari_nuclear"][index]
    cost_chadis1, cost_chadis2 = cost_GHG["vari_storage1"][index], cost_GHG["vari_storage2"][index]
    cost_OMimp = cost_GHG["vari_imp"][index]
    cost_OMcoal, cost_OMoil, cost_OMgas = cost_GHG["vari_coal"][index], cost_GHG["vari_oil"][index], \
                                          cost_GHG["vari_gas"][index]

    # GHG emission factor
    GHG_wind, GHG_solar, GHG_hydro, GHG_nuclear, GHG_stor1, GHG_stor2, GHG_imp = cost_GHG["GHG_wind"][index],\
    cost_GHG["GHG_solar"][index], cost_GHG["GHG_hydro"][index], cost_GHG["GHG_nuclear"][index],\
    cost_GHG["GHG_storage1"][index], cost_GHG["GHG_storage2"][index], cost_GHG["GHG_imp"][index]

    GHG_coal, GHG_oil, GHG_gas = cost_GHG["GHG_coal"][index], cost_GHG["GHG_oil"][index],\
                                     cost_GHG["GHG_gas"][index]

    #thermal power information

    therm_max = np.max([coal_ratio, oil_ratio, gas_ratio])

    if coal_ratio == therm_max:
        cost_invthermal = cost_invcoal
        cost_OMthermal = cost_OMcoal
        GHG_thermal = GHG_coal

    elif oil_ratio == therm_max:
        cost_invthermal = cost_invoil
        cost_OMthermal = cost_OMoil
        GHG_thermal = GHG_oil

    else:
        cost_invthermal = cost_invgas
        cost_OMthermal = cost_OMgas
        GHG_thermal = GHG_gas


    return cost_invwind, cost_invsolar, cost_invhydro, cost_invnuclear, \
           cost_invstor1, cost_invstor2, cost_invimp, cost_invthermal, \
           cost_invcoal, cost_invoil, cost_invgas, \
           cost_OMwind, cost_OMsolar, cost_OMhydro, cost_OMnuclear, \
           cost_chadis1, cost_chadis2, cost_OMimp, \
           cost_OMthermal, cost_OMcoal, cost_OMoil, cost_OMgas, \
           GHG_wind, GHG_solar, GHG_hydro, GHG_nuclear, GHG_stor1, GHG_stor2,\
           GHG_thermal, GHG_imp, GHG_coal, GHG_oil, GHG_gas, \
           duration1, duration2

