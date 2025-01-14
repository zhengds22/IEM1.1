import os
from pandas import read_csv
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib, PIL


def supply_pro0331(supply_mix2020):
    #run selected countries
    supply_run = supply_mix2020[supply_mix2020["Status"]==1]
    supply_run.insert(0,"Index", supply_run.index)
    # thermal share
    #thermal = supply_run.iloc[:, 6:9]
    thermal = supply_run[["Coal2020", "Oil2020", "Gas2020"]]
    # data integration
    generation_mix = pd.DataFrame(thermal.values, columns=("Coal", "Oil", "Gas"))
    generation_mix["Nuclear"] = supply_run["Nuclear2020"].values
    generation_mix["Hydro"] = supply_run["Hydro2020"].values
    generation_mix.insert(0,"Country",supply_run["Country"].values)
    generation_mix.insert(0, "Index", supply_run["Index"].values)

    return generation_mix



