import numpy as np, pandas as pd
import numpy as np
from scipy.interpolate import interp1d


def accum_sw(cap_solarx, cap_windx, cur_capsw):
    # acculm
    cap_windmax = np.max(cap_windx)
    cap_solarmax = np.max(cap_solarx)
    capw_now = np.empty(shape=(cur_capsw.shape[0],))#计算现役装机量
    capw_accum = np.empty(shape=(cur_capsw.shape[0],))#计算累积装机量
    for a in range(cur_capsw.shape[0]):
        if a <= 24:
            capw_now[a] = np.sum(cur_capsw[0:a+1])
            capw_accum[a] = np.sum(cur_capsw[0:a+1])
        else:
            capw_now[a] = np.sum(cur_capsw[a-24:a+1])
            capw_accum[a] = np.sum(cur_capsw[0:a+1])

    capw_now1 = capw_now/np.max(capw_now)*cap_windmax
    capw_accum1 = capw_accum/np.max(capw_now)*cap_windmax
    comp_capw = interp1d(capw_now1, capw_accum1, fill_value='extrapolate')
    capw_accumi = comp_capw(cap_windx)

    caps_now1 = capw_now/np.max(capw_now)*cap_solarmax
    caps_accum1 = capw_accum/np.max(capw_now)*cap_solarmax
    comp_caps = interp1d(caps_now1, caps_accum1, fill_value='extrapolate')
    caps_accumi = comp_caps(cap_solarx)

    return capw_accumi, caps_accumi
