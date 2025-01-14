
def get_ratio(ratio_sw, index, scenario):
    #wind/solar ratio
    ratio_sw["R_sw"]=ratio_sw["wind"+scenario]/ratio_sw["solar"+scenario]
    ratio_windsolar = ratio_sw["R_sw"][index]
    return ratio_windsolar


