from activity.pitzer_v1 import *
from activity.helper_functions import *
from activity.density_script import *

def pitzer_activity_coefficient(df, target_ion, T = 298.15):
    #print("DENSITY : ", str(density(df, T)))
    if "+" in target_ion:
        return cation_activity_coefficient(df, target_ion, T)
    
    elif "-" in target_ion:
        return anion_activity_coefficient(df, target_ion, T)