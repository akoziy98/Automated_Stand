import spd_analyze
from graphics_generation import graphics_generation
import os
#import generate_pdf

import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 22})

'''
#Params for gated SPD
params = {}
params['name'] = 'Uliana'
params['date'] = '09.06.2021'
params['dir_name'] = '!SPDStandResults//Uliana//21_06_09__19_16'
params['type'] = 'gated'
'''


#Params for freerun SPD
'''
params = {}
params['name'] = 'AGAFIA'
params['date'] = '05.06.2021'
params['dir_name'] = '!SPDStandResults//AGAFIA//21_06_08__18_54'
params['type'] = 'freerun'
'''

'''
#Params for Erazm
params = {}
params['name'] = 'Erazm'
params['date'] = '23.06.2021'
params['dir_name'] = '!SPDStandResults//Erazm//21_06_23__18_43'
params['type'] = 'gated'
params['criterion'] = 'QE'
'''

'''
#Params for Evlampy
params = {}
params['name'] = 'EVLAMPY'
params['date'] = '04.06.2021'
params['dir_name'] = '!SPDStandResults//EVLAMPY//21_06_04__18_38'
params['type'] = 'gated'
params['criterion'] = 'QE'
'''

'''
#Params for Polikarp
params = {}
params['name'] = 'Polikarp'
params['date'] = '28.06.2021'
params['dir_name'] = '!SPDStandResults//Polikarp//21_06_28__15_27'
params['type'] = 'gated'
'''

'''
grph = graphics_generation(params)
grph.get_slices()
grph.plot_dcr()
grph.plot_ap()
grph.choose_optimum()
grph.plot_TR()
grph.plot_heatmap()
'''

'''
#Params for Hariton
params = {}
params['name'] = 'Hariton'
params['date'] = '02.07.2021'
params['dir_name'] = '!SPDStandResults//Hariton//21_07_02__16_49'
params['type'] = 'gated'
'''

'''
#Params for Timofei
params = {}
params['name'] = 'Timofei'
params['date'] = '14.07.2021'
params['dir_name'] = '!SPDStandResults//Timofei//21_07_14__20_26'
params['type'] = 'gated'
'''

'''
#Params for Bronislav
params = {}
params['name'] = 'Bronislav'
params['date'] = '20.07.2021'
params['dir_name'] = '!SPDStandResults//BRONISLAV//21_07_20__21_03'
params['type'] = 'gated'
'''

#Params for Rodion
params = {}
params['name'] = 'Rodion'
params['date'] = '21.07.2021'
params['dir_name'] = '!SPDStandResults//Rodion//21_07_21__18_14'
params['type'] = 'gated'

grph = graphics_generation(params)
grph.generate_pdf()


