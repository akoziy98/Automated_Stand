import spd_analyze
from graphics_generation import graphics_generation
import os
#import generate_pdf

import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 22})


params = {}
params['name'] = 'Uliana'
params['date'] = '05.06.2021'
#params['optimal_params'] = [['Criteria', 'PDE, %', 'DCR, Hz', 'SNR', 'AP, %', 'DT, mus', 'TR, ps']]
#params['settings'] = [['Criteria','Vg', 'Vb', 'T'], ['SNR', '3ff', '64.55', '-55']]
params['temp_cart'] = 'temp_cart.png'
params['timeres'] = 'timeres.png'
params['dcr(pde)'] = 'dcr(pde).png'
params['ap(pde)'] = 'ap(pde).png'
params['dir_name'] = '!SPDStandResults//Uliana//21_05_13__18_16'

grph = graphics_generation(params)
grph.get_slices()
grph.plot_dcr()
grph.plot_ap()
grph.choose_optimum()
grph.plot_TR()
grph.plot_heatmap()

grph.generate_pdf()


colors = ['b', 'g', 'r', 'magenta', 'k', 'cyan']

