import numpy as np
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 22})

def Ap_analyze(filename, type = 'gated'):

    with open (filename) as fn:
        for line in fn:
            line = line.strip().split(';')
            cnt = [int(el) for el in line]

    Nmax = len(cnt)
    tst_sum = sum(cnt)
    #print(Nmax)

    if type == 'freerun':
        main_peak = sum(cnt[2498:2503])
        dcr = sum(cnt[4500: 4990]) / (4990 - 4500 + 1)
        # print(dcr)
        first_bin = 0
        iterator = 2510
        while first_bin == 0 and iterator < Nmax - 1:
            if cnt[iterator] != 0:
                first_bin = iterator
            iterator += 1

        deltat = 1.0000000000001598e-07
        #plt.plot(1e6 * deltat * np.arange(2600, 3450), cnt[2600:3450])
        #plt.plot(cnt[4500:4990])
        #plt.show()

        DT = (first_bin - 2499) * deltat
        cnt_new = [el - dcr for el in cnt]
        count_ap = sum(cnt_new[first_bin:3450])
        # This AP is for experimental finding
        AP = count_ap / (main_peak - dcr)
        if DT <= 0:
            DT = 1e-6
        if AP <= 0:
            AP = 1e-4
        if AP == 1:
            AP = 1e-6
        # This AP corresponds to first-order model. See the AP-DT paper for description
        AP = 1 - 1 / (1 + AP)
    elif type == 'gated':
        main_peak = sum(cnt[2498:2503])
        dcr = sum(cnt[4500: 4990]) / (4990 - 4500 + 1)
        #print(dcr)
        first_bin = 0
        iterator = 2510
        while first_bin == 0 and iterator < Nmax  -1:
            if cnt[iterator] != 0:
                first_bin = iterator
            iterator += 1

        deltat = 1.0000000000001598e-08
        DT = (first_bin - 2499) * deltat
        cnt_new = [el - dcr for el in cnt]
        count_ap = sum(cnt_new[first_bin:4900])
        #This AP is for experimental finding
        AP = count_ap / (main_peak - dcr)
        if DT <= 0:
            DT = 1e-6
        if AP <= 0:
            AP = 1e-4
        if AP == 1:
            AP = 1e-6
        #This AP corresponds to first-order model. See the AP-DT paper for description
        AP = 1 - 1 / (1 + AP)

    return DT, AP