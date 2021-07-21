def pde(AP = 0, DT = 0, R = 0, R0 = 0):
    import numpy as np
    import math
    from numpy import exp, sqrt, mean, std, diag, inf, log
    import scipy
    # from scipy.stats import norm
    import matplotlib.pyplot as plt
    import afterpulse as aft
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.rcParams.update({'font.size': 22})


    h = 6.626 * 10 ** (-34)
    c = 2.998 * 10 ** 8
    lam = 1550 * 10 ** (-9)
    E_ph = h * c / lam

    f = 100 * 10 ** 3  # pulse repetition frequncy
    T = 1 / f

    #aft_filename = 'AP11--00000.csv'
    tau = DT
    p_ap = AP
    ftau = 1 / tau
    #print('tau = ', tau, 'p-ap = ', p_ap, 'ftau = ', ftau)

    MU = 0.1
    mu = MU * 0.87312663

    #R_err = std(R_list) / sqrt(len(R_list))
    # print('R' + str(i + 1), '=', R, '+/-', R_err, 'kHz')
    R_hlp = ftau - ftau * ((1 - (1 - p_ap ** 2) / (2 * p_ap) * (1 - sqrt(1 - (4 * p_ap * R / ftau) / (1 + p_ap)))) /
                           (1 - (1 - p_ap ** 2) / (2 * p_ap) * (
                                       1 - sqrt(1 - (4 * p_ap * R0 / ftau) / (1 + p_ap))) * (
                                        1 - R / ftau * tau / T)))
    R0_sig = (f * R_hlp / (f - R * tau / T + R_hlp * (tau / T - tau // T)))
    #R_sig = (R_hlp)
    #R_sig_err = (R_err)
    #print('Rsig = ', R_sig, 'R0sig = ', R0_sig)
    pdet = R0_sig / f
    return (- log(1 - pdet) / mu)

