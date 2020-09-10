import numpy as np
import math


# Define NAM

def nam_method(x, P, T, E, area, SpinOff):
    c = []
    q = []
    qofmin = 0.4
    beta = 0.1
    pmm = 10
    carea = 1.0
    qofmin = 0.4
    beta = 0.1
    pmm = 10
    deltat = 24
    carea = 1.0
    # csnow = 0
    # Set initial states
    States = np.array([0, 0, 0.9 * x[1], 0, 0, 0, 0, 0.1])
    ## Set initial state conditions
    snow = States[0]
    u = States[1]
    l = States[2]
    if1 = States[3]
    if2 = States[4]
    of1 = States[5]
    of2 = States[6]
    bf = States[7]

    ## Set parameters
    umax = x[0]
    lmax = x[1]
    cqof = x[2]
    ckif = x[3] / deltat
    ck12 = x[4] / deltat
    tof = x[5]
    tif = x[6]
    tg = x[7]
    ckbf = x[8] / deltat
    qs = 0
    csnow = x[9]
    snowtemp = x[10]
    # Extra.area = 440/(3.6*24)
    lfrac = l / lmax
    fact = area
    # SimRR = list()
    Qsim = np.zeros(len(P))
    # f = open("result.txt", "w")
    spin = 0

    for t in range(len(P)):
        # Set boundary conditions
        prec = P[t]
        evap = E[t]
        temp = T[t]
        # --------------- Snow storage and snow melt --------------------
        if temp < snowtemp:
            snow = snow + prec
        else:
            qs = csnow * temp
            if snow < qs:
                qs = snow
                snow = 0
            else:
                snow = snow - qs
        # ---------------- Evapotranspiration module --------------------
        if temp < 0:
            u1 = u
        else:
            u1 = u + prec + qs
        if u1 > evap:
            eau = evap
            eal = 0
        else:
            eau = u1
            eal = (evap - eau) * lfrac

        u2 = min(u1 - eau, umax)

        if (lfrac > tif):
            qif = (lfrac - tif) / (1 - tif) * u2 / ckif
        else:
            qif = 0

        u3 = u1 - eau - qif

        if u3 > umax:
            pn = u3 - umax
            u = umax
        else:
            pn = 0
            u = u3
        # -------------------- Net precipitation ------------------------
        n = int(pn / pmm) + 1
        pnlst = pn - (n - 1) * pmm
        eal = eal / n

        qofsum = 0
        gsum = 0
        # ---------------------------------------------------------------
        for i in range(1, n + 1, 1):
            pn = pmm
            if i == n:
                pn = pnlst
            # ------------------ Overland flow --------------------------

            if lfrac > tof:
                qof = cqof * (lfrac - tof) / (1 - tof) * pn
            else:
                qof = 0

            qofsum = qofsum + qof

            # -------------------- Recharge -----------------------------

            if lfrac > tg:
                g = (lfrac - tg) / (1 - tg) * (pn - qof)
            else:
                g = 0

            gsum = gsum + g

            # --- Lower zone storage ---

            dl = pn - qof - g
            l = l + dl - eal

            if l > lmax:
                gsum = gsum + l - lmax
                l = lmax

            lfrac = l / lmax

        qof = qofsum
        g = gsum
        eal = n * eal
        # ------------------------ Baseflow -----------------------------

        c = math.exp(-1. / ckbf)
        bf = bf * c + g * carea * (1 - c)

        # ------------------------ Interflow ----------------------------

        c = math.exp(-1. / ck12)
        if1 = if1 * c + qif * (1 - c)
        if2 = if2 * c + if1 * (1 - c)
        # ------ Overland flow routing and overland flow component ------

        of = 0.5 * (of1 + of2) / deltat

        if of > qofmin:
            ckqof = ck12 * (of / qofmin) ** (-beta)  ## Elementwise
        else:
            ckqof = ck12

        c = math.exp(-1. / ckqof)

        # ---------------------------------------------------------------

        of1 = of1 * c + qof * (1 - c)
        of2 = of2 * c + of1 * (1 - c)
        # Update state variables
        States[0] = snow
        States[1] = u
        States[2] = l
        States[3] = if1
        States[4] = if2
        States[5] = of1
        States[6] = of2
        States[7] = bf

        # Update simulated value
        if t >= SpinOff:
            Qsim[t] = (fact * (if2 + of2 + bf))
        # s = str(q)
        # f.write(s + '\n')



    return Qsim
