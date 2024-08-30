import math
from math import acos, sqrt, cos, pi 
import numpy as np
import correctionlib
import matplotlib.pyplot as plt
import awkward as ak
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from coffea import hist, processor
import numba
from pprint import pprint
from numba import njit
from numba.typed import List, Dict
import correctionlib, rich
import correctionlib.convert
from hist import Hist
from matplotlib.colors import LogNorm
from ttdilepsolve_class import ttdilepsolve
from ttdilepsolve_class import sign
from ttdilepsolve_class import quad
from ttdilepsolve_class import sqr

filename = "root://xrootd-cms.infn.it///store/mc/RunIISummer20UL18NanoAODv9/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/250000/32D0D1A3-74EB-8146-9F3F-B392AB168FDD.root"
events = NanoEventsFactory.from_root(filename, schemaclass=NanoAODSchema, entry_stop=1).events()
WP=0.2770
bJets = events.Jet[(events.Jet.btagDeepFlavB > WP) & (ak.num(events.Jet)>=4) ]
nbJets =ak.num(bJets)
bJets_pt = events.Jet.pt[(events.Jet.btagDeepFlavB > WP) & (ak.num(events.Jet) >= 4)]
bJets_eta = events.Jet.eta[(events.Jet.btagDeepFlavB > WP) & (ak.num(events.Jet) >= 4)]
bJets_phi = events.Jet.phi[(events.Jet.btagDeepFlavB > WP) & (ak.num(events.Jet) >= 4)]
bJets_mass = events.Jet.mass[(events.Jet.btagDeepFlavB > WP) & (ak.num(events.Jet) >= 4)]

# Convert to Awkward Arrays
bJets_pt = ak.Array(bJets_pt)
bJets_eta = ak.Array(bJets_eta)
bJets_phi = ak.Array(bJets_phi)
bJets_mass = ak.Array(bJets_mass)

# Calculate jet px, py, pz, and E
bJets_px = bJets_pt * np.cos(bJets_phi)
bJets_py = bJets_pt * np.sin(bJets_phi)
bJets_pz = bJets_pt * np.sinh(bJets_eta)
bJets_E = np.sqrt(bJets_px**2 + bJets_py**2 + bJets_pz**2 + bJets_mass**2)


bJets_info = ak.zip({"E": bJets_E, "px": bJets_px, "py": bJets_py, "pz": bJets_pz})

'''
b_E = bJets_info["E"][0]
b_px = bJets_info["px"][0]
b_py = bJets_info["py"][0]
b_pz = bJets_info["pz"][0]

# Concatenating the information for the first jet into a single array
b = ak.concatenate([b_E, b_px, b_py, b_pz], axis=-1)

# Extracting the information for the second jet
bb_E = bJets_info["E"][1]
bb_px = bJets_info["px"][1]
bb_py = bJets_info["py"][1]
bb_pz = bJets_info["pz"][1]

# Concatenating the information for the second jet into a single array
bb = ak.concatenate([bb_E, bb_px, bb_py, bb_pz], axis=-1)
'''


print(nbJets)
print(bJets_px)
print(bJets_info)


def solve(ETmiss, b, bb, lp, lm, mWp, mWm, mt, mtb, pnux, pnuy, pnuz, pnubx, pnuby, pnubz, cd_diff, cubic_single_root_cmplx):
    cubic_single_root_cmplx[0] = 0
    pnuxt = [] 
    
    radic = np.sqrt(b[0]**2 - b[1]**2 - b[2]**2 - b[3]**2)
    mb = np.sqrt(max(0, radic))
    
    radic = np.sqrt(bb[0]**2 - bb[1]**2 - bb[2]**2 - bb[3]**2)
    mbb = np.sqrt(max(0, radic))
    
    radic = np.sqrt(lp[0]**2 - lp[1]**2 - lp[2]**2 - lp[3]**2)
    mlp = np.sqrt(max(0, radic))
    
    radic = np.sqrt(lm[0]**2 - lm[1]**2 - lm[2]**2 - lm[3]**2)
    mlm = np.sqrt(max(0, radic))
    
    a1 = (b[0] + lp[0]) * (mWp**2 - mlp**2) - lp[0] * (mt**2 - mb**2 - mlp**2) + 2 * b[0] * lp[0]**2 - 2 * lp[0] * (b[1]*lp[1] + b[2]*lp[2] + b[3]*lp[3])
    a2 = 2 * (b[0] * lp[1] - lp[0] * b[1])
    a3 = 2 * (b[0] * lp[2] - lp[0] * b[2])
    a4 = 2 * (b[0] * lp[3] - lp[0] * b[3])
    
    c22 = ((mWp**2 - mlp**2) * a4)**2 - 4 * (lp[0]**2 - lp[3]**2) * a1**2 - 4 * (mWp**2 - mlp**2) * lp[3] * a1 * a4
    c21 = -8 * (lp[0]**2 - lp[3]**2) * a1 * a2 + 4 * (mWp**2 - mlp**2) * (lp[1] * a4**2 - lp[3] * a2 * a4) - 8 * lp[1] * lp[3] * a1 * a4
    c11 = -8 * (lp[0]**2 - lp[3]**2) * a1 * a3 + 4 * (mWp**2 - mlp**2) * (lp[2] * a4**2 - lp[3] * a3 * a4) - 8 * lp[2] * lp[3] * a1 * a4
    c20 = -4 * (lp[0]**2 - lp[1]**2) * a4**2 - 4 * (lp[0]**2 - lp[3]**2) * a2**2 - 8 * lp[1] * lp[3] * a2 * a4
    c10 = -8 * (lp[0]**2 - lp[3]**2) * a2 * a3 + 8 * lp[1] * lp[2] * a4**2 - 8 * lp[1] * lp[3] * a3 * a4 - 8 * lp[2] * lp[3] * a2 * a4
    c00 = -4 * (lp[0]**2 - lp[2]**2) * a4**2 - 4 * (lp[0]**2 - lp[3]**2) * a3**2 - 8 * lp[2] * lp[3] * a3 * a4
    
    b1 = (bb[0] + lm[0]) * (mWm**2 - mlm**2) - lm[0] * (mtb**2 - mbb**2 - mlm**2) + 2 * bb[0] * lm[0]**2 - 2 * lm[0] * (bb[1]*lm[1] + bb[2]*lm[2] + bb[3]*lm[3])
    b2 = 2 * (bb[0] * lm[1] - lm[0] * bb[1])
    b3 = 2 * (bb[0] * lm[2] - lm[0] * bb[2])
    b4 = 2 * (bb[0] * lm[3] - lm[0] * bb[3])
    
    dp22 = ((mWm**2 - mlm**2) * b4)**2 - 4 * (lm[0]**2 - lm[3]**2) * b1**2 - 4 * (mWm**2 - mlm**2) * lm[3] * b1 * b4
    dp21 = -8 * (lm[0]**2 - lm[3]**2) * b1 * b2 + 4 * (mWm**2 - mlm**2) * (lm[1] * b4**2 - lm[3] * b2 * b4) - 8 * lm[1] * lm[3] * b1 * b4
    dp11 = -8 * (lm[0]**2 - lm[3]**2) * b1 * b3 + 4 * (mWm**2 - mlm**2) * (lm[2] * b4**2 - lm[3] * b3 * b4) - 8 * lm[2] * lm[3] * b1 * b4
    dp20 = -4 * (lm[0]**2 - lm[1]**2) * b4**2 - 4 * (lm[0]**2 - lm[3]**2) * b2**2 - 8 * lm[1] * lm[3] * b2 * b4
    dp10 = -8 * (lm[0]**2 - lm[3]**2) * b2 * b3 + 8 * lm[1] * lm[2] * b4**2 - 8 * lm[1] * lm[3] * b3 * b4 - 8 * lm[2] * lm[3] * b2 * b4
    dp00 = -4 * (lm[0]**2 - lm[2]**2) * b4**2 - 4 * (lm[0]**2 - lm[3]**2) * b3**2 - 8 * lm[2] * lm[3] * b3 * b4
    
    d22 = dp22 + ETmiss[0]**2 * dp20 + ETmiss[1]**2 * dp00 + ETmiss[0] * ETmiss[1] * dp10 + ETmiss[0] * dp21 + ETmiss[1] * dp11
    d20 = dp20
    d00 = dp00
    d10 = dp10
    d21 = -dp21 - 2 * ETmiss[0] * dp20 - ETmiss[1] * dp10
    d11 = -dp11 - 2 * ETmiss[1] * dp00 - ETmiss[0] * dp10
    
    polx = [0] * 5
    polx[0] = c00**2 * d22**2 + c11 * d22 * (c11 * d00 - c00 * d11) + c00 * c22 * (d11**2 - 2 * d00 * d22) + c22 * d00 * (c22 * d00 - c11 * d11)
    polx[1] = c00 * d21 * (2 * c00 * d22 - c11 * d11) + c00 * d11 * (2 * c22 * d10 + c21 * d11) + c22 * d00 * (2 * c21 * d00 - c11 * d10) - c00 * d22 * (c11 * d10 + c10 * d11) - 2 * c00 * d00 * (c22 * d21 + c21 * d22) - d00 * d11 * (c11 * c21 + c10 * c22) + c11 * d00 * (c11 * d21 + 2 * c10 * d22)
    polx[2] = c00**2 * (2 * d22 * d20 + d21**2) - c00 * d21 * (c11 * d10 + c10 * d11) + c11 * d20 * (c11 * d00 - c00 * d11) + c00 * d10 * (c22 * d10 - c10 * d22) + c00 * d11 * (2 * c21 * d10 + c20 * d11) + d00**2 * (2 * c22 * c20 + c21**2) - 2 * c00 * d00 * (c22 * d20 + c21 * d21 + c20 * d22) + c10 * d00 * (2 * c11 * d21 + c10 * d22) - d00 * d10 * (c11 * c21 + c10 * c22) - d00 * d11 * (c11 * c20 + c10 * c21)
    polx[3] = c00 * d21 * (2 * c00 * d20 - c10 * d10) - c00 * d20 * (c11 * d10 + c10 * d11) + c00 * d10 * (c21 * d10 + 2 * c20 * d11) - 2 * c00 * d00 * (c21 * d20 + c20 * d21) + c10 * d00 * (2 * c11 * d20 + c10 * d21) + c20 * d00 * (2 * c21 * d00 - c10 * d11) - d00 * d10 * (c11 * c20 + c10 * c21)
    polx[4] = c00**2 * d20**2 + c10 * d20 * (c10 * d00 - c00 * d10) + c20 * d10 * (c00 * d10 - c10 * d00) + c20 * d00 * (c20 * d00 - 2 * c00 * d20)
    
    pnuxt = []
    quartic(polx, pnuxt, cubic_single_root_cmplx)
    c0 = c00
    print(pnuxt)
    for i in range(len(pnuxt)):
        print("L130")
        c1 = c10 * pnuxt[i] + c11
        c2 = c20 * pnuxt[i]**2 + c21 * pnuxt[i] + c22
        d1 = d10 * pnuxt[i] + d11
        d2 = d20 * pnuxt[i]**2 + d21 * pnuxt[i] + d22
        denom = c1 * d0 - c0 * d1
        cd_diff.append(denom)
        if abs(denom) < epsilon:
            continue
        lpbz_diff = lp[0] * b[3] - b[0] * lp[3]
        lmbbz_diff = lm[0] * bb[3] - bb[0] * lm[3]
        thispnuy = (c0 * d2 - c2 * d0) / denom
        thispnubx = ETmiss[0] - pnuxt[i]
        thispnuby = ETmiss[1] - thispnuy
        thispnuz = (-a1 - a2 * pnuxt[i] - a3 * thispnuy) / a4
        thispnubz = (-b1 - b2 * thispnubx - b3 * thispnuby) / b4
        pnux.append(pnuxt[i])
        pnuy.append(thispnuy)
        pnuz.append(thispnuz)
        pnubx.append(thispnubx)
        pnuby.append(thispnuby)
        pnubz.append(thispnubz)
        print("Exiting solve function...")

    return


########################
def quartic(polx, pnuxt, cubic_single_root_cmplx):
    polxt = [x for x in polx]
    cubic_single_root_cmplx[0] = 0
    if polx[4] == 0:
        cubic(polx, pnux)
    else:
        for i in range(len(polxt)):
            polxt[i] = polx[i] / polx[-1]  # normalize to coefficient of highest order (=pnux^4)
        
        if polxt[0] == 0:
            pnuxt.append(0)
            for i in range(4):
                polxt[i] = polxt[i + 1]
            cubic(polxt, pnuxt)
        else:
            e = polxt[2] - 3 * polxt[3] * polxt[3] / 8
            f = polxt[1] + polxt[3] * polxt[3] * polxt[3] / 8 - polxt[2] * polxt[3] / 2
            g = polxt[0] - 3 * polxt[3] * polxt[3] * polxt[3] * polxt[3] / 256 + polxt[3] * polxt[3] * polxt[2] / 16 - polxt[3] * polxt[1] / 4
            
            if g == 0:
                pnuxt.append(-polxt[3] / 4)
                polxt2 = [0, 0, 0, 0]
                polxt2[0] = f
                polxt2[1] = e
                polxt2[2] = 0
                polxt2[3] = 1
                cubic(polxt2, pnux)
                for i in range(1, len(pnux)):
                    pnuxt[i] -= polxt[3] / 4
            elif f == 0:
                polxt2 = [0, 0, 0]
                polxt2[0] = g
                polxt2[1] = e
                polxt2[2] = 1
                polxt3 = []
                quadratic(polxt2, polxt3)
                for x in polxt3:
                    if x >= 0:
                        pnuxt.append(sqrt(x - polxt[3] / 4))
                        pnuxt.append(-sqrt(x - polxt[3] / 4))
            else:
                polxt2 = [0, 0, 0, 0]
                polxt2[0] = -f * f
                polxt2[1] = e * e - 4 * g
                polxt2[2] = 2 * e
                polxt2[3] = 1
                polxt3 = []
                cubic(polxt2, polxt3)
                if len(polxt3) == 1 and polxt3[0] < 0:
                    cubic_single_root_cmplx[0] += 1
                    return
                h = 0
                for x in polxt3:
                    if x > 0:
                        h = sqrt(x)
                j = (e + h * h - f / h) / 2
                polxt4 = [0, 0, 0]
                polxt4[0] = j
                polxt4[1] = h
                polxt4[2] = 1
                polxt5 = []
                quadratic(polxt4, polxt5)
                for x in polxt5:
                    pnuxt.append(x - polxt[3] / 4)
                polxt4[0] = g / j
                polxt4[1] = -h
                polxt4[2] = 1
                polxt6 = []
                quadratic(polxt4, polxt6)
                for x in polxt6:
                    pnuxt.append(x - polxt[3] / 4)
    return


######################################


def cubic(polx, pnux):
    if polx[3] == 0:
        quadratic(polx, pnux)
    else:
        q = (polx[2] * polx[2] - 3 * polx[1]) / 9
        r = (2 * polx[2] * polx[2] * polx[2] - 9 * polx[1] * polx[2] + 27 * polx[0]) / 54
        if q == 0:
            pnux.append(-polx[2] / 3)
        elif q * q * q > r * r:
            theta = acos(r / sqrt(q * q * q))
            pnux.append(-2 * sqrt(q) * cos(theta / 3) - polx[2] / 3)
            pnux.append(-2 * sqrt(q) * cos((theta + 2 * pi) / 3) - polx[2] / 3)
            pnux.append(-2 * sqrt(q) * cos((theta + 4 * pi) / 3) - polx[2] / 3)
        else:
            radicant = -r + sqrt(r * r - q * q * q)
            powthrd = pow(fabs(radicant), 1 / 3)
            a = sign(radicant) * powthrd
            b = q / a if a != 0 else 0
            pnux.append(a + b - polx[2] / 3)


##############################



def quadratic(polx, pnux):
    if polx[2] == 0:
        pnux.append(-polx[0] / polx[1])
    elif polx[1] * polx[1] == 4 * polx[2] * polx[0]:
        pnux.append(-0.5 * polx[1] / polx[2])
    elif polx[1] * polx[1] > 4 * polx[2] * polx[0]:
        q = -0.5 * (polx[1] + sign(polx[1]) * sqrt(polx[1] * polx[1] - 4 * polx[2] * polx[0]))
        pnux.append(q / polx[2])
        pnux.append(polx[0] / q)
##################################



def algebraic_pz(b, lp, mWp, mt, mb, mlp, pnux, pnuy, pnuz):
    mblp = sqrt((b[0] + lp[0]) ** 2 - (b[1] + lp[1]) ** 2 - (b[2] + lp[2]) ** 2 - (b[3] + lp[3]) ** 2)

    a1 = [(sqr(mWp) - sqr(mlp)) * lp[3] * 0.5 / (sqr(lp[0]) - sqr(lp[3])),
          lp[1] * lp[3] / (sqr(lp[0]) - sqr(lp[3])),
          lp[2] * lp[3] / (sqr(lp[0]) - sqr(lp[3]))]

    a2 = [(quad(mWp) + quad(mlp) - 2 * sqr(mWp) * sqr(mlp)) * 0.25 / (sqr(lp[0]) - sqr(lp[3])),
          -(sqr(lp[0]) - sqr(lp[1])) / (sqr(lp[0]) - sqr(lp[3])),
          -(sqr(lp[0]) - sqr(lp[2])) / (sqr(lp[0]) - sqr(lp[3])),
          (sqr(mWp) - sqr(mlp)) * lp[1] / (sqr(lp[0]) - sqr(lp[3])),
          (sqr(mWp) - sqr(mlp)) * lp[2] / (sqr(lp[0]) - sqr(lp[3])),
          2 * lp[1] * lp[2] / (sqr(lp[0]) - sqr(lp[3]))]

    b1 = [(sqr(mt) - sqr(mblp)) * (b[3] + lp[3]) * 0.5 / (sqr(b[0] + lp[0]) - sqr(b[3] + lp[3])),
          (b[1] + lp[1]) * (b[3] + lp[3]) / (sqr(b[0] + lp[0]) - sqr(b[3] + lp[3])),
          (b[2] + lp[2]) * (b[3] + lp[3]) / (sqr(b[0] + lp[0]) - sqr(b[3] + lp[3]))]

    b2 = [
    (quad(mt) + quad(mblp) - 2 * sqr(mt) * sqr(mblp)) * 0.25 / (sqr(b[0] + lp[0]) - sqr(b[3] + lp[3])),
    -(sqr(b[0] + lp[0]) - sqr(b[1] + lp[1])) / (sqr(b[0] + lp[0]) - sqr(b[3] + lp[3])),
    -(sqr(b[0] + lp[0]) - sqr(b[2] + lp[2])) / (sqr(b[0] + lp[0]) - sqr(b[3] + lp[3])),
    (sqr(mt) - sqr(mblp)) * (b[1] + lp[1]) / (sqr(b[0] + lp[0]) - sqr(b[3] + lp[3])),
    (sqr(mt) - sqr(mblp)) * (b[2] + lp[2]) / (sqr(b[0] + lp[0]) - sqr(b[3] + lp[3])),
    2 * (b[1] + lp[1]) * (b[2] + lp[2]) / (sqr(b[0] + lp[0]) - sqr(b[3] + lp[3]))
    ]
    a1val = evalterm(a1, pnux, pnuy)
    a2val = evalterm(a2, pnux, pnuy)
    b1val = evalterm(b1, pnux, pnuy)
    b2val = evalterm(b2, pnux, pnuy)

    pnuz_a = []
    pnuz_b = []

    radicant = a1val * a1val + a2val
    if radicant >= 0:
        pnuz_a.append(a1val + sqrt(radicant))
        pnuz_a.append(a1val - sqrt(radicant))
    elif fabs(radicant) < epsilon:
        pnuz_a.append(a1val)

    radicant = b1val * b1val + b2val
    if radicant >= 0:
        pnuz_b.append(b1val + sqrt(radicant))
        pnuz_b.append(b1val - sqrt(radicant))
    elif fabs(radicant) < epsilon:
        pnuz_b.append(b1val)

    if len(pnuz_a) == 0 or len(pnuz_b) == 0:
        return -1  # error

    pnuzchi_min = fabs(pnuz_a[0] - pnuz_b[0])
    a_min_ind = 0
    b_min_ind = 0
    for j in range(len(pnuz_a)):
        for k in range(len(pnuz_b)):
            pnuzchi = fabs(pnuz_a[j] - pnuz_b[k])
            if pnuzchi < pnuzchi_min:
                pnuzchi_min = pnuzchi
                a_min_ind = j
                b_min_ind = k

    if pnuzchi_min < sqrt(epsilon):
        pnuz[0] = 0.5 * (pnuz_a[a_min_ind] + pnuz_b[b_min_ind])
        return 0  # success
    else:
        return -1  # error



################################

def evalterm1(a1, pnux, pnuy):
    return a1[0] + a1[1] * pnux + a1[2] * pnuy

###############################3


def evalterm2(a2, pnux, pnuy):
    return a2[0] + a2[1] * pnux + a2[2] * pnuy + a2[3] * pnux * pnux + a2[4] * pnux * pnuy + a2[5] * pnuy * pnuy
    
    
    ############################

ETmiss = [20, 30]  # Missing transverse energy
b = [100, 10, 20, 30]  # Values for the b vector
bb = [120, 15, 25, 35]  # Values for the bb vector
lp = [80, 5, 15, 25]  # Values for the lp vector
lm = [90, 6, 16, 26]  # Values for the lm vector
mWp = 80  # Mass of W+ boson
mWm = 80  # Mass of W- boson
mt = 170  # Top quark mass
mtb = 170  # Anti-top quark mass
pnux = []  # Empty list for pnux values
pnuy = []  # Empty list for pnuy values
pnuz = []  # Empty list for pnuz values
pnubx = []  # Empty list for pnubx values
pnuby = []  # Empty list for pnuby values
pnubz = []  # Empty list for pnubz values
cd_diff = []  # Empty list for cd_diff values
cubic_single_root_cmplx = np.zeros(1)  # Initialize cubic_single_root_cmplx array


solve(ETmiss, b, bb, lp, lm, mWp, mWm, mt, mtb, pnux, pnuy, pnuz, pnubx, pnuby, pnubz, cd_diff, cubic_single_root_cmplx)


print("pnux:", pnux)
print("pnuy:", pnuy)
print("pnuz:", pnuz)
print("pnubx:", pnubx)
print("pnuby:", pnuby)
print("pnubz:", pnubz)
