from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pprint import pprint
import numba
from numba import njit
from numba.typed import List, Dict
import correctionlib, rich
import correctionlib.convert
from hist import Hist


filename = "root://xrootd-cms.infn.it///store/mc/RunIISummer20UL18NanoAODv9/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/250000/32D0D1A3-74EB-8146-9F3F-B392AB168FDD.root"
events = NanoEventsFactory.from_root(filename, schemaclass=NanoAODSchema, entry_stop=100000).events()

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


print(nbJets)
print(bJets_px)
print(bJets_info)

class TLorentzVector:
    def __init__(self, px=0.0, py=0.0, pz=0.0, E=0.0):
        self.px = px
        self.py = py
        self.pz = pz
        self.E = E
    
    def set(self, px, py, pz, E):
        self.px = px
        self.py = py
        self.pz = pz
        self.E = E
    
    def __repr__(self):
        return f"TLorentzVector(px={self.px}, py={self.py}, pz={self.pz}, E={self.E})"

class Particle:
    
    def __init__(self, mass=0.0):
        self.px = self.py = self.pz = self.E = self.m = mass
        self.pAbs = self.pt = 0.0
        self.p = [0.0, 0.0, 0.0, 0.0]

    def sintheta(self):
        R = math.sqrt(self.px**2 + self.py**2 + self.pz**2)
        theta = math.acos(self.py / R)
        return theta

    def p4(self, momx, momy, momz, energy):
        # components of 4-momenta
        self.px = self.p[0] = momx
        self.py = self.p[1] = momy
        self.pz = self.p[2] = momz
        self.E  = self.p[3] = energy
        # transverse momentum and the magnitude of the space momentum
        self.pt = math.sqrt(momx**2 + momy**2)
        self.pAbs = math.sqrt(momx**2 + momy**2 + momz**2)

    def setMass(self, mass):
        self.m = mass

    def print(self):
        print("\n({},\t{},\t{},\t{})  {}".format(self.p[0], self.p[1], self.p[2], self.p[3], self.sintheta()))

    def boost(self, parent):
        # beta and gamma values
        betax = (-1) * parent.px / parent.E
        betay = (-1) * parent.py / parent.E
        betaz = (-1) * parent.pz / parent.E
        beta2 = betax**2 + betay**2 + betaz**2
        gamma = 1.0 / math.sqrt(1.0 - beta2)
        dot = betax * self.px + betay * self.py + betaz * self.pz
        prod = gamma * (gamma * dot / (1.0 + gamma) + self.E)
        
        pX = self.px + betax * prod
        pY = self.py + betay * prod
        pZ = self.pz + betaz * prod
        e = gamma * (self.E + dot)
        
        self.p4(pX, pY, pZ, e)

def combination_indices(comb):
    if comb in {1, 7}:
        a, b = 3, 4
    elif comb in {2, 8}:
        a, b = 2, 4
    elif comb in {3, 9}:
        a, b = 2, 3
    elif comb in {4, 10}:
        a, b = 1, 4
    elif comb in {5, 11}:
        a, b = 1, 3
    elif comb in {6, 12}:
        a, b = 1, 2
    a -= 1  # Start from zero
    b -= 1  # Start from zero
    return a, b

def true_Higgs_indices(jetIsFromHiggs):
    size = len(jetIsFromHiggs)
    a = b = comb = -1
    for i in range(size):
        if jetIsFromHiggs[i] == 1:
            a = i
        elif jetIsFromHiggs[i] == -1:
            b = i

    if (a == 2 and b == 3) or (a == 3 and b == 2):
        comb = 0
    elif (a == 1 and b == 3) or (a == 3 and b == 1):
        comb = 1
    elif (a == 1 and b == 2) or (a == 2 and b == 1):
        comb = 2
    elif (a == 0 and b == 3) or (a == 3 and b == 0):
        comb = 3
    elif (a == 0 and b == 2) or (a == 2 and b == 0):
        comb = 4
    elif (a == 0 and b == 1) or (a == 1 and b == 0):
        comb = 5
    return a, b, comb

def main():
    fdilepsol = TTDilepSolve()
    # Variable declarations
    max_weight = 0.0
    max_higgs_mass = 0.0
    max_ttbar_mass = 0.0
    w_mass = 0
    t_mass = 0
    counter = 1
    # Arrays for 4-momentum components
    lp = [0.0] * 4
    lm = [0.0] * 4
    b = [0.0] * 4
    bb = [0.0] * 4
    b1 = [0.0] * 4
    b2 = [0.0] * 4
    b3 = [0.0] * 4
    b4 = [0.0] * 4
    bh1 = [0.0] * 4
    bh2 = [0.0] * 4

    nu = [0.0] * 4
    nub = [0.0] * 4
    l = [0.0] * 4
    ll = [0.0] * 4

    numb_m = 0
    numb = 0

    ETmiss = [0.0] * 2

    # Variables for jet properties
    jHiggs1_Pt = 0.0
    jHiggs1_Eta = 0.0
    jHiggs1_Phi = 0.0
    jHiggs1_E = 0.0

    jHiggs2_Pt = 0.0
    jHiggs2_Eta = 0.0
    jHiggs2_Phi = 0.0
    jHiggs2_E = 0.0

    max_weight_perCombination = [0.0] * 12
    sum_of_weights = 0.0
    number_of_solutions = 0

    # Variables used for event reading
    max_w = 0.0
    massbquark = 4.8
    massmuon = 0.106  # 0.00051; //0.10566

    a = 0
    b_jHiggs1 = 0
    b_jHiggs2 = 0
    numb_m_t = 0
    b_jHiggs1_t = 0
    b_jHiggs2_t = 0
    
    j0 = TLorentzVector()
    j1 = TLorentzVector()
    j2 = TLorentzVector()
    j3 = TLorentzVector()
    lep1 = TLorentzVector()
    lep2 = TLorentzVector()
    higgs1 = TLorentzVector()
    higgs2 = TLorentzVector()
    top = TLorentzVector()
    topbar = TLorentzVector()

    N_btags_Medium = 0
    mbb = 0.0
    
    #take input file
    filename = "root://xrootd-cms.infn.it///store/mc/RunIISummer20UL18NanoAODv9/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/250000/32D0D1A3-74EB-8146-9F3F-B392AB168FDD.root"
    events = NanoEventsFactory.from_root(filename, schemaclass=NanoAODSchema, entry_stop=100000).events()
    
    eventNumber = events.event
    Lumi = events.luminosityBlock
    Run = events.run
    electron_pt = events.Electron.pt
    electron_eta = events.Electron.eta
    electron_phi = events.Electron.phi
    electron_m = events.Electron.mass
    electron_Q = events.Electron.charge

    muon_pt = events.Muon.pt
    muon_eta = events.Muon.eta
    muon_phi = events.Muon.phi
    muon_m = events.Muon.energy
    muon_Q = events.Muon.charge
    lepPt = ak.concatenate([electron_pt, muon_pt], axis=1)
    lepEta = ak.concatenate([electron_eta, muon_eta], axis=1)
    lepPhi = ak.concatenate([electron_phi, muon_phi], axis=1)
    lepMass = ak.concatenate([electron_mass, muon_mass], axis=1)
    lepQ = ak.concatenate([electron_Q, muon_Q], axis=1)
    
    # Calculate energy
    def calculate_energy(pt, eta, phi, mass):
        px = pt * np.cos(phi)
        py = pt * np.sin(phi)
        pz = pt * np.sinh(eta)
        energy = np.sqrt(px**2 + py**2 + pz**2 + mass**2)
        return energy

    electron_E = calculate_energy(electron_pt, electron_eta, electron_phi, electron_m)
    muon_E = calculate_energy(muon_pt, muon_eta, muon_phi, muon_m)

# Concatenate Electron and Muon attributes
    lepPt = ak.concatenate([electron_pt, muon_pt], axis=1)
    lepEta = ak.concatenate([electron_eta, muon_eta], axis=1)
    lepPhi = ak.concatenate([electron_phi, muon_phi], axis=1)
    lepMass = ak.concatenate([electron_m, muon_m], axis=1)
    lepE = ak.concatenate([electron_E, muon_E], axis=1)
    lepQ = ak.concatenate([electron_Q, muon_Q], axis=1)
    njets = events.njets
    jetPt = events.jetPt
    jetEta = events.jetEta
    jetE = events.jetE
    jetPhi = events.jetPhi
    jetBTag_DeepJet = events.jetBTag_DeepJet
    jetIsFromHiggs = events.jetIsFromHiggs
    met = events.met
    metPhi = events.metPhi

    for event in events:
    # Check the number of jets
        if len(event.Jet) < 4:
           continue

    # Count medium b-tags
        N_btags_Medium = ak.sum(event.Jet.btagDeepFlavB > WP)
    
    # Uncomment if you want to skip events with less than 3 medium b-tags
    # if N_btags_Medium < 3:
    #     continue

        print(f"Processing event: {counter}")
        counter += 1
  
    # Missing transverse energy
        ETmiss = [event.MET.pt, event.MET.phi]

    # Set up Lorentz vectors for jets and leptons
        jets = ak.zip({
            "pt": event.Jet.pt,
            "eta": event.Jet.eta,
            "phi": event.Jet.phi,
            "mass": event.Jet.mass
         }, with_name="PtEtaPhiMLorentzVector")
    
        leptons = ak.zip({
            "pt": event.Jet.pt,
            "eta": event.Jet.eta,
            "phi": event.Jet.phi,
            "mass": event.Jet.mass
         }, with_name="PtEtaPhiMLorentzVector")
           
    }, with_name="PtEtaPhiMLorentzVector")

    lep1 = electrons[0] if len(electrons) > 0 else None
    lep2 = muons[0] if len(muons) > 0 else None

    if lep1 is None or lep2 is None:
        continue

    # Setup b-quark and lepton arrays
    b1 = [jets[0].energy, jets[0].px, jets[0].py, jets[0].pz]
    b2 = [jets[1].energy, jets[1].px, jets[1].py, jets[1].pz]
    b3 = [jets[2].energy, jets[2].px, jets[2].py, jets[2].pz]
    b4 = [jets[3].energy, jets[3].px, jets[3].py, jets[3].pz]

    l = [lep1.energy, lep1.px, lep1.py, lep1.pz]
    ll = [lep2.energy, lep2.px, lep2.py, lep2.pz]

    # Initialize variables for finding the maximum weight solution per event
    max_weight = -999.
    max_higgs_mass = -999.
    max_ttbar_mass = -999.
    w_mass = -1
    t_mass = -1
    numb_m = -1
    max_weight_perCombination = np.full(12, -999., dtype=float)
    sum_of_weights = 0.
    number_of_solutions = 0


    # Assign leptons based on charge
        if data["lepQ"][jentry][0] > 0:
           for index in range (4):
               lp=l
               lm=ll
        if data["lepQ"][jentry][0] < 0:
           for index in range (4):
               lp=ll
               lm=l  

    # Combinatorics
    for numb in range(1, 13):
        if numb == 1:
            for index in range (4):
            	b=b1
            	bb=b2
            	bh1=b3
            	bh2=b4
        if numb == 2:
            for index in range (4):
            	b=b1
            	bb=b3
            	bh1=b2
            	bh2=b4
        if numb == 3:
            for index in range (4):
            	b=b1
            	bb=b4
            	bh1=b3
            	bh2=b2
        if numb == 4:
            for index in range (4):
            	b=b2
            	bb=b3
            	bh1=b1
            	bh2=b4
        if numb == 5:
            for index in range (4):
            	b=b2
            	bb=b4
            	bh1=b1
            	bh2=b3
        if numb == 6:
            for index in range (4):
            	b=b3
            	bb=b4
            	bh1=b1
            	bh2=b2
        if numb == 7:
            for index in range (4):
            	bb=b1
            	b=b2
            	bh1=b3
            	bh2=b4
        if numb == 8:
            for index in range (4):
            	b=b1
            	bb=b2
            	bh1=b3
            	bh2=b4
        if numb == 9:
            for index in range (4):
            	b=b1
            	bb=b2
            	bh1=b3
            	bh2=b4
        if numb == 10:
            for index in range (4):
            	b=b1
            	bb=b2
            	bh1=b3
            	bh2=b4
        if numb == 11:
            for index in range (4):
            	b=b1
            	bb=b2
            	bh1=b3
            	bh2=b4
        if numb == 12:
            for index in range (4):
            	b=b1
            	bb=b2
            	bh1=b3
            	bh2=b4
        if numb == 13:
            for index in range (4):
            	b=b1
            	bb=b2
            	bh1=b3
            	bh2=b4
        # Additional logic to calculate weights, masses, etc.
        # This part depends on your specific analysis requirements

    counter += 1

