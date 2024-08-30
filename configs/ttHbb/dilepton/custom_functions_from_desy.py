from collections.abc import Iterable
import awkward as ak
import numpy as np
import correctionlib
from pocket_coffea.lib.triggers import get_trigger_mask
from pocket_coffea.lib.cut_functions import get_JetVetoMap_Mask
from pocket_coffea.lib.cut_definition import Cut

import copy
import importlib
import gzip
import cloudpickle
from coffea.jetmet_tools import  CorrectedMETFactory
from pocket_coffea.lib.deltaR_matching import get_matching_pairs_indices, object_matching


def dileptonic(events, params, year, processor_params, sample, isMC, **kwargs):
    # testing(processor_params, events, year)
    lepIds = abs(ak.mask(events.LeptonGood.pdgId, events.nLeptonGood == 2))

    is_em = ak.where(lepIds[:,0] + lepIds[:,1] == 24, True, False) #11 + 13 = 24
    is_em = is_em & get_trigger_mask(events, processor_params.HLT_triggers, year, isMC, 
                                     primaryDatasets=["MuonEG", "SingleMuon", "SingleEle"])
    
    is_mm = ak.where(lepIds[:,0] + lepIds[:,1] == 26, True, False) #13 + 13 = 26
    is_mm = is_mm & get_trigger_mask(events, processor_params.HLT_triggers, year, isMC, 
                                     primaryDatasets=["DoubleMuon", "SingleMuon"])
    
    is_ee = ak.where(lepIds[:,0] + lepIds[:,1] == 22, True, False) #11 + 11 = 22
    is_ee = is_ee & get_trigger_mask(events, processor_params.HLT_triggers, year, isMC, 
                                     primaryDatasets=["DoubleEle", "SingleEle"])
    
    if year in ["2022_preEE", "2022_postEE", "2023_preBPIX", "2023_postBPIX"]:
        mask_jetVetoMap = get_JetVetoMap_Mask(events, params, year, processor_params, sample, isMC, **kwargs)
    else:
        mask_jetVetoMap = True
    mask = (
        (events.nLeptonGood == 2)
        & (ak.firsts(events.LeptonGood.pt) >= params["pt_leading_lepton"])
        & (ak.mask(events.LeptonGood.pt>=params["pt_subleading_lepton"], ak.num(events.LeptonGood.pt)>=2)[:,1])
        & (ak.sum(events.LeptonGood.charge, axis=1) == 0)
        & (events.nJetGood >= params["njet"])
        & (events.nBJetGood >= params["nbjet"])
        & (events.PuppiMET.pt > params["met"])
        & mask_jetVetoMap
        & (is_em  
            | (is_mm & (
                (events.ll.mass > params["dy_window"]["stop"]) 
                | ((events.ll.mass > params["m_ee_mumu_min"])
                    & (events.ll.mass < params["dy_window"]["start"]))
            ))
            | (is_ee & (
                (events.ll.mass > params["dy_window"]["stop"]) 
                | ((events.ll.mass > params["m_ee_mumu_min"])
                    & (events.ll.mass < params["dy_window"]["start"]))
            ))
        )
    )
    return ak.where(ak.is_none(mask), False, mask)
    
    
    
def custom_dilepton(events, params, year, processor_params, **kwargs):
	if params["selec"] == "ee":
		SF = ((events.nMuonGood == 0) & (events.nElectronGood == 2))
		OS = events.ll.charge == 0
		NB = events.nBJetGood >= params["nbjet"]
		# ~ NB = events["BJetGood_L"] >= params["nbjet"]

		mask = (
			(events.nLeptonGood == 2)
			& OS & SF & NB
		)
		return ak.where(ak.is_none(mask), False, mask)
	elif params["selec"] == "em":
		SF = ((events.nMuonGood == 1) & (events.nElectronGood == 1))
		OS = events.ll.charge == 0
		NB = events.nBJetGood >= params["nbjet"]
		# ~ NB = events["BJetGood_L"] >= params["nbjet"]

		mask = (
			(events.nLeptonGood == 2)
			& OS & SF & NB
		)
		return ak.where(ak.is_none(mask), False, mask)
	elif params["selec"] == "mm":
		SF = ((events.nMuonGood == 2) & (events.nElectronGood == 0))
		OS = events.ll.charge == 0
		NB = events.nBJetGood >= params["nbjet"]
		# ~ NB = events["BJetGood_L"] >= params["nbjet"]

		mask = (
			(events.nLeptonGood == 2)
			& OS & SF & NB
		)
		return ak.where(ak.is_none(mask), False, mask)
	else:
		raise Exception("selection name not valid") 


def custom_get_dilepton(selec,nbjet = 1):
    return Cut(
        name="dilepton_cut",
        params={"selec": selec,
				"nbjet": nbjet},
        function= custom_dilepton,
    )
    



def get_ttB_id(ttBid):
    return Cut(
        name="ttB_cut",
        params={"ttBId": ttBid},
        function= ttB_masks,
    )


def ttB_masks(events, params, processor_params, year, isMC, **kwargs):
	genTtbarId = events["genTtbarId"]
	if params["ttBId"] == "ttB":
		return (((abs(genTtbarId) % 100) == 51)
			| ((abs(genTtbarId) % 100) == 52)
			| ((abs(genTtbarId) % 100) == 53)
			| ((abs(genTtbarId) % 100) == 54)
			| ((abs(genTtbarId) % 100) == 55))
	elif params["ttBId"] == "ttC":
		return (((abs(genTtbarId) % 100) == 41)
		    | ((abs(genTtbarId) % 100) == 42)
			| ((abs(genTtbarId) % 100) == 43)
			| ((abs(genTtbarId) % 100) == 44)
			| ((abs(genTtbarId) % 100) == 45))
	elif params["ttBId"] == "ttLF":
		return ~(((abs(genTtbarId) % 100) == 41)
			| ((abs(genTtbarId) % 100) == 42)
			| ((abs(genTtbarId) % 100) == 43)
			| ((abs(genTtbarId) % 100) == 44)
			| ((abs(genTtbarId) % 100) == 45)
			| ((abs(genTtbarId) % 100) == 51)
			| ((abs(genTtbarId) % 100) == 52)
			| ((abs(genTtbarId) % 100) == 53)
			| ((abs(genTtbarId) % 100) == 54)
			| ((abs(genTtbarId) % 100) == 55))


btagger = {
"btagDeepFlavB"         : {"L": 0.0614, "M": 0.3196, "T": 0.7300, "XT": 0.8184, "XXT": 0.9542},
"btagPNetB"             : {"L": 0.0499, "M": 0.2605, "T": 0.6915, "XT": 0.8033, "XXT": 0.9664},
"btagRobustParTAK4B"    : {"L": 0.0897, "M": 0.4510, "T": 0.8604, "XT": 0.9234, "XXT": 0.9893},
}

def custom_btagging(Jet, btag, wp):
    return Jet[Jet[btag] > btagger[btag][wp]]
    
def custom_nBtagMin(events, params, year, processor_params, **kwargs):
    '''Mask for min N jets with minpt and passing btagging.
    The btag params will come from the processor, not from the parameters
    '''
    if params["coll"] == "BJetGood":
        # No need to apply the btaggin on the jet
        # Assume that the collection of clean bjets has been created
        if params["minpt"] > 0.0:
            return (
                ak.sum((events.BJetGood.pt >= params["minpt"]), axis=1) >= params["N"]
            )
        else:
            return events.nBJetGood >= params["N"]
    else:
        if params["minpt"] > 0.0:
            return (
                ak.sum(
                    (
                        events[params["coll"]][params["btagger"]]
                        > btagger[params["btagger"]][params["wp"]]
                    )
                    & (events[params["coll"]].pt >= params["minpt"]),
                    axis=1,
                )
                >= params["N"]
            )
        else:
            return (
                ak.sum(
                    (
                        events[params["coll"]][params["btagger"]]
                        > btagger[params["btagger"]][params["wp"]]
                    ),
                    axis=1,
                )
                >= params["N"]
            )
            

def custom_get_nBtagMin(N, minpt=0, coll="BJetGood", wp="M", name=None , btag = "btagDeepFlavB"):
    if name == None:
        name = f"n{coll}_btagMin{N}_pt{minpt}"
    return Cut(
        name=name, params={"N": N, "coll": coll, "minpt": minpt, "wp": wp, "btagger" : btag}, function=custom_nBtagMin
    )
