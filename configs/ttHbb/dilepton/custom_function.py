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
    
    
    


