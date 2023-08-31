import numpy as np
import awkward as ak
from pocket_coffea.lib.cut_definition import Cut
from pocket_coffea.lib.triggers import get_trigger_mask


def trigger_mask(events, params, **kwargs):
    mask = np.zeros(len(events), dtype='bool')
    for trigger in params["triggers"]:
        mask = mask | events.HLT[trigger]
    assert (params["category"] in ["pass", "fail"]), "The allowed categories for the trigger selection are 'pass' and 'fail'"
    if params["category"] == "fail":
        mask = ~mask
    return mask




###########################
## Functions to count objects
def count_objects_gt(events, params, **kwargs):
    '''
    Count the number of objects in `params["object"]` and
    keep only events with larger (>) amount than `params["value"]`.
    '''
    mask = ak.num(events[params["object"]], axis=1) > params["value"]
    return ak.where(ak.is_none(mask), False, mask)


def count_objects_lt(events, params, year, sample):
    '''
    Count the number of objects in `params["object"]` and
    keep only events with smaller (<) amount than `params["value"]`.
    '''
    mask = ak.num(events[params["object"]], axis=1) < params["value"]
    return ak.where(ak.is_none(mask), False, mask)


def count_objects_eq(events, params, year, sample):
    '''
    Count the number of objects in `params["object"]` and
    keep only events with same (==) amount than `params["value"]`.
    '''
    mask = ak.num(events[params["object"]], axis=1) == params["value"]
    return ak.where(ak.is_none(mask), False, mask)


def ele_leading_f(events, **kargs):

    mask = ((abs(ak.firsts(events.LeptonGood.pdgId)) == 11) )
    return ak.where(ak.is_none(mask), False, mask)

def mu_leading_f(events, **kargs):

    mask = ((abs(ak.firsts(events.LeptonGood.pdgId)) == 13))
    return ak.where(ak.is_none(mask), False, mask)
    
##################################
# Min number of objects with pt cut


def min_nObj(events, params, **kwargs):
    if f"n{params['coll']}" in events.fields:
        return events[f"n{params['coll']}"] >= params["N"]
    return ak.num(events[params['coll']]) >= params["N"]


def min_nObj_minPt(events, params, **kwargs):
    return ak.sum(events[params["coll"]].pt >= params["minpt"], axis=1) >= params["N"]


def less_nObj(events, params, **kwargs):
    if f"n{params['coll']}" in events.fields:
        return events[f"n{params['coll']}"] < params["N"]
    return ak.num(events[params['coll']]) < params["N"]


def eq_nObj(events, params, **kwargs):
    if f"n{params['coll']}" in events.fields:
        return events[f"n{params['coll']}"] == params["N"]
    return ak.num(events[params['coll']]) == params["N"]


def eq_nObj_minPt(events, params, **kwargs):
    return ak.sum(events[params["coll"]].pt >= params["minpt"], axis=1) == params["N"]


def get_nObj_min(N, minpt=None, coll="JetGood", name=None):
    '''
    Factory function which creates a cut for minimum number of objects.
    Optionally a minimum pT is requested.

    :param N: request >= N objects
    :param coll: collection to use
    :param minpt: minimum pT
    :param name: name for the cut, by defaul it is built as n{coll}_min{N}_pt{minpt}

    :returns: a Cut object
    '''
    if name == None:
        if minpt:
            name = f"n{coll}_min{N}_pt{minpt}"
        else:
            name = f"n{coll}_min{N}"
    if minpt:
        return Cut(
            name=name,
            params={"N": N, "coll": coll, "minpt": minpt},
            function=min_nObj_minPt,
        )
    else:
        return Cut(name=name, params={"N": N, "coll": coll}, function=min_nObj)


def get_nObj_eq(N, minpt=None, coll="JetGood", name=None):
    '''
    Factory function which creates a cut for == number of objects.
    Optionally a minimum pT is requested.

    :param N: request == N objects
    :param coll: collection to use
    :param minpt: minimum pT
    :param name: name for the cut, by defaul it is built as n{coll}_eq{N}_pt{minpt}

    :returns: a Cut object
    '''
    if name == None:
        if minpt:
            name = f"n{coll}_eq{N}_pt{minpt}"
        else:
            name = f"n{coll}_eq{N}"
    if minpt:
        return Cut(
            name=name,
            params={"N": N, "coll": coll, "minpt": minpt},
            function=eq_nObj_minPt,
        )
    else:
        return Cut(name=name, params={"N": N, "coll": coll}, function=eq_nObj)


def get_nObj_less(N, coll="JetGood", name=None):
    '''
    Factory function which creates a cut for < number of objects.

    :param N: request < N objects
    :param coll: collection to use
    :param name: name for the cut, by defaul it is built as n{coll}_less{N}
    :returns: a Cut object
    '''
    if name == None:
        name = f"n{coll}_less{N}"
    return Cut(name=name, params={"N": N, "coll": coll}, function=less_nObj)


def dilepton_triggerSF(events, params, year, sample, **kwargs):
    MET = events[params["METbranch"][year]]

    # Masks for same-flavor (SF) and opposite-sign (OS)
    SF = ((events.nMuonGood == 2) & (events.nElectronGood == 0)) | (
        (events.nMuonGood == 0) & (events.nElectronGood == 2)
    )
    OS = events.ll.charge == 0
    # SFOS = SF & OS
    not_SF = (events.nMuonGood == 1) & (events.nElectronGood == 1)

    mask = (
        (events.nLeptonGood == 2)
        & (ak.firsts(events.LeptonGood.pt) > params["pt_leading_lepton"])
        & (MET.pt > params["met"])
        & OS
        & (events.ll.mass > params["mll"])  # Opposite sign
        &
        # If same-flavour we exclude a mll mass interval
        (
            (
                SF
                & (
                    (events.ll.mass < params["mll_SFOS"]["low"])
                    | (events.ll.mass > params["mll_SFOS"]["high"])
                )
            )
            | not_SF
        )
    )

    # Pad None values with False
    return ak.where(ak.is_none(mask), False, mask)

