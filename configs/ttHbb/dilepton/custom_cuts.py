import numpy as np
import awkward as ak
from pocket_coffea.lib.cut_definition import Cut
import custom_cut_functions as cuts_f

def get_trigger_passfail(triggers, category):
    return Cut(
        name=f"{'_'.join(triggers)}_{category}",
        params={"triggers": triggers, "category": category},
        function=cuts_f.trigger_mask
    )

def get_ht_above(minht, name=None):
    if name == None:
        name = f"minht{minht}"
    return Cut(
        name=name,
        params={"minht": minht},
        function=cuts_f.ht_above
    )

def get_ht_below(maxht, name=None):
    if name == None:
        name = f"maxht{maxht}"
    return Cut(
        name=name,
        params={"maxht": maxht},
        function=cuts_f.ht_below
    )

dilepton_presel_triggerSF= Cut(
    name="dilepton_triggerSF",
    params={
        "METbranch": {
            '2016_PreVFP': "MET",
            '2016_PostVFP': "MET",
            '2017': "MET",
            '2018': "MET",
        },
        "pt_leading_lepton": 25,
        "mll": 20,
        "mll_SFOS": {'low': 76, 'high': 106},
#	"met": 0,
        "njet": 2,
        "nbjet": 1,
        "pt_leading_lepton": 25,
        "met": 40,
        "mll": 20,
        "mll_SFOS": {'low': 76, 'high': 106},

    },
    function=cuts_f.dilepton_triggerSF,
)
