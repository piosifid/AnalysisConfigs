import numpy as np
import awkward as ak
from pocket_coffea.lib.cut_definition import Cut
import custom_cut_functions as cuts_f
from collections.abc import Iterable
from pocket_coffea.lib.cut_definition import Cut

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

dilepton_triggerSF_presel= Cut(
    name="dilepton_triggerSF",
    params={
        
        "pt_leading_lepton": 25,
        "mll": 20,
        "mll_SFOS": {'low': 76, 'high': 106},
 	"met": 40,
        "njet": 4,
        "nbjet": 2,
 
    },
    function=cuts_f.dilepton_triggerSF,
)
