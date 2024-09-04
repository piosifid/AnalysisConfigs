import awkward as ak
from pocket_coffea.workflows.tthbb_base_processor import ttHbbBaseProcessor
from pocket_coffea.lib.deltaR_matching import metric_eta, metric_phi
from pocket_coffea.lib.deltaR_matching import object_matching
from pocket_coffea.lib.parton_provenance import get_partons_provenance_ttHbb, get_partons_provenance_ttbb4F, get_partons_provenance_tt5F
from pocket_coffea.lib.objects import get_dilepton

from pocket_coffea.lib.objects import (
    btagging,
)



class ttHbb_Run3Test(ttHbbBaseProcessor):
    def __init__(self, cfg) -> None:
        super().__init__(cfg=cfg)
        self.isRun3 = True if self.params["run_period"]=='Run3' else False

    def apply_object_preselection(self, variation):
        super().apply_object_preselection(variation=variation)
        self.events["ll"] = get_dilepton(
            self.events.ElectronGood, self.events.MuonGood
        )
        
    
	
   
