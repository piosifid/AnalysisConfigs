import awkward as ak
from pocket_coffea.workflows.tthbb_base_processor import ttHbbBaseProcessor
from pocket_coffea.lib.deltaR_matching import metric_eta, metric_phi
from pocket_coffea.lib.deltaR_matching import object_matching
from pocket_coffea.lib.parton_provenance import get_partons_provenance_ttHbb, get_partons_provenance_ttbb4F, get_partons_provenance_tt5F
from pocket_coffea.lib.objects import get_dilepton

from pocket_coffea.lib.objects import (
    btagging,
)

from custom_functions_from_desy import custom_btagging

class ttHbb_Run3Test(ttHbbBaseProcessor):
    def __init__(self, cfg) -> None:
        super().__init__(cfg=cfg)
        self.isRun3 = True if self.params["run_period"]=='Run3' else False

    def apply_object_preselection(self, variation):
        super().apply_object_preselection(variation=variation)
        self.events["ll"] = get_dilepton(
            self.events.ElectronGood, self.events.MuonGood
        )
        self.events["BJetGood_L"] = btagging(
            self.events["JetGood"], self.params.btagging.working_point[self._year], wp="L"
        )
        self.events["BJetGood_M"] = btagging(
            self.events["JetGood"], self.params.btagging.working_point[self._year], wp="M"
        )
        self.events["BJetGood_T"] = btagging(
            self.events["JetGood"], self.params.btagging.working_point[self._year], wp="T"
        )
        self.events["BJetGood_DF"] = custom_btagging(
            self.events["JetGood"], btag="btagDeepFlavB", wp="L"
        )
        self.events["BJetGood_PN"] = custom_btagging(
            self.events["JetGood"], btag="btagPNetB", wp="L"
        )
        self.events["BJetGood_R"] = custom_btagging(
            self.events["JetGood"], btag="btagRobustParTAK4B", wp="L"
        )
    
	
    def define_common_variables_after_presel(self, variation):
        super().define_common_variables_before_presel(variation=variation)
        # ~ print(self.events["BJetGood_L"]["btagDeepFlavB"][1][0:50])
        # ~ print(ak.to_list(self.events["BJetGood_M"]["btagDeepFlavB"][:50]))
        # ~ print(ak.to_list(self.events["BJetGood_L"]["btagDeepFlavB"][:50]))
        # ~ print(self.events["BJetGood_L"])

        # Compute deltaR(b, b) of all possible b-jet pairs.
        # We require deltaR > 0 to exclude the deltaR between the jets with themselves
        deltaR = ak.flatten(self.events["BJetGood"].metric_table(self.events["BJetGood"]), axis=2)
        deltaEta = ak.flatten(self.events["BJetGood"].metric_table(self.events["BJetGood"], metric=metric_eta), axis=2)
        deltaPhi = ak.flatten(self.events["BJetGood"].metric_table(self.events["BJetGood"], metric=metric_phi), axis=2)
        deltaR = deltaR[deltaR > 0.]
        deltaEta = deltaEta[deltaEta > 0.]
        deltaPhi = deltaPhi[deltaPhi > 0.]

        # Get the deltaR with no possibility of repetition of identical b-jet pairs

        # Get all the possible combinations of b-jet pairs
        pairs = ak.argcombinations(self.events["BJetGood"], 2, axis=1)
        b1 = self.events["BJetGood"][pairs.slot0]
        b2 = self.events["BJetGood"][pairs.slot1]

        # Compute deltaR between the pairs
        deltaR_unique = b1.delta_r(b2)
        idx_pairs_sorted = ak.argsort(deltaR_unique, axis=1)
        pairs_sorted = pairs[idx_pairs_sorted]

        # Compute the minimum deltaR(b, b), deltaEta(b, b), deltaPhi(b, b) and the invariant mass of the closest b-jet pair
        self.events["deltaRbb_min"] = ak.min(deltaR, axis=1)
        self.events["deltaEtabb_min"] = ak.min(deltaEta, axis=1)
        self.events["deltaPhibb_min"] = ak.min(deltaPhi, axis=1)
        self.events["mbb"] = (self.events["BJetGood"][pairs_sorted.slot0] + self.events["BJetGood"][pairs_sorted.slot1]).mass
        
        
        
        # calculation for DeepFlavB btagger
        deltaR = ak.flatten(self.events["BJetGood_DF"].metric_table(self.events["BJetGood_DF"]), axis=2)
        deltaEta = ak.flatten(self.events["BJetGood_DF"].metric_table(self.events["BJetGood_DF"], metric=metric_eta), axis=2)
        deltaPhi = ak.flatten(self.events["BJetGood_DF"].metric_table(self.events["BJetGood_DF"], metric=metric_phi), axis=2)
        deltaR = deltaR[deltaR > 0.]
        deltaEta = deltaEta[deltaEta > 0.]
        deltaPhi = deltaPhi[deltaPhi > 0.]
        pairs = ak.argcombinations(self.events["BJetGood_DF"], 2, axis=1)
        b1 = self.events["BJetGood_DF"][pairs.slot0]
        b2 = self.events["BJetGood_DF"][pairs.slot1]
        deltaR_unique = b1.delta_r(b2)
        idx_pairs_sorted = ak.argsort(deltaR_unique, axis=1)
        pairs_sorted = pairs[idx_pairs_sorted]
        self.events["deltaRbb_min_DF"] = ak.min(deltaR, axis=1)
        self.events["deltaEtabb_min_DF"] = ak.min(deltaEta, axis=1)
        self.events["deltaPhibb_min_DF"] = ak.min(deltaPhi, axis=1)
        self.events["mbb_DF"] = (self.events["BJetGood_DF"][pairs_sorted.slot0] + self.events["BJetGood_DF"][pairs_sorted.slot1]).mass
        
        
        # calculation for PNetB btagger
        deltaR = ak.flatten(self.events["BJetGood_PN"].metric_table(self.events["BJetGood_PN"]), axis=2)
        deltaEta = ak.flatten(self.events["BJetGood_PN"].metric_table(self.events["BJetGood_PN"], metric=metric_eta), axis=2)
        deltaPhi = ak.flatten(self.events["BJetGood_PN"].metric_table(self.events["BJetGood_PN"], metric=metric_phi), axis=2)
        deltaR = deltaR[deltaR > 0.]
        deltaEta = deltaEta[deltaEta > 0.]
        deltaPhi = deltaPhi[deltaPhi > 0.]
        pairs = ak.argcombinations(self.events["BJetGood_PN"], 2, axis=1)
        b1 = self.events["BJetGood_PN"][pairs.slot0]
        b2 = self.events["BJetGood_PN"][pairs.slot1]
        deltaR_unique = b1.delta_r(b2)
        idx_pairs_sorted = ak.argsort(deltaR_unique, axis=1)
        pairs_sorted = pairs[idx_pairs_sorted]
        self.events["deltaRbb_min_PN"] = ak.min(deltaR, axis=1)
        self.events["deltaEtabb_min_PN"] = ak.min(deltaEta, axis=1)
        self.events["deltaPhibb_min_PN"] = ak.min(deltaPhi, axis=1)
        self.events["mbb_PN"] = (self.events["BJetGood_PN"][pairs_sorted.slot0] + self.events["BJetGood_PN"][pairs_sorted.slot1]).mass
        
        
        # calculation for RobustParTAK4B btagger
        deltaR = ak.flatten(self.events["BJetGood_R"].metric_table(self.events["BJetGood_R"]), axis=2)
        deltaEta = ak.flatten(self.events["BJetGood_R"].metric_table(self.events["BJetGood_R"], metric=metric_eta), axis=2)
        deltaPhi = ak.flatten(self.events["BJetGood_R"].metric_table(self.events["BJetGood_R"], metric=metric_phi), axis=2)
        deltaR = deltaR[deltaR > 0.]
        deltaEta = deltaEta[deltaEta > 0.]
        deltaPhi = deltaPhi[deltaPhi > 0.]
        pairs = ak.argcombinations(self.events["BJetGood_R"], 2, axis=1)
        b1 = self.events["BJetGood_R"][pairs.slot0]
        b2 = self.events["BJetGood_R"][pairs.slot1]
        deltaR_unique = b1.delta_r(b2)
        idx_pairs_sorted = ak.argsort(deltaR_unique, axis=1)
        pairs_sorted = pairs[idx_pairs_sorted]
        self.events["deltaRbb_min_R"] = ak.min(deltaR, axis=1)
        self.events["deltaEtabb_min_R"] = ak.min(deltaEta, axis=1)
        self.events["deltaPhibb_min_R"] = ak.min(deltaPhi, axis=1)
        self.events["mbb_R"] = (self.events["BJetGood_R"][pairs_sorted.slot0] + self.events["BJetGood_R"][pairs_sorted.slot1]).mass
