from pocket_coffea.workflows.tthbb_base_processor import ttHbbBaseProcessor
from pocket_coffea.lib.objects import (
  #  jet_correction,
    lepton_selection,
 #   jet_selection,
    btagging,
    get_dilepton,
)
from pocket_coffea.lib.hist_manager import Axis


class dileptonTriggerProcessor(ttHbbBaseProcessor):
    def __init__(self, cfg) -> None:
        super().__init__(cfg=cfg)

        self.output_format["trigger_efficiency"] = {cat: {} for cat in self._categories}

    def apply_object_preselection(self, variation):
        super().apply_object_preselection(variation)
        self.events["ll"] = get_dilepton(
                self.events.ElectronGood, self.events.MuonGood
            )
