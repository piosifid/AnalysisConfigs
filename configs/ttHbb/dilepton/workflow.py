from pocket_coffea.workflows.tthbb_base_processor import ttHbbBaseProcessor
from pocket_coffea.lib.objects import (
    jet_correction,
    lepton_selection,
    jet_selection,
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


    def define_custom_axes_extra(self):
        # Add 'era' axis for data samples
        if not self._isMC:
            # Extract the eras from the lumi dictionary to build the categorical axis "era"
            eras = list(filter(lambda x : 'tot' not in x, sorted(set(self.cfg.parameters.lumi.picobarns[self._year].keys()))))
            self.custom_axes += [
                Axis(
                    coll="metadata",
                    field="era",
                    name="era",
                    bins=eras,
                    type="strcat",
                    growth=False,
                    label="Era",
                )
            ]
