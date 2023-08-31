from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.cut_functions import get_nObj_min, get_nObj_eq, get_HLTsel
from pocket_coffea.parameters.histograms import *
from pocket_coffea.parameters.cuts import passthrough
from math import pi

import workflow
from workflow import dileptonTriggerProcessor

import custom_cut_functions
import custom_cuts
from custom_cut_functions import *
from custom_cuts import *
import os
localdir = os.path.dirname(os.path.abspath(__file__))

# Loading default parameters
from pocket_coffea.parameters import defaults
default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir+"/params")

# adding object preselection
parameters = defaults.merge_parameters_from_files(default_parameters,
                                                  f"{localdir}/params/object_preselection_semileptonic.yaml",
                                                  f"{localdir}/params/triggers.yaml",
                                                  f"{localdir}/params/plotting_style.yaml",
                                                  update=True)
cfg = Configurator(
    parameters = parameters,
    datasets = {
        "jsons": [f"{localdir}/datasets/backgrounds_MC_ttbar.json",
                  f"{localdir}/datasets/DATA_MET.json",
                    ],
        "filter" : {
            "samples": ["TTToSemiLeptonic",
                        "TTTo2L2Nu",
                        "DATA_MET"],
            "samples_exclude" : [],
            "year": ['2018']
        },
    },

    workflow = dileptonTriggerProcessor,

    skim = [get_HLTsel(primaryDatasets=["MET"])],
    preselections = [dilepton_triggerSF_presel, 
                       get_nObj_eq(2, 15., "ElectronGood"),
                       get_nObj_min(1, 25., "ElectronGood"),],
    categories = {
        "DoubleEle_pass" : [
            get_HLTsel(primaryDatasets=["DoubleEle"]),
        ],
	"DoubleEle_fail" : [
            get_HLTsel(primaryDatasets=["DoubleEle"], invert=True),
        ],
        "inclusive" : [passthrough],
    },

    weights = {
        "common": {
            "inclusive": ["genWeight","lumi","XS",
                          "pileup",
                          "sf_ele_reco", "sf_ele_id",
                          "sf_mu_id","sf_mu_iso",
                          ],
            "bycategory" : {
            }
        },
        "bysample": {
        }
    },

    variations = {
        "weights": {
            "common": {
                "inclusive": [ "pileup" ],
                "bycategory" : {
                }
            },
            "bysample": {
            }
        },
        "shape": {
            "common":{
                "inclusive": [ "JES_Total", "JER" ]
            }
        }
    },
    
    variables = {
        **muon_hists(coll="MuonGood"),
        **muon_hists(coll="MuonGood", pos=0),
        "MuonGood_pt" : HistConf(
            [
                Axis(coll="MuonGood", field="pt", type="variable",
                     bins=[20, 30, 40, 60, 80, 100, 200],
                     label="Muon $p_{T}$ [GeV]",
                     lim=(20,200))
            ]
        ),
        "ElectronGood_pt" : HistConf(
            [
                Axis(coll="ElectronGood", field="pt", type="variable",
                     bins=[20, 30, 40, 60, 80, 100, 200],
                     label="Electron $p_{T}$ [GeV]",
                     lim=(20,200))
            ]
        ),
        "ElectronGood_pt_1" : HistConf(
            [
                Axis(coll="ElectronGood", field="pt", pos=0, type="variable",
                     bins=[25, 30, 40, 60, 80, 100, 200],
                     label="leading electron $p_{T}$ [GeV]",
                     lim=(20,200))
            ]
        ),
        "MuonGood_pt_1" : HistConf(
            [
                Axis(coll="MuonGood", field="pt", pos=0, type="variable",
                     bins=[25, 30, 40, 60, 80, 100, 200],
                     label="leading muon $p_{T}$ [GeV]",
                     lim=(20,200))
            ]
        ),
        "ElectronGood_pt_2" : HistConf(
            [
                Axis(coll="ElectronGood", field="pt", pos=1, type="variable",
                     bins=[20, 30, 40, 60, 80, 100, 200],
                     label=" subleading electron $p_{T}$ [GeV]",
                     lim=(20,200))
            ]
        ),
        "MuonGood_pt_2" : HistConf(
            [
                Axis(coll="MuonGood", field="pt", pos=1, type="variable",
                     bins=[20, 30, 40, 60, 80, 100, 200],
                     label="subleading muon $p_{T}$ [GeV]",
                     lim=(20,200))
            ]
        ),
        "ElectronGood_etaSC_1" : HistConf(
            [
                Axis(coll="ElectronGood", field="etaSC", pos=0, type="variable",
                     bins=[-2.5, -2.0, -1.5660, -1.4442, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4442, 1.5660, 2.0, 2.5],
                     label="Electron Supercluster $\eta$",
                     lim=(-2.5,2.5))
            ]
        ),
        "ElectronGood_eta_1" : HistConf(
            [
                Axis(coll="ElectronGood", field="eta", pos=0, type="variable",
                     bins=[-2.5, -2.0, -1.5660, -1.4442, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4442, 1.5660, 2.0, 2.5],
                     label="leading electron $\eta$",
                     lim=(-2.5,2.5))
            ]
        ),
        "ElectronGood_eta_2" : HistConf(
            [
                Axis(coll="ElectronGood", field="eta", pos=1, type="variable",
                     bins=[-2.5, -2.0, -1.5660, -1.4442, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4442, 1.5660, 2.0, 2.5],
                     label="subleading electron $\eta$",
                     lim=(-2.5,2.5))
            ]
        ),
        "MuonGood_eta_1" : HistConf(
            [
                Axis(coll="MuonGood", field="eta", pos=0, type="variable",
                     bins=[-2.5, -2.0, -1.5660, -1.4442, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4442, 1.5660, 2.0, 2.5],
                     label="leading muon $\eta$",
                     lim=(-2.5,2.5))
            ]
        ),
        "MuonGood_eta_2" : HistConf(
            [
                Axis(coll="MuonGood", field="eta", pos=2, type="variable",
                     bins=[-2.5, -2.0, -1.5660, -1.4442, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4442, 1.5660, 2.0, 2.5],
                     label="subleading muon $\eta$",
                     lim=(-2.5,2.5))
            ]
        ),
        "ElectronGood_phi_1" : HistConf(
            [
                Axis(coll="ElectronGood", field="phi", pos=0,
                     bins=12, start=-pi, stop=pi,
                     label="Electron $\phi$"),
            ]
        ),
        **jet_hists(coll="JetGood"),
        **count_hist(name="nMuons", coll="MuonGood",bins=3, start=0, stop=3),
        **count_hist(name="nElectrons", coll="ElectronGood",bins=3, start=0, stop=3),
        **count_hist(name="nLeptons", coll="LeptonGood",bins=3, start=0, stop=5),
        **count_hist(name="nJets", coll="JetGood",bins=8, start=0, stop=8),
        **count_hist(name="nBJets", coll="BJetGood",bins=6, start=4, stop=10),
        "ht" : HistConf(
            [
                Axis(coll="events", field="JetGood_Ht", bins=[40, 60, 80, 100, 125, 150, 175, 200], label="$H_T$ [GeV]", lim=(20,200))
            ]
        ),
        "MET" : HistConf(
            [
                Axis(coll="MET", field="pt", bins=[40, 60, 80, 100, 125, 150, 175, 200], label="$MET$ [GeV]", lim=(20,200))
            ]
        ),
        "muon_leading_ele_pt_leading" : HistConf(
            [
                Axis(coll="MuonGood", field="pt", pos=0, type="variable",
                     bins=[20, 50, 80,120, 200],
                     label="leading muon $p_{T}$ [GeV]",
                     lim=(10,200)),
                Axis(coll="ElectronGood", field="pt", pos=0, type="variable",
                     bins=[20, 50, 80,120, 200],
                     label="leading electron $p_{T}$ [GeV]",
                     lim=(0,200)),
            ]
        ),
        "electron_eta_pt_leading" : HistConf(
            [
                Axis(coll="ElectronGood", field="pt", pos=0, type="variable",
                     bins=[20, 50, 80, 120, 200],
                     label="leading electron $p_{T}$ [GeV]",
                     lim=(20,200)),
                Axis(coll="ElectronGood", field="eta", pos=0, type="variable",
                     bins=[-2.5, -1.5, -0.9, -0.4, 0.0, 0.4, 0.9, 1.5, 2.5 ],
                     label="leading electron $\eta$",
                     lim=(-2.5,2.5)),
            ]
        ),
        "ele_eta_pt_leading" : HistConf(
            [
                Axis(coll="ElectronGood", field="pt", pos=0, type="variable",
                     bins=[20, 50, 80,120, 200],
                     label="leading electron $p_{T}$ [GeV]",
                     lim=(20,200)),
                Axis(coll="ElectronGood", field="eta", pos=0, type="variable",
                     bins=[-2.5, -1.5, -0.9, -0.4, 0.0, 0.4, 0.9, 1.5, 2.5 ],
                     label="leading electron $\eta$",
                     lim=(-2.5,2.5)),
            ]
        ),
        "electron_eta_pt_subleading" : HistConf(
            [
                Axis(coll="ElectronGood", field="pt", pos=1, type="variable",
                     bins=[20, 50, 80, 120, 200],
                     label="subleading electron $p_{T}$ [GeV]",
                     lim=(20,200)),
                Axis(coll="ElectronGood", field="eta", pos=1, type="variable",
                     bins=[-2.5, -1.5, -0.9, -0.4, 0.0, 0.4, 0.9, 1.5, 2.5 ],
                     label="subleading electron $\eta$",
                     lim=(-2.5,2.5)),
            ]
        ),
        "electron_pt_both" : HistConf(
            [
                Axis(name="whattt", coll="ElectronGood", field="pt", pos=0, type="variable",
                     bins=[20, 50, 90, 200],
                     label="leading electron $p_{T}$ [GeV]",
                     lim=(20,200)),
                Axis(coll="ElectronGood", field="pt", pos=1, type="variable",
                     bins=[20, 50, 90, 200],
                     label="subleading electron $p_{T}$ [GeV]",
                     lim=(20,200)),
            ]
        ),
        "muon_pt_both" : HistConf(
          [
                Axis(coll="MuonGood", field="pt", pos=0, type="variable",
                     bins=[20, 50, 80,120, 200],
                     label="leading muon $p_{T}$ [GeV]",
                     lim=(20,200)),
                Axis(name="whattt", coll="MuonGood", field="pt", pos=1, type="variable",
                     bins=[20, 50, 90, 200],
                     label="subleading muon $p_{T}$ [GeV]",
                     lim=(20,200)),
            ]
        ),
        "electron_eta_both" : HistConf(
            [
                Axis(name="leading electron h", coll="ElectronGood", field="eta", pos=0, type="variable",
                     bins=[-2.5, -1.2, 0, 1.2, 2.5 ],
                     label="leading electron $\eta$",
                     lim=(-2.5,2.5)),
                Axis(coll="ElectronGood", field="eta", pos=1, type="variable",
                     bins=[-2.5, -1.2, 0, 1.2, 2.5 ],
                     label="subleading electron $\eta$",
                     lim=(-2.5,2.5)),
            ]
        ),
        "muon_eta_both" : HistConf(
            [
                Axis(name="leading muon h", coll="MuonGood", field="eta", pos=0, type="variable",
                     bins=[-2.5, -1.2, 0, 1.2, 2.5 ],
                     label="leading muon $\eta$",
                     lim=(-2.5,2.5)),
                Axis(coll="MuonGood", field="eta", pos=1, type="variable",
                     bins=[-2.5, -1.2, 0, 1.2, 2.5 ],
                     label="subleading muon $\eta$",
                     lim=(-2.5,2.5)),
            ]
        ),
        "muon_eta_pt_leading" : HistConf(
            [
                Axis(coll="MuonGood", field="pt", pos=0, type="variable",
                     bins=[20, 50, 80,120, 200],
                     label="leading muon $p_{T}$ [GeV]",
                     lim=(20,200)),
                Axis(coll="MuonGood", field="eta", pos=0, type="variable",
                     bins=[-2.5, -1.5, -0.9, -0.4, 0.0, 0.4, 0.9, 1.5, 2.5 ],
                     label="leading muon $\eta$",
                     lim=(-2.5,2.5)),
            ]
        ),
        "muon_eta_pt_subleading" : HistConf(
            [
                Axis(coll="MuonGood", field="pt", pos=1, type="variable",
                     bins=[20, 50, 80,120, 200],
                     label="subleading muon $p_{T}$ [GeV]",
                     lim=(20,200)),
                Axis(coll="MuonGood", field="eta", pos=1, type="variable",
                     bins=[-2.5, -1.5, -0.9, -0.4, 0.0, 0.4, 0.9, 1.5, 2.5 ],
                     label="subleading muon $\eta$",
                     lim=(-2.5,2.5)),
            ]
        ),
        "electron_phi_pt_leading" : HistConf(
            [
                Axis(coll="ElectronGood", field="pt", pos=0, type="variable",
                     bins=[30, 35, 40, 50, 60, 70, 80, 90, 100, 130, 200],
                     label="Electron $p_{T}$ [GeV]",
                     lim=(0,200)),
                Axis(coll="ElectronGood", field="phi", pos=0,
                     bins=12, start=-pi, stop=pi,
                     label="Electron $\phi$",
                     lim=(-pi,pi)),
            ]
        ),
    },
)

run_options = {
        "executor"       : "dask/lxplus",
        "env"            : "singularity",
        "workers"        : 1,
        "scaleout"       : 100,
        "worker_image"   : "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/pocketcoffea:lxplus-cc7-latest",
        "queue"          : "microcentury",
        "walltime"       : "02:00:00",
        "mem_per_worker" : "4GB", # GB
        "disk_per_worker" : "1GB", # GB
        "exclusive"      : False,
        "chunk"          : 200000,
        "retries"        : 50,
        "treereduction"  : 20,
        "adapt"          : False,
        
    }

if "dask/lxplus"  in run_options["executor"]:
    import cloudpickle
    cloudpickle.register_pickle_by_value(workflow)
    cloudpickle.register_pickle_by_value(custom_cut_functions)
    cloudpickle.register_pickle_by_value(custom_cuts)
