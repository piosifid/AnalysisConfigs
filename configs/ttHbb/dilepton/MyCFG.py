from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.cut_functions import get_nObj_eq, get_nObj_min, get_nObj_less, get_HLTsel, get_nBtagMin, get_nElectron, get_nMuon
from pocket_coffea.parameters.histograms import *
#from params.custom_histograms import *
from pocket_coffea.parameters.cuts import passthrough

import workflow
from workflow import dileptonTriggerProcessor
from math import pi
import custom_cut_functions
import custom_cuts
from custom_cut_functions import *
from custom_cuts import *
# ~ from params.binning import bins
# ~ from params.axis_settings import axis_settings
import os
localdir = os.path.dirname(os.path.abspath(__file__))

# Loading default parameters
from pocket_coffea.parameters import defaults
default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir+"/params")



year = "2023_postBPix"
parameters = defaults.merge_parameters_from_files(default_parameters,
                                                  f"{localdir}/params/object_preselection.yaml",
                                                  f"{localdir}/params/triggers.yaml",
                                         #         f"{localdir}/params/lumi.yaml",
                                                  update=True)
parameters["run_period"] = "Run3"

cfg = Configurator(
    parameters = parameters,
    datasets = {
        "jsons":[
			f"{localdir}/datasets/Run3_MC_ttBkg.json",
		#	f"{localdir}/datasets/Run3_MC_ttBkg_2.json",
			f"{localdir}/datasets/Run3_DATA_DoubleEle.json",
		#	f"{localdir}/datasets/Run3_DATA_DoubleMuon.json",
		#	f"{localdir}/datasets/Run3_DATA_MuonEG.json",
		#	f"{localdir}/datasets/backgrounds_MC_ttbar.json",
			
		],
        "filter" : {
            "samples": [
				"DATA_EGamma",
		#		"DATA_DoubleMuon",
		#		"TTToLNu2Q",
		#		"DATA_MuonEG",
				"TTTo2L2Nu",
		#		"TTToSemiLeptonic",
				
            ],
            "samples_exclude" : [],
            "year": [year]
            # ~ "year": ['2022_preEE','2022_postEE','2023_preBPix','2023_postBPix']
        }
        },
        
        
    workflow = dileptonTriggerProcessor,
    #workflow_options = {},
     # Skimming and categorization
    skim = [
    # get_HLTsel(primaryDatasets=["DoubleEle", "SingleEle"])
             ],
             
    preselections = [dilepton_triggerSF_presel, 
                     get_nObj_eq(2, 15., "ElectronGood"),
                     get_nObj_min(1, 25., "ElectronGood"), 
                     get_nObj_min(2, 15., "JetGood"), 
                    
                     
    ],
    
    categories = {
        "DE" : [   get_nObj_min(2, 15., "ElectronGood"), get_nObj_min(1, 25., "ElectronGood"),
                   get_HLTsel(primaryDatasets=["DoubleEle", "SingleEle"])
        ],
        "inclusive" : [passthrough], 
        
    },
       weights = {
        "common": {
            "inclusive": [ "genWeight",
                           "lumi",
                           "XS",
                           "pileup",
				  "sf_ele_reco",
				  "sf_ele_id", 
				#   "sf_ele_trigger",
				  "sf_mu_id",
				  "sf_mu_iso", 

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
                "inclusive": [  ],
                "bycategory" : {
                }
            },
            "bysample": {
            }
        },
        "shape": {
            "common":{
                "inclusive": [  ]
            }
        }
    },
    
     variables = {
        "ElectronGood_pt" : HistConf(
            [
                Axis(coll="ElectronGood", field="pt", type="variable",
                     bins=[25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250],
                     label="Electron $p_{T}$ [GeV]",
                     lim=(20,250))
            ]
        ),
        "ElectronGood_eta" : HistConf(
            [
                Axis(coll="ElectronGood", field="eta", type="variable",
                     bins=[-2.5, -2.0, -1.5660, -1.4442, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4442, 1.5660, 2.0, 2.5],
                     label="electron $\eta$",
                     lim=(-2.5,2.5))
            ]
        ),
        "ElectronGood_phi" : HistConf(
            [
                Axis(coll="ElectronGood", field="phi",
                     bins=12, start=-pi, stop=pi,
                     label="Electron $\phi$"),
            ]
        ),
        "ElectronGood_pt_1" : HistConf(
            [
                Axis(coll="ElectronGood", field="pt", pos=0, type="variable",
                     bins=[25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250],
                     label="Leading Electron $p_{T}$ [GeV]",
                     lim=(20,200))
            ]
        ),
        "ElectronGood_pt_2" : HistConf(
            [
                Axis(coll="ElectronGood", field="pt", pos=1, type="variable",
                     bins=[25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250],
                     label="Subleading Electron $p_{T}$ [GeV]",
                     lim=(20,200))
            ]
        ),
        "ElectronGood_eta_1" : HistConf(
            [
                Axis(coll="ElectronGood", field="eta", pos=0, type="variable",
                     bins=[-2.5, -2.0, -1.5660, -1.4442, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4442, 1.5660, 2.0, 2.5],
                     label="Leading Electron $\eta$",
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
        "ElectronGood_eta_2" : HistConf(
            [
                Axis(coll="ElectronGood", field="eta", pos=1, type="variable",
                     bins=[-2.5, -2.0, -1.5660, -1.4442, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4442, 1.5660, 2.0, 2.5],
                     label="Subleading Electron $\eta$",
                     lim=(-2.5,2.5))
            ]
        ),
        "ElectronGood_phi_2" : HistConf(
            [
                Axis(coll="ElectronGood", field="phi", pos=1,
                     bins=12, start=-pi, stop=pi,
                     label="Subleading Electron $\phi$"),
            ]
        ),
        "ht" : HistConf(
            [
                Axis(coll="events", field="JetGood_Ht", bins=[40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500], label="$H_T$ [GeV]", lim=(20,500))
            ]
        ),
        "MET" : HistConf(
            [
                Axis(coll="MET", field="pt", bins=[0, 20, 40, 60, 80, 100, 125, 150, 175, 200], label="$MET$ [GeV]", lim=(0,200))
            ]
        ),
        **jet_hists(coll="JetGood"),
        **count_hist(name="nJets", coll="JetGood",bins=10, start=2, stop=12),  
        **count_hist(name="nMuons", coll="MuonGood",bins=3, start=0, stop=3),
        **count_hist(name="nElectrons", coll="ElectronGood",bins=3, start=0, stop=3),
        **count_hist(name="nLeptons", coll="LeptonGood",bins=3, start=0, stop=5),
       
       },
       
       
)
 
 
 
run_options = {
        "executor"       : "dask/lxplus",
        "env"            : "singularity",
        "workers"        : 250,
        "scaleout"       : 100,
        "worker_image"   : "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/pocketcoffea:lxplus-el9-latest",
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
    
    
    
if "dask"  in run_options["executor"]:
    import cloudpickle
    cloudpickle.register_pickle_by_value(workflow)
    cloudpickle.register_pickle_by_value(custom_cut_functions)
    cloudpickle.register_pickle_by_value(custom_cuts)    

