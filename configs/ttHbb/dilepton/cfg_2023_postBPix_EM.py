from pocket_coffea.utils.configurator import Configurator
from pocket_coffea.lib.cut_functions import get_nObj_eq, get_nObj_min, get_nObj_less, get_HLTsel, get_nBtagMin, get_nElectron, get_nMuon, get_nPVgood, goldenJson, eventFlags
from pocket_coffea.parameters.histograms import *
#from params.custom_histograms import *
from pocket_coffea.parameters.cuts import passthrough
from pocket_coffea.lib.weights.common import common_weights
import workflow_from_desy
from workflow_from_desy import ttHbb_Run3Test
from math import pi
import custom_functions_from_desy
import custom_cuts_from_desy
from custom_functions_from_desy import *
from custom_cuts_from_desy import *
# ~ from params.binning import bins
# ~ from params.axis_settings import axis_settings
import os
localdir = os.path.dirname(os.path.abspath(__file__))

# Loading default parameters
from pocket_coffea.parameters import defaults
default_parameters = defaults.get_default_parameters()
defaults.register_configuration_dir("config_dir", localdir+"/params")



year = "2023_preBPix"
parameters = defaults.merge_parameters_from_files(default_parameters,
                                                  f"{localdir}/params/object_preselection.yaml",
                                                  f"{localdir}/params/triggers.yaml",
                                               #   f"{localdir}/params/lumi.yaml",
                                                  update=True)
parameters["run_period"] = "Run3"

cfg = Configurator(
    parameters = parameters,
    datasets = {
        "jsons":[
			f"{localdir}/datasets/Run3_MC_ttBkg.json",
		#	f"{localdir}/datasets/Run3_MC_ttBkg_2.json",
		#	f"{localdir}/datasets/Run3_DATA_DoubleEle.json",
		#	f"{localdir}/datasets/Run3_DATA_DoubleMuon.json",
			f"{localdir}/datasets/Run3_DATA_MuonEG.json",
		#	f"{localdir}/datasets/backgrounds_MC_ttbar.json",
			
		],
        "filter" : {
            "samples": [
		#		"DATA_EGamma",
		#		"DATA_DoubleMuon",
		#		"TTToLNu2Q",
	  		        "DATA_MuonEG",
				"TTTo2L2Nu",
		#		"TTToSemiLeptonic",
				
            ],
            "samples_exclude" : [],
            "year": [year],
           
            # ~ "year": ['2022_preEE','2022_postEE','2023_preBPix','2023_postBPix']
        }
        },
        
        
    workflow = ttHbb_Run3Test,
    #workflow_options = {},
     # Skimming and categorization
    skim = [ get_nPVgood(1), eventFlags, goldenJson, get_nObj_min(2, 15., "Jet"),
             ],
             
    preselections = [dileptonic_presel 
                      
    ],
    
    categories = {
         "DATA_EGamma"  : {
                get_nObj_eq(1, 15., "ElectronGood"),
                get_nObj_eq(1, 15., "MuonGood"),
                get_nObj_eq(2, 15., "LeptonGood"),
                get_nObj_min(1, 25., "LeptonGood"),
                get_nBtagMin(N=2, coll="BJetGood_M",wp="M"),
                
                get_HLTsel(primaryDatasets=["MuonEG"]), #i.e. HLT_eleXX_eleXX OR HLT_eleXX
                # ~ get_HLTsel(primaryDatasets=["DoubleEle"]), #i.e. HLT_eleXX_eleXX OR HLT_eleXX
                
            },
        "inclusive" : [passthrough], 
        
    },
       weights_classes = common_weights,
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
        "LeptonGood_pt" : HistConf(
            [
                Axis(coll="LeptonGood", field="pt", type="variable",
                bins=list(range(30, 310, 10)),  # Generate bins from 30 to 300 with a step of 10
                label="Lepton $p_{T}$ [GeV]",
                lim=(20,300))
            ]
        ),
        "LeptonGood_eta" : HistConf(
            [
                Axis(coll="LeptonGood", field="eta", type="variable",
                     bins=[-2.5, -2.0, -1.5660, -1.4442, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4442, 1.5660, 2.0, 2.5],
                     label="lepton $\eta$",
                     lim=(-2.5,2.5))
            ]
        ),
        "LeptonGood_phi" : HistConf(
            [
                Axis(coll="LeptonGood", field="phi",
                     bins=12, start=-pi, stop=pi,
                     label="Lepton $\phi$"),
            ]
        ),
        "ElectronGood_pt_1" : HistConf(
            [
                Axis(coll="ElectronGood", field="pt", pos=0, type="variable",
                bins=list(range(30, 310, 10)),  # Generate bins from 30 to 300 with a step of 10
                label="Leading Electron $p_{T}$ [GeV]",
                lim=(20,300))
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
        "ElectronGood_phi_1" : HistConf(
            [
                Axis(coll="ElectronGood", field="phi", pos=0,
                     bins=12, start=-pi, stop=pi,
                     label="Electron $\phi$"),
            ]
        ),
         "MuonGood_pt_2" : HistConf(
            [
                Axis(coll="MuonGood", field="pt", pos=0, type="variable",
                bins=list(range(30, 310, 10)),  # Generate bins from 30 to 300 with a step of 10
                label="Leading Muon $p_{T}$ [GeV]",
                lim=(20,300))
            ]
        ),
        "MuonGood_eta_2" : HistConf(
            [
                Axis(coll="MuonGood", field="eta", pos=0, type="variable",
                     bins=[-2.5, -2.0, -1.5660, -1.4442, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4442, 1.5660, 2.0, 2.5],
                     label="leading muon $\eta$",
                     lim=(-2.5,2.5))
            ]
        ),
        "MuonGood_phi_2" : HistConf(
            [
                Axis(coll="MuonGood", field="phi", pos=0,
                     bins=12, start=-pi, stop=pi,
                     label="Muon $\phi$"),
            ]
        ),
                "MET" : HistConf(
            [
                Axis(coll="MET", field="pt", bins=[40, 60, 80, 100, 125, 150, 175, 200], label="$MET$ [GeV]", lim=(20,200))
            ]
        ),        "ht" : HistConf(
            [
                Axis(coll="events", field="JetGood_Ht", bins=[40, 60, 80, 100, 125, 150, 175, 200], label="$H_T$ [GeV]", lim=(20,200))
            ]
        ),
         **count_hist(name="nJets", coll="JetGood",bins=10, start=2, stop=12),  
         **count_hist(name="nMuons", coll="MuonGood",bins=3, start=0, stop=3),
         **count_hist(name="nElectrons", coll="ElectronGood",bins=3, start=0, stop=3),
         **count_hist(name="nLeptons", coll="LeptonGood",bins=3, start=0, stop=5),
       
       },
       
       
)
       
 
 
run_options = {
        "executor"       : "dask/lxplus",
        "env"            : "singularity",
        "workers"        : 350,
        "scaleout"       : 100,
        "worker_image"   : "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/pocketcoffea:lxplus-el9-latest",
        "queue"          : "longlunch",
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
    cloudpickle.register_pickle_by_value(workflow_from_desy)
    cloudpickle.register_pickle_by_value(custom_functions_from_desy)
    cloudpickle.register_pickle_by_value(custom_cuts_from_desy)    

