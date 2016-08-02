import FWCore.ParameterSet.Config as cms

import HLTrigger.HLTfilters.hltHighLevel_cfi
hltHM = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
hltHM.HLTPaths = ['HLT_PAPixelTracks_Multiplicity100_v*',
                  'HLT_PAPixelTracks_Multiplicity130_v*',
                  'HLT_PAPixelTracks_Multiplicity160_v*'
                  #'HLT_PAPixelTracks_Multiplicity190_v*',
                  #'HLT_PAPixelTracks_Multiplicity220_v*'
]

hltHM.andOr = cms.bool(True)
hltHM.throw = cms.bool(False)
				
ana = cms.EDAnalyzer('SC_MixedHarmonics_53X',
                                                  vertexSrc = cms.InputTag('offlinePrimaryVertices'),
                		  		  trackSrc = cms.InputTag('generalTracks'),
                                                  towerSrc = cms.InputTag("towerMaker"),
                                                  offlineDCA = cms.untracked.double(3.0),
                                                  offlineptErr = cms.untracked.double(0.1),
				  		  offlinenhits = cms.untracked.double(0),
						  offlineChi2 = cms.untracked.double(999.9),
						  useCentrality = cms.untracked.bool(False),
				                  useEtaGap = cms.untracked.bool(False),
                                                  reverseBeam = cms.untracked.bool(False),
						  doEffCorrection = cms.untracked.bool(False),
						  Nmin = cms.untracked.int32(1),
                                                  Nmax = cms.untracked.int32(10000),
                                                  eff = cms.untracked.int32(1),
                                                  n1 = cms.untracked.int32(2),
                                                  n2 = cms.untracked.int32(2),
                                                  n3 = cms.untracked.int32(2),
                                                  n4 = cms.untracked.int32(2),
                                                  vzLow = cms.untracked.double(0.0),
                                                  vzHigh = cms.untracked.double(15.0),
                                                  ptLow = cms.untracked.double(0.3),
						  ptHigh = cms.untracked.double(3.0),						 
						  etaTracker = cms.untracked.double(2.4),
     					          gapValue = cms.untracked.double(2.0),
						  etaLowHF = cms.untracked.double(4.4),
                                                  etaHighHF = cms.untracked.double(5.0),
                                                  etaBins = cms.untracked.vdouble(-2.4,-2.3,-2.2,-2.1,-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,
                                                                                  -1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,
                                                                                  -0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,
                                                                                  1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4),
                                                  dEtaBins = cms.untracked.vdouble(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
                                                                                   1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,
                                                                                   2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.8),
						  ptBins = cms.untracked.vdouble(
                                                                      0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
                                                                      0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
                                                                      1.0, 1.05, 1.1, 1.15, 1.2,
                                                                      1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                                                                      2.5, 3.0, 4.0, 5.0, 7.5, 10.0, 12.0, 15.0,
                                                                      20.0, 25.0, 30.0, 45.0, 60.0, 90.0, 120.0, 
                                                                      180.0, 300.0, 500.0
                                                                      ),
						  centBins = cms.untracked.vdouble(4487.37, 4370.52, 4279.22, 4199.03, 4113.36, 4030.29, 3947.37, 3870.82, 3797.94, 3721.48, 3648.81, 3575.99, 3506.18, 3440.74, 3374.5, 3310.49, 3249.72, 3190.49, 3127.43, 3066.91, 3012.16, 2954.08, 2897.16, 2840.3, 2786.54, 2735.06, 2682.83, 2631.95, 2580.71, 2529.93, 2483.34, 2436.59, 2389.05, 2343.58, 2300.27, 2256.49, 2210.35, 2167.14, 2128.09, 2086.24, 2044.85, 2002.72, 1962.42, 1925.23, 1889.2, 1851.68, 1815.58, 1778.47, 1743.48, 1706.47, 1671.08, 1636.7, 1604.94, 1571.63, 1539.86, 1508.37, 1477.12, 1445.73, 1417.7, 1387.98, 1359.02, 1330.3, 1301.45, 1274.07, 1246.54, 1219.36, 1191.97, 1165.77, 1140.4, 1114.92, 1091.98, 1067.94, 1043.67, 1019.66, 995.39, 970.466, 947.786, 924.75, 902.723, 879.824, 859.262, 838.212, 817.18, 796.627, 776.494, 757.142, 737.504, 719.604, 701.142, 684.043, 665.89, 648.427, 630.224, 612.877, 596.435, 580.397, 565.396, 550.272, 535.204, 520.48, 505.854, 491.648, 477.531, 463.192, 449.773, 436.806, 423.944, 410.4, 397.962, 386.135, 374.47, 362.499, 351.17, 339.635, 328.402, 317.875, 307.348, 296.957, 287.002, 276.94, 267.822, 258.796, 249.366, 239.974, 231.563, 223.362, 214.902, 206.818, 199.417, 191.609, 184.184, 177.042, 169.839, 163.579, 157.186, 151.136, 145.165, 139.213, 133.218, 127.748, 122.445, 117.458, 112.715, 108.179, 103.713, 99.2518, 94.8864, 90.7892, 86.692, 82.819, 79.0331, 75.4791, 71.8774, 68.5738, 65.5363, 62.6369, 59.7441, 57.0627, 54.3838, 51.7242, 49.1577, 46.7914, 44.4615, 42.3374, 40.2863, 38.2674, 36.3979, 34.4769, 32.7274, 30.9911, 29.3998, 27.7739, 26.2442, 24.795, 23.3496, 21.8717, 20.5263, 19.2405, 18.08, 16.9542, 15.882, 14.8344, 13.8014, 12.7824, 11.8165, 10.8308, 9.94351, 9.08363, 8.20773, 7.40535, 6.57059, 5.81859, 5.0626, 4.32634, 3.57026, 2.83467, 2.09189, 1.36834, 0.673038, 0.0)


)

