import FWCore.ParameterSet.Config as cms

from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
OuterTrackerMonitorTTCluster = DQMEDAnalyzer('OuterTrackerMonitorTTCluster',

    TopFolderName = cms.string('SiOuterTracker'),
    TTClusters    = cms.InputTag("TTClustersFromPhase2TrackerDigis", "ClusterInclusive"),

# Number of clusters per layer
    TH1TTCluster_Barrel = cms.PSet(
        Nbinsx = cms.int32(7),
        xmax = cms.double(7.5),
        xmin = cms.double(0.5)
        ),

# Number of clusters per disc
    TH1TTCluster_ECDiscs = cms.PSet(
        Nbinsx = cms.int32(6),
        xmax = cms.double(6.5),
        xmin = cms.double(0.5)
        ),

# Number of clusters per EC ring
    TH1TTCluster_ECRings = cms.PSet(
        Nbinsx = cms.int32(16),
        xmin = cms.double(0.5),
        xmax = cms.double(16.5)
        ),

# Cluster eta distribution
    TH1TTCluster_Eta = cms.PSet(
        Nbinsx = cms.int32(45),
        xmax = cms.double(5.0),
        xmin = cms.double(-5.0)
        ),

# Cluster phi distribution
    TH1TTCluster_Phi = cms.PSet(
        Nbinsx = cms.int32(60),
        xmax = cms.double(3.5),
        xmin = cms.double(-3.5)
        ),

# Cluster R distribution
    TH1TTCluster_R = cms.PSet(
        Nbinsx = cms.int32(45),
        xmax = cms.double(120),
        xmin = cms.double(0)
        ),

# Cluster nearest distance distribution
    TH1TTCluster_Nearest = cms.PSet(
        Nbinsx = cms.int32(1000),
        xmax = cms.double(10),
        xmin = cms.double(0)
        ),

# Cluster PS nearest distance distribution
    TH1TTCluster_Nearest_PS = cms.PSet(
        Nbinsx = cms.int32(1000),
        xmax = cms.double(5),
        xmin = cms.double(0)
        ),

# Cluster 2S nearest distance distribution
    TH1TTCluster_Nearest_2S = cms.PSet(
        Nbinsx = cms.int32(1111),
        xmax = cms.double(5),
        xmin = cms.double(0)
        ),

# Cluster Width vs. I/O sensor
    TH2TTCluster_Width = cms.PSet(
        Nbinsx = cms.int32(7),
        xmax = cms.double(6.5),
        xmin = cms.double(-0.5),
        Nbinsy = cms.int32(2),
        ymax = cms.double(1.5),
        ymin = cms.double(-0.5)
        ),

# TTCluster forward/backward endcap y vs x
    TH2TTCluster_Position = cms.PSet(
        Nbinsx = cms.int32(960),
        xmax = cms.double(120),
        xmin = cms.double(-120),
        Nbinsy = cms.int32(960),
        ymax = cms.double(120),
        ymin = cms.double(-120)
        ),

# TTCluster forward/backward endcap local y vs x
    TH2TTCluster_Local_Position = cms.PSet(
        Nbinsx = cms.int32(1000),
        xmax = cms.double(5),
        xmin = cms.double(-5),
        Nbinsy = cms.int32(1000),
        ymax = cms.double(5),
        ymin = cms.double(-5)
        ),

# TTCluster 2S local x
    TH2TTCluster_Local_2S_Position = cms.PSet(
        Nbinsx = cms.int32(1000),
        xmax = cms.double(5),
        xmin = cms.double(-5),
        ),

#TTCluster #rho vs z
    TH2TTCluster_RZ = cms.PSet(
        Nbinsx = cms.int32(900),
        xmax = cms.double(300),
        xmin = cms.double(-300),
        Nbinsy = cms.int32(900),
        ymax = cms.double(120),
        ymin = cms.double(0)
        ),
#TTCluster Pixel y position vs width
    TTCluster_YvWidth_Pixel = cms.PSet(
        Nbinsx = cms.int32(1000),
        xmax = cms.double(5),
        xmin = cms.double(-5),
        Nbinsy = cms.int32(7),
        ymax = cms.double(6.5),
        ymin = cms.double(-0.5),
        ),
)
