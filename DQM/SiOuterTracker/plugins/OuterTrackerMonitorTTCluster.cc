// -*- C++ -*-
//
// Package:    SiOuterTracker
// Class:      SiOuterTracker
//
/**\class SiOuterTracker OuterTrackerMonitorTTCluster.cc
 DQM/SiOuterTracker/plugins/OuterTrackerMonitorTTCluster.cc

 Description: [one line class summary]

 Implementation:
 [Notes on implementation]
 */
//
// Original Author:  Isabelle Helena J De Bruyn
//         Created:  Mon, 10 Feb 2014 13:57:08 GMT
//

// system include files
#include <memory>
#include <numeric>
#include <vector>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"

class OuterTrackerMonitorTTCluster : public DQMEDAnalyzer {
public:
  explicit OuterTrackerMonitorTTCluster(const edm::ParameterSet &);
  ~OuterTrackerMonitorTTCluster() override;
  void analyze(const edm::Event &, const edm::EventSetup &) override;
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  void dqmBeginRun(const edm::Run &iRun, const edm::EventSetup &iSetup) override;
  // TTCluster stacks
  MonitorElement *Cluster_IMem_Barrel = nullptr;
  MonitorElement *Cluster_IMem_Endcap_Disc = nullptr;
  MonitorElement *Cluster_IMem_Endcap_Ring = nullptr;
  MonitorElement *Cluster_IMem_Endcap_Ring_Fw[5] = {nullptr, nullptr, nullptr, nullptr, nullptr};
  MonitorElement *Cluster_IMem_Endcap_Ring_Bw[5] = {nullptr, nullptr, nullptr, nullptr, nullptr};
  MonitorElement *Cluster_OMem_Barrel = nullptr;
  MonitorElement *Cluster_OMem_Endcap_Disc = nullptr;
  MonitorElement *Cluster_OMem_Endcap_Ring = nullptr;
  MonitorElement *Cluster_OMem_Endcap_Ring_Fw[5] = {nullptr, nullptr, nullptr, nullptr, nullptr};
  MonitorElement *Cluster_OMem_Endcap_Ring_Bw[5] = {nullptr, nullptr, nullptr, nullptr, nullptr};
  MonitorElement *Cluster_W = nullptr;
  MonitorElement *Cluster_IMem_W = nullptr;
  MonitorElement *Cluster_OMem_W = nullptr;
  MonitorElement *Cluster_Phi = nullptr;
  MonitorElement *Cluster_R = nullptr;
  MonitorElement *Cluster_Eta = nullptr;
  MonitorElement *Cluster_Nearest = nullptr;
  MonitorElement *Cluster_Nearest_PS = nullptr;
  MonitorElement *Cluster_Nearest_2S = nullptr;

  MonitorElement *Cluster_Barrel_XY = nullptr;
  MonitorElement *Cluster_Endcap_Fw_XY = nullptr;
  MonitorElement *Cluster_Endcap_Bw_XY = nullptr;
  MonitorElement *Cluster_Barrel_Local_PS_XY = nullptr;
  MonitorElement *Cluster_Endcap_Fw_Local_PS_XY = nullptr;
  MonitorElement *Cluster_Endcap_Bw_Local_PS_XY = nullptr;
  MonitorElement *Cluster_Barrel_Local_2S_XY = nullptr;
  MonitorElement *Cluster_Endcap_Fw_Local_2S_XY = nullptr;
  MonitorElement *Cluster_Endcap_Bw_Local_2S_XY = nullptr;
  MonitorElement *Cluster_Barrel_Local_Pixel_XY = nullptr;
  MonitorElement *Cluster_Endcap_Fw_Local_Pixel_XY = nullptr;
  MonitorElement *Cluster_Endcap_Bw_Local_Pixel_XY = nullptr;
  MonitorElement *Cluster_Barrel_Local_Strip_XY = nullptr;
  MonitorElement *Cluster_Endcap_Fw_Local_Strip_XY = nullptr;
  MonitorElement *Cluster_Endcap_Bw_Local_Strip_XY = nullptr;
  MonitorElement *Cluster_RZ = nullptr;
  MonitorElement *Cluster_YvWidth_Pixel = nullptr;

private:
  edm::ParameterSet conf_;
  edm::EDGetTokenT<edmNew::DetSetVector<TTCluster<Ref_Phase2TrackerDigi_>>> tagTTClustersToken_;
  std::string topFolderName_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topoToken_;
  const TrackerGeometry *tkGeom_ = nullptr;
  const TrackerTopology *tTopo_ = nullptr;
};

//
// constructors and destructor
//
OuterTrackerMonitorTTCluster::OuterTrackerMonitorTTCluster(const edm::ParameterSet &iConfig)
    : conf_(iConfig),
      geomToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord, edm::Transition::BeginRun>()),
      topoToken_(esConsumes<TrackerTopology, TrackerTopologyRcd, edm::Transition::BeginRun>()) {
  topFolderName_ = conf_.getParameter<std::string>("TopFolderName");
  tagTTClustersToken_ = consumes<edmNew::DetSetVector<TTCluster<Ref_Phase2TrackerDigi_>>>(
      conf_.getParameter<edm::InputTag>("TTClusters"));
}

OuterTrackerMonitorTTCluster::~OuterTrackerMonitorTTCluster() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//
void OuterTrackerMonitorTTCluster::dqmBeginRun(const edm::Run &iRun, const edm::EventSetup &iSetup) {
  tkGeom_ = &(iSetup.getData(geomToken_));
  tTopo_ = &(iSetup.getData(topoToken_));
}

// ------------ method called for each event  ------------
void OuterTrackerMonitorTTCluster::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
  /// Track Trigger Clusters
  edm::Handle<edmNew::DetSetVector<TTCluster<Ref_Phase2TrackerDigi_>>> Phase2TrackerDigiTTClusterHandle;
  iEvent.getByToken(tagTTClustersToken_, Phase2TrackerDigiTTClusterHandle);

  /// Loop over the input Clusters
  typename edmNew::DetSetVector<TTCluster<Ref_Phase2TrackerDigi_>>::const_iterator inputIter;
  typename edmNew::DetSet<TTCluster<Ref_Phase2TrackerDigi_>>::const_iterator contentIter;
  typename edmNew::DetSet<TTCluster<Ref_Phase2TrackerDigi_>>::const_iterator contentIter2;

  // Adding protection
  if (!Phase2TrackerDigiTTClusterHandle.isValid())
    return;

  for (inputIter = Phase2TrackerDigiTTClusterHandle->begin(); inputIter != Phase2TrackerDigiTTClusterHandle->end();
       ++inputIter) {
    for (contentIter = inputIter->begin(); contentIter != inputIter->end(); ++contentIter) {
      // Make reference cluster
      edm::Ref<edmNew::DetSetVector<TTCluster<Ref_Phase2TrackerDigi_>>, TTCluster<Ref_Phase2TrackerDigi_>> tempCluRef =
          edmNew::makeRefTo(Phase2TrackerDigiTTClusterHandle, contentIter);

      DetId detIdClu = tkGeom_->idToDet(tempCluRef->getDetId())->geographicalId();
      unsigned int memberClu = tempCluRef->getStackMember();
      unsigned int widClu = tempCluRef->findWidth();

      MeasurementPoint mp = tempCluRef->findAverageLocalCoordinates();
      const GeomDet *theGeomDet = tkGeom_->idToDet(detIdClu);
      LocalPoint posCluLocal = theGeomDet->topology().localPosition(mp);
      Global3DPoint posClu = theGeomDet->surface().toGlobal(theGeomDet->topology().localPosition(mp));

      // PS (=0) or 2S (=1)
      unsigned int typeClu = 0;
      if (detIdClu.subdetId() == static_cast<int>(StripSubdetector::TOB))
      {
        unsigned int barrel_layer = tTopo_->layer(detIdClu);
        if (barrel_layer==4 || barrel_layer==5 || barrel_layer==6) typeClu=1;
      }
      else if (detIdClu.subdetId() == static_cast<int>(StripSubdetector::TID))
      { 
        unsigned int endcap_disk = tTopo_->layer(detIdClu);
        unsigned int endcap_ring = tTopo_->tidRing(detIdClu);
        if (endcap_disk<=2) {
          if (endcap_ring>10) typeClu=1;
        }
        else if (endcap_disk>2) {
          if (endcap_ring>7) typeClu=1;
        }
      }

      double r = posClu.perp();
      double z = posClu.z();

      Cluster_W->Fill(widClu, memberClu);
      if (memberClu == 0) {
        Cluster_IMem_W->Fill(widClu);
      } else {
        Cluster_OMem_W->Fill(widClu);
      }

      Cluster_Eta->Fill(posClu.eta());
      Cluster_Phi->Fill(posClu.phi());
      Cluster_R->Fill(r);
      Cluster_RZ->Fill(z, r);

      // Iterate over other clusters to calculate closest distance
      double d_nearest = 9999.0;
      for (contentIter2 = inputIter->begin(); contentIter2 != inputIter->end(); ++contentIter2) {
        if (contentIter==contentIter2) continue;

        // Make reference cluster 2
        edm::Ref<edmNew::DetSetVector<TTCluster<Ref_Phase2TrackerDigi_>>, TTCluster<Ref_Phase2TrackerDigi_>> tempCluRef2 =
            edmNew::makeRefTo(Phase2TrackerDigiTTClusterHandle, contentIter2);

        unsigned int memberClu2 = tempCluRef2->getStackMember();
        DetId detIdClu2 = tkGeom_->idToDet(tempCluRef2->getDetId())->geographicalId();
        MeasurementPoint mp2 = tempCluRef2->findAverageLocalCoordinates();
        const GeomDet *theGeomDet2 = tkGeom_->idToDet(detIdClu2);
        Global3DPoint posClu2 = theGeomDet2->surface().toGlobal(theGeomDet2->topology().localPosition(mp2));


        if (detIdClu==detIdClu2 && memberClu==memberClu2) {
          double d = abs(posClu.x()-posClu2.x()) + abs(posClu.y()-posClu2.y()) + abs(posClu.z()-posClu2.z());
          if (d < d_nearest) d_nearest = d;
        }
      }

      Cluster_Nearest->Fill(d_nearest);
      if (typeClu)
        Cluster_Nearest_2S->Fill(d_nearest);
      else
        Cluster_Nearest_PS->Fill(d_nearest);
        
      if (detIdClu.subdetId() == static_cast<int>(StripSubdetector::TOB))  // Phase 2 Outer Tracker Barrel
      {
        if (memberClu == 0)
          Cluster_IMem_Barrel->Fill(tTopo_->layer(detIdClu));
        else
          Cluster_OMem_Barrel->Fill(tTopo_->layer(detIdClu));

        Cluster_Barrel_XY->Fill(posClu.x(), posClu.y());

        if (typeClu) {
          Cluster_Barrel_Local_2S_XY->Fill(posCluLocal.x());
        }
        else {
          Cluster_Barrel_Local_PS_XY->Fill(posCluLocal.x(), posCluLocal.y());
          if (memberClu == 0) {
            Cluster_Barrel_Local_Pixel_XY->Fill(posCluLocal.x(), posCluLocal.y());
            Cluster_YvWidth_Pixel->Fill(posCluLocal.y(), widClu);
          }
          else {
            Cluster_Barrel_Local_Strip_XY->Fill(posCluLocal.x());
          }
        }
      }                                                                         // end if isBarrel
      else if (detIdClu.subdetId() == static_cast<int>(StripSubdetector::TID))  // Phase 2 Outer Tracker Endcap
      {
        if (memberClu == 0) {
          Cluster_IMem_Endcap_Disc->Fill(tTopo_->layer(detIdClu));  // returns wheel
          Cluster_IMem_Endcap_Ring->Fill(tTopo_->tidRing(detIdClu));
        } else {
          Cluster_OMem_Endcap_Disc->Fill(tTopo_->layer(detIdClu));  // returns wheel
          Cluster_OMem_Endcap_Ring->Fill(tTopo_->tidRing(detIdClu));
        }

        if (posClu.z() > 0) {
          Cluster_Endcap_Fw_XY->Fill(posClu.x(), posClu.y());

          if (typeClu) {
            Cluster_Endcap_Fw_Local_2S_XY->Fill(posCluLocal.x());
          }
          else {
            Cluster_Endcap_Fw_Local_PS_XY->Fill(posCluLocal.x(), posCluLocal.y());
            if (memberClu == 0) {
              Cluster_Endcap_Fw_Local_Pixel_XY->Fill(posCluLocal.x(), posCluLocal.y());
              Cluster_YvWidth_Pixel->Fill(posCluLocal.y(), widClu);
            }
            else{
              Cluster_Endcap_Fw_Local_Strip_XY->Fill(posCluLocal.x());
            }
          }

          if (memberClu == 0)
            Cluster_IMem_Endcap_Ring_Fw[tTopo_->layer(detIdClu) - 1]->Fill(tTopo_->tidRing(detIdClu));
          else
            Cluster_OMem_Endcap_Ring_Fw[tTopo_->layer(detIdClu) - 1]->Fill(tTopo_->tidRing(detIdClu));
        } else {
          Cluster_Endcap_Bw_XY->Fill(posClu.x(), posClu.y());

          if (typeClu) {
            Cluster_Endcap_Bw_Local_2S_XY->Fill(posCluLocal.x());
          }
          else {
            Cluster_Endcap_Bw_Local_PS_XY->Fill(posCluLocal.x(), posCluLocal.y());
            if (memberClu == 0) {
              Cluster_Endcap_Bw_Local_Pixel_XY->Fill(posCluLocal.x(), posCluLocal.y());
              Cluster_YvWidth_Pixel->Fill(posCluLocal.y(), widClu);
            }
            else {
              Cluster_Endcap_Bw_Local_Strip_XY->Fill(posCluLocal.x());
            }
          }

          if (memberClu == 0)
            Cluster_IMem_Endcap_Ring_Bw[tTopo_->layer(detIdClu) - 1]->Fill(tTopo_->tidRing(detIdClu));
          else
            Cluster_OMem_Endcap_Ring_Bw[tTopo_->layer(detIdClu) - 1]->Fill(tTopo_->tidRing(detIdClu));
        }
      }  // end if isEndcap
    }    // end loop contentIter
  }      // end loop inputIter
}  // end of method

// ------------ method called once each job just before starting event loop
// ------------
void OuterTrackerMonitorTTCluster::bookHistograms(DQMStore::IBooker &iBooker,
                                                  edm::Run const &run,
                                                  edm::EventSetup const &es) {
  std::string HistoName;
  const int numDiscs = 5;

  iBooker.setCurrentFolder(topFolderName_ + "/Clusters/NClusters");

  // NClusters
  edm::ParameterSet psTTCluster_Barrel = conf_.getParameter<edm::ParameterSet>("TH1TTCluster_Barrel");
  HistoName = "NClusters_IMem_Barrel";
  Cluster_IMem_Barrel = iBooker.book1D(HistoName,
                                       HistoName,
                                       psTTCluster_Barrel.getParameter<int32_t>("Nbinsx"),
                                       psTTCluster_Barrel.getParameter<double>("xmin"),
                                       psTTCluster_Barrel.getParameter<double>("xmax"));
  Cluster_IMem_Barrel->setAxisTitle("Barrel Layer", 1);
  Cluster_IMem_Barrel->setAxisTitle("# L1 Clusters", 2);

  HistoName = "NClusters_OMem_Barrel";
  Cluster_OMem_Barrel = iBooker.book1D(HistoName,
                                       HistoName,
                                       psTTCluster_Barrel.getParameter<int32_t>("Nbinsx"),
                                       psTTCluster_Barrel.getParameter<double>("xmin"),
                                       psTTCluster_Barrel.getParameter<double>("xmax"));
  Cluster_OMem_Barrel->setAxisTitle("Barrel Layer", 1);
  Cluster_OMem_Barrel->setAxisTitle("# L1 Clusters", 2);

  edm::ParameterSet psTTCluster_ECDisc = conf_.getParameter<edm::ParameterSet>("TH1TTCluster_ECDiscs");
  HistoName = "NClusters_IMem_Endcap_Disc";
  Cluster_IMem_Endcap_Disc = iBooker.book1D(HistoName,
                                            HistoName,
                                            psTTCluster_ECDisc.getParameter<int32_t>("Nbinsx"),
                                            psTTCluster_ECDisc.getParameter<double>("xmin"),
                                            psTTCluster_ECDisc.getParameter<double>("xmax"));
  Cluster_IMem_Endcap_Disc->setAxisTitle("Endcap Disc", 1);
  Cluster_IMem_Endcap_Disc->setAxisTitle("# L1 Clusters", 2);

  HistoName = "NClusters_OMem_Endcap_Disc";
  Cluster_OMem_Endcap_Disc = iBooker.book1D(HistoName,
                                            HistoName,
                                            psTTCluster_ECDisc.getParameter<int32_t>("Nbinsx"),
                                            psTTCluster_ECDisc.getParameter<double>("xmin"),
                                            psTTCluster_ECDisc.getParameter<double>("xmax"));
  Cluster_OMem_Endcap_Disc->setAxisTitle("Endcap Disc", 1);
  Cluster_OMem_Endcap_Disc->setAxisTitle("# L1 Clusters", 2);

  edm::ParameterSet psTTCluster_ECRing = conf_.getParameter<edm::ParameterSet>("TH1TTCluster_ECRings");
  HistoName = "NClusters_IMem_Endcap_Ring";
  Cluster_IMem_Endcap_Ring = iBooker.book1D(HistoName,
                                            HistoName,
                                            psTTCluster_ECRing.getParameter<int32_t>("Nbinsx"),
                                            psTTCluster_ECRing.getParameter<double>("xmin"),
                                            psTTCluster_ECRing.getParameter<double>("xmax"));
  Cluster_IMem_Endcap_Ring->setAxisTitle("Endcap Ring", 1);
  Cluster_IMem_Endcap_Ring->setAxisTitle("# L1 Clusters", 2);

  HistoName = "NClusters_OMem_Endcap_Ring";
  Cluster_OMem_Endcap_Ring = iBooker.book1D(HistoName,
                                            HistoName,
                                            psTTCluster_ECRing.getParameter<int32_t>("Nbinsx"),
                                            psTTCluster_ECRing.getParameter<double>("xmin"),
                                            psTTCluster_ECRing.getParameter<double>("xmax"));
  Cluster_OMem_Endcap_Ring->setAxisTitle("Endcap Ring", 1);
  Cluster_OMem_Endcap_Ring->setAxisTitle("# L1 Clusters", 2);

  for (int i = 0; i < numDiscs; i++) {
    HistoName = "NClusters_IMem_Disc+" + std::to_string(i + 1);
    Cluster_IMem_Endcap_Ring_Fw[i] = iBooker.book1D(HistoName,
                                                    HistoName,
                                                    psTTCluster_ECRing.getParameter<int32_t>("Nbinsx"),
                                                    psTTCluster_ECRing.getParameter<double>("xmin"),
                                                    psTTCluster_ECRing.getParameter<double>("xmax"));
    Cluster_IMem_Endcap_Ring_Fw[i]->setAxisTitle("Endcap Ring", 1);
    Cluster_IMem_Endcap_Ring_Fw[i]->setAxisTitle("# L1 Clusters ", 2);
  }

  for (int i = 0; i < numDiscs; i++) {
    HistoName = "NClusters_IMem_Disc-" + std::to_string(i + 1);
    Cluster_IMem_Endcap_Ring_Bw[i] = iBooker.book1D(HistoName,
                                                    HistoName,
                                                    psTTCluster_ECRing.getParameter<int32_t>("Nbinsx"),
                                                    psTTCluster_ECRing.getParameter<double>("xmin"),
                                                    psTTCluster_ECRing.getParameter<double>("xmax"));
    Cluster_IMem_Endcap_Ring_Bw[i]->setAxisTitle("Endcap Ring", 1);
    Cluster_IMem_Endcap_Ring_Bw[i]->setAxisTitle("# L1 Clusters ", 2);
  }

  for (int i = 0; i < numDiscs; i++) {
    HistoName = "NClusters_OMem_Disc+" + std::to_string(i + 1);
    Cluster_OMem_Endcap_Ring_Fw[i] = iBooker.book1D(HistoName,
                                                    HistoName,
                                                    psTTCluster_ECRing.getParameter<int32_t>("Nbinsx"),
                                                    psTTCluster_ECRing.getParameter<double>("xmin"),
                                                    psTTCluster_ECRing.getParameter<double>("xmax"));
    Cluster_OMem_Endcap_Ring_Fw[i]->setAxisTitle("Endcap Ring", 1);
    Cluster_OMem_Endcap_Ring_Fw[i]->setAxisTitle("# L1 Clusters ", 2);
  }

  for (int i = 0; i < numDiscs; i++) {
    HistoName = "NClusters_OMem_Disc-" + std::to_string(i + 1);
    Cluster_OMem_Endcap_Ring_Bw[i] = iBooker.book1D(HistoName,
                                                    HistoName,
                                                    psTTCluster_ECRing.getParameter<int32_t>("Nbinsx"),
                                                    psTTCluster_ECRing.getParameter<double>("xmin"),
                                                    psTTCluster_ECRing.getParameter<double>("xmax"));
    Cluster_OMem_Endcap_Ring_Bw[i]->setAxisTitle("Endcap Ring", 1);
    Cluster_OMem_Endcap_Ring_Bw[i]->setAxisTitle("# L1 Clusters ", 2);
  }

  iBooker.setCurrentFolder(topFolderName_ + "/Clusters");

  // Cluster Width
  edm::ParameterSet psTTClusterWidth = conf_.getParameter<edm::ParameterSet>("TH2TTCluster_Width");
  HistoName = "Cluster_W";
  Cluster_W = iBooker.book2D(HistoName,
                             HistoName,
                             psTTClusterWidth.getParameter<int32_t>("Nbinsx"),
                             psTTClusterWidth.getParameter<double>("xmin"),
                             psTTClusterWidth.getParameter<double>("xmax"),
                             psTTClusterWidth.getParameter<int32_t>("Nbinsy"),
                             psTTClusterWidth.getParameter<double>("ymin"),
                             psTTClusterWidth.getParameter<double>("ymax"));
  Cluster_W->setAxisTitle("L1 Cluster Width", 1);
  Cluster_W->setAxisTitle("Stack Member", 2);

  HistoName = "Cluster_IMem_W";
  Cluster_IMem_W = iBooker.book1D(HistoName,
                             HistoName,
                             psTTClusterWidth.getParameter<int32_t>("Nbinsx"),
                             psTTClusterWidth.getParameter<double>("xmin"),
                             psTTClusterWidth.getParameter<double>("xmax"));
  Cluster_IMem_W->setAxisTitle("L1 Cluster Width", 1);
  Cluster_IMem_W->setAxisTitle("# L1 Clusters", 2);

  HistoName = "Cluster_OMem_W";
  Cluster_OMem_W = iBooker.book1D(HistoName,
                             HistoName,
                             psTTClusterWidth.getParameter<int32_t>("Nbinsx"),
                             psTTClusterWidth.getParameter<double>("xmin"),
                             psTTClusterWidth.getParameter<double>("xmax"));
  Cluster_OMem_W->setAxisTitle("L1 Cluster Width", 1);
  Cluster_OMem_W->setAxisTitle("# L1 Clusters", 2);

  // Cluster eta distribution
  edm::ParameterSet psTTClusterEta = conf_.getParameter<edm::ParameterSet>("TH1TTCluster_Eta");
  HistoName = "Cluster_Eta";
  Cluster_Eta = iBooker.book1D(HistoName,
                               HistoName,
                               psTTClusterEta.getParameter<int32_t>("Nbinsx"),
                               psTTClusterEta.getParameter<double>("xmin"),
                               psTTClusterEta.getParameter<double>("xmax"));
  Cluster_Eta->setAxisTitle("#eta", 1);
  Cluster_Eta->setAxisTitle("# L1 Clusters ", 2);

  // Cluster phi distribution
  edm::ParameterSet psTTClusterPhi = conf_.getParameter<edm::ParameterSet>("TH1TTCluster_Phi");
  HistoName = "Cluster_Phi";
  Cluster_Phi = iBooker.book1D(HistoName,
                               HistoName,
                               psTTClusterPhi.getParameter<int32_t>("Nbinsx"),
                               psTTClusterPhi.getParameter<double>("xmin"),
                               psTTClusterPhi.getParameter<double>("xmax"));
  Cluster_Phi->setAxisTitle("#phi", 1);
  Cluster_Phi->setAxisTitle("# L1 Clusters", 2);

  // Cluster R distribution
  edm::ParameterSet psTTClusterR = conf_.getParameter<edm::ParameterSet>("TH1TTCluster_R");
  HistoName = "Cluster_R";
  Cluster_R = iBooker.book1D(HistoName,
                             HistoName,
                             psTTClusterR.getParameter<int32_t>("Nbinsx"),
                             psTTClusterR.getParameter<double>("xmin"),
                             psTTClusterR.getParameter<double>("xmax"));
  Cluster_R->setAxisTitle("R [cm]", 1);
  Cluster_R->setAxisTitle("# L1 Clusters", 2);

  // Cluster nearest distance distribution
  edm::ParameterSet psTTClusterNearest = conf_.getParameter<edm::ParameterSet>("TH1TTCluster_Nearest");
  HistoName = "Cluster_Nearest";
  Cluster_Nearest = iBooker.book1D(HistoName,
                             HistoName,
                             psTTClusterNearest.getParameter<int32_t>("Nbinsx"),
                             psTTClusterNearest.getParameter<double>("xmin"),
                             psTTClusterNearest.getParameter<double>("xmax"));
  Cluster_Nearest->setAxisTitle("Cluster Distance [cm]", 1);
  Cluster_Nearest->setAxisTitle("# L1 Clusters", 2);

  // PS Cluster nearest distance distribution
  edm::ParameterSet psTTClusterNearest_PS = conf_.getParameter<edm::ParameterSet>("TH1TTCluster_Nearest_PS");
  HistoName = "Cluster_Nearest_PS";
  Cluster_Nearest_PS = iBooker.book1D(HistoName,
                             HistoName,
                             psTTClusterNearest_PS.getParameter<int32_t>("Nbinsx"),
                             psTTClusterNearest_PS.getParameter<double>("xmin"),
                             psTTClusterNearest_PS.getParameter<double>("xmax"));
  Cluster_Nearest_PS->setAxisTitle("PS Cluster Distance [cm]", 1);
  Cluster_Nearest_PS->setAxisTitle("# L1 Clusters", 2);

  // 2S Cluster nearest distance distribution
  edm::ParameterSet psTTClusterNearest_2S = conf_.getParameter<edm::ParameterSet>("TH1TTCluster_Nearest_2S");
  HistoName = "Cluster_Nearest_2S";
  Cluster_Nearest_2S = iBooker.book1D(HistoName,
                             HistoName,
                             psTTClusterNearest_2S.getParameter<int32_t>("Nbinsx"),
                             psTTClusterNearest_2S.getParameter<double>("xmin"),
                             psTTClusterNearest_2S.getParameter<double>("xmax"));
  Cluster_Nearest_2S->setAxisTitle("2S Cluster Distance [cm]", 1);
  Cluster_Nearest_2S->setAxisTitle("# L1 Clusters", 2);

  iBooker.setCurrentFolder(topFolderName_ + "/Clusters/Position");

  // Position plots
  edm::ParameterSet psTTCluster_Barrel_XY = conf_.getParameter<edm::ParameterSet>("TH2TTCluster_Position");
  HistoName = "Cluster_Barrel_XY";
  Cluster_Barrel_XY = iBooker.book2D(HistoName,
                                     HistoName,
                                     psTTCluster_Barrel_XY.getParameter<int32_t>("Nbinsx"),
                                     psTTCluster_Barrel_XY.getParameter<double>("xmin"),
                                     psTTCluster_Barrel_XY.getParameter<double>("xmax"),
                                     psTTCluster_Barrel_XY.getParameter<int32_t>("Nbinsy"),
                                     psTTCluster_Barrel_XY.getParameter<double>("ymin"),
                                     psTTCluster_Barrel_XY.getParameter<double>("ymax"));
  Cluster_Barrel_XY->setAxisTitle("L1 Cluster Barrel position x [cm]", 1);
  Cluster_Barrel_XY->setAxisTitle("L1 Cluster Barrel position y [cm]", 2);

  edm::ParameterSet psTTCluster_Endcap_Fw_XY = conf_.getParameter<edm::ParameterSet>("TH2TTCluster_Position");
  HistoName = "Cluster_Endcap_Fw_XY";
  Cluster_Endcap_Fw_XY = iBooker.book2D(HistoName,
                                        HistoName,
                                        psTTCluster_Endcap_Fw_XY.getParameter<int32_t>("Nbinsx"),
                                        psTTCluster_Endcap_Fw_XY.getParameter<double>("xmin"),
                                        psTTCluster_Endcap_Fw_XY.getParameter<double>("xmax"),
                                        psTTCluster_Endcap_Fw_XY.getParameter<int32_t>("Nbinsy"),
                                        psTTCluster_Endcap_Fw_XY.getParameter<double>("ymin"),
                                        psTTCluster_Endcap_Fw_XY.getParameter<double>("ymax"));
  Cluster_Endcap_Fw_XY->setAxisTitle("L1 Cluster Forward Endcap position x [cm]", 1);
  Cluster_Endcap_Fw_XY->setAxisTitle("L1 Cluster Forward Endcap position y [cm]", 2);

  edm::ParameterSet psTTCluster_Endcap_Bw_XY = conf_.getParameter<edm::ParameterSet>("TH2TTCluster_Position");
  HistoName = "Cluster_Endcap_Bw_XY";
  Cluster_Endcap_Bw_XY = iBooker.book2D(HistoName,
                                        HistoName,
                                        psTTCluster_Endcap_Bw_XY.getParameter<int32_t>("Nbinsx"),
                                        psTTCluster_Endcap_Bw_XY.getParameter<double>("xmin"),
                                        psTTCluster_Endcap_Bw_XY.getParameter<double>("xmax"),
                                        psTTCluster_Endcap_Bw_XY.getParameter<int32_t>("Nbinsy"),
                                        psTTCluster_Endcap_Bw_XY.getParameter<double>("ymin"),
                                        psTTCluster_Endcap_Bw_XY.getParameter<double>("ymax"));
  Cluster_Endcap_Bw_XY->setAxisTitle("L1 Cluster Backward Endcap position x [cm]", 1);
  Cluster_Endcap_Bw_XY->setAxisTitle("L1 Cluster Backward Endcap position y [cm]", 2);

  edm::ParameterSet psTTCluster_Barrel_Local_PS_XY = conf_.getParameter<edm::ParameterSet>("TH2TTCluster_Local_Position");
  HistoName = "Cluster_Barrel_Local_PS_XY";
  Cluster_Barrel_Local_PS_XY = iBooker.book2D(HistoName,
                                     HistoName,
                                     psTTCluster_Barrel_Local_PS_XY.getParameter<int32_t>("Nbinsx"),
                                     psTTCluster_Barrel_Local_PS_XY.getParameter<double>("xmin"),
                                     psTTCluster_Barrel_Local_PS_XY.getParameter<double>("xmax"),
                                     psTTCluster_Barrel_Local_PS_XY.getParameter<int32_t>("Nbinsy"),
                                     psTTCluster_Barrel_Local_PS_XY.getParameter<double>("ymin"),
                                     psTTCluster_Barrel_Local_PS_XY.getParameter<double>("ymax"));
  Cluster_Barrel_Local_PS_XY->setAxisTitle("L1 PS Cluster Barrel local position x [cm]", 1);
  Cluster_Barrel_Local_PS_XY->setAxisTitle("L1 PS Cluster Barrel local position y [cm]", 2);

  edm::ParameterSet psTTCluster_Barrel_Local_2S_XY = conf_.getParameter<edm::ParameterSet>("TH2TTCluster_Local_2S_Position");
  HistoName = "Cluster_Barrel_Local_2S_XY";
  Cluster_Barrel_Local_2S_XY = iBooker.book1D(HistoName,
                                     HistoName,
                                     psTTCluster_Barrel_Local_2S_XY.getParameter<int32_t>("Nbinsx"),
                                     psTTCluster_Barrel_Local_2S_XY.getParameter<double>("xmin"),
                                     psTTCluster_Barrel_Local_2S_XY.getParameter<double>("xmax"));
  Cluster_Barrel_Local_2S_XY->setAxisTitle("L1 2S Cluster Barrel local position x [cm]", 1);
  Cluster_Barrel_Local_2S_XY->setAxisTitle("# L1 2S Clusters", 2);

  edm::ParameterSet psTTCluster_Barrel_Local_Pixel_XY = conf_.getParameter<edm::ParameterSet>("TH2TTCluster_Local_Position");
  HistoName = "Cluster_Barrel_Local_Pixel_XY";
  Cluster_Barrel_Local_Pixel_XY = iBooker.book2D(HistoName,
                                     HistoName,
                                     psTTCluster_Barrel_Local_Pixel_XY.getParameter<int32_t>("Nbinsx"),
                                     psTTCluster_Barrel_Local_Pixel_XY.getParameter<double>("xmin"),
                                     psTTCluster_Barrel_Local_Pixel_XY.getParameter<double>("xmax"),
                                     psTTCluster_Barrel_Local_Pixel_XY.getParameter<int32_t>("Nbinsy"),
                                     psTTCluster_Barrel_Local_Pixel_XY.getParameter<double>("ymin"),
                                     psTTCluster_Barrel_Local_Pixel_XY.getParameter<double>("ymax"));
  Cluster_Barrel_Local_Pixel_XY->setAxisTitle("L1 Pixel Cluster Barrel local position x [cm]", 1);
  Cluster_Barrel_Local_Pixel_XY->setAxisTitle("L1 Pixel Cluster Barrel local position y [cm]", 2);

  edm::ParameterSet psTTCluster_Barrel_Local_Strip_XY = conf_.getParameter<edm::ParameterSet>("TH2TTCluster_Local_2S_Position");
  HistoName = "Cluster_Barrel_Local_Strip_XY";
  Cluster_Barrel_Local_Strip_XY = iBooker.book1D(HistoName,
                                     HistoName,
                                     psTTCluster_Barrel_Local_Strip_XY.getParameter<int32_t>("Nbinsx"),
                                     psTTCluster_Barrel_Local_Strip_XY.getParameter<double>("xmin"),
                                     psTTCluster_Barrel_Local_Strip_XY.getParameter<double>("xmax"));
  Cluster_Barrel_Local_Strip_XY->setAxisTitle("L1 Strip Cluster Barrel local position x [cm]", 1);
  Cluster_Barrel_Local_Strip_XY->setAxisTitle("# L1 Strip Clusters", 2);

  edm::ParameterSet psTTCluster_Endcap_Fw_Local_PS_XY = conf_.getParameter<edm::ParameterSet>("TH2TTCluster_Local_Position");
  HistoName = "Cluster_Endcap_Fw_Local_PS_XY";
  Cluster_Endcap_Fw_Local_PS_XY = iBooker.book2D(HistoName,
                                        HistoName,
                                        psTTCluster_Endcap_Fw_Local_PS_XY.getParameter<int32_t>("Nbinsx"),
                                        psTTCluster_Endcap_Fw_Local_PS_XY.getParameter<double>("xmin"),
                                        psTTCluster_Endcap_Fw_Local_PS_XY.getParameter<double>("xmax"),
                                        psTTCluster_Endcap_Fw_Local_PS_XY.getParameter<int32_t>("Nbinsy"),
                                        psTTCluster_Endcap_Fw_Local_PS_XY.getParameter<double>("ymin"),
                                        psTTCluster_Endcap_Fw_Local_PS_XY.getParameter<double>("ymax"));
  Cluster_Endcap_Fw_Local_PS_XY->setAxisTitle("L1 PS Cluster Forward Endcap local position x [cm]", 1);
  Cluster_Endcap_Fw_Local_PS_XY->setAxisTitle("L1 PS Cluster Forward Endcap local position y [cm]", 2);

  edm::ParameterSet psTTCluster_Endcap_Fw_Local_2S_XY = conf_.getParameter<edm::ParameterSet>("TH2TTCluster_Local_2S_Position");
  HistoName = "Cluster_Endcap_Fw_Local_2S_XY";
  Cluster_Endcap_Fw_Local_2S_XY = iBooker.book1D(HistoName,
                                        HistoName,
                                        psTTCluster_Endcap_Fw_Local_2S_XY.getParameter<int32_t>("Nbinsx"),
                                        psTTCluster_Endcap_Fw_Local_2S_XY.getParameter<double>("xmin"),
                                        psTTCluster_Endcap_Fw_Local_2S_XY.getParameter<double>("xmax"));
  Cluster_Endcap_Fw_Local_2S_XY->setAxisTitle("L1 2S Cluster Forward Endcap local position x [cm]", 1);
  Cluster_Endcap_Fw_Local_2S_XY->setAxisTitle("# L1 2S Clusters", 2);

  edm::ParameterSet psTTCluster_Endcap_Fw_Local_Pixel_XY = conf_.getParameter<edm::ParameterSet>("TH2TTCluster_Local_Position");
  HistoName = "Cluster_Endcap_Fw_Local_Pixel_XY";
  Cluster_Endcap_Fw_Local_Pixel_XY = iBooker.book2D(HistoName,
                                        HistoName,
                                        psTTCluster_Endcap_Fw_Local_Pixel_XY.getParameter<int32_t>("Nbinsx"),
                                        psTTCluster_Endcap_Fw_Local_Pixel_XY.getParameter<double>("xmin"),
                                        psTTCluster_Endcap_Fw_Local_Pixel_XY.getParameter<double>("xmax"),
                                        psTTCluster_Endcap_Fw_Local_Pixel_XY.getParameter<int32_t>("Nbinsy"),
                                        psTTCluster_Endcap_Fw_Local_Pixel_XY.getParameter<double>("ymin"),
                                        psTTCluster_Endcap_Fw_Local_Pixel_XY.getParameter<double>("ymax"));
  Cluster_Endcap_Fw_Local_Pixel_XY->setAxisTitle("L1 Pixel Cluster Forward Endcap local position x [cm]", 1);
  Cluster_Endcap_Fw_Local_Pixel_XY->setAxisTitle("L1 Pixel Cluster Forward Endcap local position y [cm]", 2);

  edm::ParameterSet psTTCluster_Endcap_Fw_Local_Strip_XY = conf_.getParameter<edm::ParameterSet>("TH2TTCluster_Local_2S_Position");
  HistoName = "Cluster_Endcap_Fw_Local_Strip_XY";
  Cluster_Endcap_Fw_Local_Strip_XY = iBooker.book1D(HistoName,
                                        HistoName,
                                        psTTCluster_Endcap_Fw_Local_Strip_XY.getParameter<int32_t>("Nbinsx"),
                                        psTTCluster_Endcap_Fw_Local_Strip_XY.getParameter<double>("xmin"),
                                        psTTCluster_Endcap_Fw_Local_Strip_XY.getParameter<double>("xmax"));
  Cluster_Endcap_Fw_Local_Strip_XY->setAxisTitle("L1 Strip Cluster Forward Endcap local position x [cm]", 1);
  Cluster_Endcap_Fw_Local_Strip_XY->setAxisTitle("# L1 Strip Clusters", 2);

  edm::ParameterSet psTTCluster_Endcap_Bw_Local_PS_XY = conf_.getParameter<edm::ParameterSet>("TH2TTCluster_Local_Position");
  HistoName = "Cluster_Endcap_Bw_Local_PS_XY";
  Cluster_Endcap_Bw_Local_PS_XY = iBooker.book2D(HistoName,
                                        HistoName,
                                        psTTCluster_Endcap_Bw_Local_PS_XY.getParameter<int32_t>("Nbinsx"),
                                        psTTCluster_Endcap_Bw_Local_PS_XY.getParameter<double>("xmin"),
                                        psTTCluster_Endcap_Bw_Local_PS_XY.getParameter<double>("xmax"),
                                        psTTCluster_Endcap_Bw_Local_PS_XY.getParameter<int32_t>("Nbinsy"),
                                        psTTCluster_Endcap_Bw_Local_PS_XY.getParameter<double>("ymin"),
                                        psTTCluster_Endcap_Bw_Local_PS_XY.getParameter<double>("ymax"));
  Cluster_Endcap_Bw_Local_PS_XY->setAxisTitle("L1 PS Cluster Backward Endcap local position x [cm]", 1);
  Cluster_Endcap_Bw_Local_PS_XY->setAxisTitle("L1 PS Cluster Backward Endcap local position y [cm]", 2);

  edm::ParameterSet psTTCluster_Endcap_Bw_Local_2S_XY = conf_.getParameter<edm::ParameterSet>("TH2TTCluster_Local_2S_Position");
  HistoName = "Cluster_Endcap_Bw_Local_2S_XY";
  Cluster_Endcap_Bw_Local_2S_XY = iBooker.book1D(HistoName,
                                        HistoName,
                                        psTTCluster_Endcap_Bw_Local_2S_XY.getParameter<int32_t>("Nbinsx"),
                                        psTTCluster_Endcap_Bw_Local_2S_XY.getParameter<double>("xmin"),
                                        psTTCluster_Endcap_Bw_Local_2S_XY.getParameter<double>("xmax"));
  Cluster_Endcap_Bw_Local_2S_XY->setAxisTitle("L1 2S Cluster Backward Endcap local position x [cm]", 1);
  Cluster_Endcap_Bw_Local_2S_XY->setAxisTitle("# L1 2S Clusters", 2);

  edm::ParameterSet psTTCluster_Endcap_Bw_Local_Pixel_XY = conf_.getParameter<edm::ParameterSet>("TH2TTCluster_Local_Position");
  HistoName = "Cluster_Endcap_Bw_Local_Pixel_XY";
  Cluster_Endcap_Bw_Local_Pixel_XY = iBooker.book2D(HistoName,
                                        HistoName,
                                        psTTCluster_Endcap_Bw_Local_Pixel_XY.getParameter<int32_t>("Nbinsx"),
                                        psTTCluster_Endcap_Bw_Local_Pixel_XY.getParameter<double>("xmin"),
                                        psTTCluster_Endcap_Bw_Local_Pixel_XY.getParameter<double>("xmax"),
                                        psTTCluster_Endcap_Bw_Local_Pixel_XY.getParameter<int32_t>("Nbinsy"),
                                        psTTCluster_Endcap_Bw_Local_Pixel_XY.getParameter<double>("ymin"),
                                        psTTCluster_Endcap_Bw_Local_Pixel_XY.getParameter<double>("ymax"));
  Cluster_Endcap_Bw_Local_Pixel_XY->setAxisTitle("L1 Pixel Cluster Backward Endcap local position x [cm]", 1);
  Cluster_Endcap_Bw_Local_Pixel_XY->setAxisTitle("L1 Pixel Cluster Backward Endcap local position y [cm]", 2);

  edm::ParameterSet psTTCluster_Endcap_Bw_Local_Strip_XY = conf_.getParameter<edm::ParameterSet>("TH2TTCluster_Local_2S_Position");
  HistoName = "Cluster_Endcap_Bw_Local_Strip_XY";
  Cluster_Endcap_Bw_Local_Strip_XY = iBooker.book1D(HistoName,
                                        HistoName,
                                        psTTCluster_Endcap_Bw_Local_Strip_XY.getParameter<int32_t>("Nbinsx"),
                                        psTTCluster_Endcap_Bw_Local_Strip_XY.getParameter<double>("xmin"),
                                        psTTCluster_Endcap_Bw_Local_Strip_XY.getParameter<double>("xmax"));
  Cluster_Endcap_Bw_Local_Strip_XY->setAxisTitle("L1 Strip Cluster Backward Endcap local position x [cm]", 1);
  Cluster_Endcap_Bw_Local_Strip_XY->setAxisTitle("# L1 Strip Clusters", 2);

  // TTCluster #rho vs. z
  edm::ParameterSet psTTCluster_RZ = conf_.getParameter<edm::ParameterSet>("TH2TTCluster_RZ");
  HistoName = "Cluster_RZ";
  Cluster_RZ = iBooker.book2D(HistoName,
                              HistoName,
                              psTTCluster_RZ.getParameter<int32_t>("Nbinsx"),
                              psTTCluster_RZ.getParameter<double>("xmin"),
                              psTTCluster_RZ.getParameter<double>("xmax"),
                              psTTCluster_RZ.getParameter<int32_t>("Nbinsy"),
                              psTTCluster_RZ.getParameter<double>("ymin"),
                              psTTCluster_RZ.getParameter<double>("ymax"));
  Cluster_RZ->setAxisTitle("L1 Cluster position z [cm]", 1);
  Cluster_RZ->setAxisTitle("L1 Cluster position #rho [cm]", 2);

  // TTCluster Pixel y pos vs. width
  edm::ParameterSet psTTCluster_YvWidth_Pixel = conf_.getParameter<edm::ParameterSet>("TTCluster_YvWidth_Pixel");
  HistoName = "Cluster_YvWidth_Pixel";
  Cluster_YvWidth_Pixel = iBooker.book2D(HistoName,
                              HistoName,
                              psTTCluster_YvWidth_Pixel.getParameter<int32_t>("Nbinsx"),
                              psTTCluster_YvWidth_Pixel.getParameter<double>("xmin"),
                              psTTCluster_YvWidth_Pixel.getParameter<double>("xmax"),
                              psTTCluster_YvWidth_Pixel.getParameter<int32_t>("Nbinsy"),
                              psTTCluster_YvWidth_Pixel.getParameter<double>("ymin"),
                              psTTCluster_YvWidth_Pixel.getParameter<double>("ymax"));
  Cluster_YvWidth_Pixel->setAxisTitle("L1 Cluster position y [cm]", 1);
  Cluster_YvWidth_Pixel->setAxisTitle("L1 Cluster width", 2);

}  // end of method

DEFINE_FWK_MODULE(OuterTrackerMonitorTTCluster);
