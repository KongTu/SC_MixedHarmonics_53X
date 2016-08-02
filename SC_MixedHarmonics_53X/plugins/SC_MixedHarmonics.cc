// -*- C++ -*-
//
// Package:    SC_MixedHarmonics/SC_MixedHarmonics
// Class:      SC_MixedHarmonics
// 
/**\class SC_MixedHarmonics SC_MixedHarmonics.cc SC_MixedHarmonics/SC_MixedHarmonics/plugins/SC_MixedHarmonics.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Zhoudunming Tu
//         Created:  Mon, 01 Aug 2016 09:01:02 GMT
//
//


#include "SC_MixedHarmonics/SC_MixedHarmonics/interface/SC_MixedHarmonicsBase.h"


SC_MixedHarmonics::SC_MixedHarmonics(const edm::ParameterSet& iConfig)
{

  trackName_  =  iConfig.getParameter<edm::InputTag>("trackName");
  vertexName_ =  iConfig.getParameter<edm::InputTag>("vertexName");
  towerName_ =  iConfig.getParameter<edm::InputTag>("towerName");

  trackSrc_ = consumes<reco::TrackCollection>(trackName_);
  vertexSrc_ = consumes<reco::VertexCollection>(vertexName_);
  towerSrc_ = consumes<CaloTowerCollection>(towerName_);

  Nmin_ = iConfig.getUntrackedParameter<int>("Nmin");
  Nmax_ = iConfig.getUntrackedParameter<int>("Nmax");

  n1_ = iConfig.getUntrackedParameter<int>("n1");
  n2_ = iConfig.getUntrackedParameter<int>("n2");
  n3_ = iConfig.getUntrackedParameter<int>("n3");
  n4_ = iConfig.getUntrackedParameter<int>("n4");
  
  useCentrality_ = iConfig.getUntrackedParameter<bool>("useCentrality");
  reverseBeam_ = iConfig.getUntrackedParameter<bool>("reverseBeam");
  doEffCorrection_ = iConfig.getUntrackedParameter<bool>("doEffCorrection");
  useEtaGap_ = iConfig.getUntrackedParameter<bool>("useEtaGap");

  eff_ = iConfig.getUntrackedParameter<int>("eff");

  etaTracker_ = iConfig.getUntrackedParameter<double>("etaTracker");

  gapValue_ = iConfig.getUntrackedParameter<double>("gapValue");
  
  etaLowHF_ = iConfig.getUntrackedParameter<double>("etaLowHF");
  etaHighHF_ = iConfig.getUntrackedParameter<double>("etaHighHF");
  
  vzLow_ = iConfig.getUntrackedParameter<double>("vzLow");
  vzHigh_ = iConfig.getUntrackedParameter<double>("vzHigh");
  
  ptLow_ = iConfig.getUntrackedParameter<double>("ptLow");
  ptHigh_ = iConfig.getUntrackedParameter<double>("ptHigh");

  offlineptErr_ = iConfig.getUntrackedParameter<double>("offlineptErr", 0.0);
  offlineDCA_ = iConfig.getUntrackedParameter<double>("offlineDCA", 0.0);
  offlineChi2_ = iConfig.getUntrackedParameter<double>("offlineChi2", 0.0);
  offlinenhits_ = iConfig.getUntrackedParameter<double>("offlinenhits", 0.0);
  
  etaBins_ = iConfig.getUntrackedParameter<std::vector<double>>("etaBins");
  dEtaBins_ = iConfig.getUntrackedParameter<std::vector<double>>("dEtaBins");
  ptBins_ = iConfig.getUntrackedParameter<std::vector<double>>("ptBins");
  centBins_ = iConfig.getUntrackedParameter<std::vector<double>>("centBins");

}


SC_MixedHarmonics::~SC_MixedHarmonics()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SC_MixedHarmonics::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexSrc_,vertices);
  double bestvz=-999.9, bestvx=-999.9, bestvy=-999.9;
  double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
  const reco::Vertex & vtx = (*vertices)[0];
  bestvz = vtx.z(); 
  bestvx = vtx.x(); 
  bestvy = vtx.y();
  bestvzError = vtx.zError(); 
  bestvxError = vtx.xError(); 
  bestvyError = vtx.yError();

  //first selection; vertices
  if( fabs(bestvz) < vzLow_ || fabs(bestvz) > vzHigh_ ) return;

  vtxZ->Fill( bestvz );

  Handle<CaloTowerCollection> towers;
  iEvent.getByToken(towerSrc_, towers);

  Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(trackSrc_, tracks);

  int nTracks = 0;
  for(unsigned it = 0; it < tracks->size(); it++){

     const reco::Track & trk = (*tracks)[it];
  
     math::XYZPoint bestvtx(bestvx,bestvy,bestvz);

        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
        double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError); 

        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt() > 0.1 ) continue;
        if(fabs(dzvtx/dzerror) > 3.0) continue;
        if(fabs(dxyvtx/dxyerror) > 3.0) continue;
        if(trk.pt() < 0.4 || fabs(trk.eta()) > 2.4) continue;
        nTracks++;//count multiplicity

  }

  if( !useCentrality_ ) if( nTracks < Nmin_ || nTracks >= Nmax_ ) return;

  double etHFtowerSumPlus = 0.0;
  double etHFtowerSumMinus = 0.0;
  double etHFtowerSum = 0.0;
  
  if( useCentrality_ ){

    for( unsigned i = 0; i<towers->size(); ++ i){
       const CaloTower & tower = (*towers)[ i ];
       double eta = tower.eta();
       bool isHF = tower.ietaAbs() > 29;
          if(isHF && eta > 0){
            etHFtowerSumPlus += tower.pt();
          }
          if(isHF && eta < 0){
            etHFtowerSumMinus += tower.pt();
          }
    }
    etHFtowerSum=etHFtowerSumPlus + etHFtowerSumMinus;

    int bin = -1;
    for(int j=0; j<200; j++){
      if( etHFtowerSum >= centBins_[j] ){
         bin = j; break;
      }
    }

    int hiBin = bin;
    if( hiBin < Nmin_ || hiBin >= Nmax_ ) return;
    cbinHist->Fill( hiBin );

  }

  Ntrk->Fill( nTracks );

  const int NetaBins = etaBins_.size() - 1 ;

/*
The SC(m,n) = <<cos(m+n-m-n)>> - <<cos(m-m)>><<cos(n-n)>>
 */

/*
define all the ingredients for 2-, 3-, and 4-particle cumulants,
where Q_coefficient_power is used in the following names 
 */

//2-particle correlator

  TComplex Q_k1_1[NetaBins], Q_k2_1[NetaBins], Q_k1k2_2[NetaBins], Q_m1_1[NetaBins],
  Q_m2_1[NetaBins], Q_m1m2_2[NetaBins], Q_eta_0_1[NetaBins], Q_eta_0_2[NetaBins]; 

//4-particle correlator, including the ingredient for 3-particle correlator

  TComplex Q_n1_1, Q_n2_1, Q_n3_1, Q_n4_1;
  
  TComplex Q_n1n2_2, Q_n1n3_2, Q_n1n4_2, Q_n2n3_2, Q_n2n4_2, Q_n3n4_2;
  
  TComplex Q_n1n2n3_3, Q_n1n2n4_3, Q_n1n3n4_3, Q_n2n3n4_3;
  
  TComplex Q_n1n2n3n4_4;
  
  TComplex Q_0_1, Q_0_2, Q_0_3, Q_0_4;

//------------------------------------------------------------------

//Start filling Q-vectors;
  for(unsigned it = 0; it < tracks->size(); it++){

     const reco::Track & trk = (*tracks)[it];
  
     math::XYZPoint bestvtx(bestvx,bestvy,bestvz);

        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
        double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError); 
        double nhits = trk.numberOfValidHits();
        double chi2n = trk.normalizedChi2();
        double nlayers = trk.hitPattern().trackerLayersWithMeasurement();
        chi2n = chi2n/nlayers;
        double phi = trk.phi();
        double trkEta = trk.eta();

        double weight = 1.0;
        if( doEffCorrection_ ) { weight = 1.0/effTable[eff_]->GetBinContent( effTable[eff_]->FindBin(trk.eta(), trk.pt()) );}

        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt() > offlineptErr_ ) continue;
        if(fabs(dzvtx/dzerror) > offlineDCA_) continue;
        if(fabs(dxyvtx/dxyerror) > offlineDCA_) continue;
        if(chi2n > offlineChi2_ ) continue;
        if(nhits < offlinenhits_ ) continue;
        if(trk.pt() < ptLow_ || trk.pt() > ptHigh_ ) continue;
        if(fabs(trkEta) > etaTracker_ ) continue;

        for(int eta = 0; eta < NetaBins; eta++){
          if( trkEta > etaBins_[eta] && trkEta < etaBins_[eta+1] ){
          
            //for use of 2-particle:
            Q_k1_1[eta] += q_vector(n1_, 1, weight, phi);
            Q_k2_1[eta] += q_vector(-n1_, 1, weight, phi);
            Q_k1k2_2[eta] += q_vector(0.0, 2, weight, phi);

            Q_m1_1[eta] += q_vector(n2_, 1, weight, phi);
            Q_m2_1[eta] += q_vector(-n2_, 1, weight, phi);
            Q_m1m2_2[eta] += q_vector(0.0, 2, weight, phi);

            Q_eta_0_1[eta] += q_vector(0,1,weight,phi);
            Q_eta_0_2[eta] += q_vector(0,2,weight,phi);

          }
        }

        //for use of 4-particle:
        Q_n1_1 += q_vector(n1_, 1, weight, phi);
        Q_n2_1 += q_vector(n2_, 1, weight, phi);
        Q_n3_1 += q_vector(n3_, 1, weight, phi);
        Q_n4_1 += q_vector(n4_, 1, weight, phi);

        Q_n1n2_2 += q_vector(n1_+n2_, 2, weight, phi);
        Q_n1n3_2 += q_vector(n1_+n3_, 2, weight, phi);
        Q_n1n4_2 += q_vector(n1_+n4_, 2, weight, phi);
        Q_n2n3_2 += q_vector(n2_+n3_, 2, weight, phi);
        Q_n2n4_2 += q_vector(n2_+n4_, 2, weight, phi);
        Q_n3n4_2 += q_vector(n3_+n4_, 2, weight, phi);

        Q_n1n2n3_3 += q_vector(n1_+n2_+n3_, 3, weight, phi);
        Q_n1n2n4_3 += q_vector(n1_+n2_+n4_, 3, weight, phi);
        Q_n1n3n4_3 += q_vector(n1_+n3_+n4_, 3, weight, phi);
        Q_n2n3n4_3 += q_vector(n2_+n3_+n4_, 3, weight, phi);

        Q_n1n2n3n4_4 += q_vector(n1_+n2_+n3_+n4_, 4, weight, phi);
        
        Q_0_1 += q_vector(0,1,weight,phi);
        Q_0_2 += q_vector(0,2,weight,phi);
        Q_0_3 += q_vector(0,3,weight,phi);
        Q_0_4 += q_vector(0,4,weight,phi);

  }
    

/*
calculate 2-particle cumulant with a gap
 */

  for(int ieta = 0; ieta < NetaBins; ieta++){
    for(int jeta = 0; jeta < NetaBins; jeta++){

      double deltaEta = fabs( etaBins_[ieta] - etaBins_[jeta] );

      TComplex N_2_k;
      TComplex N_2_m;
      TComplex D_2;

      if( ieta == jeta ){

        N_2_k = Q_k1_1[ieta]*Q_k2_1[jeta] - Q_k1k2_2[ieta];
        N_2_m = Q_m1_1[ieta]*Q_m2_1[jeta] - Q_m1m2_2[ieta];
        D_2 = Q_eta_0_1[ieta]*Q_eta_0_1[jeta] - Q_eta_0_2[ieta];

      }
      else{
        
        N_2_k = Q_k1_1[ieta]*Q_k2_1[jeta];
        N_2_m = Q_m1_1[ieta]*Q_m2_1[jeta];
        D_2 = Q_eta_0_1[ieta]*Q_eta_0_1[jeta];
      }
      
      if( useEtaGap_ ){
        if(deltaEta < gapValue_) continue;

        c2_k->Fill(N_2_k.Re()/D_2.Re(), D_2.Re());
        c2_m->Fill(N_2_m.Re()/D_2.Re(), D_2.Re());

        c2_k_imag->Fill(N_2_k.Im()/D_2.Re(), D_2.Re());
        c2_m_imag->Fill(N_2_m.Im()/D_2.Re(), D_2.Re());
      }
      else{

        c2_k->Fill(N_2_k.Re()/D_2.Re(), D_2.Re());
        c2_m->Fill(N_2_m.Re()/D_2.Re(), D_2.Re());

        c2_k_imag->Fill(N_2_k.Im()/D_2.Re(), D_2.Re());
        c2_m_imag->Fill(N_2_m.Im()/D_2.Re(), D_2.Re());
        
      }
    }
  }

/*
calculate 3-particle cumulant
 */

  TComplex N_3, D_3;

  N_3 = Q_n1_1*Q_n2_1*Q_n3_1 - Q_n1n2_2*Q_n3_1 - Q_n2_1*Q_n1n3_2 - Q_n1_1*Q_n2n3_2 - Q_n1n2n3_3.operator*(2);
  D_3 = Q_0_1*Q_0_1*Q_0_1 - (Q_0_2*Q_0_1).operator*(3) + Q_0_3.operator*(2);

  c3->Fill( N_3.Re()/D_3.Re(), D_3.Re());
  c3_imag->Fill( N_3.Im()/D_3.Re(), D_3.Re());
/*
calculate 4-particle cumulant
 */  

  TComplex N_4, D_4;

  N_4 = Q_n1_1*Q_n2_1*Q_n3_1*Q_n4_1 - Q_n1n2_2*Q_n3_1*Q_n4_1 - Q_n2_1*Q_n1n3_2*Q_n4_1
      - Q_n1_1*Q_n2n3_2*Q_n4_1 + (Q_n1n2n3_3*Q_n4_1).operator*(2) - Q_n2_1*Q_n3_1*Q_n1n4_2
      + Q_n2n3_2*Q_n1n4_2 - Q_n1_1*Q_n3_1*Q_n2n4_2 + Q_n1n3_2*Q_n2n4_2
      + (Q_n3_1*Q_n1n2n4_3).operator*(2) - Q_n1_1*Q_n2_1*Q_n3n4_2 + Q_n1n2_2*Q_n3n4_2
      + (Q_n2_1*Q_n1n3n4_3).operator*(2) + (Q_n1_1*Q_n2n3n4_3).operator*(2) - Q_n1n2n3n4_4.operator*(6);

  D_4 = Q_0_1*Q_0_1*Q_0_1*Q_0_1 - (Q_0_1*Q_0_1*Q_0_2).operator*(6) + (Q_0_2*Q_0_2).operator*(3) + (Q_0_1*Q_0_3).operator*(8) - Q_0_4.operator*(6);
 
  c4->Fill( N_4.Re()/D_4.Re(), D_4.Re());
  c4_imag->Fill( N_4.Im()/D_4.Re(), D_4.Re());


}
// ------------ method called once each job just before starting event loop  ------------
void 
SC_MixedHarmonics::beginJob()
{
  edm::Service<TFileService> fs;
    
  TH1D::SetDefaultSumw2();

  edm::FileInPath fip1("SC_MixedHarmonics/SC_MixedHarmonics/data/Hydjet_eff_mult_v1.root");
  TFile f1(fip1.fullPath().c_str(),"READ");
  for(int i = 0; i < 5; i++){
     effTable[i] = (TH2D*)f1.Get(Form("rTotalEff3D_%d",i));
  }

  Ntrk = fs->make<TH1D>("Ntrk",";Ntrk",5000,0,5000);
  vtxZ = fs->make<TH1D>("vtxZ",";vz", 400,-20,20);
  cbinHist = fs->make<TH1D>("cbinHist",";cbin",200,0,200);

  c2_k = fs->make<TH1D>("c2_k",";c2", 20000,-1,1);
  c2_m = fs->make<TH1D>("c2_m",";c2", 20000,-1,1);
  c3 = fs->make<TH1D>("c3",";c3", 20000,-1,1);
  c4 = fs->make<TH1D>("c4",";c4", 20000,-1,1);

  c2_k_imag = fs->make<TH1D>("c2_k_imag",";c2", 20000,-1,1);
  c2_m_imag = fs->make<TH1D>("c2_m_imag",";c2", 20000,-1,1);
  c3_imag = fs->make<TH1D>("c3_imag",";c3", 20000,-1,1);
  c4_imag = fs->make<TH1D>("c4_imag",";c4", 20000,-1,1);
}

TComplex 
SC_MixedHarmonics::q_vector(double n, double p, double w, double phi) 
{
  double term1 = pow(w,p);
  TComplex e(1, n*phi, 1);
  return term1*e;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SC_MixedHarmonics::endJob() 
{
}
void 
SC_MixedHarmonics::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
SC_MixedHarmonics::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
SC_MixedHarmonics::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
SC_MixedHarmonics::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SC_MixedHarmonics::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
}

//define this as a plug-in
DEFINE_FWK_MODULE(SC_MixedHarmonics);
