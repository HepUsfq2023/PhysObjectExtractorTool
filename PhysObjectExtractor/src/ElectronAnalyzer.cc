// -*- C++ -*-
//
// Package:    Electron/ElectronAnalyzer
// Class:      ElectronAnalyzer
// PROYECTO FINAL

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/ValueMap.h"

//class to extract electron information
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

//class to extract MET information
#include "DataFormats/PatCandidates/interface/MET.h"

//class to extract jet information
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

//Transient track for impact parameter
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include<vector>
#include "Math/Vector4D.h"
#include "TRandom3.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


//using reco::TrackCollection;

class ElectronAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ElectronAnalyzer(const edm::ParameterSet&);
      ~ElectronAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup& ) override;
      virtual void endJob() override;

      edm::EDGetTokenT<pat::ElectronCollection> electronToken_, electronToken2_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      
      //------------To extract Met's Information----------
      edm::EDGetTokenT<pat::METCollection> metToken_;
      edm::EDGetTokenT<pat::METCollection> rawToken_;
      
      //declare the input tag for PFJetCollection
      edm::EDGetTokenT<pat::JetCollection> jetToken_;

      // ----------member data ---------------------------

      TTree *mtree;
      int numelectron; //number of electrons in the event
      std::vector<float> electron_e;
      std::vector<float> electron_pt;
      std::vector<float> electron_px;
      std::vector<float> electron_py;
      std::vector<float> electron_pz;
      std::vector<float> electron_eta;
      std::vector<float> electron_phi;
      std::vector<float> electron_ch;
      std::vector<float> electron_iso;
      std::vector<bool> electron_veto;//
      std::vector<bool> electron_isLoose;
      std::vector<bool> electron_isMedium;
      std::vector<bool> electron_isTight;
      std::vector<float> electron_dxy;
      std::vector<float> electron_dz;
      std::vector<float> electron_dxyError;
      std::vector<float> electron_dzError;
      std::vector<int> electron_ismvaLoose;
      std::vector<int> electron_ismvaTight;
      std::vector<double> electron_ip3d;	
      std::vector<double> electron_sip3d;
      float electron_tmassmin;
      float Ht;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ElectronAnalyzer::ElectronAnalyzer(const edm::ParameterSet& iConfig):
 electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
 vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
 metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
 jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets")))
{
   //now do what ever initialization is needed

   edm::Service<TFileService> fs;
   mtree = fs->make<TTree>("Events", "Events");
  
  mtree->Branch("numberelectron",&numelectron);   
  mtree->GetBranch("numberelectron")->SetTitle("number of electrons");
  mtree->Branch("electron_e",&electron_e);
  mtree->GetBranch("electron_e")->SetTitle("electron energy");
  mtree->Branch("electron_pt",&electron_pt);
  mtree->GetBranch("electron_pt")->SetTitle("electron transverse momentum");
  mtree->Branch("electron_px",&electron_px);
  mtree->GetBranch("electron_px")->SetTitle("electron momentum x-component");
  mtree->Branch("electron_py",&electron_py);
  mtree->GetBranch("electron_py")->SetTitle("electron momentum y-component");
  mtree->Branch("electron_pz",&electron_pz);
  mtree->GetBranch("electron_pz")->SetTitle("electron momentum z-component");
  mtree->Branch("electron_eta",&electron_eta);
  mtree->GetBranch("electron_eta")->SetTitle("electron pseudorapidity");
  mtree->Branch("electron_phi",&electron_phi);
  mtree->GetBranch("electron_phi")->SetTitle("electron polar angle");
  mtree->Branch("electron_ch",&electron_ch);
  mtree->GetBranch("electron_ch")->SetTitle("electron charge");
  mtree->Branch("electron_iso",&electron_iso);
  mtree->GetBranch("electron_iso")->SetTitle("electron isolation");
  mtree->Branch("electron_veto",&electron_veto);//
  mtree->GetBranch("electron_veto")->SetTitle("electron veto");//
  mtree->Branch("electron_isLoose",&electron_isLoose);
  mtree->GetBranch("electron_isLoose")->SetTitle("electron tagged loose");
  mtree->Branch("electron_isMedium",&electron_isMedium);
  mtree->GetBranch("electron_isMedium")->SetTitle("electron tagged medium");
  mtree->Branch("electron_isTight",&electron_isTight);
  mtree->GetBranch("electron_isTight")->SetTitle("electron tagged tight");
  mtree->Branch("electron_dxy",&electron_dxy);
  mtree->GetBranch("electron_dxy")->SetTitle("electron transverse plane impact parameter (mm)");
  mtree->Branch("electron_dz",&electron_dz);
  mtree->GetBranch("electron_dz")->SetTitle("electron longitudinal impact parameter (mm)");
  mtree->Branch("electron_dxyError",&electron_dxyError);
  mtree->GetBranch("electron_dxyError")->SetTitle("electron transverse impact parameter uncertainty (mm)");
  mtree->Branch("electron_dzError",&electron_dzError);
  mtree->GetBranch("electron_dzError")->SetTitle("electron longitudinal impact parameter uncertainty (mm)");
  mtree->Branch("electron_ismvaLoose",&electron_ismvaLoose);
  mtree->GetBranch("electron_ismvaLoose")->SetTitle("electron mva Loose");
  mtree->Branch("electron_ismvaTight",&electron_ismvaTight);
  mtree->GetBranch("electron_ismvaTight")->SetTitle("electron mva Tight");
  mtree->Branch("electron_ip3d",&electron_ip3d);
  mtree->GetBranch("electron_ip3d")->SetTitle("electron impact parameter in 3d");
  mtree->Branch("electron_sip3d",&electron_sip3d);
  mtree->GetBranch("electron_sip3d")->SetTitle("electron significance on impact parameter in 3d");
  mtree->Branch("electron_tmassmin",&electron_tmassmin);
  mtree->GetBranch("electron_tmassmin")->SetTitle("electron mt min");
  mtree->Branch("Ht",&Ht);
  mtree->GetBranch("Ht")->SetTitle("Ht");
}

//Destructor
ElectronAnalyzer::~ElectronAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ElectronAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
   using namespace edm;

   Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_, electrons);

   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   const reco::Vertex &primaryVertex = vertices->front();
   math::XYZPoint pv(vertices->begin()->position());

   numelectron = 0;
   electron_e.clear();
   electron_pt.clear();
   electron_px.clear();
   electron_py.clear();
   electron_pz.clear();
   electron_eta.clear();
   electron_phi.clear();
   electron_ch.clear();
   electron_iso.clear();
   electron_veto.clear();
   electron_isLoose.clear();
   electron_isMedium.clear();
   electron_isTight.clear();
   electron_dxy.clear();
   electron_dz.clear();
   electron_dxyError.clear();
   electron_dzError.clear();
   electron_ismvaLoose.clear();
   electron_ismvaTight.clear();
   electron_ip3d.clear();
   electron_sip3d.clear();
   Ht = 0;
   
   int cutpt1 = 17;
   int cutpt2 = 12;
   int index1 = -1;
   int index2 = -1;
   int i1 = -1;
   int i2 = -1;
   int cpt1 = 8;
   int cpt2 = 8;
   int numjets = 0;
   

    for (const pat::Electron &el : *electrons)
    {
      electron_e.push_back(el.energy());
      electron_pt.push_back(el.pt());
      electron_px.push_back(el.px());
      electron_py.push_back(el.py());
      electron_pz.push_back(el.pz());
      electron_eta.push_back(el.eta());
      electron_phi.push_back(el.phi());
      electron_ch.push_back(el.charge());
      electron_iso.push_back(el.ecalPFClusterIso());
      electron_veto.push_back(el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-veto"));//
      electron_isLoose.push_back(el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-loose"));
      electron_isMedium.push_back(el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-medium"));
      electron_isTight.push_back(el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-tight"));
      electron_dxy.push_back(el.gsfTrack()->dxy(pv));
      electron_dz.push_back(el.gsfTrack()->dz(pv));
      electron_dxyError.push_back(el.gsfTrack()->d0Error());
      electron_dzError.push_back(el.gsfTrack()->dzError());
      electron_ismvaLoose.push_back(el.electronID("mvaEleID-Spring15-25ns-nonTrig-V1-wp90"));
      electron_ismvaTight.push_back(el.electronID("mvaEleID-Spring15-25ns-nonTrig-V1-wp80"));

      //get impact parameter in 3D
      // https://github.com/cms-sw/cmssw/blob/CMSSW_7_6_X/PhysicsTools/PatAlgos/plugins/PATElectronProducer.cc
      // This is needed by the IPTools methods from the tracking group
      edm::ESHandle<TransientTrackBuilder> trackBuilder;
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder);      
      reco::TransientTrack tt = trackBuilder->build(el.gsfTrack());
      std::pair<bool,Measurement1D> ip3dpv = IPTools::absoluteImpactParameter3D(tt, primaryVertex);
      electron_ip3d.push_back(ip3dpv.second.value());
      electron_sip3d.push_back(ip3dpv.second.significance());
      //std::cout<<"ip3d vanilla = "<<el.ip3d()<<"\t ip3d from iptools = "<<ip3dpv.second.value()<<std::endl;
      
      
      //-------------------Preselección 1 de Electrones de Interés-------------------------------------------------
      if(el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-loose")){
	 if(electron_pt.at(numelectron)>cutpt1){
	   cutpt1 = electron_pt.at(numelectron);
	   index1 = numelectron;
	 
	 }
	 
	 else if(electron_pt.at(numelectron)>cutpt2 && index1 != index2){
	   cutpt2 = electron_pt.at(numelectron);
	   index2 = numelectron;
	
	 }

       }
       
      if(index1 == -1 && index2 == -1){
        if(electron_pt.at(numelectron)>cpt1){
	   cpt1 = electron_pt.at(numelectron);
	   i1 = numelectron;
	 
	 }
	 
	 else if(electron_pt.at(numelectron)>cpt2 && i1 != i2){
	   cpt2 = electron_pt.at(numelectron);
	   i2 = numelectron;
	
	 }
    
       }
      
      numelectron++;
    }
    
    //-------------------- Momento transversal y ángulo phi para el MET-----------
    Handle<pat::METCollection> mets;
    iEvent.getByToken(metToken_, mets);
    
    const pat::MET &met = mets->front();
    
    float met_pt = met.pt();   
    float met_phi = met.phi();
    
    //----------------Cálculo con el primer criterio de selección-------------------
    if(index1 != -1 && index2 != -1){
    //---------------Cálculo de Delta Phi 1-----------------------------------------
    float dp1 = std::fmod(electron_pt.at(index1) - met_phi, 2.0 * M_PI);
    if (dp1 < -M_PI) {
      dp1 += 2.0 * M_PI;
    }
    else if (dp1 > M_PI) {
      dp1 -= 2.0 * M_PI;
    }
    
    //--------------Cálculo de mt 1-------------------------------------------------
    float electron_tmass1 = sqrt(2.0 * electron_pt.at(index1) * met_pt * (1.0 - cos(dp1)));
    
    //---------------Cálculo de Delta Phi 2-----------------------------------------
    float dp2 = std::fmod(electron_pt.at(index2) - met_phi, 2.0 * M_PI);
    if (dp2 < -M_PI) {
      dp2 += 2.0 * M_PI;
    }
    else if (dp2 > M_PI) {
      dp2 -= 2.0 * M_PI;
    }
    
    //--------------Cálculo de mt 2-------------------------------------------------
    float electron_tmass2 = sqrt(2.0 * electron_pt.at(index2) * met_pt * (1.0 - cos(dp2)));
    
    //--------------Cálculo de mt min-----------------------------------------------
    
    if(electron_tmass1 > electron_tmass2){
      electron_tmassmin = electron_tmass2;
    
    }
    
    if(electron_tmass1 < electron_tmass2){
      electron_tmassmin = electron_tmass1;
    
    }
    
    }
    
    
    
    //---------------Cálculo con el segundo criterio de selección-------------------
    
    if(index1 == -1 && index2 == -1 && i1 != -1 && i2 != -1){
    //---------------Definimos los Jets---------------------------------------------
    Handle<pat::JetCollection> jets;
    iEvent.getByToken(jetToken_, jets);
    
    if(jets.isValid()){
    for (const pat::Jet &jet : *jets){
      pat::Jet uncorrJet = jet.correctedJet(0);
      if(uncorrJet.pt() > 40 && uncorrJet.eta()< 2.4){
        Ht = Ht + uncorrJet.pt();
        numjets++;
      }
    }
    
    }
    
    if(Ht> 300 && numjets>=2){
      float dp1 = std::fmod(electron_pt.at(i1) - met_phi, 2.0 * M_PI);
    if (dp1 < -M_PI) {
      dp1 += 2.0 * M_PI;
    }
    else if (dp1 > M_PI) {
      dp1 -= 2.0 * M_PI;
    }
    
    //--------------Cálculo de mt 1-------------------------------------------------
    float electron_tmass1 = sqrt(2.0 * electron_pt.at(i1) * met_pt * (1.0 - cos(dp1)));
    
    //---------------Cálculo de Delta Phi 2-----------------------------------------
    float dp2 = std::fmod(electron_pt.at(i2) - met_phi, 2.0 * M_PI);
    if (dp2 < -M_PI) {
      dp2 += 2.0 * M_PI;
    }
    else if (dp2 > M_PI) {
      dp2 -= 2.0 * M_PI;
    }
    
    //--------------Cálculo de mt 2-------------------------------------------------
    float electron_tmass2 = sqrt(2.0 * electron_pt.at(i2) * met_pt * (1.0 - cos(dp2)));
    
    //--------------Cálculo de mt min-----------------------------------------------
    
    if(electron_tmass1 > electron_tmass2){
      electron_tmassmin = electron_tmass2;
    
    }
    
    if(electron_tmass1 < electron_tmass2){
      electron_tmassmin = electron_tmass1;
    
    }
    
    }
    
    }
    
    else{
        electron_tmassmin = 0;
    }

  mtree->Fill();
  return;

}


// ------------ method called once each job just before starting event loop  ------------
void
ElectronAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
ElectronAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronAnalyzer);
