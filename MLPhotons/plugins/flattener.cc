// -*- C++ -*-
//
// Package:    ML_Photons/flattener
// Class:      flattener
// 
/**\class flattener flattener.cc ML_Photons/flattener/plugins/flattener.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Steven Clark
//         Created:  Tue, 30 Mar 2021 17:29:55 GMT
//
//

#include "flattener.h"
//

#include "FWCore/Utilities/interface/InputTag.h"

using namespace edm;
using namespace fastjet;
using namespace fastjet::contrib;

class flattener : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns> {
   public:
      explicit flattener( const edm::ParameterSet& );
      ~flattener() override;

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze( const edm::Event & event, const edm::EventSetup & ) override;
      virtual void endJob() override;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void clearVars();

      // ----------member data ---------------------------
      const edm::InputTag triggerResultsTag;
      const edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken;

      //Gen Particles
      const edm::InputTag genpartTag;
      const edm::EDGetTokenT<vector<reco::GenParticle>> genpartToken;

      //Jets
      const edm::InputTag patjetTag;
      const edm::EDGetTokenT<vector<pat::Jet>> patjetToken;

      //PF Candidates
      const edm::InputTag pfcandTag;
      const edm::EDGetTokenT<vector<pat::PackedCandidate>> pfcandToken;

      //METs
      const edm::InputTag metTag;
      const edm::EDGetTokenT<vector<pat::MET>> metToken;

      //Muons
      const edm::InputTag muonTag;
      const edm::EDGetTokenT<vector<pat::Muon>> muonToken;

      //Electrons
      const edm::InputTag electronTag;
      const edm::EDGetTokenT<vector<pat::Electron>> electronToken;
      
      //PAT Photons
      const edm::InputTag patPhoTag;
      const edm::EDGetTokenT<vector<pat::Photon>> patPhoToken;

      //Primary Vertex
      const edm::InputTag pvtxTag;
      const edm::EDGetTokenT<vector<reco::Vertex>> pvtxToken;

      //Primary Vertex
      const edm::InputTag svtxTag;
      const edm::EDGetTokenT<vector<reco::VertexCompositePtrCandidate>> svtxToken;

      //Trigger file
      std::string trigfile;
      
      //bTagger Name
      std::string btagName;
      
      //RUCLU
      const edm::InputTag ruclu_etaTag;
      const edm::InputTag ruclu_phiTag;
      const edm::InputTag ruclu_energyTag;
      const edm::InputTag ruclu_r1Tag;
      const edm::InputTag ruclu_r2Tag;
      const edm::InputTag ruclu_r3Tag;
      const edm::InputTag monophoTag;
      const edm::InputTag diphoTag;
      const edm::InputTag hadronTag;
      const edm::InputTag moeTag;

      const edm::EDGetTokenT<std::vector<float>> ruclu_etaToken;
      const edm::EDGetTokenT<std::vector<float>> ruclu_phiToken;
      const edm::EDGetTokenT<std::vector<float>> ruclu_energyToken;
      const edm::EDGetTokenT<std::vector<float>> ruclu_r1Token;
      const edm::EDGetTokenT<std::vector<float>> ruclu_r2Token;
      const edm::EDGetTokenT<std::vector<float>> ruclu_r3Token;
      const edm::EDGetTokenT<std::vector<float>> monophoToken;
      const edm::EDGetTokenT<std::vector<float>> diphoToken; //dcom
      const edm::EDGetTokenT<std::vector<float>> hadronToken;
      const edm::EDGetTokenT<std::vector<float>> moeToken;

      TTree* tree;

      long run;
      long id;
      int lumiSec;
      double wgt;
      bool isMC;

      std::vector<int> trigs;

      std::vector<float> jet_pt;
      std::vector<float> jet_eta;
      std::vector<float> jet_phi;
      std::vector<float> jet_energy;
      std::vector<float> jet_mass;
      std::vector<float> jet_btag;

      std::vector<float> met_pt;
      std::vector<float> met_eta;
      std::vector<float> met_phi;
      std::vector<float> met_energy;
      std::vector<float> met_mass;

      std::vector<float> muon_pt;
      std::vector<float> muon_eta;
      std::vector<float> muon_phi;
      std::vector<float> muon_energy;
      std::vector<float> muon_mass;

      std::vector<float> electron_pt;
      std::vector<float> electron_eta;
      std::vector<float> electron_phi;
      std::vector<float> electron_energy;
      std::vector<float> electron_mass;

      std::vector<float> patpho_pt;
      std::vector<float> patpho_eta;
      std::vector<float> patpho_phi;
      std::vector<float> patpho_energy;
      std::vector<float> patpho_mass;
      std::vector<float> patpho_r9;
      std::vector<float> patpho_sieie;
      std::vector<float> patpho_hoe;
      std::vector<float> patpho_chargedHadronIso;
      std::vector<float> patpho_neutralHadronIso;
      std::vector<float> patpho_photonIso;
      std::vector<int> patpho_cutBased;
      std::vector<float> patpho_mvaID;
      std::vector<bool> patpho_hasPixelSeed;
      std::vector<bool> patpho_passElectronVeto;

      std::vector<float> genpart_pt;
      std::vector<float> genpart_eta;
      std::vector<float> genpart_phi;
      std::vector<float> genpart_energy;
      std::vector<float> genpart_mass;
      std::vector<int> genpart_pdgid;
      std::vector<int> genpart_status;
      std::vector<int> genpart_motherpdgid;

      std::vector<float> pvtx_x;
      std::vector<float> pvtx_y;
      std::vector<float> pvtx_z;
      std::vector<float> pvtx_chi2;
      std::vector<float> pvtx_ndof;
      std::vector<int> pvtx_size;

      std::vector<float> svtx_x;
      std::vector<float> svtx_y;
      std::vector<float> svtx_z;
      std::vector<float> svtx_chi2;
      std::vector<float> svtx_ndof;

      std::vector<float> ruclu_eta;
      std::vector<float> ruclu_phi;
      std::vector<float> ruclu_energy;
      std::vector<float> ruclu_r1;
      std::vector<float> ruclu_r2;
      std::vector<float> ruclu_r3;
      std::vector<float> ruclu_pfIso;
      std::vector<float> ruclu_dipho;
      std::vector<float> ruclu_monopho;
      std::vector<float> ruclu_hadron;
      std::vector<float> ruclu_moe;
      std::vector<bool> ruclu_isPhoton;
      std::vector<bool> ruclu_hasPixelSeed;
      std::vector<bool> ruclu_passElectronVeto;
  };

flattener::flattener(const edm::ParameterSet& iConfig):
  triggerResultsTag (iConfig.getParameter<edm::InputTag>("TriggerInputTag_HLT")),
  triggerResultsToken (consumes<edm::TriggerResults>(triggerResultsTag)),

  genpartTag (iConfig.getParameter<edm::InputTag>("genpartInputTag")),
  genpartToken (consumes<vector<reco::GenParticle>>(genpartTag)),

  patjetTag (iConfig.getParameter<edm::InputTag>("patjetInputTag")),
  patjetToken (consumes<vector<pat::Jet>>(patjetTag)),

  pfcandTag (iConfig.getParameter<edm::InputTag>("pfcandInputTag")),
  pfcandToken (consumes<vector<pat::PackedCandidate>>(pfcandTag)),

  metTag (iConfig.getParameter<edm::InputTag>("metInputTag")),
  metToken (consumes<vector<pat::MET>>(metTag)),

  muonTag (iConfig.getParameter<edm::InputTag>("muonInputTag")),
  muonToken (consumes<vector<pat::Muon>>(muonTag)),

  electronTag (iConfig.getParameter<edm::InputTag>("electronInputTag")),
  electronToken (consumes<vector<pat::Electron>>(electronTag)),

  patPhoTag (iConfig.getParameter<edm::InputTag>("patPhoInputTag")),
  patPhoToken (consumes<vector<pat::Photon>>(patPhoTag)),

  pvtxTag (iConfig.getParameter<edm::InputTag>("pvtxInputTag")),
  pvtxToken (consumes<vector<reco::Vertex>>(pvtxTag)),

  svtxTag (iConfig.getParameter<edm::InputTag>("svtxInputTag")),
  svtxToken (consumes<vector<reco::VertexCompositePtrCandidate>>(svtxTag)),

  //Trigger file
  trigfile(iConfig.getParameter<std::string>("tfile_path")),
  
  //Trigger file
  btagName(iConfig.getParameter<std::string>("btag_name")),

  //RUCLU Stuff
  ruclu_etaTag (iConfig.getParameter<edm::InputTag>("ruclu_etaTag")),
  ruclu_phiTag (iConfig.getParameter<edm::InputTag>("ruclu_phiTag")),
  ruclu_energyTag (iConfig.getParameter<edm::InputTag>("ruclu_energyTag")),
  ruclu_r1Tag (iConfig.getParameter<edm::InputTag>("ruclu_r1Tag")),
  ruclu_r2Tag (iConfig.getParameter<edm::InputTag>("ruclu_r2Tag")),
  ruclu_r3Tag (iConfig.getParameter<edm::InputTag>("ruclu_r3Tag")),
  monophoTag (iConfig.getParameter<edm::InputTag>("monophoInputTag")),
  diphoTag (iConfig.getParameter<edm::InputTag>("diphoInputTag")), //dcom
  hadronTag (iConfig.getParameter<edm::InputTag>("hadronInputTag")),
  moeTag (iConfig.getParameter<edm::InputTag>("moeInputTag")),

  ruclu_etaToken (consumes<std::vector<float>>(ruclu_etaTag)),
  ruclu_phiToken (consumes<std::vector<float>>(ruclu_phiTag)),
  ruclu_energyToken (consumes<std::vector<float>>(ruclu_energyTag)),
  ruclu_r1Token (consumes<std::vector<float>>(ruclu_r1Tag)),
  ruclu_r2Token (consumes<std::vector<float>>(ruclu_r2Tag)),
  ruclu_r3Token (consumes<std::vector<float>>(ruclu_r3Tag)),
  monophoToken (consumes<std::vector<float>>(monophoTag)),
  diphoToken (consumes<std::vector<float>>(diphoTag)), //dcom
  hadronToken (consumes<std::vector<float>>(hadronTag)),
  moeToken (consumes<std::vector<float>>(moeTag)),

  wgt (iConfig.getParameter<double>("weightInput"))

{

  isMC = iConfig.getParameter<bool>("isMC");
  
  // Access the TFileService
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree", "tree");

  tree->Branch("lumiSec", &lumiSec, "lumiSec/i");
  tree->Branch("run", &run, "run/i");
  tree->Branch("id", &id, "id/i");

  tree->Branch("wgt", &wgt );

  tree->Branch("triggers", &trigs);

  tree->Branch("jet_pt", &jet_pt );
  tree->Branch("jet_eta", &jet_eta );
  tree->Branch("jet_phi", &jet_phi );
  tree->Branch("jet_energy", &jet_energy );
  tree->Branch("jet_mass", &jet_mass );
  tree->Branch("jet_btag", &jet_btag );

  tree->Branch("met_pt", &met_pt );
  tree->Branch("met_eta", &met_eta );
  tree->Branch("met_phi", &met_phi );
  tree->Branch("met_energy", &met_energy );
  tree->Branch("met_mass", &met_mass );

  tree->Branch("muon_pt", &muon_pt );
  tree->Branch("muon_eta", &muon_eta );
  tree->Branch("muon_phi", &muon_phi );
  tree->Branch("muon_energy", &muon_energy );
  tree->Branch("muon_mass", &muon_mass );

  tree->Branch("electron_pt", &electron_pt );
  tree->Branch("electron_eta", &electron_eta );
  tree->Branch("electron_phi", &electron_phi );
  tree->Branch("electron_energy", &electron_energy );
  tree->Branch("electron_mass", &electron_mass );

  tree->Branch("patpho_pt", &patpho_pt );
  tree->Branch("patpho_eta", &patpho_eta );
  tree->Branch("patpho_phi", &patpho_phi );
  tree->Branch("patpho_energy", &patpho_energy );
  tree->Branch("patpho_mass", &patpho_mass );
  tree->Branch("patpho_r9", &patpho_r9 );
  tree->Branch("patpho_sieie", &patpho_sieie );
  tree->Branch("patpho_hoe", &patpho_hoe );
  tree->Branch("patpho_chargedHadronIso", &patpho_chargedHadronIso );
  tree->Branch("patpho_neutralHadronIso", &patpho_neutralHadronIso );
  tree->Branch("patpho_photonIso", &patpho_photonIso );
  tree->Branch("patpho_cutBased", &patpho_cutBased );
  tree->Branch("patpho_mvaID", &patpho_mvaID );
  tree->Branch("patpho_hasPixelSeed", &patpho_hasPixelSeed);
  tree->Branch("patpho_passElectronVeto", &patpho_passElectronVeto);
  
  if (isMC){
    tree->Branch("genpart_pt", &genpart_pt);
    tree->Branch("genpart_eta", &genpart_eta);
    tree->Branch("genpart_phi", &genpart_phi);
    tree->Branch("genpart_energy", &genpart_energy);
    tree->Branch("genpart_mass", &genpart_mass);
    tree->Branch("genpart_pdgid", &genpart_pdgid);
    tree->Branch("genpart_status", &genpart_status);
    tree->Branch("genpart_motherpdgid", &genpart_motherpdgid);
  }

  tree->Branch("pvtx_x", &pvtx_x );
  tree->Branch("pvtx_y", &pvtx_y );
  tree->Branch("pvtx_z", &pvtx_z );
  tree->Branch("pvtx_chi2", &pvtx_chi2 );
  tree->Branch("pvtx_ndof", &pvtx_ndof );
  tree->Branch("pvtx_size", &pvtx_size );

  tree->Branch("svtx_x", &svtx_x );
  tree->Branch("svtx_y", &svtx_y );
  tree->Branch("svtx_z", &svtx_z );
  tree->Branch("svtx_chi2", &svtx_chi2 );
  tree->Branch("svtx_ndof", &svtx_ndof );

  tree->Branch("ruclu_eta", &ruclu_eta );
  tree->Branch("ruclu_phi", &ruclu_phi );
  tree->Branch("ruclu_energy", &ruclu_energy );
  tree->Branch("ruclu_r1", &ruclu_r1 );
  tree->Branch("ruclu_r2", &ruclu_r2 );
  tree->Branch("ruclu_r3", &ruclu_r3 );
  tree->Branch("ruclu_pfIso", &ruclu_pfIso );
  tree->Branch("ruclu_dipho", &ruclu_dipho ); //dcom
  tree->Branch("ruclu_monopho", &ruclu_monopho );
  tree->Branch("ruclu_hadron", &ruclu_hadron );
  tree->Branch("ruclu_moe", &ruclu_moe );
  tree->Branch("ruclu_isPhoton", &ruclu_isPhoton );
  tree->Branch("ruclu_hasPixelSeed", &ruclu_hasPixelSeed );
  tree->Branch("ruclu_passElectronVeto", &ruclu_passElectronVeto );

}

flattener::~flattener()
{
}

// member functions

// ------------ method called to produce the data  ------------
void
flattener::analyze( const edm::Event & iEvent, const edm::EventSetup & iSetup)
{
  run = iEvent.eventAuxiliary().run();
  lumiSec = iEvent.eventAuxiliary().luminosityBlock();
  id = iEvent.eventAuxiliary().id().event();
  
  //Triggering
  Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerResultsToken, triggerResults); // Get this event's trigger info

  //Get triggers from file
  std::ifstream myfile(trigfile.c_str());
  if(!myfile){std::cout << "Can't open file" << std::endl;}
  std::string str;
  std::vector<std::string> datatriggernames;
  while(std::getline(myfile, str)){
    datatriggernames.push_back(str);
  }
  myfile.close();
  //

  int tcount = 0;

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerResults); // get all trigger names
  for( unsigned int dtn = 0; dtn < datatriggernames.size(); dtn += 1 ){ //Loop through trigger name array
    for (unsigned int i = 0, n = triggerResults->size(); i < n; ++i) // loop aver all names
    {
      if (names.triggerName(i).find( datatriggernames.at(dtn) )!=std::string::npos) // if passes a specific trigger:
      {
        if (triggerResults->accept(i)) {
          trigs.push_back(1);
          tcount += 1;
          }
        else{trigs.push_back(0);}
      }
    }
  }

  if (isMC){
    //GenParticles
    Handle<vector<reco::GenParticle>> genpart;
    iEvent.getByToken(genpartToken, genpart);
    for (auto gp_iter = genpart->begin(); gp_iter != genpart->end(); ++gp_iter){
      genpart_pt.push_back(gp_iter->pt()); 
      genpart_eta.push_back(gp_iter->eta()); 
      genpart_phi.push_back(gp_iter->phi()); 
      genpart_energy.push_back(gp_iter->energy()); 
      genpart_mass.push_back(gp_iter->mass()); 
      genpart_pdgid.push_back(gp_iter->pdgId());
      genpart_status.push_back(gp_iter->status());
      if (gp_iter->mother() != nullptr){
        genpart_motherpdgid.push_back(gp_iter->mother()->pdgId());
      }
      else{
        genpart_motherpdgid.push_back(0);
      }
    }
  }

  //Jets
  Handle<vector<pat::Jet>> patJet;
  iEvent.getByToken(patjetToken, patJet); 
  for (auto jet_iter = patJet->begin(); jet_iter != patJet->end(); ++jet_iter){
    jet_pt.push_back(jet_iter->pt());
    jet_eta.push_back(jet_iter->eta());
    jet_phi.push_back(jet_iter->phi());
    jet_energy.push_back(jet_iter->energy());
    jet_mass.push_back(jet_iter->mass());
    jet_btag.push_back(jet_iter->bDiscriminator(btagName.c_str()));
  }

  //METs
  Handle<vector<pat::MET>> met;
  iEvent.getByToken(metToken, met); 
  for (auto met_iter = met->begin(); met_iter != met->end(); ++met_iter){
    met_pt.push_back(met_iter->pt());
    met_eta.push_back(met_iter->eta());
    met_phi.push_back(met_iter->phi());
    met_energy.push_back(met_iter->energy());
    met_mass.push_back(met_iter->mass());
  }
  
  //Muons
  Handle<vector<pat::Muon>> muon;
  iEvent.getByToken(muonToken, muon);
  for (auto muon_iter = muon->begin(); muon_iter != muon->end(); ++muon_iter){
    muon_pt.push_back(muon_iter->pt());
    muon_eta.push_back(muon_iter->eta());
    muon_phi.push_back(muon_iter->phi());
    muon_energy.push_back(muon_iter->energy());
    muon_mass.push_back(muon_iter->mass());
  }

  //Electrons
  Handle<vector<pat::Electron>> electron;
  iEvent.getByToken(electronToken, electron); 
  for (auto electron_iter = electron->begin(); electron_iter != electron->end(); ++electron_iter){
    electron_pt.push_back(electron_iter->pt());
    electron_eta.push_back(electron_iter->eta());
    electron_phi.push_back(electron_iter->phi());
    electron_energy.push_back(electron_iter->energy());
    electron_mass.push_back(electron_iter->mass());
  }

  //PAT Photons
  Handle<vector<pat::Photon>> patpho;
  iEvent.getByToken(patPhoToken, patpho); 
  for (auto patpho_iter = patpho->begin(); patpho_iter != patpho->end(); ++patpho_iter){
    patpho_pt.push_back(patpho_iter->pt());
    patpho_eta.push_back(patpho_iter->eta());
    patpho_phi.push_back(patpho_iter->phi());
    patpho_energy.push_back(patpho_iter->energy());
    patpho_mass.push_back(patpho_iter->mass());
    patpho_r9.push_back(patpho_iter->r9());
    patpho_sieie.push_back(patpho_iter->sigmaIetaIeta());
    patpho_hoe.push_back(patpho_iter->hadTowOverEm());
    patpho_chargedHadronIso.push_back(patpho_iter->chargedHadronIso());
    patpho_neutralHadronIso.push_back(patpho_iter->neutralHadronIso());
    patpho_photonIso.push_back(patpho_iter->photonIso());
    patpho_hasPixelSeed.push_back(patpho_iter->hasPixelSeed());
    patpho_passElectronVeto.push_back(patpho_iter->passElectronVeto());

    if (patpho_iter->photonID("cutBasedPhotonID-Fall17-94X-V2-loose") == 1){
      patpho_cutBased.push_back(1);
    }
    else if (patpho_iter->photonID("cutBasedPhotonID-Fall17-94X-V2-medium") == 1){
      patpho_cutBased.push_back(2);
    }
    else if (patpho_iter->photonID("cutBasedPhotonID-Fall17-94X-V2-tight") == 1){
      patpho_cutBased.push_back(3);
    }
    else {
      patpho_cutBased.push_back(0);
    }
    patpho_mvaID.push_back(patpho_iter->userFloat("PhotonMVAEstimatorRunIIFall17v2Values"));
  }

  //Primary Vertex
  Handle<vector<reco::Vertex>> pvtx;
  iEvent.getByToken(pvtxToken, pvtx);
  pvtx_size.push_back(pvtx->size()); // @size?
  for (auto pvtx_iter = pvtx->begin(); pvtx_iter != pvtx->end(); ++pvtx_iter){
    pvtx_x.push_back(pvtx_iter->x());
    pvtx_y.push_back(pvtx_iter->y());
    pvtx_z.push_back(pvtx_iter->z());
    pvtx_chi2.push_back(pvtx_iter->chi2());
    pvtx_ndof.push_back(pvtx_iter->ndof());
  }

  //Secondary Vertex
  Handle<vector<reco::VertexCompositePtrCandidate>> svtx;
  iEvent.getByToken(svtxToken, svtx);
  for (auto svtx_iter = svtx->begin(); svtx_iter != svtx->end(); ++svtx_iter){
    svtx_x.push_back(svtx_iter->vx());
    svtx_y.push_back(svtx_iter->vy());
    svtx_z.push_back(svtx_iter->vz());
    //svtx_chi2.push_back(svtx_iter->chi2_()); //TODO: Fix this
    //svtx_ndof.push_back(svtx_iter->ndof_());
  }

  //RUCLU
  Handle<vector<float> > ru_etas;
  Handle<vector<float> > ru_phis;
  Handle<vector<float> > ru_energys;
  Handle<vector<float> > ru_r1s;
  Handle<vector<float> > ru_r2s;
  Handle<vector<float> > ru_r3s;
  Handle<vector<float> > monophos;
  Handle<vector<float> > diphos; // dcom
  Handle<vector<float> > hadrons;
  Handle<vector<float> > moe;

  iEvent.getByToken(ruclu_etaToken, ru_etas);
  iEvent.getByToken(ruclu_phiToken, ru_phis);
  iEvent.getByToken(ruclu_energyToken, ru_energys);
  iEvent.getByToken(ruclu_r1Token, ru_r1s);
  iEvent.getByToken(ruclu_r2Token, ru_r2s);
  iEvent.getByToken(ruclu_r3Token, ru_r3s);
  iEvent.getByToken(monophoToken, monophos);
  iEvent.getByToken(diphoToken, diphos); // dcom
  iEvent.getByToken(hadronToken, hadrons);
  iEvent.getByToken(moeToken, moe);

  Handle<vector<pat::PackedCandidate>> pfcand;
  iEvent.getByToken(pfcandToken, pfcand);

  if(diphos.isValid() == 1){
    for (unsigned int i=0; i<ru_etas.product()->size(); ++i) {
      ruclu_eta.push_back(ru_etas.product()->at(i));
      ruclu_phi.push_back(ru_phis.product()->at(i));
      ruclu_energy.push_back(ru_energys.product()->at(i));
      ruclu_r1.push_back(ru_r1s.product()->at(i));
      ruclu_r2.push_back(ru_r2s.product()->at(i));
      ruclu_r3.push_back(ru_r3s.product()->at(i));
      ruclu_dipho.push_back(diphos.product()->at(i)); // dcom
      ruclu_monopho.push_back(monophos.product()->at(i));
      ruclu_hadron.push_back(hadrons.product()->at(i));
      ruclu_moe.push_back(moe.product()->at(i));

      // Calculate pfIso
      float pfCandE = 0;
      for (auto pfcand_iter = pfcand->begin(); pfcand_iter != pfcand->end(); ++pfcand_iter){

        // Calculate dR
        float pi = 3.14159265358979323846;
        float dPhi = pfcand_iter->phi() - ru_phis.product()->at(i);
        if (dPhi > pi) dPhi -= 2*pi;
        if (dPhi <= -pi) dPhi += 2*pi;

        float dR = sqrt(pow(pfcand_iter->eta() - ru_etas.product()->at(i), 2) + pow(dPhi, 2));
        if (dR < 0.3){
          pfCandE += pfcand_iter->energy();
        }
      }
      if (pfCandE == 0) pfCandE = 1;
      ruclu_pfIso.push_back(ru_energys.product()->at(i) / pfCandE);

      // Calculate isPhoton
      bool isPhoton = false;
      for (auto patpho_iter = patpho->begin(); patpho_iter != patpho->end(); ++patpho_iter){
        float dR = sqrt(pow(patpho_iter->eta() - ru_etas.product()->at(i), 2) + pow(patpho_iter->phi() - ru_phis.product()->at(i), 2));
        if (dR < 0.15){
          isPhoton = true;
          ruclu_isPhoton.push_back(true);
          ruclu_hasPixelSeed.push_back(patpho_iter->hasPixelSeed());
          ruclu_passElectronVeto.push_back(patpho_iter->passElectronVeto());
          break;
        }
      }
      if (!isPhoton){
        ruclu_isPhoton.push_back(false);
        ruclu_hasPixelSeed.push_back(false);
        ruclu_passElectronVeto.push_back(false);
      }

    }
  }
  //else{
  //  ruclu_eta.push_back(-999.);
  //  ruclu_phi.push_back(-999.);
  //  ruclu_energy.push_back(-999.);
  //  ruclu_r1.push_back(-999.);
  //  ruclu_r2.push_back(-999.);
  //  ruclu_r3.push_back(-999.);
  //  ruclu_dipho.push_back(-999.);
  //  ruclu_dipho.push_back(-999.);
  //  ruclu_monopho.push_back(-999.);
  //  ruclu_hadron.push_back(-999.);
  //  ruclu_moe.push_back(-999.);
  //}

  tree->Fill();
  clearVars();

}

void 
flattener::clearVars(){

  trigs.clear();

  genpart_pt.clear();
  genpart_eta.clear();
  genpart_phi.clear();
  genpart_energy.clear();
  genpart_mass.clear();
  genpart_pdgid.clear();
  genpart_status.clear();
  genpart_motherpdgid.clear();

  jet_pt.clear();
  jet_eta.clear();
  jet_phi.clear();
  jet_energy.clear();
  jet_mass.clear();
  jet_btag.clear();

  met_pt.clear();
  met_eta.clear();
  met_phi.clear();
  met_energy.clear();
  met_mass.clear();

  muon_pt.clear();
  muon_eta.clear();
  muon_phi.clear();
  muon_energy.clear();
  muon_mass.clear();

  electron_pt.clear();
  electron_eta.clear();
  electron_phi.clear();
  electron_energy.clear();
  electron_mass.clear();

  patpho_pt.clear();
  patpho_eta.clear();
  patpho_phi.clear();
  patpho_energy.clear();
  patpho_mass.clear();
  patpho_r9.clear();
  patpho_sieie.clear();
  patpho_hoe.clear();
  patpho_chargedHadronIso.clear();
  patpho_neutralHadronIso.clear();
  patpho_photonIso.clear();
  patpho_cutBased.clear();
  patpho_mvaID.clear();
  patpho_hasPixelSeed.clear();
  patpho_passElectronVeto.clear();

  pvtx_x.clear();
  pvtx_y.clear();
  pvtx_z.clear();
  pvtx_chi2.clear();
  pvtx_ndof.clear();
  pvtx_size.clear();

  svtx_x.clear();
  svtx_y.clear();
  svtx_z.clear();
  svtx_chi2.clear();
  svtx_ndof.clear();

  ruclu_eta.clear();
  ruclu_phi.clear();
  ruclu_energy.clear();
  ruclu_r1.clear();
  ruclu_r2.clear();
  ruclu_r3.clear();
  ruclu_pfIso.clear();
  ruclu_dipho.clear();
  ruclu_monopho.clear();
  ruclu_hadron.clear();
  ruclu_moe.clear();
  ruclu_isPhoton.clear();
  ruclu_hasPixelSeed.clear();
  ruclu_passElectronVeto.clear();
}

void
flattener::beginJob()
{
}

void
flattener::endJob() {
}

void 
flattener::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
}

void 
flattener::endRun(edm::Run const&, edm::EventSetup const&)
{
}

void
flattener::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(flattener);
