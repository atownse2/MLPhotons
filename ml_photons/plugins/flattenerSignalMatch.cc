// -*- C++ -*-
//
// Package:    ML_Photons/flattenerSignalMatch
// Class:      flattenerSignalMatch
// 
/**\class flattenerSignalMatch flattenerSignalMatch.cc ML_Photons/flattenerSignalMatch/plugins/flattenerSignalMatch.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Steven Clark
//         Created:  Tue, 30 Mar 2021 17:29:55 GMT
//
//

#include "flattenerSignalMatch.h"
//

#include "FWCore/Utilities/interface/InputTag.h"

using namespace edm;
using namespace fastjet;
using namespace fastjet::contrib;

class flattenerSignalMatch : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns> {
   public:
      explicit flattenerSignalMatch( const edm::ParameterSet& );
      ~flattenerSignalMatch() override;

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

      int isDirectPrompt;

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

      //Gen a
      std::vector<float> a_pt;
      std::vector<float> a_eta;
      std::vector<float> a_phi;
      std::vector<float> a_energy;
      std::vector<float> a_mass;

      //Gen Photons
      std::vector<float> genpho_pt;
      std::vector<float> genpho_eta;
      std::vector<float> genpho_phi;
      std::vector<float> genpho_energy;
      std::vector<float> genpho_mass;

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
      std::vector<float> patpho_wp90;

      std::vector<float> pvtx_x;
      std::vector<float> pvtx_y;
      std::vector<float> pvtx_z;
      std::vector<float> pvtx_chi2;
      std::vector<float> pvtx_ndof;

      std::vector<float> svtx_x;
      std::vector<float> svtx_y;
      std::vector<float> svtx_z;
      std::vector<float> svtx_chi2;
      std::vector<float> svtx_ndof;

      std::vector<float> ruclu_etas;
      std::vector<float> ruclu_phis;
      std::vector<float> ruclu_energys;
      std::vector<float> ruclu_r1s;
      std::vector<float> ruclu_r2s;
      std::vector<float> ruclu_r3s;
      std::vector<float> dipho_scores;
      std::vector<float> monopho_scores;
      std::vector<float> hadron_scores;
      std::vector<float> moes;

  };

flattenerSignalMatch::flattenerSignalMatch(const edm::ParameterSet& iConfig):
  triggerResultsTag (iConfig.getParameter<edm::InputTag>("TriggerInputTag_HLT")),
  triggerResultsToken (consumes<edm::TriggerResults>(triggerResultsTag)),

  genpartTag (iConfig.getParameter<edm::InputTag>("genpartInputTag")),
  genpartToken (consumes<vector<reco::GenParticle>>(genpartTag)),

  patjetTag (iConfig.getParameter<edm::InputTag>("patjetInputTag")),
  patjetToken (consumes<vector<pat::Jet>>(patjetTag)),

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
  
  // Access the TFileService
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree", "tree");

  tree->Branch("lumiSec", &lumiSec, "lumiSec/i");
  tree->Branch("run", &run, "run/i");
  tree->Branch("id", &id, "id/i");

  tree->Branch("isDirectPrompt", &isDirectPrompt, "isDirectPrompt/i");

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

  //Generator-level as
  tree->Branch("a_pt", &a_pt );
  tree->Branch("a_eta", &a_eta );
  tree->Branch("a_phi", &a_phi );
  tree->Branch("a_energy", &a_energy );
  tree->Branch("a_mass", &a_mass );

  //Generator-level photons
  tree->Branch("genpho_pt", &genpho_pt );
  tree->Branch("genpho_eta", &genpho_eta );
  tree->Branch("genpho_phi", &genpho_phi );
  tree->Branch("genpho_energy", &genpho_energy );
  tree->Branch("genpho_mass", &genpho_mass );

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
  tree->Branch("patpho_wp90", &patpho_wp90);

  tree->Branch("pvtx_x", &pvtx_x );
  tree->Branch("pvtx_y", &pvtx_y );
  tree->Branch("pvtx_z", &pvtx_z );
  tree->Branch("pvtx_chi2", &pvtx_chi2 );
  tree->Branch("pvtx_ndof", &pvtx_ndof );

  tree->Branch("svtx_x", &svtx_x );
  tree->Branch("svtx_y", &svtx_y );
  tree->Branch("svtx_z", &svtx_z );
  tree->Branch("svtx_chi2", &svtx_chi2 );
  tree->Branch("svtx_ndof", &svtx_ndof );

  tree->Branch("ruclu_eta", &ruclu_etas );
  tree->Branch("ruclu_phi", &ruclu_phis );
  tree->Branch("ruclu_energy", &ruclu_energys );
  tree->Branch("ruclu_r1", &ruclu_r1s );
  tree->Branch("ruclu_r2", &ruclu_r2s );
  tree->Branch("ruclu_r3", &ruclu_r3s );
  tree->Branch("diphoScores", &dipho_scores ); //dcom
  tree->Branch("monophoScores", &monopho_scores );
  tree->Branch("hadronScores", &hadron_scores );
  tree->Branch("moe", &moes );

}

flattenerSignalMatch::~flattenerSignalMatch()
{
}

// member functions

// ------------ method called to produce the data  ------------
void
flattenerSignalMatch::analyze( const edm::Event & iEvent, const edm::EventSetup & iSetup)
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
  

  //Jets
  Handle<vector<pat::Jet>> patJet;
  iEvent.getByToken(patjetToken, patJet); // Get this event's trigger info
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
  iEvent.getByToken(metToken, met); // Get this event's trigger info
  for (auto met_iter = met->begin(); met_iter != met->end(); ++met_iter){
    met_pt.push_back(met_iter->pt());
    met_eta.push_back(met_iter->eta());
    met_phi.push_back(met_iter->phi());
    met_energy.push_back(met_iter->energy());
    met_mass.push_back(met_iter->mass());
  }
    
  //Loop through gen particles
  Handle<vector<reco::GenParticle>> genpart;
  iEvent.getByToken(genpartToken, genpart);
  for (auto gp_iter = genpart->begin(); gp_iter != genpart->end(); ++gp_iter){
    if (gp_iter->pdgId() == 90000054 && gp_iter->numberOfDaughters()==2){
      if( gp_iter->daughter(0)->pdgId() == 22 && gp_iter->daughter(1)->pdgId() == 22 ){
        a_pt.push_back(gp_iter->pt());
        a_eta.push_back(gp_iter->eta());
        a_phi.push_back(gp_iter->phi());
        a_energy.push_back(gp_iter->energy());
        a_mass.push_back(gp_iter->mass());

        genpho_pt.push_back(gp_iter->daughter(0)->pt());
        genpho_eta.push_back(gp_iter->daughter(0)->eta());
        genpho_phi.push_back(gp_iter->daughter(0)->phi());
        genpho_energy.push_back(gp_iter->daughter(0)->energy());
        genpho_mass.push_back(gp_iter->daughter(0)->mass());

        genpho_pt.push_back(gp_iter->daughter(1)->pt());
        genpho_eta.push_back(gp_iter->daughter(1)->eta());
        genpho_phi.push_back(gp_iter->daughter(1)->phi());
        genpho_energy.push_back(gp_iter->daughter(1)->energy());
        genpho_mass.push_back(gp_iter->daughter(1)->mass());

        }
      }
    } //genpart loop
  
  //Muons
  Handle<vector<pat::Muon>> muon;
  iEvent.getByToken(muonToken, muon); // Get this event's trigger info
  for (auto muon_iter = muon->begin(); muon_iter != muon->end(); ++muon_iter){
    muon_pt.push_back(muon_iter->pt());
    muon_eta.push_back(muon_iter->eta());
    muon_phi.push_back(muon_iter->phi());
    muon_energy.push_back(muon_iter->energy());
    muon_mass.push_back(muon_iter->mass());
  }

  //Electrons
  Handle<vector<pat::Electron>> electron;
  iEvent.getByToken(electronToken, electron); // Get this event's trigger info
  for (auto electron_iter = electron->begin(); electron_iter != electron->end(); ++electron_iter){
     electron_pt.push_back(electron_iter->pt());
     electron_eta.push_back(electron_iter->eta());
     electron_phi.push_back(electron_iter->phi());
     electron_energy.push_back(electron_iter->energy());
     electron_mass.push_back(electron_iter->mass());
   }
  
  //PAT Photons
  Handle<vector<pat::Photon>> patpho;
  iEvent.getByToken(patPhoToken, patpho); // Get this event's trigger info
  for (auto patpho_iter = patpho->begin(); patpho_iter != patpho->end(); ++patpho_iter){
     patpho_pt.push_back(patpho_iter->pt());
     patpho_eta.push_back(patpho_iter->eta());
     patpho_phi.push_back(patpho_iter->phi());
     patpho_energy.push_back(patpho_iter->energy());
     patpho_mass.push_back(patpho_iter->mass());
     //patpho_wp90.push_back(patpho_iter->userFloat("PhotonMVAEstimatorRunIIFall17v1Values"));
     patpho_wp90.push_back(patpho_iter->photonID("mvaPhoID-RunIIFall17-v1-wp90")); //https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2
   }

  //Primary Vertex
  Handle<vector<reco::Vertex>> pvtx;
  iEvent.getByToken(pvtxToken, pvtx); // Get this event's trigger info
  for (auto pvtx_iter = pvtx->begin(); pvtx_iter != pvtx->end(); ++pvtx_iter){
    pvtx_x.push_back(pvtx_iter->x());
    pvtx_y.push_back(pvtx_iter->y());
    pvtx_z.push_back(pvtx_iter->z());
    pvtx_chi2.push_back(pvtx_iter->chi2());
    pvtx_ndof.push_back(pvtx_iter->ndof());
  }

  //Secondary Vertex
  Handle<vector<reco::VertexCompositePtrCandidate>> svtx;
  iEvent.getByToken(svtxToken, svtx); // Get this event's trigger info
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

  if(diphos.isValid() == 1){
    for (unsigned int i=0; i<ru_etas.product()->size(); ++i) {
      ruclu_etas.push_back(ru_etas.product()->at(i));
      ruclu_phis.push_back(ru_phis.product()->at(i));
      ruclu_energys.push_back(ru_energys.product()->at(i));
      ruclu_r1s.push_back(ru_r1s.product()->at(i));
      ruclu_r2s.push_back(ru_r2s.product()->at(i));
      ruclu_r3s.push_back(ru_r3s.product()->at(i));
      dipho_scores.push_back(diphos.product()->at(i)); // dcom
      monopho_scores.push_back(monophos.product()->at(i));
      hadron_scores.push_back(hadrons.product()->at(i));
      moes.push_back(moe.product()->at(i));
    }
  }

  tree->Fill();
  clearVars();

}

void 
flattenerSignalMatch::clearVars(){

  trigs.clear();

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

  a_pt.clear();
  a_eta.clear();
  a_phi.clear();
  a_energy.clear();
  a_mass.clear();

  genpho_pt.clear();
  genpho_eta.clear();
  genpho_phi.clear();
  genpho_energy.clear();
  genpho_mass.clear();

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
  patpho_wp90.clear();

  pvtx_x.clear();
  pvtx_y.clear();
  pvtx_z.clear();
  pvtx_chi2.clear();
  pvtx_ndof.clear();

  svtx_x.clear();
  svtx_y.clear();
  svtx_z.clear();
  svtx_chi2.clear();
  svtx_ndof.clear();

  ruclu_etas.clear();
  ruclu_phis.clear();
  ruclu_energys.clear();
  ruclu_r1s.clear();
  ruclu_r2s.clear();
  ruclu_r3s.clear();
  dipho_scores.clear();
  monopho_scores.clear();
  hadron_scores.clear();
  moes.clear();
}

void
flattenerSignalMatch::beginJob()
{
}

void
flattenerSignalMatch::endJob() {
}

void 
flattenerSignalMatch::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
}

void 
flattenerSignalMatch::endRun(edm::Run const&, edm::EventSetup const&)
{
}

void
flattenerSignalMatch::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(flattenerSignalMatch);