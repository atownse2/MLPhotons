// -*- C++ -*-
//
// Package:    ML_Photons/ml_photons
// Class:      ml_photons
// 
/**\class ml_photons ml_photons.cc ML_Photons/ml_photons/plugins/ml_photons.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Steven Clark
//         Created:  Tue, 30 Mar 2021 17:29:55 GMT
//
//

#include "ml_photons.h"
//
using namespace cms::Ort;
using namespace edm;
using namespace fastjet;
using namespace fastjet::contrib;

//
// class declaration
//

//class ml_photons : public edm::stream::EDProducer<edm::GlobalCache<ONNXRuntime>> {
class ml_photons : public edm::stream::EDProducer<> {
   public:
      explicit ml_photons( const edm::ParameterSet& );
      //explicit ml_photons(const edm::ParameterSet&, const ONNXRuntime *);
      ~ml_photons() override;

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      //static std::unique_ptr<ONNXRuntime> initializeGlobalCache(const edm::ParameterSet &);
      //static void globalEndJob(const ONNXRuntime *);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      // ----------member data ---------------------------
      const double MATCH_DR;
      EDGetTokenT<std::vector<reco::CaloCluster>> token_clusters;
      EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> token_HEE;
      EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> token_HEB;
      EDGetTokenT<std::vector<reco::Vertex>> token_vtx;
      EDGetTokenT<edm::TriggerResults> triggerResultsToken;
      ONNXRuntime ort_class;
      //ONNXRuntime ort_regress;
  };

ml_photons::ml_photons(const edm::ParameterSet& iConfig):
//ml_photons::ml_photons(const edm::ParameterSet& iConfig, const ONNXRuntime *cache):
  MATCH_DR(iConfig.getParameter<double>("MATCH_DeltaR") ),
  token_clusters(consumes<std::vector<reco::CaloCluster>>(iConfig.getParameter<edm::InputTag>("CluInputTag"))),
  token_HEE(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >>(iConfig.getParameter<edm::InputTag>("HEEInputTag"))),
  token_HEB(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >>(iConfig.getParameter<edm::InputTag>("HEBInputTag"))),
  token_vtx(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("VtxInputTag"))),
  triggerResultsToken(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerInputTag_HLT"))),
  ort_class(iConfig.getParameter<edm::FileInPath>("classifier_path").fullPath())
  //ort_regress(iConfig.getParameter<edm::FileInPath>("regressor_path").fullPath())
{
}


ml_photons::~ml_photons()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
//std::unique_ptr<ONNXRuntime> ml_photons::initializeGlobalCache(const edm::ParameterSet &iConfig) {
    //return std::make_unique<ONNXRuntime>(iConfig.getParameter<edm::FileInPath>("classifier_path").fullPath());
//}

//
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

void
ml_photons::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  //desc.add<edm::FileInPath>("classifier_path",
  //                            edm::FileInPath("ML_Photons/ml_photons/plugins/classifier.onnx"));
 
  //desc.add<edm::FileInPath>("regressor_path",
  //                            edm::FileInPath("ML_Photons/ml_photons/plugins/regressor.onnx"));

  desc.setUnknown();
  //descriptions.addDefault(desc);
  descriptions.addWithDefaultLabel(desc);
}

//void ml_photons::globalEndJob() {}
//void ml_photons::globalEndJob(const ONNXRuntime *cache) {}

// ------------ method called to produce the data  ------------
void
ml_photons::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

/////////////////////////////////////////////////////////////////////////////
  int evt_id = iEvent.eventAuxiliary().id().event();
  int evt_run = iEvent.eventAuxiliary().id().run();
  int evt_lumi = iEvent.eventAuxiliary().id().luminosityBlock();


  Handle<std::vector<reco::CaloCluster>> CLS_pho;
  iEvent.getByToken(token_clusters, CLS_pho);
  Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> HEB;
  iEvent.getByToken(token_HEB, HEB);
  Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerResultsToken, triggerResults); // Get this event's trigger info
  Handle<std::vector<reco::Vertex>> vtx;
  iEvent.getByToken(token_vtx, vtx);


// TODO: Add triggering here

//
// TODO: Add vertexing here

//
  //////////////////////////////////////////////////////
  //BEGIN WITH CLUSTERING
  //////////////////////////////////////////////////////
  edm::ESHandle<CaloGeometry> pG; 
  iSetup.get<CaloGeometryRecord>().get(pG);
  const CaloGeometry cG = *pG;
  const CaloSubdetectorGeometry* EBgeom = cG.getSubdetectorGeometry(DetId::Ecal, 1);
	
  std::vector<cluster> Clusters;
  float maxE = 0.;
  int minEta = 999; int maxEta = -999; int minPhi = 999; int maxPhi = -999; //FIXME take out minmax stuff
  //bool fail_kTowerRecovered = false;
  for ( edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>::const_iterator cBegin = HEB->begin(), cEnd = HEB->end(), ic = cBegin; ic != cEnd; ++ic )
  { // Loop over the crystals in the EcalRecHitas collection
    if (ic->energy() > 0.) // let's do some zero suppression to save space
    {
      //if(ic->checkFlag(EcalRecHit::kTowerRecovered) == 1){
      //  fail_kTowerRecovered = true;
      //  continue;
      //}

      EBDetId ebDID(ic->detid()); //Need this to get the ieta/iphi
      auto cell = EBgeom->getGeometry(ebDID);
      float eta, phi, E;
      bool ncrack;
      eta = cell->getPosition().eta(); phi = cell->getPosition().phi(); E = ic->energy();
      ncrack = ebDID.isNextToBoundary(ebDID);
      cluster tC(eta, phi, ebDID.ieta(), ebDID.iphi(), E, ncrack);
      if (E > maxE) maxE = E;
      if (ebDID.ieta() < minEta) minEta = ebDID.ieta();
      if (ebDID.ieta() > maxEta) maxEta = ebDID.ieta();
      if (ebDID.iphi() < minPhi) minPhi = ebDID.iphi();
      if (ebDID.iphi() > maxPhi) maxPhi = ebDID.iphi();
      Clusters.push_back(tC);
      }
  }
  
  unsigned int sizebefore = 1;
  unsigned int sizeafter = 0;

  while (sizebefore != sizeafter)
  {
    sort(Clusters.begin(), Clusters.end(), &comparator);
    sizebefore = Clusters.size();

    //std::cout<<"----------------->  "<<sizebefore<<std::endl;
    Clusters = DoPairings(Clusters, MATCH_DR);
    sizeafter = Clusters.size();
  } 
  
  // Now we've got all the clusters

/////////////////////////////////////////////////////////////////////////////

  static const int isize = 30;
  std::vector<std::string> input_names;
  for(int i=0; i<isize*isize; i++){
    input_names.push_back("input");
  }
  std::vector<std::string> output_names ; //Can be blank

  std::vector<std::vector<float>> input(isize*isize, std::vector<float>(isize*isize,0.0)); //For some reason, this has to be of size 900x900 :/
  std::vector<float> outputs(3, 0); // init to zeros

  std::ofstream output_class ;
  output_class.open("classifier_scores.csv", std::ios_base::app);

  for (auto C : Clusters)
  {
    C.makeImage();
    std::vector<std::vector<float>> img = C.image;

    //Perform the Classification
    outputs = ort_class.run(input_names, img, output_names, 1)[0];

    //Compute softmax of Classifier output
    float denom = 0.0;
    for(unsigned int ii=0; ii< outputs.size(); ii++){
      denom += exp(outputs.at(ii));
    }

    output_class << evt_id << ", ";

    for(unsigned int ii=0; ii< outputs.size(); ii++){
      output_class << exp(outputs.at(ii)) / denom << ", ";
    }
    output_class << std::endl;
  }

}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
ml_photons::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
ml_photons::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
ml_photons::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
ml_photons::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
ml_photons::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
ml_photons::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

//define this as a plug-in
DEFINE_FWK_MODULE(ml_photons);
