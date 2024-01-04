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

#include "FWCore/Utilities/interface/InputTag.h"

using namespace cms::Ort;
using namespace edm;
using namespace fastjet;
using namespace fastjet::contrib;

class ml_photons : public edm::stream::EDProducer<> {
   public:
      explicit ml_photons( const edm::ParameterSet& );
      ~ml_photons() override;

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      // ----------member data ---------------------------
      const double MATCH_DR = 0.15; //Cluster Size
      EDGetTokenT<std::vector<reco::CaloCluster>> token_clusters;
      EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> token_HEE;
      EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> token_HEB;
      ONNXRuntime ort_class;
      ONNXRuntime ort_regress;
      std::string cname_;

  };

ml_photons::ml_photons(const edm::ParameterSet& iConfig):
  token_clusters(consumes<std::vector<reco::CaloCluster>>(iConfig.getParameter<edm::InputTag>("CluInputTag"))),
  token_HEE(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >>(iConfig.getParameter<edm::InputTag>("HEEInputTag"))),
  token_HEB(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >>(iConfig.getParameter<edm::InputTag>("HEBInputTag"))),
  ort_class(iConfig.getParameter<std::string>("classifier_path")),
  ort_regress(iConfig.getParameter<std::string>("regressor_path")),
  cname_(iConfig.getParameter<std::string>("cluster_name"))
{

  produces<std::vector<float>> (cname_+"Eta");
  produces<std::vector<float>> (cname_+"Phi");
  produces<std::vector<float>> (cname_+"E");
  produces<std::vector<float>> (cname_+"MoE");
  produces<std::vector<float>> (cname_+"Monopho");
  produces<std::vector<float>> (cname_+"Dipho");
  produces<std::vector<float>> (cname_+"Hadron");
  produces<std::vector<float>> (cname_+"R1");
  produces<std::vector<float>> (cname_+"R2");
  produces<std::vector<float>> (cname_+"R3");

}


ml_photons::~ml_photons()
{
}

// member functions
void
ml_photons::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addWithDefaultLabel(desc);
}


// ------------ method called to produce the data  ------------
void
ml_photons::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  Handle<std::vector<reco::CaloCluster>> CLS_pho;
  iEvent.getByToken(token_clusters, CLS_pho);
  Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> HEB;
  iEvent.getByToken(token_HEB, HEB);

////////////////////////////////////////////////////////
  //BEGIN WITH CLUSTERING
////////////////////////////////////////////////////////
  edm::ESHandle<CaloGeometry> pG; 
  iSetup.get<CaloGeometryRecord>().get(pG);
  const CaloGeometry cG = *pG;
  const CaloSubdetectorGeometry* EBgeom = cG.getSubdetectorGeometry(DetId::Ecal, 1);
	
  std::vector<cluster> Clusters;
  float maxE = 0.;
  int minEta = 999; int maxEta = -999; int minPhi = 999; int maxPhi = -999; //FIXME take out minmax stuff
  for ( edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>::const_iterator cBegin = HEB->begin(), cEnd = HEB->end(), ic = cBegin; ic != cEnd; ++ic )
  { // Loop over the crystals in the EcalRecHits collection
    if (ic->energy() > 0.) // let's do some zero suppression to save space
    {

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

    Clusters = DoPairings(Clusters, MATCH_DR);
    sizeafter = Clusters.size();
  } 
  
/////////////////////////////////////////////////////////////////////////////
// Now we've got all the clusters

  //Outputs
  std::unique_ptr<std::vector<float>> cluster_Eta_( new std::vector<float> );
  std::unique_ptr<std::vector<float>> cluster_Phi_( new std::vector<float> );
  std::unique_ptr<std::vector<float>> cluster_E_( new std::vector<float> );
  std::unique_ptr<std::vector<float>> cluster_MoE_( new std::vector<float> );
  std::unique_ptr<std::vector<float>> cluster_mono_( new std::vector<float> );
  std::unique_ptr<std::vector<float>> cluster_di_( new std::vector<float> );
  std::unique_ptr<std::vector<float>> cluster_had_( new std::vector<float> );
  std::unique_ptr<std::vector<float>> cluster_r1_( new std::vector<float> );
  std::unique_ptr<std::vector<float>> cluster_r2_( new std::vector<float> );
  std::unique_ptr<std::vector<float>> cluster_r3_( new std::vector<float> );

  //Setting up data types for onnx runtime
  static const int isize = 30;
  std::vector<std::string> class_input_names;
  for(int i=0; i<isize*isize; i++){
    class_input_names.push_back("img");
  }

  std::vector<std::string> regress_input_names;
  regress_input_names.emplace_back("img");
  regress_input_names.emplace_back("eta");

  std::vector<std::string> class_output_names ; //Initialize vectors to hold ML output. Can be blank
  std::vector<std::string> regress_output_names ; 

  std::vector<std::vector<float>> input(isize*isize, std::vector<float>(isize*isize,0.0)); //Possible fix: Classifier takes 900x900 vector. not sure why... :/
  std::vector<float> class_outputs(3, 0); // init to zeros
  std::vector<float> regress_outputs(1, 0); // init to zeros

  for (auto C : Clusters)
  {
    //make Image from clusters
    C.makeImage();

    //Create inputs for classifier and regressor
    std::vector<std::vector<float>> img = C.image;
    std::vector<float> r_img = img.at(0);
    float eta = C.vec.Eta();
    std::vector<float> eta_v = {eta};

    //Normalize regressor image
    float img_sum = std::accumulate(r_img.begin(), r_img.end(), 0.0);
    for (unsigned int i=0; i<r_img.size(); i++){
      r_img.at(i) = r_img.at(i) / img_sum;
    }

    //Set up multiple inputs for regressor
    FloatArrays regress_data_;
    regress_data_.emplace_back(r_img);
    regress_data_.emplace_back(eta_v);

    //Perform the Classification
    class_outputs = ort_class.run(class_input_names, img, {}, class_output_names, 1)[0];

    //Perform the Regression
    regress_outputs = ort_regress.run(regress_input_names, regress_data_, {}, regress_output_names, 1)[0];

    //Compute softmax of Classifier output
    float denom = 0.0;
    for(unsigned int ii=0; ii< class_outputs.size(); ii++){
      denom += exp(class_outputs.at(ii));
    }

    cluster_Eta_->emplace_back( C.vec.Eta() );
    cluster_Phi_->emplace_back( C.vec.Phi() );
    cluster_E_->emplace_back( C.getTotalE() );

    cluster_r1_->emplace_back( C.compute_En( 1. ) /  C.compute_En( 0. ));
    cluster_r2_->emplace_back( C.compute_En( 2. ) /  C.compute_En( 0. ));
    cluster_r3_->emplace_back( C.compute_En( 3. ) /  C.compute_En( 0. ));

    cluster_MoE_->emplace_back( regress_outputs.at(0) );

    cluster_mono_->emplace_back(exp( class_outputs.at(0) ) / denom);
    cluster_di_->emplace_back(exp( class_outputs.at(1) ) / denom);
    cluster_had_->emplace_back(exp( class_outputs.at(2) ) / denom);
  }

  iEvent.put(std::move(cluster_Eta_), cname_+"Eta");
  iEvent.put(std::move(cluster_Phi_), cname_+"Phi");
  iEvent.put(std::move(cluster_E_), cname_+"E");
  iEvent.put(std::move(cluster_MoE_), cname_+"MoE");
  iEvent.put(std::move(cluster_mono_), cname_+"Monopho");
  iEvent.put(std::move(cluster_di_), cname_+"Dipho");
  iEvent.put(std::move(cluster_had_), cname_+"Hadron");
  iEvent.put(std::move(cluster_r1_), cname_+"R1");
  iEvent.put(std::move(cluster_r2_), cname_+"R2");
  iEvent.put(std::move(cluster_r3_), cname_+"R3");

}

void
ml_photons::beginStream(edm::StreamID)
{
}

void
ml_photons::endStream() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ml_photons);
