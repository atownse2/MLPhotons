// -*- C++ -*-
//
// Package:    MLPhotons/
// Class:      MLPhotonProducer
// 
/**\class MLPhotonProducer MLPhotonProducer.cc MLPhotonProducer/MLPhotonProducer/plugins/MLPhotonProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Steven Clark
//         Created:  Tue, 30 Mar 2021 17:29:55 GMT
//
//

#include <iostream>
#include <vector>
#include <memory>

// Framework
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// #include "Geometry/Records/interface/CaloGeometryRecord.h"
// #include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
// #include "Geometry/CaloTopology/interface/CaloTopology.h"
// #include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

// #include "DataFormats/VertexReco/interface/Vertex.h"
// #include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "RecoEgamma/EgammaMLPhotonProducers/interface/MLPhotonProducer.h"

MLPhotonProducer::MLPhotonProducer(const edm::ParameterSet& iConfig):
  token_clusters(consumes<std::vector<reco::CaloCluster>>(iConfig.getParameter<edm::InputTag>("CluInputTag"))),
  token_HEE(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >>(iConfig.getParameter<edm::InputTag>("HEEInputTag"))),
  token_HEB(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >>(iConfig.getParameter<edm::InputTag>("HEBInputTag"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("VtxInputTag"))),
  ort_class(iConfig.getParameter<std::string>("classifier_path")),
  ort_regress(iConfig.getParameter<std::string>("regressor_path")),
  MLPhotonCollection_(iConfig.getParameter<std::string>("cluster_name"))
{
  produces<reco::MLPhotonCollection>(MLPhotonCollection_);
}

MLPhotonProducer::~MLPhotonProducer(){}

void MLPhotonProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<std::vector<reco::CaloCluster>> CLS_pho;
  iEvent.getByToken(token_clusters, CLS_pho);
  edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> HEB;
  iEvent.getByToken(token_HEB, HEB);

////////////////////////////////////////////////////////
  //BEGIN WITH CLUSTERING
////////////////////////////////////////////////////////
  edm::ESHandle<CaloGeometry> pG; 
  iSetup.get<CaloGeometryRecord>().get(pG);
  const CaloGeometry cG = *pG;
  const CaloSubdetectorGeometry* EBgeom = cG.getSubdetectorGeometry(DetId::Ecal, 1);
	
  std::vector<Cluster> clusters;

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
      Cluster tC(eta, phi, ebDID.ieta(), ebDID.iphi(), E, ncrack);
      if (E > maxE) maxE = E;
      if (ebDID.ieta() < minEta) minEta = ebDID.ieta();
      if (ebDID.ieta() > maxEta) maxEta = ebDID.ieta();
      if (ebDID.iphi() < minPhi) minPhi = ebDID.iphi();
      if (ebDID.iphi() > maxPhi) maxPhi = ebDID.iphi();
      clusters.push_back(tC);
      }
  }
  
  unsigned int sizebefore = 1;
  unsigned int sizeafter = 0;

  while (sizebefore != sizeafter)
  {
    sort(clusters.begin(), clusters.end(), [](const auto& a, const auto& b) {return a.vec.Mag2() > b.vec.Mag2();});
    sizebefore = clusters.size();
    clusters = DoPairings(clusters, MATCH_DR);
    sizeafter = clusters.size();
  } 
  
  /////////////////////////////////////////////////////////////////////////////
  // Now we've got all the clusters

  //Outputs vector of MLPhotons
  // std::unique_ptr<std::vector<MLPhoton>> mlphotons(new std::vector<MLPhoton>);
  auto mlphotons = std::make_unique<reco::MLPhotonCollection>();

  // Get the primary vertex
  edm::Handle<reco::VertexCollection> vtxs;
  iEvent.getByToken(vtxToken_, vtxs);
  const reco::Vertex &pvtx = vtxs->front();

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

  for (auto C : clusters)
  {
    //make Image from clusters
    C.makeImage();

    //Create inputs for classifier and regressor
    std::vector<std::vector<float>> img = C.image;
    std::vector<float> r_img = img.at(0);
    std::vector<float> eta_v = {static_cast<float>(C.Eta())};

    //Normalize regressor image
    float img_sum = std::accumulate(r_img.begin(), r_img.end(), 0.0);
    for (unsigned int i=0; i<r_img.size(); i++){
      r_img.at(i) = r_img.at(i) / img_sum;
    }

    //Set up multiple inputs for regressor
    cms::Ort::FloatArrays regress_data_;
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

    float moe = regress_outputs.at(0); // mass/energy
    float energy = C.getTotalE();
    float eta = C.Eta();
    float phi = C.Phi();

    // Get LorentzVector
    math::PtEtaPhiMLorentzVector p4 = calculateLorentzVector(moe, energy, eta, phi, pvtx.z());

    reco::MLPhoton mlpho(p4, pvtx.position());

    mlpho.set_moe(moe);

    mlpho.set_diphoton_score(exp( class_outputs.at(0) ) / denom);
    mlpho.set_monophoton_score(exp( class_outputs.at(1) ) / denom);
    mlpho.set_hadron_score(exp( class_outputs.at(2) ) / denom);

    mlpho.set_r1( C.compute_En( 1. ) /  C.compute_En( 0. ));
    mlpho.set_r2( C.compute_En( 2. ) /  C.compute_En( 0. ));
    mlpho.set_r3( C.compute_En( 3. ) /  C.compute_En( 0. ));

    // mlphotons->emplace_back(mlpho);
    mlphotons->push_back(mlpho);
  }

  iEvent.put(std::move(mlphotons), MLPhotonCollection_);
}

bool MLPhotonProducer::comparator(const Cluster& C1, const Cluster& C2){
	return C1.vec.Mag2() > C2.vec.Mag2();
}

std::vector<Cluster> MLPhotonProducer::DoPairings(std::vector<Cluster> inC, float R){
	std::vector<Cluster> outC;
	while (inC.size() > 0) {
		//std::cout<<"in size > 0 loop"<<std::endl;
		Cluster C1(inC.at(0));

		if (inC.size() == 1) {
			//std::cout<<"only one element left"<<std::endl;
			outC.push_back(C1);	
			inC.erase(inC.begin()+0);
		} else {
			//std::cout<<"looking for pair:"<<std::endl;
			std::pair<int, float> BestPair = FindNearest(C1, inC);
			//std::cout<<"Nearest index and dR: "<<BestPair.first<<", "<<BestPair.second<<std::endl;
			if (BestPair.second > R){
				Cluster C2(inC.at(BestPair.first));
				outC.push_back(C1);
				outC.push_back(C2);
				inC.erase(inC.begin()+BestPair.first);
				inC.erase(inC.begin()+0);
			} else {
				Cluster C2(inC.at(BestPair.first));
				C1.combine(C2);
				outC.push_back(C1);
				inC.erase(inC.begin()+BestPair.first);
				inC.erase(inC.begin()+0);
			}
		}
	}
	return outC;
}

std::pair<int, float> MLPhotonProducer::FindNearest(Cluster C1, std::vector<Cluster> inC){
	float minR = 999.;
	int	index = 1;
	for (unsigned int i = 1; i < inC.size(); ++i)
	{ 
		float DR = C1.deltaR(inC.at(i));
		if (DR < minR)
		{
			minR = DR;
			index = i;
		}
	}
	return std::make_pair(index, minR);
}

math::PtEtaPhiMLorentzVector MLPhotonProducer::calculateLorentzVector(float moe, float energy, float eta, float phi, float zpv){
	float R = 129.0;
	float theta = 2.0 * std::atan(std::exp(-1.0 * std::abs(eta)));
	int sign = eta >= 0 ? 1 : -1;
	float z = R / std::tan(theta) * sign;
	float zp = std::abs(zpv) + std::abs(z);
	float thetaprime = std::atan(R / zp);
	float etaprime = std::log(std::tan(thetaprime / 2.0));

	if ((zpv < 0 && z >= 0) || (zpv >= 0 && z >= 0 && z > zpv) || (zpv < 0 && z < 0 && z >= zpv)) {
		etaprime *= -1.0;
	}

	float pt = energy * std::sin(std::atan(std::exp(-1 * etaprime)) * 2);
	math::PtEtaPhiMLorentzVector p4(pt, eta, phi, moe * energy);
	return p4;
}


//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MLPhotonProducer);
