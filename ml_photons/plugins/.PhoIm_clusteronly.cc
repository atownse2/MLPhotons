/**\class PhoIm_clusteronly PhoIm_clusteronly.cc phoimproducer/PhoIm/plugins/PhoIm_clusteronly.cc

## Description: [one line class summary]
	Produces inputs for the image photon taggers.
## Implementation:

*/
//
// Original Author:  Marc Osherson
//			Created:  Wed, 30 Jan 2019 17:20:33 GMT
//			Latest Update: 23-7-2019
//
//			Updated For Cluster: 5/27/2020, Steven Clark


#include "PhoIm.h"
//#include "Match_Cluster.h"
using namespace edm;
using namespace fastjet;
using namespace fastjet::contrib;

class PhoIm_clusteronly : public edm::stream::EDProducer<> {
   public:
      explicit PhoIm_clusteronly(const edm::ParameterSet&);
      ~PhoIm_clusteronly();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


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
};

PhoIm_clusteronly::PhoIm_clusteronly(const edm::ParameterSet& iConfig):
  MATCH_DR(iConfig.getParameter<double>("MATCH_DeltaR") ),
  token_clusters(consumes<std::vector<reco::CaloCluster>>(iConfig.getParameter<edm::InputTag>("CluInputTag"))),
  token_HEE(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >>(iConfig.getParameter<edm::InputTag>("HEEInputTag"))),
  token_HEB(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >>(iConfig.getParameter<edm::InputTag>("HEBInputTag"))),
  token_vtx(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("VtxInputTag"))),
  triggerResultsToken(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerInputTag_HLT")))

{    
}

PhoIm_clusteronly::~PhoIm_clusteronly(){}

void PhoIm_clusteronly::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  int evt_id = iEvent.eventAuxiliary().id().event();
  int evt_run = iEvent.eventAuxiliary().id().run();
  int evt_lumi = iEvent.eventAuxiliary().id().luminosityBlock();

  //Only processing 1/10 of data
  if (evt_id % 10 == 0){

  Handle<std::vector<reco::CaloCluster>> CLS_pho;
  iEvent.getByToken(token_clusters, CLS_pho);
	Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> HEB;
	iEvent.getByToken(token_HEB, HEB);
  Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerResultsToken, triggerResults); // Get this event's trigger info
  Handle<std::vector<reco::Vertex>> vtx;
  iEvent.getByToken(token_vtx, vtx);

  // Triggering:  =====+
  int trig_low = 0;
  int trig_high = 0;
  int trig_single = 0;

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerResults); // get all trigger names
  for (unsigned int i = 0, n = triggerResults->size(); i < n; ++i) // loop aver all names
  {
    if (names.triggerName(i).find("HLT_DoublePhoton60")!=std::string::npos) // if passes a specific trigger:
    {
      if (triggerResults->accept(i)) {trig_low = 1;}
    }
    if (names.triggerName(i).find("HLT_DoublePhoton85")!=std::string::npos) // if passes a specific trigger:
    {
      if (triggerResults->accept(i)) {trig_high = 1;}
    }
    if (names.triggerName(i).find("HLT_Photon175")!=std::string::npos) // if passes a specific trigger:
    {
      if (triggerResults->accept(i)) {trig_single = 1;}
    }
  }

  int TrigBit = ( ( (trig_low << 1) + trig_high ) << 1 ) + trig_single;
  //int TrigBit = (trig_low << 1 ) + trig_single;

  //if(TrigBit > 0) //Letting in everything
  //{

  // Vertexing
  float vtx_x = 0.;
  float vtx_y = 0.;
  float vtx_z = 0.;
  for ( std::vector<reco::Vertex>::const_iterator vBegin = vtx->begin(), vEnd = vtx->end(), iv = vBegin; iv != vEnd; ++iv ){
    if ( not iv->isFake() and iv->ndof() > 4 and abs(iv->position().z()) <= 24 and iv->position().Rho() <= 2 ){
      vtx_x = iv->position().x();
      vtx_y = iv->position().y();
      vtx_z = iv->position().z();
      break;
    }
  }
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
  bool fail_kTowerRecovered = false;
  for ( edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>::const_iterator cBegin = HEB->begin(), cEnd = HEB->end(), ic = cBegin; ic != cEnd; ++ic )
  { // Loop over the crystals in the EcalRecHitas collection
    if (ic->energy() > 0.) // let's do some zero suppression to save space
    {
      if(ic->checkFlag(EcalRecHit::kTowerRecovered) == 1){
        fail_kTowerRecovered = true;
        continue;
      }

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

  ///////////////////////////////////////////////////////////////////
  //LOOP THROUGH CLUSTERS, MATCH GEN PARTICLES
  ///////////////////////////////////////////////////////////////////
  std::ofstream output_c;
  output_c.open("ntuple.csv", std::ios_base::app);
  output_c << TrigBit << "," << evt_id << "," << evt_lumi << "," << evt_run << ",";
  for (auto C : Clusters)
  {
      //////////////////////////////////////////
      //OUTPUT EVERYTHING
      //////////////////////////////////////////

      //std::ofstream output_all;
      //output_all.open("clusters.csv", std::ios_base::app);
      output_c << vtx_x<<","<< vtx_y<<"," << vtx_z<<","<< fail_kTowerRecovered << ",";
      C.PRINT_Eta_Phi(output_c);
    }
  output_c<<std::endl;
  //} //Trigger
  }
}


void PhoIm_clusteronly::beginStream(edm::StreamID)
{}
void PhoIm_clusteronly::endStream()
{}

void PhoIm_clusteronly::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
DEFINE_FWK_MODULE(PhoIm_clusteronly);
