// system include files
#include <memory>
#include <iostream>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <array>
#include <string>

// CMSSW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TH2.h"
#include "TCanvas.h"

#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/Nsubjettiness.hh"

#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGenericDetId.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenStatusFlags.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "FWCore/Framework/interface/makeRefToBaseProdFrom.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"


using namespace edm;

class cluster
	{
		public:

    static const int isize=30;

		TVector3 vec;
		std::vector<int> ietas;
		std::vector<int> iphis;
		std::vector<int> ncracks;
		std::vector<float> Es;
    std::vector<int> xcoords;
    std::vector<int> ycoords;
    std::vector<std::vector<float>> image;
    std::string simage;
    //std::vector<std::vector<float>> image(isize*isize, std::vector<float>(isize*isize,0.0));
		cluster(const float &Eta, const float &Phi, const int &iEta, const int &iPhi, const float &E, const bool &nCrack);
		cluster(const cluster & clu);

		~cluster();

		void combine(cluster C);
		void PRINT();
	  void PRINT(std::ofstream& outfile);
	  void PRINT_Eta_Phi(std::ofstream& outfile);
	  float getTotalE(void);
    void makeImage();
    float compute_En(float);
    std::string MakeString(int, int, float);

	};

	cluster::cluster(const float &Eta, const float &Phi, const int &iEta, const int &iPhi, const float &E, const bool &nCrack)
	{

		vec.SetPtEtaPhi(1.0, Eta, Phi);
		vec.SetMag(E);
		Es.push_back(E);	
		ietas.push_back(iEta);
		iphis.push_back(iPhi);
    ncracks.push_back(nCrack);
	}
	cluster::cluster(const cluster & C)
	{
		vec.SetPtEtaPhi(1.0, C.vec.Eta(), C.vec.Phi());
		vec.SetMag(C.vec.Mag());
		Es.insert( Es.end(), C.Es.begin(), C.Es.end() );
		ietas.insert( ietas.end(), C.ietas.begin(), C.ietas.end() );
		iphis.insert( iphis.end(), C.iphis.begin(), C.iphis.end() );
		ncracks.insert( ncracks.end(), C.ncracks.begin(), C.ncracks.end() );
	}
	cluster::~cluster(void)
	{

	}
	void cluster::combine(cluster C)
	{
		vec+=C.vec;
		ietas.insert( ietas.end(), C.ietas.begin(), C.ietas.end() );
		iphis.insert( iphis.end(), C.iphis.begin(), C.iphis.end() );
		ncracks.insert( ncracks.end(), C.ncracks.begin(), C.ncracks.end() );
		Es.insert( Es.end(), C.Es.begin(), C.Es.end() );
	}
	void cluster::PRINT(void)
	{
		//std::cout<<vec.Eta()<<","<<vec.Phi()<<":";
		for (unsigned int i = 0; i < Es.size(); ++i)
		{
			std::cout<<Es[i]<<","<<ietas[i]<<","<<iphis[i]<< "," << ncracks[i] << ":";
		}
		std::cout<<std::endl;
	}
	void cluster::PRINT(std::ofstream& outfile)
	{
		//outfile <<vec.Eta()<<","<<vec.Phi()<<",";
		for (unsigned int i = 0; i < Es.size(); ++i)
		{
			outfile<<" "<<ietas[i]<<" "<<iphis[i]<<" "<< ncracks[i] <<" "<<Es[i]<<' ';
		}
		outfile<<",";
		//outfile<<std::endl;
	}

	void cluster::PRINT_Eta_Phi(std::ofstream& outfile)
	{
		outfile << vec.Eta()<<","<<vec.Phi()<<",";
		for (unsigned int i = 0; i < Es.size(); ++i)
		{
			//outfile<<" "<<ietas[i]<<" "<<iphis[i]<<" "<<Es[i]<<' ';
			outfile<<" "<<ietas[i]<<" "<<iphis[i]<<" "<< ncracks[i] <<" "<<Es[i]<<' ';
		}
		outfile<<",";
		//outfile<<std::endl;
	}

	float cluster::getTotalE(void)
	{
    float totalE = 0;
    totalE = std::accumulate(Es.begin(), Es.end(), 0.0);
    return totalE;
  }
	
  std::string cluster::MakeString(int xx, int yy, float ee)
	{
      std::string sx(std::to_string(xx));
      std::string sy(std::to_string(yy));
      std::string se(std::to_string(ee));
      std::string sp(" ");
      std::string sc(", ");

      std::string myString ("");
      myString.append(sx);
      myString.append(sp);
      myString.append(sy);
      myString.append(sp);
      myString.append(se);
      myString.append(sp);
      myString.append(sc);

    return myString;
  }

	void cluster::makeImage()
	{
    int iphi_max = *std::max_element(iphis.begin(), iphis.end());
    int iphi_min = *std::min_element(iphis.begin(), iphis.end());

    if(iphi_max - iphi_min > isize && iphi_max > 340 && iphi_min < 20){  //This loop fixes wrapped events
      for(unsigned int ii=0; ii<iphis.size(); ii++){
        if(iphis[ii] < isize){
          iphis[ii] = 360 + iphis[ii] + 1;
        }
      }
    }

    //Add minimum + half isize to each element of ieta, iphi
    int ieta_min = *std::min_element(ietas.begin(), ietas.end());
    iphi_min = *std::min_element(iphis.begin(), iphis.end());

    for(int& d : ietas){d += ieta_min + (isize/2) ;}
    for(int& d : iphis){d += iphi_min + (isize/2) ;}

    int cxval = ietas[0];
    int cyval = iphis[0];

    //std::vector<int> xcoords;
    //std::vector<int> ycoords;

    for(unsigned int ii=0; ii<ietas.size(); ii++){
      xcoords.push_back(ietas[ii] - cxval + (isize / 2));
      ycoords.push_back(iphis[ii] - cyval + (isize / 2));
    }

    int y_min = *std::min_element(ycoords.begin(), ycoords.end());
    while(y_min < 0){
      for(int& d : ycoords){d += abs(y_min);}
      y_min = *std::min_element(ycoords.begin(), ycoords.end());
    }

    int y_max = *std::max_element(ycoords.begin(), ycoords.end());
    if(y_max > isize){
      for(int& d : ycoords){d -= y_max - isize;}
    }

    ///////////////////////////
    std::vector<std::vector<float>> img(isize, std::vector<float>(isize)); //initialize vector of size isize x isize to all zeros
    for(unsigned int ii=0; ii<xcoords.size(); ii++){
      int x = xcoords[ii];
      int y = ycoords[ii];
      float e = Es[ii];
      if(x >= isize){x=isize-1;}
      if(y >= isize){y=isize-1;}
      if(x < 0){x=0;}
      if(y < 0){y=0;}
      img[x][y] += e;
      
      std::string sbit = MakeString(x, y, e);
      simage.append(sbit);
    }

    ////////////////////////
    //Flatten Vector
    ///////////////////////

    std::vector<float>flatImg;//new 1D vector
    // using nested for loops to construct a new 1D vector
    for(unsigned int i=0; i<img.size(); i++){
     for(unsigned int j=0; j<img[i].size(); j++){
       flatImg.push_back(img[i][j]);
     }
    }

    for(unsigned int i=0; i<isize*isize; i++){
      image.push_back(flatImg);
    }

    return;
  }


  float cluster::compute_En(float n)
  {
    float totalE = 0.;
    float numerator = 0.;
    totalE = std::accumulate(Es.begin(), Es.end(), 0.0);
    for(unsigned int ii=0; ii < xcoords.size(); ii++){
      int xx = xcoords.at(ii);
      int yy = ycoords.at(ii);
      numerator += pow( (xx * xx + yy * yy), n/2 ) ;
    }

    return numerator / totalE;
  }


std::pair<int, float> FindNearest(cluster C1, std::vector<cluster> inC)
{
	float minR = 999.;
	int	index = 1;
	for (unsigned int i = 1; i < inC.size(); ++i)
	{ 
		float DR = C1.vec.DeltaR(inC.at(i).vec);
		if (DR < minR)
		{
			minR = DR;
			index = i;
		}
	}
	return std::make_pair(index, minR);
}

std::vector<cluster> DoPairings(std::vector<cluster> inC, float R)
{
	std::vector<cluster> outC;
	while (inC.size() > 0)
	{
		//std::cout<<"in size > 0 loop"<<std::endl;
		cluster C1(inC.at(0));
		//std::cout<<"just removed first elemenet"<<std::endl;
		if (inC.size() == 1)
		{
			//std::cout<<"only one element left"<<std::endl;
			outC.push_back(C1);	
			inC.erase(inC.begin()+0);
		}
		else
		{
			//std::cout<<"looking for pair:"<<std::endl;
			std::pair<int, float> BestPair = FindNearest(C1, inC);
			//std::cout<<"Nearest index and dR: "<<BestPair.first<<", "<<BestPair.second<<std::endl;
			if (BestPair.second > R)
			{
				cluster C2(inC.at(BestPair.first));
				outC.push_back(C1);
				outC.push_back(C2);
				inC.erase(inC.begin()+BestPair.first);
				inC.erase(inC.begin()+0);
			}
			else
			{
				cluster C2(inC.at(BestPair.first));
				C1.combine(C2);
				outC.push_back(C1);
				inC.erase(inC.begin()+BestPair.first);
				inC.erase(inC.begin()+0);
			}
		}
	}
	return outC;
}

bool comparator(const cluster& C1, const cluster& C2)
{
	return C1.vec.Mag() > C2.vec.Mag();
}

