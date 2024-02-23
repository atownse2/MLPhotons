// Description: Class to store ML photon information
//

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

// #include "DataFormats/Math/interface/LorentzVector.h"

namespace reco {

	class MLPhoton : public RecoCandidate {
		
		public:

			// Default constructor
			MLPhoton() : RecoCandidate() {}

			MLPhoton(const PolarLorentzVector& p4, const Point vtx);

			~MLPhoton();

			// Getter methods
			float moe() const {return moe_;}

			float diphoton_score() const {return diphoton_score_;}
			float monophoton_score() const {return monophoton_score_;}
			float hadron_score() const {return hadron_score_;}

			float r1() const {return r1_;}
			float r2() const {return r2_;}
			float r3() const {return r3_;}

			// Setter methods
			void set_moe(float moe) {moe_ = moe;}

			void set_diphoton_score(float score) {diphoton_score_ = score;}
			void set_monophoton_score(float score) {monophoton_score_ = score;}
			void set_hadron_score(float score) {hadron_score_ = score;}

			void set_r1(float r1) {r1_ = r1;}
			void set_r2(float r2) {r2_ = r2;}
			void set_r3(float r3) {r3_ = r3;}


		private:

			bool overlap(const Candidate&) const override;		
			
			PolarLorentzVector p4_;

			float moe_;
			float energy_;
			float eta_;
			float phi_;

			float diphoton_score_;
			float monophoton_score_;
			float hadron_score_;

			float r1_;
			float r2_;
			float r3_;

			// float pfIso;

	};

}