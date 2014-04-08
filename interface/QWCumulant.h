#include <TComplex.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom3.h>
#include <TFile.h>
#include "QWConst.h"
//
// constants, enums and typedefs
//

//#define QW_DEBUG 1
//#define QW_PEREVENT 1

#define PR(x) cout << "!!QW!! " << __LINE__ << " DEBUG OUTPUT " << (#x) << " = " << (x) << endl;
#define PRD(x) cout << "!!QW!! " << __LINE__ << " DEBUG OUTPUT " << (#x) << endl;
//
// class declaration
//

const int NMAX_TRK = 5000;
typedef struct QWEvent_ {
	int     Cent;
	int     Mult;
	double  vz;
	double  Pt[NMAX_TRK];
	double  Eta[NMAX_TRK];
	double  Phi[NMAX_TRK];
	int     Charge[NMAX_TRK];
	double	rEff[NMAX_TRK];
	double	rFak[NMAX_TRK];
} QWEvent;

class QWCumuIO {
public:
	QWCumuIO() : M(0), S1(0), S2(0), S3(0), S4(0) {;};
	void 		AddParticle(double, double ww = 1.);
	void		FinishEvent();

	double		M;
	double		S1;
	double		S2;
	double		S3;
	double		S4;
	double		Q2[NOrder];
	double		Q4[NOrder];
	TComplex	Q[NOrder];
	TComplex	Qw1[NOrder];
	TComplex	Qw2[NOrder];
	TComplex	Qw3[NOrder];
	TComplex	Qstar[NOrder];
	TComplex	Qstarw1[NOrder];
	TComplex	Qstarw2[NOrder];
	TComplex	Qstarw3[NOrder];
};

///////////////// Class ////////////////////////////

class QWQFlow {
public:
	explicit QWQFlow(TFileDirectory);
	void 	Finalize(int, QWCumuIO* , QWCumuIO* , QWCumuIO*);
//	void	Fill(int cent, double pt, double eta);
private:
	TH1D	* hC2[N_vn];
	TH1D	* hC4[N_vn];
	TH1D	* hC6[N_vn];
	TH1D	* hC8[N_vn];
	TH1D	* hC2p[N_vn];
	TH1D	* hC4p[N_vn];
	TH1D	* hC6p[N_vn];
	TH1D	* hC8p[N_vn];

	TH1D	* hC2w[N_vn];
	TH1D	* hC4w[N_vn];

	TH1D	* hC2x4[N_vn];
	TH1D	* hC2x2r[N_vn];
	TH1D	* hC2x4r[N_vn];
	TH1D	* hC4x2r[N_vn];
	TH1D	* hC4x4r[N_vn];

	TH1D	* hC2x6[N_vn];
	TH1D	* hC4x6[N_vn];
	TH1D	* hC2x8[N_vn];
	TH1D	* hC4x8[N_vn];
	TH1D	* hC6x8[N_vn];

	TH1D	* hReQ[N_vn];
	TH1D	* hImQ[N_vn];

	TH1D	* hReQ2[N_vn];
	TH1D	* hImQ2[N_vn];

	TH1D	* hReQ31[N_vn];
	TH1D	* hImQ31[N_vn];
	TH1D	* hReQ32[N_vn];
	TH1D	* hImQ32[N_vn];

	TH1D	* hW22;
	TH1D	* hW42;
	TH1D	* hW62;
	TH1D	* hW82;

	TH1D	* hW2x4;
	TH1D	* hW2x2r;
	TH1D	* hW2x4r;
	TH1D	* hW4x2r;
	TH1D	* hW4x4r;

	TH1D	* hW2x6;
	TH1D	* hW4x6;
	TH1D	* hW2x8;
	TH1D	* hW4x8;
	TH1D	* hW6x8;

	TH1D	* hW;
	TH1D	* hW2;
	TH1D	* hW3;
	TH1D	* hW4;
	TH1D	* hW5;
	TH1D	* hW6;
	TH1D	* hW7;
	TH1D	* hW8;

	TH1D	* hM2;
	TH1D	* hM4;

	TH1D	* hCNT;

};

//////////////////////////////////////////////////////////////////////
class QWCumulant : public edm::EDAnalyzer {
	public:
		explicit QWCumulant(const edm::ParameterSet&);
		~QWCumulant();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() ;
		virtual void analyze(const edm::Event&, const edm::EventSetup&);
		virtual void endJob() ;

		virtual void beginRun(edm::Run const&, edm::EventSetup const&);
		virtual void endRun(edm::Run const&, edm::EventSetup const&);
		virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
		virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

	/////////////////////////////////////////////
		double GenPhi();
		int getNoffCent(const edm::Event&, const edm::EventSetup&, int& Noff);
		//void Cumu(QWCumuIO*);
		void dumpCumu(QWCumuIO*);
		TRandom3 * gRandom;
		//QWCumuIO * InclusiveQ();
		QWCumuIO * CumuDBG();
		// ----------member data ---------------------------
		edm::InputTag tracks_; //used to select what tracks to read from configuration file
		edm::InputTag centrality_;	// centrality
		edm::InputTag vertexSrc_;
		edm::InputTag correctHist_;

		edm::InputTag fweight_;
		edm::InputTag facceptance_;
	/////////////////////////////////////////////
		double 	minvz_, maxvz_;
		double 	dzdzerror_;
		//double 	d0d0error_;
		double 	chi2_;
		double 	pterrorpt_;
		double 	rfpmineta_, rfpmaxeta_;
		double 	poimineta_, poimaxeta_;
		double 	rfpptmin_, rfpptmax_;
		double 	poiptmin_, poiptmax_;
		int 	charge_;
		bool	simu_;
		double	simuv2_;
		double 	simuv3_;

		bool	bFak;
		bool	bEff;
		bool	bacc;
		bool	bPhiEta;
		bool	bCentNoff;
//		bool	bHLTHM2013;
		int 	Noffmin_;
		int 	Noffmax_;

		QWEvent * t;
		TFile	* fEffFak;
		TFile	* facc;
	/////////////////////////////////////////////
		TH1I * hMult[nCentBins];
		TH1I * hNoff[nCentBins];
		TH1D * hPt[nCentBins];
		TH2D * hPhiEta[nCentBins][nPtBins][2];

		TH2D	* hdNdPtdEta[20];
		TH2D	* hdNdPtdEtaPt[20];

		TH2D * hEff_cbin[20];
		TH2D * hFak_cbin[20];

		TH2D * hacc[nCentBins][nPtBins][2];

		QWQFlow * QRFP;
		QWQFlow * QPt[nPtBins];
		QWQFlow * QEta[nEtaBins];
};



