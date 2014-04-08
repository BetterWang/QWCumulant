// -*- C++ -*-
//
// Package:    QWCumulant
// Class:      QWCumulant
// 
/**\class QWCumulant QWCumulant.cc QWAna/QWCumulant/src/QWCumulant.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Quan Wang
//         Created:  Tue Oct 16 16:33:30 EDT 2012
// $Id: QWCumulant.cc,v 1.1 2013/01/15 15:56:58 qwang Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TH1.h"
#include "TH2.h"
#include "TComplex.h"


#include "QWAna/QWCumulant/interface/QWCumulant.h"


using namespace std;

//#ifdef QW_DEBUG
//
// constructors and destructor
//
QWCumulant::QWCumulant(const edm::ParameterSet& iConfig)
	:
		tracks_(iConfig.getUntrackedParameter<edm::InputTag>("tracks_"))
	,	centrality_(iConfig.getParameter<edm::InputTag>("centrality_"))
	,	vertexSrc_(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc_"))
	,	bacc(false)
{
	//now do what ever initialization is needed
	minvz_ = iConfig.getUntrackedParameter<double>("minvz_", -15.);
	maxvz_ = iConfig.getUntrackedParameter<double>("maxvz_", 15.);
	dzdzerror_ = iConfig.getUntrackedParameter<double>("dzdzerror_", 10.);
	//d0d0error_ = iConfig.getUntrackedParameter<double>("d0d0error_", 3.);
	chi2_ = iConfig.getUntrackedParameter<double>("chi2_", 40);
	pterrorpt_ = iConfig.getUntrackedParameter<double>("pterrorpt_", 0.05);
	rfpmineta_ = iConfig.getUntrackedParameter<double>("rfpmineta_", -2.);
	rfpmaxeta_ = iConfig.getUntrackedParameter<double>("rfpmaxeta_", 2.);
	poimineta_ = iConfig.getUntrackedParameter<double>("poimineta_", rfpmineta_);
	poimaxeta_ = iConfig.getUntrackedParameter<double>("poimaxeta_", rfpmaxeta_);
	fweight_ = iConfig.getUntrackedParameter<edm::InputTag>("fweight_", string("NA"));
	facceptance_ = iConfig.getUntrackedParameter<edm::InputTag>("facceptance_", string("NA"));
	rfpptmin_ = iConfig.getUntrackedParameter<double>("rfpptmin_", 0.1);
	rfpptmax_ = iConfig.getUntrackedParameter<double>("rfpptmax_", 100);
	poiptmin_ = iConfig.getUntrackedParameter<double>("poiptmin_", rfpptmin_);
	poiptmax_ = iConfig.getUntrackedParameter<double>("poiptmax_", rfpptmax_);
	charge_ = iConfig.getUntrackedParameter<int>("charge_", 0);
	simu_ = iConfig.getUntrackedParameter<bool>("simu_", false);
	simuv2_ = iConfig.getUntrackedParameter<double>("simuv2_", 0.);
	simuv3_ = iConfig.getUntrackedParameter<double>("simuv3_", 0.);
	bFak = iConfig.getUntrackedParameter<bool>("bFak_", false);
	bEff = iConfig.getUntrackedParameter<bool>("bEff_", false);
	bPhiEta = iConfig.getUntrackedParameter<bool>("bPhiEta_", false);
	bCentNoff = iConfig.getUntrackedParameter<bool>("bCentNoff_", false);
//	bHLTHM2013 = iConfig.getUntrackedParameter<bool>("bHLTHM2013_", false);
	Noffmin_ = iConfig.getUntrackedParameter<int>("Noffmin_", 0);
	Noffmax_ = iConfig.getUntrackedParameter<int>("Noffmax_", 10000);

	if ( simu_ ) gRandom = new TRandom3;

	string streff = fweight_.label();
	if ( streff == string("NA") ) {
		bFak = false;
		bEff = false;
		fEffFak = 0;
	} else {
		fEffFak = new TFile(streff.c_str());
		if ( !fEffFak->IsOpen() ) {
			bFak = false;
			bEff = false;
		} else {
			cout << "!!! Using particle weight " << streff << endl;
			if ( bFak ) cout << "!!! Apply Fak correction" << endl;
			if ( bEff ) cout << "!!! Apply Eff correction" << endl;
			for ( int i = 0; i < 20; i++ ) {
				hEff_cbin[i] = (TH2D*) fEffFak->Get("rTotalEff3D");
				hFak_cbin[i] = (TH2D*) fEffFak->Get(Form("rFak_cbin%i", i));
			}
		}
	}
	string stracc = facceptance_.label();
	if ( stracc == string("NA") ) {
		bacc = false;
		facc = 0;
	} else {
		facc = new TFile(stracc.c_str());
		if ( !facc->IsOpen() ) {
			bacc = false;
		} else {
			cout << "!!! Using acceptance weight " << stracc << endl;
			bacc = true;
			for ( int cent = 0; cent < nCentBins; cent++ ) {
				for ( int ipt = 0; ipt < nPtBins; ipt++ ) {
					hacc[cent][ipt][0] = (TH2D*) facc->Get(Form("hPhiEta_%i_%i_0", cent, ipt));
					hacc[cent][ipt][1] = (TH2D*) facc->Get(Form("hPhiEta_%i_%i_1", cent, ipt));
				}
			}
		}
	}


	//
	//cout << __LINE__ << "\t" << tracks_.label().c_str() << "\t|" << tracks_.instance() << "\t|" << tracks_.process() << endl;
	//
	t = new QWEvent;
	memset(t, 0, sizeof(QWEvent));
	//
	//
	edm::Service<TFileService> fs;
	for ( int cent = 0; cent < nCentBins; cent++ ) {
		hMult[cent]	= fs->make<TH1I>(Form("hMult_%i", cent),"", 5000, 0, 5000);
		hNoff[cent]	= fs->make<TH1I>(Form("hNoff_%i", cent),"", 5000, 0, 5000);
		hPt[cent]	= fs->make<TH1D>(Form("hPt_%i", cent), "", 20000, 0, 100);
		if ( bPhiEta ) {
			for ( int i = 0; i < nPtBins; i++ ) {
				hPhiEta[cent][i][0] = fs->make<TH2D>(Form("hPhiEta_%i_%i_0", cent, i), "", 512, -Pi, Pi, 480, -2.4, 2.4);
				hPhiEta[cent][i][1] = fs->make<TH2D>(Form("hPhiEta_%i_%i_1", cent, i), "", 512, -Pi, Pi, 480, -2.4, 2.4);
			}
		}
	}

	for ( int cent = 0; cent < 20; cent++ ) {
		hdNdPtdEta[cent] = fs->make<TH2D>(Form("hdNdPtdEta_%i", cent), Form("hdNdPtdEta_%i", cent), nEtaBins, etabins, 38, fakpt );
		hdNdPtdEtaPt[cent] = fs->make<TH2D>(Form("hdNdPtdEtaPt_%i", cent), Form("hdNdPtdEta_%i", cent), nEtaBins, etabins, 38, fakpt );
	}

	TFileDirectory rfp = fs->mkdir("RFP");
	QRFP = new QWQFlow(rfp);

	for ( int i = 0; i < nPtBins; i++ ) {
		TFileDirectory dir = fs->mkdir(Form("POI_pt_%i_%i", int(10*ptbins[i]), int(10*ptbins[i+1]) ));
		QPt[i] = new QWQFlow(dir);
	}

	for ( int i = 0; i < nEtaBins; i++ ) {
		TFileDirectory dir = fs->mkdir(Form("POI_eta_%i_%i", int(10*etabins[i]), int(10*etabins[i+1]) ));
		QEta[i] = new QWQFlow(dir);
	}


}


QWCumulant::~QWCumulant()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

int
QWCumulant::getNoffCent(const edm::Event& iEvent, const edm::EventSetup& iSetup, int& Noff)
{
	// very hard coded Noff track centrality cut
	using namespace edm;
	using namespace reco;
//	int Noff = 0;

	Handle<VertexCollection> vertexCollection;
	iEvent.getByLabel(vertexSrc_, vertexCollection);
	const VertexCollection * recoVertices = vertexCollection.product();

	int primaryvtx = 0;
	math::XYZPoint v1( (*recoVertices)[primaryvtx].position().x(), (*recoVertices)[primaryvtx].position().y(), (*recoVertices)[primaryvtx].position().z() );
	double vxError = (*recoVertices)[primaryvtx].xError();
	double vyError = (*recoVertices)[primaryvtx].yError();
	double vzError = (*recoVertices)[primaryvtx].zError();


	Handle<TrackCollection> tracks;
	iEvent.getByLabel(tracks_,tracks);
	for(TrackCollection::const_iterator itTrack = tracks->begin();
		itTrack != tracks->end();                      
		++itTrack) {

		if ( !itTrack->quality(reco::TrackBase::highPurity) ) continue;
		if ( itTrack->charge() == 0 ) continue;
		if ( itTrack->pt() < 0.4 ) continue;

		double d0 = -1.* itTrack->dxy(v1);
		double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
		double dz=itTrack->dz(v1);
		double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
		if ( fabs(itTrack->eta()) > 2.4 ) continue;
		if ( fabs( dz/dzerror ) > 3. ) continue;
		if ( fabs( d0/derror ) > 3. ) continue;
		if ( itTrack->ptError()/itTrack->pt() > pterrorpt_ ) continue;
//		bool b_pix = itTrack->numberOfValidHits() < 7;
//		if ( b_pix ) {
//			if ( fabs( dz/dzerror ) > dzdzerror_ ) continue;
//			if ( itTrack->normalizedChi2() > chi2_ ) continue;
//		} else {
//			// full track
//			if ( fabs( dz/dzerror ) > 3. ) continue;
//			if ( fabs( d0/derror ) > 3. ) continue;
//			if ( itTrack->ptError()/itTrack->pt() > pterrorpt_ ) continue;
//			if ( itTrack->numberOfValidHits() < 12 ) continue;
//		}

		Noff++;
	}

	int cent = nCentNoff-1;
	while ( CentNoffCut[cent] <= Noff ) cent--;
	hNoff[cent]->Fill(Noff);
	return cent;
}
// ------------ method called for each event  ------------
	void
QWCumulant::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace reco;

	// vertex
	Handle<VertexCollection> vertexCollection;
	iEvent.getByLabel(vertexSrc_, vertexCollection);
	const VertexCollection * recoVertices = vertexCollection.product();

	int primaryvtx = 0;
	math::XYZPoint v1( (*recoVertices)[primaryvtx].position().x(), (*recoVertices)[primaryvtx].position().y(), (*recoVertices)[primaryvtx].position().z() );
	double vxError = (*recoVertices)[primaryvtx].xError();
	double vyError = (*recoVertices)[primaryvtx].yError();
	double vzError = (*recoVertices)[primaryvtx].zError();

//	for ( unsigned int i = 0; i < recoVertices->size(); i++ ) {
//		size_t daughter = (*recoVertices)[i].tracksSize();
//		cout << "i = " << i << "\tnTracks = " << daughter <<"\t vz = " << (*recoVertices)[i].position().z() << endl;
//		//cout << "i = " << i << "\ttrkSize = " << "\t vz = " << (*recoVertices)[i].position().z() << endl;
//	}
	double vz = (*recoVertices)[primaryvtx].z();
	if (vz < minvz_ || vz > maxvz_) {
		return;
	}
	
	// centrality
	int bin = 0;
	int cbin = 0;
	int Noff = 0;

	if ( bCentNoff ) {
		cbin = getNoffCent( iEvent, iSetup, Noff);
		if ( (Noff < Noffmin_) or (Noff >= Noffmax_) ) {
			return;
		}
	} else {
		edm::Handle<int> ch;
		iEvent.getByLabel(centrality_,ch);
		bin = *(ch.product());
		while ( centbins[cbin+1] < bin*2.5+0.1 ) cbin++;
	}
	bin = cbin;

	// track
	Handle<TrackCollection> tracks;
	iEvent.getByLabel(tracks_,tracks);
	t->Mult = 0;
	t->Cent = cbin;
	t->vz = vz;
	//cout << __LINE__ << "\t" << cbin << endl;
	QWCumuIO Qinc;
	QWCumuIO Qpt[nPtBins];
	QWCumuIO Qptq[nPtBins];
	QWCumuIO Qeta[nEtaBins];
	QWCumuIO Qetaq[nEtaBins];

	for(TrackCollection::const_iterator itTrack = tracks->begin();
			itTrack != tracks->end();                      
			++itTrack) {
		if ( itTrack->charge() == 0 ) continue;

		double d0 = -1.* itTrack->dxy(v1);
		double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
		double dz=itTrack->dz(v1);
		double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);

		if ( fabs(itTrack->eta()) > 2.4 ) continue;
		if ( fabs( dz/dzerror ) > 3. ) continue;
		if ( fabs( d0/derror ) > 3. ) continue;
		if ( itTrack->ptError()/itTrack->pt() > pterrorpt_ ) continue;
		/*
		bool b_pix = itTrack->numberOfValidHits() < 7;
		if ( b_pix ) {
			if ( fabs( dz/dzerror ) > dzdzerror_ ) continue;
			if ( itTrack->normalizedChi2() > chi2_ ) continue;
		} else {
			// full track
			if ( fabs( dz/dzerror ) > 3. ) continue;
			if ( fabs( d0/derror ) > 3. ) continue;
			if ( itTrack->ptError()/itTrack->pt() > pterrorpt_ ) continue;
			if ( itTrack->numberOfValidHits() < 12 ) continue;
		}
		*/

		t->Charge[t->Mult] = itTrack->charge();
		if ( (charge_ == 1) && (t->Charge[t->Mult]<0) ) continue;
		if ( (charge_ == -1) && (t->Charge[t->Mult]>0) ) continue;

		t->Pt[t->Mult] = itTrack->pt();
		if ( t->Pt[t->Mult] >= ptbins[nPtBins] || t->Pt[t->Mult] <= ptbins[0] ) continue;
		t->Eta[t->Mult] = itTrack->eta();
		if ( (t->Eta[t->Mult] < etabins[0]) || (t->Eta[t->Mult] >= etabins[nEtaBins]) ) continue;
		bool bRFP = true;
		bool bPOIpt = true;
		bool bPOIeta = true;

		if ( bEff ) {
			t->rEff[t->Mult] = hEff_cbin[bin]->GetBinContent( hEff_cbin[bin]->FindBin(t->Eta[t->Mult], t->Pt[t->Mult] ) );
		} else {
			t->rEff[t->Mult] = 1.;
		}
		if ( bFak ) {
			t->rFak[t->Mult] = hFak_cbin[bin]->GetBinContent( hFak_cbin[bin]->FindBin(t->Eta[t->Mult], t->Pt[t->Mult] ) );
		} else {
			t->rFak[t->Mult] = 0.;
		}
		if ( t->rEff[t->Mult] <= 0.1 or TMath::IsNaN(t->rEff[t->Mult]) ) continue;
		double weight = (1.-t->rFak[t->Mult])/t->rEff[t->Mult];

		double phi = itTrack->phi();
		if ( simu_ ) {
			phi = GenPhi();
		}

		double wacc = 1.;
		int ipt=0;
		while ( t->Pt[t->Mult] > ptbins[ipt+1] ) ipt++;
		if ( bacc ) {
			wacc = 1./hacc[t->Cent][ipt][t->Charge[t->Mult]>0]->GetBinContent(hacc[t->Cent][ipt][t->Charge[t->Mult]>0]->FindBin(phi, t->Eta[t->Mult]));
		}
		if ( bPhiEta ) hPhiEta[t->Cent][ipt][t->Charge[t->Mult]>0]->Fill(phi, t->Eta[t->Mult], wacc);

		weight *= wacc;

		if ( (t->Pt[t->Mult] < rfpptmin_) || (t->Pt[t->Mult] > rfpptmax_) || itTrack->eta() < rfpmineta_ || itTrack->eta() > rfpmaxeta_ ) bRFP = false;
		if ( (itTrack->eta() < poimineta_) || (itTrack->eta() > poimaxeta_) ) bPOIpt = false;
		if ( (t->Pt[t->Mult] < poiptmin_) || (t->Pt[t->Mult] > poiptmax_) ) bPOIeta = false;

		if (bRFP) {
			Qinc.AddParticle(phi, weight);
		}
		if (bPOIpt) {
			int npt=0;
			while ( t->Pt[t->Mult] > ptbins[npt+1] ) npt++;
			Qpt[npt].AddParticle(phi, weight);
			if ( bRFP ) Qptq[npt].AddParticle(phi, weight);
		}
		if (bPOIeta) {
			int neta=0;
			while ( t->Eta[t->Mult] > etabins[neta+1] ) neta++;
			Qeta[neta].AddParticle(phi, weight);
			if ( bRFP ) Qetaq[neta].AddParticle(phi, weight);
		}

		hdNdPtdEta[bin]->Fill(t->Eta[t->Mult], t->Pt[t->Mult]);
		hdNdPtdEtaPt[bin]->Fill(t->Eta[t->Mult], t->Pt[t->Mult], t->Pt[t->Mult]);

		t->Phi[t->Mult] = phi;
		hPt[t->Cent]->Fill(t->Pt[t->Mult]);

		t->Mult++;
	}
	hMult[t->Cent]->Fill(t->Mult);
	Qinc.FinishEvent();
	for ( int i = 0; i < nPtBins; i++ ) {
		Qptq[i].FinishEvent();
		Qpt[i].FinishEvent();
	}
	for ( int i = 0; i < nEtaBins; i++ ) {
		Qetaq[i].FinishEvent();
		Qeta[i].FinishEvent();
	}

	QRFP->Finalize(t->Cent, &Qinc, &Qinc, &Qinc);
	for ( int i = 0; i < nPtBins; i++ ) {
		QPt[i]->Finalize(t->Cent, &Qinc, &Qpt[i], &Qptq[i]);
	}
	for ( int i = 0; i < nEtaBins; i++ ) {
		QEta[i]->Finalize(t->Cent, &Qinc, &Qeta[i], &Qetaq[i]);
	}
}




// ------------------ QWCumuIO -------------------------
void
QWCumuIO::AddParticle(double phi, double ww)
{
	M++;
	S1 += ww;
	S2 += ww*ww;
	S3 += ww*ww*ww;
	S4 += ww*ww*ww*ww;
	for ( int i = 1; i < NOrder; i++ ) {
		Q[i] += TComplex(cos(i*phi), sin(i*phi));
		Qw1[i] += TComplex(cos(i*phi), sin(i*phi))*ww;
		Qw2[i] += TComplex(cos(i*phi), sin(i*phi))*ww*ww;
		Qw3[i] += TComplex(cos(i*phi), sin(i*phi))*ww*ww*ww;
	}
}

void
QWCumuIO::FinishEvent()
{
	for ( int i = 1; i < NOrder; i++ ) {
		Qstar[i] = TComplex::Conjugate(Q[i]);
		Qstarw1[i] = TComplex::Conjugate(Qw1[i]);
		Qstarw2[i] = TComplex::Conjugate(Qw2[i]);
		Qstarw3[i] = TComplex::Conjugate(Qw3[i]);

		Q2[i] = TMath::Power(TComplex::Abs(Q[i]), 2);
		Q4[i] = TMath::Power(Q2[i], 2);
	}
}


// ------------  generate phi -----------------------------------
double
QWCumulant::GenPhi()
{
	double max = 1+2*simuv2_+2*simuv3_;
	double phi = gRandom->Uniform(-Pi, Pi);
	while ( gRandom->Uniform(max) > 1 + 2*simuv2_*cos(2*phi) + 2*simuv3_*cos(3*phi) ) phi = gRandom->Uniform(-Pi, Pi);
	return phi;
}

// ------------ method called once each job just before starting event loop  ------------
	void 
QWCumulant::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
QWCumulant::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
	void 
QWCumulant::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
	void 
QWCumulant::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
	void 
QWCumulant::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
	void 
QWCumulant::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
QWCumulant::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);

	//Specify that only 'tracks' is allowed
	//To use, remove the default given above and uncomment below
	//ParameterSetDescription desc;
	//desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
	//descriptions.addDefault(desc);
}

//////////////////////////////////////////
//void 
//QWCumulant::Cumu(QWCumuIO *q)
//{
//	for ( unsigned int i = 0; i < q->M; i++ ) {
//		for( int n = 1; n <= 8; n++ ) {
//			q->Q[n] += q->w[i]*TComplex(cos( (n)*q->Phi[i]), sin( (n)*q->Phi[i]) );
//		}
//	}
//	for( int n = 1; n < 9; n++ ) {
//		q->Qstar[n] = TComplex::Conjugate(q->Q[n]);
//		q->Q2[n] = TMath::Power(TComplex::Abs(q->Q[n]), 2);
//		q->Q4[n] = TMath::Power(q->Q2[n], 2);
//	}
//
//}
//
//////////////////////////////////////////
void
QWCumulant::dumpCumu(QWCumuIO *q)
{
	/*
	cout << __LINE__ << "\tDump q " << endl;
	cout << "M = " << q->M << endl;
	for ( int i = 1; i <= 8; i++ ) {
		cout << "Q[" << i << "] = " << q->Q[i] 
			<< " \t" << "Q*[" << i << "] = " << q->Qstar[i]
			<< " \t" << "Q2[" << i << "] = " << q->Q2[i]
			<< " \t" << "Q4[" << i << "] = " << q->Q4[i] << endl;
	}
	*/
	return;
}

/////////////////// QWQFlow ///////////////////////
QWQFlow::QWQFlow(TFileDirectory dir)
{
	for ( int i = 1; i < N_vn; i++ ) {
		hC2[i] = dir.make<TH1D>(Form("hC2_%i", i), Form("hC2_%i", i), nCentBins, 0, nCentBins);
		hC4[i] = dir.make<TH1D>(Form("hC4_%i", i), Form("hC4_%i", i), nCentBins, 0, nCentBins);
		hC6[i] = dir.make<TH1D>(Form("hC6_%i", i), Form("hC6_%i", i), nCentBins, 0, nCentBins);
		hC8[i] = dir.make<TH1D>(Form("hC8_%i", i), Form("hC8_%i", i), nCentBins, 0, nCentBins);
		hC2p[i] = dir.make<TH1D>(Form("hC2p_%i", i), Form("hC2_%i", i), nCentBins, 0, nCentBins);
		hC4p[i] = dir.make<TH1D>(Form("hC4p_%i", i), Form("hC4_%i", i), nCentBins, 0, nCentBins);
		hC6p[i] = dir.make<TH1D>(Form("hC6p_%i", i), Form("hC6_%i", i), nCentBins, 0, nCentBins);
		hC8p[i] = dir.make<TH1D>(Form("hC8p_%i", i), Form("hC8_%i", i), nCentBins, 0, nCentBins);

		hC2w[i] = dir.make<TH1D>(Form("hC2w_%i", i), Form("hC2w_%i", i), nCentBins, 0, nCentBins);
		hC4w[i] = dir.make<TH1D>(Form("hC4w_%i", i), Form("hC4w_%i", i), nCentBins, 0, nCentBins);

		hC2x4[i] = dir.make<TH1D>(Form("hC2x4_%i", i), Form("hC2x4_%i", i), nCentBins, 0, nCentBins);
		hC2x2r[i] = dir.make<TH1D>(Form("hC2x2r_%i", i), Form("hC2x2r_%i", i), nCentBins, 0, nCentBins);
		hC2x4r[i] = dir.make<TH1D>(Form("hC2x4r_%i", i), Form("hC2x4r_%i", i), nCentBins, 0, nCentBins);
		hC4x2r[i] = dir.make<TH1D>(Form("hC4x2r_%i", i), Form("hC4x2r_%i", i), nCentBins, 0, nCentBins);
		hC4x4r[i] = dir.make<TH1D>(Form("hC4x4r_%i", i), Form("hC4x4r_%i", i), nCentBins, 0, nCentBins);

		hC2x6[i] = dir.make<TH1D>(Form("hC2x6_%i", i), Form("hC2x6_%i", i), nCentBins, 0, nCentBins);
		hC4x6[i] = dir.make<TH1D>(Form("hC4x6_%i", i), Form("hC4x6_%i", i), nCentBins, 0, nCentBins);
		hC2x8[i] = dir.make<TH1D>(Form("hC2x8_%i", i), Form("hC2x8_%i", i), nCentBins, 0, nCentBins);
		hC4x8[i] = dir.make<TH1D>(Form("hC4x8_%i", i), Form("hC4x8_%i", i), nCentBins, 0, nCentBins);
		hC6x8[i] = dir.make<TH1D>(Form("hC6x8_%i", i), Form("hC6x8_%i", i), nCentBins, 0, nCentBins);

		hReQ[i] = dir.make<TH1D>(Form("hReQ_%i", i), Form("hReQ_%i", i), nCentBins, 0, nCentBins);
		hImQ[i] = dir.make<TH1D>(Form("hImQ_%i", i), Form("hImQ_%i", i), nCentBins, 0, nCentBins);
		hReQ2[i] = dir.make<TH1D>(Form("hReQ2_%i", i), Form("hReQ2_%i", i), nCentBins, 0, nCentBins);
		hImQ2[i] = dir.make<TH1D>(Form("hImQ2_%i", i), Form("hImQ2_%i", i), nCentBins, 0, nCentBins);
		hReQ31[i] = dir.make<TH1D>(Form("hReQ31_%i", i), Form("hReQ31_%i", i), nCentBins, 0, nCentBins);
		hImQ31[i] = dir.make<TH1D>(Form("hImQ31_%i", i), Form("hImQ31_%i", i), nCentBins, 0, nCentBins);
		hReQ32[i] = dir.make<TH1D>(Form("hReQ32_%i", i), Form("hReQ32_%i", i), nCentBins, 0, nCentBins);
		hImQ32[i] = dir.make<TH1D>(Form("hImQ32_%i", i), Form("hImQ32_%i", i), nCentBins, 0, nCentBins);
	}

	hW2x4 = dir.make<TH1D>("hW2x4", "hW2x4", nCentBins, 0, nCentBins);
	hW2x2r = dir.make<TH1D>("hW2x2r", "hW2x2r", nCentBins, 0, nCentBins);
	hW2x4r = dir.make<TH1D>("hW2x4r", "hW2x4r", nCentBins, 0, nCentBins);
	hW4x2r = dir.make<TH1D>("hW4x2r", "hW4x2r", nCentBins, 0, nCentBins);
	hW4x4r = dir.make<TH1D>("hW4x4r", "hW4x4r", nCentBins, 0, nCentBins);

	hW2x6 = dir.make<TH1D>("hW2x6", "hW2x6", nCentBins, 0, nCentBins);
	hW4x6 = dir.make<TH1D>("hW4x6", "hW4x6", nCentBins, 0, nCentBins);
	hW2x8 = dir.make<TH1D>("hW2x8", "hW2x8", nCentBins, 0, nCentBins);
	hW4x8 = dir.make<TH1D>("hW4x8", "hW4x8", nCentBins, 0, nCentBins);
	hW6x8 = dir.make<TH1D>("hW6x8", "hW6x8", nCentBins, 0, nCentBins);

	hW22 = dir.make<TH1D>("hW22", "hW22", nCentBins, 0, nCentBins);
	hW42 = dir.make<TH1D>("hW42", "hW42", nCentBins, 0, nCentBins);
	hW62 = dir.make<TH1D>("hW62", "hW62", nCentBins, 0, nCentBins);
	hW82 = dir.make<TH1D>("hW82", "hW82", nCentBins, 0, nCentBins);

	hW = dir.make<TH1D>("hW", "hW", nCentBins, 0, nCentBins);
	hW2 = dir.make<TH1D>("hW2", "hW2", nCentBins, 0, nCentBins);
	hW3 = dir.make<TH1D>("hW3", "hW3", nCentBins, 0, nCentBins);
	hW4 = dir.make<TH1D>("hW4", "hW4", nCentBins, 0, nCentBins);
	hW5 = dir.make<TH1D>("hW5", "hW5", nCentBins, 0, nCentBins);
	hW6 = dir.make<TH1D>("hW6", "hW6", nCentBins, 0, nCentBins);
	hW7 = dir.make<TH1D>("hW7", "hW7", nCentBins, 0, nCentBins);
	hW8 = dir.make<TH1D>("hW8", "hW8", nCentBins, 0, nCentBins);

	hM2 = dir.make<TH1D>("hM2", "hM2", nCentBins, 0, nCentBins);
	hM4 = dir.make<TH1D>("hM4", "hM4", nCentBins, 0, nCentBins);

	hCNT = dir.make<TH1D>("hCNT", "hCNT", nCentBins, 0, nCentBins);

}

/*
void
QWQFlow::Fill(int cent, double pt, double eta)
{
	hdNdPtdEta[cent]->Fill(eta, pt);
	hdNdPtdEtaPt[cent]->Fill(eta, pt, pt);
	return;
}
*/

void
QWQFlow::Finalize(int cent, QWCumuIO* R, QWCumuIO* p, QWCumuIO* q)
{
	double m = R->M;
	double mp = p->M;
	double mq = q->M;

	double w2 = mp*m-mq;
	double w3 = (mp*m-2*mq)*(m-1);
	double w4 = (mp*m-3*mq)*(m-1)*(m-2);
	double w5 = m*(m-1)*(m-2)*(m-3)*(m-4);
	double w6 = m*(m-1)*(m-2)*(m-3)*(m-4)*(m-5);
	double w7 = m*(m-1)*(m-2)*(m-3)*(m-4)*(m-5)*(m-6);
	double w8 = m*(m-1)*(m-2)*(m-3)*(m-4)*(m-5)*(m-6)*(m-7);
	double rw2 = m*(m-1);
	double rw4 = m*(m-1)*(m-2)*(m-3);
	double rw5 = m*(m-1)*(m-2)*(m-3)*(m-4);
	double rw6 = m*(m-1)*(m-2)*(m-3)*(m-4)*(m-5);
	double rw7 = m*(m-1)*(m-2)*(m-3)*(m-4)*(m-5)*(m-6);
	double rw8 = m*(m-1)*(m-2)*(m-3)*(m-4)*(m-5)*(m-6)*(m-7);

	if ( (w2 == 0) || (w4 == 0) || (w6==0) || (w8==0)) return;
	//if ( (w2 == 0) || (w4 == 0) ) return;
	hW->Fill(cent, mp);
	hW2->Fill(cent, w2);
	hW3->Fill(cent, w3);
	hW4->Fill(cent, w4);
	hW5->Fill(cent, w5);
	hW6->Fill(cent, w6);
	hW7->Fill(cent, w7);
	hW8->Fill(cent, w8);

	hCNT->Fill(cent);

	hW22->Fill(cent, TMath::Power(w2, 2));
	hW42->Fill(cent, TMath::Power(w4, 2));
	hW62->Fill(cent, TMath::Power(w6, 2));
	hW82->Fill(cent, TMath::Power(w8, 2));

	hW2x4->Fill(cent, w2*w4);
	hW2x2r->Fill(cent, w2*rw2);
	hW4x2r->Fill(cent, w4*rw2);
	hW2x4r->Fill(cent, w2*rw4);
	hW4x4r->Fill(cent, w4*rw4);

	hW2x6->Fill(cent, w2*w6);
	hW4x6->Fill(cent, w4*w6);
	hW2x8->Fill(cent, w2*w8);
	hW4x8->Fill(cent, w4*w8);
	hW6x8->Fill(cent, w6*w8);

	if ( R == p ) {
		hM2->Fill(cent, R->S1 * R->S1 - R->S2 );
		hM4->Fill(cent, pow(R->S1, 4) - 6 * R->S2 * R->S1 * R->S1 + 8 * R->S3 * R->S1 + 3 * R->S2 * R->S2 - 6 * R->S4 );
	} else {
		hM2->Fill(cent, mp * R->S1 - q->S1 );
		hM4->Fill(cent, mp * ( pow(R->S1, 3) - 3 * R->S1 * R->S2 + 2 * R->S3 )
				- 3 * ( q->S1 * (R->S1 * R->S1 - R->S2 ) + 2 * ( q->S3 - q->S2 * R->S1 ) ) );
	}

	for ( int n = 1; n < N_vn; n++ ) {
		TComplex c2 = p->Q[n]*R->Qstar[n] - q->M;
		hC2[n]->Fill(cent, c2.Re());
		hC2p[n]->Fill(cent, c2.Re()*c2.Re()/w2);

		TComplex c2w;
		if ( R == p ) c2w = R->Qw1[n]*R->Qstarw1[n] - R->S2;
		else c2w = p->Q[n]*R->Qstarw1[n] - q->S1;
		hC2w[n]->Fill(cent, c2w.Re());

		TComplex c4 = (p->Q[n] * R->Q[n] * R->Qstar[n] * R->Qstar[n])
			- (q->Q[2*n] * R->Qstar[n] * R->Qstar[n])
			- (p->Q[n] * R->Q[n] * R->Qstar[2*n])
			- (2. * m * p->Q[n] * R->Qstar[n])
			- (2. * mq * R->Q2[n])
			+ (7. * q->Q[n] * R->Qstar[n] )
			- (R->Q[n] * q->Qstar[n])
			+ (q->Q[2*n] * R->Qstar[2*n])
			+ (2. * p->Q[n] * R->Qstar[n])
			+ (2. * mq * m)
			- (6. * mq);
		//cout << __LINE__ << "\tc4 = " << c << endl;
		hC4[n]->Fill(cent, c4.Re());
		hC4p[n]->Fill(cent, c4.Re()*c4.Re()/w4);

		TComplex c6 = R->Q[n] * R->Qstar[n] * R->Q[n] * R->Qstar[n] * R->Q[n] * R->Qstar[n]
			+ 9. * R->Q[2*n] * R->Qstar[2*n] * R->Q[n] * R->Qstar[n]
			- 6. * R->Q[2*n] * R->Q[n] * R->Qstar[n] * R->Qstar[n] * R->Qstar[n]
			+ 4. * R->Q[3*n] * R->Qstar[n] * R->Qstar[n] * R->Qstar[n]
			- 12.* R->Q[3*n] * R->Qstar[2*n] * R->Qstar[n]
			+ 18. * (m-4) * R->Q[2*n] * R->Qstar[n] * R->Qstar[n]
			+ 4. * R->Q[3*n] * R->Qstar[3*n]
			- 9. * (m-4) * ( R->Q[n] * R->Qstar[n] * R->Q[n] * R->Qstar[n] + R->Q[2*n] * R->Qstar[2*n] )
			+ 18. * (m-2) * (m-5) * R->Q[n] * R->Qstar[n]
			- 6. * m * (m-4) * (m-5);
      		hC6[n]->Fill(cent, c6.Re());
      		hC6p[n]->Fill(cent, c6.Re()*c6.Re()/w6);

	      	TComplex c8 = R->Q[n] * R->Qstar[n] * R->Q[n] * R->Qstar[n] * R->Q[n] * R->Qstar[n] * R->Q[n] * R->Qstar[n]
			- 12. * R->Q[2*n] * R->Q[n] * R->Q[n] * R->Qstar[n] * R->Qstar[n] * R->Qstar[n] * R->Qstar[n]
			+ 6. * R->Q[2*n] * R->Q[2*n] * R->Qstar[n] * R->Qstar[n] * R->Qstar[n] * R->Qstar[n]
			+ 16. * R->Q[3*n] * R->Q[n] * R->Qstar[n] * R->Qstar[n] * R->Qstar[n] * R->Qstar[n]
			- 96. * R->Q[3*n] * R->Q[n] * R->Qstar[2*n] * R->Qstar[n] * R->Qstar[n]
			- 12. * R->Q[4*n] * R->Qstar[n] * R->Qstar[n] * R->Qstar[n] * R->Qstar[n]
			- 36. * R->Q[2*n] * R->Q[2*n] * R->Qstar[2*n] * R->Qstar[n] * R->Qstar[n]
			+ 96. * (m-6) * R->Q[2*n] * R->Q[n] * R->Qstar[n] * R->Qstar[n] * R->Qstar[n]
			+ 72. * R->Q[4*n] * R->Qstar[2*n] * R->Qstar[n] * R->Qstar[n]
			+ 48. * R->Q[3*n] * R->Q[n] * R->Qstar[2*n] * R->Qstar[2*n]
			- 64. * (m-6) * R->Q[3*n] * R->Qstar[n] * R->Qstar[n] * R->Qstar[n]
			+ 192. * (m-6) * R->Q[3*n] * R->Qstar[2*n] * R->Qstar[n]
			- 96. * R->Q[4*n] * R->Qstar[3*n] * R->Qstar[n]
			- 36. * R->Q[4*n] * R->Qstar[2*n] * R->Qstar[2*n]
			- 144. * (m-7) * (m-4) * R->Q[2*n] * R->Qstar[n] * R->Qstar[n]
			+ 36. * R->Q[4*n] * R->Qstar[4*n]
			+ 64. * R->Q[3*n] * R->Qstar[3*n] * R->Q[n] * R->Qstar[n]
			- 64. * (m-6) * R->Q[3*n] * R->Qstar[3*n]
			+ 9. * R->Q[2*n] * R->Qstar[2*n] * R->Q[2*n] * R->Qstar[2*n]
			+ 36. * R->Q[n] * R->Qstar[n] * R->Q[n] * R->Qstar[n] * R->Q[2*n] * R->Qstar[2*n]
			- 144. * (m-6) * R->Q[2*n] * R->Qstar[2*n] * R->Q[n] * R->Qstar[n]
			+ 72. * (m-7) * (m-4) * ( R->Q[2*n] * R->Qstar[2*n] + R->Q[n] * R->Qstar[n] * R->Q[n] * R->Qstar[n] )
			- 16. * (m-6) * R->Q[n] * R->Qstar[n] * R->Q[n] * R->Qstar[n] * R->Q[n] * R->Qstar[n]
			- 96. * (m-7) * (m-6) * (m-2) * R->Q[n] * R->Qstar[n]
			+ 24. * m * (m-7) * (m-6) * (m-5);
		hC8[n]->Fill(cent, c8.Re());
		hC8p[n]->Fill(cent, c8.Re()*c8.Re()/w8);

		TComplex c4w;
		if ( R == p ) {
			c4w = R->Qw1[n]*R->Qw1[n]*R->Qstarw1[n]*R->Qstarw1[n]
				+ R->Qw2[2*n] * R->Qstarw2[2*n]
				- 2. * ( R->Qw2[2*n] * R->Qstarw1[n] * R->Qstarw1[n] )
				+ 8. * ( R->Qw3[n] * R->Qstarw1[n] )
				- 4. * R->S2 * R->Qw1[n] * R->Qstarw1[n]
				- 6. * R->S4
				+ 2. * R->S2 * R->S2;
		} else {
			c4w = p->Q[n] * R->Qw1[n] * R->Qstarw1[n] * R->Qstarw1[n]
				- q->Qw1[2*n] * R->Qstarw1[n] * R->Qstarw1[n]
				- p->Q[n] * R->Qw1[n] * R->Qstarw2[2*n]
				- 2. * R->S2 * p->Q[n] * R->Qstarw1[n]
				- 2. * q->S1 * R->Qw1[n] * R->Qstarw1[n]
				+ 7. * q->Qw2[n] * R->Qstarw1[n]
				- R->Qw1[n] * q->Qstarw2[n]
				+ q->Qstarw1[2*n] * R->Qstarw2[2*n]
				+ 2. * p->Q[n] * R->Qstarw3[n]
				+ 2. * q->S1 * R->S2
				- 6. * q->S3;
		}
		hC4w[n]->Fill(cent, c4w.Re());

		double m2 = R->Q2[n] - R->M;
		double m4 = R->Q4[n] + pow(R->Q2[2*n], 2)
			- 2 * (R->Q[2*n] * R->Qstar[n] * R->Qstar[n]).Re()
			- 4 * (R->M -2) * R->Q2[n] + 2 * R->M * (R->M -3);
		hC2x4[n]->Fill(cent, c2.Re()*c4.Re());
		hC2x2r[n]->Fill(cent, m2*c2.Re());
		hC2x4r[n]->Fill(cent, m4*c2.Re());
		hC4x2r[n]->Fill(cent, m2*c4.Re());
		hC4x4r[n]->Fill(cent, m4*c4.Re());

		hC2x6[n]->Fill(cent, c2.Re()*c6.Re());
		hC4x6[n]->Fill(cent, c4.Re()*c6.Re());
		hC2x8[n]->Fill(cent, c2.Re()*c8.Re());
		hC4x8[n]->Fill(cent, c4.Re()*c8.Re());
		hC6x8[n]->Fill(cent, c6.Re()*c8.Re());
	}

	for ( int n = 1; n < N_vn; n++ ) {
		hReQ[n]->Fill(cent, p->Q[n].Re());
		hImQ[n]->Fill(cent, p->Q[n].Im());

		TComplex c = p->Q[n]*R->Q[n] - q->Q[2*n];
		hReQ2[n]->Fill(cent, c.Re());
		hImQ2[n]->Fill(cent, c.Im());

		c = p->Q[n] * ( R->Q2[n] - m ) - ( q->Q[2*n] * R->Qstar[n] + mq * R->Q[n] - 2. * q->Q[n] );
		hReQ31[n]->Fill(cent, c.Re());
		hImQ31[n]->Fill(cent, c.Im());
		
		c = p->Q[n] * R->Qstar[n] * R->Qstar[n] - p->Q[n] * R->Qstar[2*n] - ( 2 * mq * R->Qstar[n] - 2. * q->Qstar[n]);
		hReQ32[n]->Fill(cent, c.Re());
		hImQ32[n]->Fill(cent, c.Im());
	}
}


//define this as a plug-in
DEFINE_FWK_MODULE(QWCumulant);
