#include <TMath.h>

//#define QW_PbPb_PT 1
//#define QW_CENT10 1

const double Pi = TMath::Pi();
const double Pi2 = 2*Pi;
const int N_vn = 6+1;
const int NOrder = 4*N_vn-1;
//#ifdef QW_CENT10
//const double centbins[]={0,10,20,30,40,50,60,70,80,90,100};
//#else
const double centbins[]={0,5,10,15,20,25,30,35,40,50,60,70,80,90,100};
//#endif
//const Int_t nCentBins = sizeof(centbins)/sizeof(double)-1;
const Int_t nCentBins = 20;
/*
const Int_t centmap[] = {
	0,	// 0	2.5
	0,	// 1	5.
	1,	// 2	7.5
	1,	// 3	10.
	2,	// 4	12.5
	2,	// 5	15.
	3,	// 6	17.5
	3,	// 7	20.
	4,	// 8	22.5
	4,	// 9	25.
	5,	// 10	27.5
	5,	// 11	30
	6,	// 12	32.5
	6,	// 13	35.
	7,	// 14	37.5
	7,	// 15	40.
	8,	// 16	42.5
	8,	// 17	45.
	8,	// 18	47.5
	8,	// 19	50.
	9,	// 20	52.5
	9,	// 21	55.
	9,	// 22	57.5
	9,	// 23	60.
	10,	// 24	62.5
	10,	// 25	65.
	10,	// 26	67.5
	10,	// 27	70
	11,	// 28	72.5
	11,	// 29	75.
	11,	// 30	77.5
	11,	// 31	80.
	12,	// 32	82.5
	12,	// 33	85.
	12,	// 34	87.5
	12,	// 35	90.
	13,	// 36	92.5
	13,	// 37	95.
	13,	// 38	97.5
	13	// 39	100
};
*/
const Int_t nBin = 40;

//const double ptbins_acc[19]={
//	0.3,  0.4,  0.5,  0.6,  0.8,  1.0,  1.2,  1.6,  2.0,
//	2.5,  3.0,  3.5,  4.0,  5.0,  6.0,  8.0, 10.0, 12.0, 1000.0};
//const int nPtAcc = 19;
	
//#ifdef QW_PbPb_PT
//const double ptbins[]={
//	0.3,  0.4,  0.5,  0.6,  0.8,  1.0,  1.2,  1.6,  2.0,
//	2.5,  3.0,  3.5,  4.0,  5.0,  6.0,  8.0, 10.0, 12.0, 16.0,
//	20.0, 30.0, 40.0, 100.0}; // PbPb pt binning nPtBins = 22

//#else

const double ptbins[] = {
        0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0,
        10.0, 12.0, 16.0, 20.0, 30.0, 40.0, 100.0}; // pPb pt binning nPtBins = 17;
//const double ptbins[]={
//	0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
//       	2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 
//	4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 6.0, 7.0, 8.0, 9.0,10.0,11.0,12.0,13.0,14.0,15.0,
//       16.0,17.0,18.0,19.0,20.0,21.0,22.0,23.0,24.0,25.0,26.0,27.0,28.0,29.0,30.0,40.0,50.0,60.0,60.0,80.0,
//       90.0,100.0}; // pPb pt binning nPtBins = 17;
////#endif

const double fakpt[] = {
	0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.5, 3, 4, 6, 8, 12
}; // 39 nbins=38


/*
const double ptbins[]={
	0.4,  1.0,  2.0, 3.0,  4.0,  5.0,  6.0,  8.0, 10.0, 12.0, 14.0,
	17.0, 20.0, 30.0, 40.0, 100.0};
*/
const Int_t nPtBins = sizeof(ptbins)/sizeof(double)-1;

const double etabins[] = {
	-2.4, -2.2, -2, -1.8, -1.6, -1.4, -1.2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4
};
const Int_t nEtaBins = sizeof(etabins)/sizeof(double)-1;

//const Int_t CentNoffCut[] = {100000, 350,290,260,240,220,185,150,120,100,80,70,65,60,55,50,40,0};
const Int_t CentNoffCut[] = {100000, 350, 320, 300, 260, 240, 220, 185, 150, 120, 100, 80, 60, 50, 40, 30, 20, 10, 0};
//const Int_t CentNoffCut[] = {100000,         300, 260,      220, 185, 150, 120, 110, 90, 35, 0};
const Int_t nCentNoff = sizeof(CentNoffCut)/sizeof(Int_t);
