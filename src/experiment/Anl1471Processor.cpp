/** \file Anl1471Processor.cpp
 * \brief A class to process data from ANL1471 experiment using
 * VANDLE.
 *
 *\author S. V. Paulauskas, update by S. Z. Taylor
 *\date July 14, 2015
 */
#include <fstream>
#include <iostream>

#include <cmath>

#include "BarBuilder.hpp"
#include "DammPlotIds.hpp"
#include "DoubleBetaProcessor.hpp"
#include "DetectorDriver.hpp"
#include "GeProcessor.hpp"
#include "GetArguments.hpp"
#include "Globals.hpp"
#include "Anl1471Processor.hpp"
#include "RawEvent.hpp"
#include "TimingMapBuilder.hpp"
#include "VandleProcessor.hpp"

double Anl1471Processor::qdc_;
double Anl1471Processor::tof_;
double Anl1471Processor::vandle_;
double Anl1471Processor::beta_;
double Anl1471Processor::ge_;

double Anl1471Processor::VID_;
double Anl1471Processor::BID_;
double Anl1471Processor::GamEn_;
double Anl1471Processor::SNRBL_;
double Anl1471Processor::SNRBR_;
double Anl1471Processor::SNRVL_;
double Anl1471Processor::SNRVR_;


namespace dammIds {
    namespace experiment {
        const int DD_DEBUGGING0  = 0;
        const int DD_DEBUGGING1  = 1;
        const int DD_DEBUGGING2  = 2;
        const int DD_DEBUGGING3  = 3;
        const int DD_DEBUGGING4  = 4;
        const int DD_DEBUGGING5  = 5;
        const int DD_DEBUGGING6  = 6;
        const int DD_DEBUGGING7  = 7;
        const int DD_DEBUGGING8  = 8;
        const int DD_DEBUGGING9  = 9;
        const int DD_PROTONBETA2TDIFF_VS_BETA2EN = 13;
        const int D_ENERGY = 14;
        const int D_ENERGYBETA = 15;
        const int DD_PROTONGAMMATDIFF_VS_GAMMAEN = 16;
    }
}//namespace dammIds

using namespace std;
using namespace dammIds::experiment;

void Anl1471Processor::DeclarePlots(void) {
    DeclareHistogram2D(DD_DEBUGGING0, SC, SD, "QDC CTof- No Tape Move");
    DeclareHistogram2D(DD_DEBUGGING1, SC, SD, "QDC ToF Ungated");
    DeclareHistogram2D(DD_DEBUGGING2, SC, SC, "Cor ToF vs. Gamma E");
    DeclareHistogram1D(DD_DEBUGGING3, S7, "Vandle Multiplicity");
    DeclareHistogram1D(DD_DEBUGGING4, S7, "Beta Multiplicity");
    DeclareHistogram2D(DD_DEBUGGING5, SC, SC, "Mult2 Sym Plot Tof ");
    DeclareHistogram1D(DD_DEBUGGING6, SE, "LaBr3 RAW");
    DeclareHistogram2D(DD_PROTONBETA2TDIFF_VS_BETA2EN, SD, SA, "BetaProton Tdiff vs. Beta Energy");

    const int energyBins1  = SD;
    DeclareHistogram1D(D_ENERGY, energyBins1,
		       "Gamma singles ungated");
    DeclareHistogram1D(D_ENERGYBETA, energyBins1,
                       "Gamma singles beta gated");
    DeclareHistogram2D(DD_PROTONGAMMATDIFF_VS_GAMMAEN,
		       SD, SB, "GammaProton TDIFF vs. Gamma Energy");
}

Anl1471Processor::Anl1471Processor() : EventProcessor(OFFSET, RANGE, "Anl1471PRocessor") {
    associatedTypes.insert("vandle");
    associatedTypes.insert("beta");
    associatedTypes.insert("ge");

    char hisFileName[32];
    GetArgument(1, hisFileName, 32);
    string temp = hisFileName;
    temp = temp.substr(0, temp.find_first_of(" "));
#ifdef useroot
    stringstream rootname;
    rootname << temp << ".root";
    rootfile_ = new TFile(rootname.str().c_str(),"RECREATE");
    roottree_ = new TTree("ANL","");
    roottree_->Branch("vandle",&vandle_,"VID_:SNRVL_:SNRVR_:QDCVL_:QDCVR_");
    roottree_->Branch("beta",&beta_,"BID_:SNRBL_:SNRBR_:QDCBL_:QDCBR_");
    roottree_->Branch("ge",&ge_,"GamEn");
    qdctof_ = new TH2D("qdctof","",1000,-100,900,16000,0,16000);
    Vsize_ = new TH1D("Vsize","",40,0,40);
    Bsize_ = new TH1D("Bsize","",40,0,40);
#endif
}

//Destructor closes/writes files and class
Anl1471Processor::~Anl1471Processor() {
#ifdef useroot
    rootfile_->Write();
    rootfile_->Close();
    delete(rootfile_);
#endif
}


///We do nothing here since we're completely dependent on the resutls of others
bool Anl1471Processor::PreProcess(RawEvent &event){
    if (!EventProcessor::PreProcess(event))
        return(false);
    return(true);
}


//where evrything is done
bool Anl1471Processor::Process(RawEvent &event) {
    if (!EventProcessor::Process(event))
        return(false);
    double plotMult_ = 2;
    double plotOffset_ = 1000;

    BarMap vbars, betas;
    map<unsigned int, pair<double,double> > lrtBetas;
    vector<ChanEvent*> geEvts;
    vector<vector<AddBackEvent> > geAddback;

    if(event.GetSummary("vandle")->GetList().size() != 0)
        vbars = ((VandleProcessor*)DetectorDriver::get()->
            GetProcessor("VandleProcessor"))->GetBars();
    if(event.GetSummary("beta:double")->GetList().size() != 0) {
        betas = ((DoubleBetaProcessor*)DetectorDriver::get()->
            GetProcessor("DoubleBetaProcessor"))->GetBars();
        lrtBetas = ((DoubleBetaProcessor*)DetectorDriver::get()->
            GetProcessor("DoubleBetaProcessor"))->GetLowResBars();
    }
    if(event.GetSummary("ge")->GetList().size() != 0) {
        geEvts = ((GeProcessor*)DetectorDriver::get()->
            GetProcessor("GeProcessor"))->GetGeEvents();
        geAddback = ((GeProcessor*)DetectorDriver::get()->
            GetProcessor("GeProcessor"))->GetAddbackEvents();
    }

#ifdef useroot
    Vsize_->Fill(vbars.size());
    Bsize_->Fill(betas.size());
#endif

    //Obtain some useful logic statuses
    //    double lastProtonTime =
    //	TreeCorrelator::get()->place("logic_t1_0")->last().time;
    // bool isTapeMoving = TreeCorrelator::get()->place("TapeMove")->status();

    //int bananaNum = 2;
    //bool hasMultOne = vbars.size() == 1;
    //bool hasMultTwo = vbars.size() == 2;
    //bool isFirst = true;

    plot(DD_DEBUGGING3, vbars.size());
    plot(DD_DEBUGGING4, betas.size());
    //Begin processing for VANDLE bars
    for (BarMap::iterator it = vbars.begin(); it !=  vbars.end(); it++) {
        TimingDefs::TimingIdentifier barId = (*it).first;
        BarDetector bar = (*it).second;

        if(!bar.GetHasEvent() || bar.GetType() == "small")
            continue;

        //unsigned int barLoc = barId.first;
        TimingCalibration cal = bar.GetCalibration();

        for(BarMap::iterator itStart = betas.begin();
	    itStart != betas.end(); itStart++) {
	    unsigned int startLoc = (*itStart).first.first;
            BarDetector start = (*itStart).second;
            if(!start.GetHasEvent())
                continue;

            double tofOffset = cal.GetTofOffset(startLoc);
            double tof = bar.GetCorTimeAve() -
                start.GetCorTimeAve() + tofOffset;

            double corTof =
                ((VandleProcessor*)DetectorDriver::get()->
		 GetProcessor("VandleProcessor"))->
		CorrectTOF(tof, bar.GetFlightPath(), cal.GetZ0());

	    //	    bool notPrompt = corTof > 45.;
	    // bool inPeel = histo.BananaTest(bananaNum,
            //corTof*plotMult_+plotOffset_,
            //bar.GetQdc());
	    bool isLowStart = start.GetQdc() < 300;


#ifdef useroot
        qdctof_->Fill(tof,bar.GetQdc());
        qdc_ = bar.GetQdc();
        tof_ = tof;
        roottree_->Fill();
        qdc_ = tof_ = VID_ = BID_ = SNRVL_ = SNRVR_ = -9999;
	GamEn_ = SNRBL_ = SNRBR_ = vandle_ = beta_ = ge_ = -9999;
#endif

	    plot(DD_DEBUGGING1, tof*plotMult_+plotOffset_, bar.GetQdc());
	    //if(!isTapeMoving && !isLowStart)
	    //	plot(DD_DEBUGGING0, corTof*plotMult_+plotOffset_,bar.GetQdc());
	    //if(hasMultOne)
	    //	plot(DD_DEBUGGING4, corTof*plotMult_+plotOffset_, bar.GetQdc());

	    /*
	    ///Starting to look for 2n coincidences in VANDLE
	    BarMap::iterator itTemp = it;
	    itTemp++;

	    for (BarMap::iterator it2 = itTemp; it2 !=  vbars.end(); it2++) {
		TimingDefs::TimingIdentifier barId2 = (*it2).first;
		BarDetector bar2 = (*it2).second;

		if(!bar.GetHasEvent())
		    continue;

		unsigned int barLoc2 = barId2.first;

		bool isAdjacent = abs((int)barLoc2 - (int)barLoc) < 1;

		TimingCalibration cal2 = bar2.GetCalibration();

		double tofOffset2 = cal2.GetTofOffset(startLoc);
		double tof2 = bar2.GetCorTimeAve() -
		    start.GetCorTimeAve() + tofOffset2;

		double corTof2 =
		    ((VandleProcessor*)DetectorDriver::get()->
		     GetProcessor("VandleProcessor"))->
		    CorrectTOF(tof2, bar2.GetFlightPath(), cal2.GetZ0());

		bool inPeel2 = histo.BananaTest(bananaNum,
						corTof2*plotMult_+plotOffset_,
						bar2.GetQdc());

		if(hasMultTwo && inPeel && inPeel2 && !isAdjacent) {
		    plot(DD_DEBUGGING5, corTof*plotMult_+plotOffset_,
			 corTof2*plotMult_+plotOffset_);
		    plot(DD_DEBUGGING5, corTof2*plotMult_+plotOffset_,
			 corTof*plotMult_+plotOffset_);
		}
	    }
*/
        } // for(TimingMap::iterator itStart
    } //(BarMap::iterator itBar
    //End processing for VANDLE bars

    /*
    //------------------ Double Beta Processing --------------
    for(map<unsigned int, pair<double,double> >::iterator it = lrtBetas.begin();
    it != lrtBetas.end(); it++)
	plot(DD_PROTONBETA2TDIFF_VS_BETA2EN, it->second.second,
	(it->second.first - lastProtonTime) /
	(10e-3/Globals::get()->clockInSeconds()) );
    */
	     
    /*
    //----------------- GE Processing -------------------
    bool hasBeta = TreeCorrelator::get()->place("Beta")->status();
    double clockInSeconds = Globals::get()->clockInSeconds();
    //plot with 10 ms bins
    const double plotResolution = 10e-3 / clockInSeconds;
    
    for (vector<ChanEvent*>::iterator it1 = geEvts.begin();
	 it1 != geEvts.end(); ++it1) {
        ChanEvent *chan = *it1;
	
        double gEnergy = chan->GetCalEnergy();
        double gTime   = chan->GetCorrectedTime();
        if (gEnergy < 10.) //hard coded fix later.
            continue;
	
        plot(D_ENERGY, gEnergy);
	if(hasBeta)
	    plot(D_ENERGYBETA, gEnergy);
	plot(DD_PROTONGAMMATDIFF_VS_GAMMAEN, gEnergy ,
	     (gTime - lastProtonTime) / plotResolution) ;
    }*/

    EndProcess();
    return(true);
}
