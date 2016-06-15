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
#include "HighResTimingData.hpp"

double Anl1471Processor::qdc_;
double Anl1471Processor::tof;
//double Anl1471Processor::vandle_;
//double Anl1471Processor::beta_;
//double Anl1471Processor::ge_;

//double Anl1471Processor::VID;
//double Anl1471Processor::BID;
//double Anl1471Processor::GamEn;
//double Anl1471Processor::SNRBL;
//double Anl1471Processor::SNRBR;
//double Anl1471Processor::SNRVL;
//double Anl1471Processor::SNRVR;

static HighResTimingData::HrtRoot leftVandle;
static HighResTimingData::HrtRoot rightVandle;
static HighResTimingData::HrtRoot leftBeta;
static HighResTimingData::HrtRoot rightBeta;

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
    }
}//namespace dammIds

using namespace std;
using namespace dammIds::experiment;

void Anl1471Processor::DeclarePlots(void) {
    DeclareHistogram2D(DD_DEBUGGING0, SB, SD, "left-MaxValvsTDIFF");
    DeclareHistogram2D(DD_DEBUGGING1, SB, SD, "right-MaxValvsTDIFF");
    DeclareHistogram2D(DD_DEBUGGING2, SB, S6, "TDIFF-vandle");
    DeclareHistogram1D(DD_DEBUGGING3, S7, "Vandle Multiplicity");
    DeclareHistogram1D(DD_DEBUGGING4, S7, "Beta Multiplicity");
    DeclareHistogram2D(DD_DEBUGGING5, SC, SC, "DEBUG5");
    DeclareHistogram1D(DD_DEBUGGING6, SE, "DEBUG6");
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
    //roottree_->Branch("vandle",&vandle_,"VID:SNRVL:SNRVR:QDCVL:QDCVR");
    //roottree_->Branch("beta",&beta_,"BID:SNRBL:SNRBR:QDCBL:QDCBR");
    //roottree_->Branch("ge",&ge_,"GamEn");
    roottree_->Branch("leftV",&leftVandle,"qdc/D:time:snr:wtime:phase:abase:sbase:id/I");
    roottree_->Branch("rightV",&rightVandle,"qdc/D:time:snr:wtime:phase:abase:sbase:id/I");
    roottree_->Branch("leftB",&leftBeta,"qdc/D:time:snr:wtime:phase:abase:sbase:id/I");
    roottree_->Branch("rightB",&rightBeta,"qdc/D:time:snr:wtime:phase:abase:sbase:id/I");
    qdctof_ = new TH2D("qdctof","",1000,-100,900,16000,0,16000);
    Vsize = new TH1D("Vsize","",40,0,40);
    Bsize = new TH1D("Bsize","",40,0,40);
#endif
}

//Destructor closes/writes files and class
Anl1471Processor::~Anl1471Processor() {
#ifdef useroot
    rootfile_->Write();
    rootfile_->Close();
    delete(rootfile_);
    //delete(roottree_);
    //delete(qdctof_);
    //delete(Vsize);
    //delete(Bsize);
#endif
}



//where everything is done
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
    Vsize->Fill(vbars.size());
    Bsize->Fill(betas.size());
#endif
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
            // double corTof =
            //     ((VandleProcessor*)DetectorDriver::get()->
	    // 	 GetProcessor("VandleProcessor"))->
	    // 	CorrectTOF(tof, bar.GetFlightPath(), cal.GetZ0());
	    //	    bool isLowStart = start.GetQdc() < 300;

	    //stuff to fill root tree
	    bar.GetLeftSide().FillRootStructure(leftVandle);
	    bar.GetRightSide().FillRootStructure(rightVandle);

	    start.GetLeftSide().FillRootStructure(leftBeta);
	    start.GetRightSide().FillRootStructure(rightBeta);


	    if (barId.first == 2){
		plot(DD_DEBUGGING0,
		     bar.GetTimeDifference()*2+1000,
		     bar.GetLeftSide().GetMaximumValue());
		plot(DD_DEBUGGING1,
		     bar.GetTimeDifference()*2+1000,
		     bar.GetRightSide().GetMaximumValue());
		}
	    plot(DD_DEBUGGING2,
		 bar.GetTimeDifference()*2+1000, barId.first);

	    //VID=(*it).first.first;
	    //SNRVL=bar.GetLeftSide().GetSignalToNoiseRatio();
	    // SNRVR=bar.GetRightSide().GetSignalToNoiseRatio();
	    //QDCVL=bar.GetLeftSide().GetTraceQdc();
	    //QDCVR=bar.GetRightSide().GetTraceQdc();

#ifdef useroot
        qdctof_->Fill(tof,bar.GetQdc());
        qdc_ = bar.GetQdc();
        tof = tof;
        roottree_->Fill();
	bar.GetLeftSide().ZeroRootStructure(leftVandle);
	bar.GetRightSide().ZeroRootStructure(rightVandle);
	start.GetLeftSide().ZeroRootStructure(leftBeta);
	start.GetRightSide().ZeroRootStructure(rightBeta);
        qdc_ = tof = -9999;
	    //VID = BID = SNRVL = SNRVR = -9999;
	    //GamEn = SNRBL = SNRBR = vandle_ = beta_ = ge_ = -9999;
#endif

	    plot(DD_DEBUGGING1, tof*plotMult_+plotOffset_, bar.GetQdc());

        } // for(TimingMap::iterator itStart
    } //(BarMap::iterator itBar
    //End processing for VANDLE bars

    EndProcess();
    return(true);
}
