/** \file DoubleBetaProcessor.cpp
 *\brief A DoubleBeta processor class that can be used to analyze double
 * beta detectors.
 *\author S. V. Paulauskas
 *\date October 26, 2014
 */
#include "BarBuilder.hpp"
#include "DammPlotIds.hpp"
#include "DoubleBetaProcessor.hpp"
#include "Globals.hpp"
#include "RawEvent.hpp"
#include "TimingMapBuilder.hpp"
#include "TreeCorrelator.hpp"

namespace dammIds {
    namespace doublebeta {
        const int DD_SINGLESQDC = 0;//!< ID for the singles QDC
        const int DD_QDC  = 1; //!< ID for the Bar QDC of the double beta detector
        const int DD_TDIFF = 2;//!< ID to plot the Time Difference between ends
        const int DD_PP1 = 3;//!< ID to plot the phase-phase for your favorite bar (4)
       	const int DD_QDCTDIFF1 = 4;//!< QDC vs. TDiff for your favorite bar (4), timing signals
        const int DD_PP2 = 5;//!< ID to plost the phase-phase for your favorite bar (5)
        const int DD_QDCTDIFF2 = 6;//!< QDC vs. TDiff for your favorite bar (5), energy signals
	const int DD_PP3 = 7;//!< ID to plot the phase-phase for your favorite bar (6)
        const int DD_QDCTDIFF3 = 8;//!< QDC vs. TDiff for your favorite bar (6), timing signals
        const int DD_PP4 = 9;//!< ID to plost the phase-phase for your favorite bar (7)
        const int DD_QDCTDIFF4 = 10;//!< QDC vs. TDiff for your favorite bar (7), energy signals
	/* 
	const int DD_SNQDC1L = 11;//!< Signal to noise ratio vs qdc bar 1 left
	const int DD_SNQDC1R = 12;//!< Signal to noise ratio vs qdc bar 1 right
	const int DD_SNQDC2L = 13;//!< Signal to noise ratio vs qdc bar 2 left
	const int DD_SNQDC2R = 14;//!< Signal to noise ratio vs qdc bar 2 right
	const int DD_SNQDC3L = 15;//!< Signal to noise ratio vs qdc bar 3 left
	const int DD_SNQDC3R = 16;//!< Signal to noise ratio vs qdc bar 3 right
	const int DD_SNQDC4L = 17;//!< Signal to noise ratio vs qdc bar 4 left
	const int DD_SNQDC4R = 18;//!< Signal to noise ratio vs qdc bar 4 right
	*/

	const int DD_SNvsSumQDC1 = 11;//!< Signal to noise ratio vs qdc bar 1 summed
	const int DD_SNvsSumQDC2 = 12;//!< Signal to noise ratio vs qdc bar 2 summed
	const int DD_SNvsSumQDC3 = 13;//!< Signal to noise ratio vs qdc bar 3 summed
	const int DD_SNvsSumQDC4 = 14;//!< Signal to noise ratio vs qdc bar 4 summed
	


	/* 
	const int DD_BETAWALK1 = 19;//!< Walk Correction, QDC single side vs TDIFF (Bar 1)
	const int DD_BETAWALK2 = 20;//!< Walk Correction, QDC single side vs TDIFF (Bar 2)
	const int DD_BETAWALK3 = 21;//!< Walk Correction, QDC single side vs TDIFF (Bar 3)
	const int DD_BETAWALK4 = 22;//!< Walk Correction, QDC single side vs TDIFF (Bar 4)
	*/
	/*
	const int DD_QDCLR1 = 23;//!<Right vs Left QDC(Bar 1)
	const int DD_QDCLR2 = 24;//!<Right vs Left QDC(Bar 2)
	const int DD_QDCLR3 = 25;//!<Right vs Left QDC(Bar 3)
	const int DD_QDCLR4 = 26;//!<Right vs Left QDC(Bar 4)
	*/
}
}

using namespace std;
using namespace dammIds::doublebeta;

DoubleBetaProcessor::DoubleBetaProcessor():
    EventProcessor(OFFSET, RANGE, "DoubleBetaProcessor") {
    associatedTypes.insert("beta");
}
/*    
      #ifdef useroot
      tfile_ = new TFile("test.root","RECREATE");
      thist_ = new TH3D("interesting","",200,-100,100,2000,0,2000,120,0,60);
      #endif
      }
    
      DoubleBetaProcessor::~DoubleBetaProcessor() {
	        
#ifdef useroot
      tfile_->Write();
      tfile_->Close();
      #endif
      }
    
*/
void DoubleBetaProcessor::DeclarePlots(void) {
    DeclareHistogram2D(DD_SINGLESQDC, SD, S4, "Location vs. Singles QDC");
    DeclareHistogram2D(DD_QDC, SD, S4, "Location vs. Coincident QDC");
    DeclareHistogram2D(DD_TDIFF, SB, S3, "Location vs. Time Difference");
    DeclareHistogram2D(DD_PP1, SC, SC,"Phase vs. Phase - Bar 4 Only");
    DeclareHistogram2D(DD_QDCTDIFF1, SC, SC,"(Bar 4)TDiff vs. Coincident QDC");
    DeclareHistogram2D(DD_PP2, SC, SC,"Phase vs. Phase - Bar 5 Only");
    DeclareHistogram2D(DD_QDCTDIFF2, SC, SC,"(Bar 5)TDiff vs. Coincident QDC");
    DeclareHistogram2D(DD_PP3, SC, SC,"Phase vs. Phase - Bar 6 Only");
    DeclareHistogram2D(DD_QDCTDIFF3, SC, SC,"(Bar 6)TDiff vs. Coincident QDC");
    DeclareHistogram2D(DD_PP4, SC, SC,"Phase vs. Phase - Bar 7L Only");
    DeclareHistogram2D(DD_QDCTDIFF4, SC, SC,"(Bar 7L)TDiff vs. Coincident QDC");
    /*
    DeclareHistogram2D(DD_SNQDC1L, SC, S6,"(Bar 1L)SN ratio vs QDC)");
    DeclareHistogram2D(DD_SNQDC1R, SC, S6,"(Bar 1R)SN ratio vs QDC)");
    DeclareHistogram2D(DD_SNQDC2L, SC, S6,"(Bar 2L)SN ratio vs QDC)");
    DeclareHistogram2D(DD_SNQDC2R, SC, S6,"(Bar 2R)SN ratio vs QDC)");
    DeclareHistogram2D(DD_SNQDC3L, SC, S6,"(Bar 3L)SN ratio vs QDC)");
    DeclareHistogram2D(DD_SNQDC3R, SC, S6,"(Bar 3R)SN ratio vs QDC)");
    DeclareHistogram2D(DD_SNQDC4L, SC, S6,"(Bar 4L)SN ratio vs QDC)");
    DeclareHistogram2D(DD_SNQDC4R, SC, S6,"(Bar 4R)SN ratio vs QDC)");
    */
 
    DeclareHistogram2D(DD_SNvsSumQDC1, SC, S6,"Bar 1 SNR vs QDC Sum");
    DeclareHistogram2D(DD_SNvsSumQDC2, SC, S6,"Bar 2 SNR vs QDC Sum");
    DeclareHistogram2D(DD_SNvsSumQDC3, SC, S6,"Bar 3 SNR vs QDC Sum");
    DeclareHistogram2D(DD_SNvsSumQDC4, SC, S6,"Bar 4 SNR vs QDC Sum");


    /*
    DeclareHistogram2D(DD_BETAWALK1, SC, SE,"Beta Walk Correction, QDC(bar 1) vs. TDIFF");
    DeclareHistogram2D(DD_BETAWALK2, SC, SE,"Beta Walk Correction, QDC(bar 2) vs. TDIFF");
    DeclareHistogram2D(DD_BETAWALK3, SC, SE,"Beta Walk Correction, QDC(bar 3) vs. TDIFF");
    DeclareHistogram2D(DD_BETAWALK4, SC, SE,"Beta Walk Correction, QDC(bar 4) vs. TDIFF");
    */
    /*
    DeclareHistogram2D(DD_QDCLR1, SC, SC,"Bar 1 R vs. L QDC");
    DeclareHistogram2D(DD_QDCLR2, SC, SC,"Bar 2 R vs. L QDC");
    DeclareHistogram2D(DD_QDCLR3, SC, SC,"Bar 3 R vs. L QDC");
    DeclareHistogram2D(DD_QDCLR4, SC, SC,"Bar 4 R vs. L QDC");
    */
}

bool DoubleBetaProcessor::PreProcess(RawEvent &event) {
    if (!EventProcessor::PreProcess(event))
        return(false);
    lrtbars_.clear();
    bars_.clear();
    
    
    static const vector<ChanEvent*> & events =
        event.GetSummary("beta:double")->GetList();

    BarBuilder builder(events);
    builder.BuildBars();

    lrtbars_ = builder.GetLrtBarMap();
    bars_ = builder.GetBarMap();

    double resolution = 2;
    double offset = 1500;

    for(map<unsigned int, pair<double,double> >::iterator it = lrtbars_.begin();
	it != lrtbars_.end(); it++) {
	stringstream place;
	place << "DoubleBeta" << (*it).first;
	EventData data((*it).second.first, (*it).second.second, (*it).first);
	TreeCorrelator::get()->place(place.str())->activate(data);
    }
    
    for(BarMap::const_iterator it = bars_.begin(); it != bars_.end(); it++) {
        unsigned int barNum = (*it).first.first;
        plot(DD_QDC, (*it).second.GetLeftSide().GetTraceQdc(), barNum * 2);
        plot(DD_QDC, (*it).second.GetRightSide().GetTraceQdc(), barNum * 2 + 1);
        plot(DD_TDIFF, (*it).second.GetTimeDifference()*resolution + offset, barNum);

        if(barNum == 0) {

	    plot(DD_SNvsSumQDC1, ((*it).second.GetLeftSide().GetTraceQdc()+(*it).second.GetRightSide().GetTraceQdc())/2, 
		 ((*it).second.GetLeftSide().GetSignalToNoiseRatio()+(*it).second.GetRightSide().GetSignalToNoiseRatio())/2);


	    //plot (DD_QDCLR1,(*it).second.GetLeftSide().GetTraceQdc() , (*it).second.GetRightSide().GetTraceQdc());
	 
	  /*	  if(it->second.GetLeftSide().GetTraceQdc() > 4000){
	    for(unsigned int j = 0; j<it->second.GetRightSide().GetTrace()->size(); j++){
	      cout << j << " " << it->second.GetRightSide().GetTrace()->at(j) << endl;
	    }
	  cout << endl << endl;
	}
	  */
	    /*
	  plot(DD_SNQDC1L, (*it).second.GetLeftSide().GetTraceQdc(),
	       (*it).second.GetLeftSide().GetSignalToNoiseRatio());
	  plot(DD_SNQDC1R, (*it).second.GetRightSide().GetTraceQdc(),
	       (*it).second.GetLeftSide().GetSignalToNoiseRatio());
	    */
	  plot(DD_PP1, (*it).second.GetLeftSide().GetPhase()*resolution,
                        (*it).second.GetRightSide().GetPhase()*resolution);
            plot(DD_QDCTDIFF1, (*it).second.GetTimeDifference()*resolution+offset,
             (*it).second.GetQdc());
	    //	    if((*it).second.GetLeftSide().GetTraceQdc() > 3000.0){
	    //plot(DD_BETAWALK1,  (*it).second.GetTimeDifference()*resolution+offset,
	    //	   (*it).second.GetRightSide().GetTraceQdc());
	    }
        
	
	if(barNum == 1) {
	 
	    plot(DD_SNvsSumQDC2, ((*it).second.GetLeftSide().GetTraceQdc()+(*it).second.GetRightSide().GetTraceQdc())/2, 
		 ((*it).second.GetLeftSide().GetSignalToNoiseRatio()+(*it).second.GetRightSide().GetSignalToNoiseRatio())/2);

	    //	    plot (DD_QDCLR2,(*it).second.GetLeftSide().GetTraceQdc() , (*it).second.GetRightSide().GetTraceQdc());

	  //thist_->Fill((*it).second.GetTimeDifference(),
	  // 	       (*it).second.GetRightSide().GetTraceQdc(),
	  //	       (*it).second.GetRightSide().GetSignalToNoiseRatio());

	    /*
	  plot(DD_SNQDC2L, (*it).second.GetLeftSide().GetTraceQdc(),
	       (*it).second.GetLeftSide().GetSignalToNoiseRatio());
	  plot(DD_SNQDC2R, (*it).second.GetRightSide().GetTraceQdc(),
	       (*it).second.GetLeftSide().GetSignalToNoiseRatio());
	    */
	     plot(DD_PP2, (*it).second.GetLeftSide().GetPhase()*resolution,
                        (*it).second.GetRightSide().GetPhase()*resolution);
             plot(DD_QDCTDIFF2, (*it).second.GetTimeDifference()*resolution+offset,
             (*it).second.GetLeftSide().GetTraceQdc());
	     //if((*it).second.GetLeftSide().GetTraceQdc() > 3000.0){
	     //plot(DD_BETAWALK2,  (*it).second.GetTimeDifference()*resolution+offset,
	     //	   (*it).second.GetRightSide().GetTraceQdc());
	    }
        
	if(barNum == 2) {
	    
	    plot(DD_SNvsSumQDC3, ((*it).second.GetLeftSide().GetTraceQdc()+(*it).second.GetRightSide().GetTraceQdc())/2, 
		 ((*it).second.GetLeftSide().GetSignalToNoiseRatio()+(*it).second.GetRightSide().GetSignalToNoiseRatio())/2);

	    //	    plot (DD_QDCLR3,(*it).second.GetLeftSide().GetTraceQdc() , (*it).second.GetRightSide().GetTraceQdc());

	    /*
	    plot(DD_SNQDC3L, (*it).second.GetLeftSide().GetTraceQdc(),
		 (*it).second.GetLeftSide().GetSignalToNoiseRatio());
	    plot(DD_SNQDC3R, (*it).second.GetRightSide().GetTraceQdc(),
		 (*it).second.GetLeftSide().GetSignalToNoiseRatio());
	    */
	  plot(DD_PP3, (*it).second.GetLeftSide().GetPhase()*resolution,
                        (*it).second.GetRightSide().GetPhase()*resolution);
            plot(DD_QDCTDIFF3, (*it).second.GetTimeDifference()*resolution+offset,
             (*it).second.GetLeftSide().GetTraceQdc());
	    //if((*it).second.GetLeftSide().GetTraceQdc() > 3000.0){
	    //plot(DD_BETAWALK3,  (*it).second.GetTimeDifference()*resolution+offset,
	    //	   (*it).second.GetRightSide().GetTraceQdc());
	    }
        
	
	if(barNum == 3) {
	    
	    plot(DD_SNvsSumQDC4, ((*it).second.GetLeftSide().GetTraceQdc()+(*it).second.GetRightSide().GetTraceQdc())/2, 
		 ((*it).second.GetLeftSide().GetSignalToNoiseRatio()+(*it).second.GetRightSide().GetSignalToNoiseRatio())/2);	    


//	    plot (DD_QDCLR4,(*it).second.GetLeftSide().GetTraceQdc() , (*it).second.GetRightSide().GetTraceQdc());

	    /*
	    plot(DD_SNQDC4L, (*it).second.GetLeftSide().GetTraceQdc(),
	       (*it).second.GetLeftSide().GetSignalToNoiseRatio());
	    plot(DD_SNQDC4R, (*it).second.GetRightSide().GetTraceQdc(),
	       (*it).second.GetLeftSide().GetSignalToNoiseRatio());
	    */
	     plot(DD_PP4, (*it).second.GetLeftSide().GetPhase()*resolution,
                        (*it).second.GetRightSide().GetPhase()*resolution);
             plot(DD_QDCTDIFF4, (*it).second.GetTimeDifference()*resolution+offset,
             (*it).second.GetLeftSide().GetTraceQdc());
	     //  plot(DD_QDCTDIFF5, (*it).second.GetTimeDifference()*resolution+offset,
             //(*it).second.GetRightSide().GetTraceQdc());

	     //if((*it).second.GetLeftSide().GetTraceQdc() > 3000.0){
	     //plot(DD_BETAWALK4,  (*it).second.GetTimeDifference()*resolution+offset,
	     //	   (*it).second.GetRightSide().GetTraceQdc());
	}
    }
    return(true);
}

bool DoubleBetaProcessor::Process(RawEvent &event) {
    if (!EventProcessor::Process(event))
        return(false);
    EndProcess();
    return(true);
}
