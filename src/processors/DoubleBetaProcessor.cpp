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
      /* const int DD_BETAWALK1 = 11;//!< Walk Correction, QDC single side vs TDIFF (Bar 1)
       const int DD_BETAWALK2 = 12;//!< Walk Correction, QDC single side vs TDIFF (Bar 2)
       const int DD_BETAWALK3 = 13;//!< Walk Correction, QDC single side vs TDIFF (Bar 3)
       const int DD_BETAWALK4 = 14;//!< Walk Correction, QDC single side vs TDIFF (Bar 4)
      */
        const int DD_PP5 = 15;//!< ID to plost the phase-phase for your favorite bar (1)
        const int DD_QDCTDIFF5 = 16;//!< QDC vs. TDiff for your favorite bar (1), energy signals
       }
}

using namespace std;
using namespace dammIds::doublebeta;

DoubleBetaProcessor::DoubleBetaProcessor():
    EventProcessor(dammIds::doublebeta::OFFSET, dammIds::doublebeta::RANGE,
                   "Double Beta") {
    associatedTypes.insert("beta");
}

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
    DeclareHistogram2D(DD_PP4, SC, SC,"Phase vs. Phase - Bar 7 Only");
    DeclareHistogram2D(DD_QDCTDIFF4, SC, SC,"(Bar 7)TDiff vs. Coincident QDC");
    DeclareHistogram2D(DD_PP5, SC, SC,"Phase vs. Phase - Bar 1 Only");
    DeclareHistogram2D(DD_QDCTDIFF5, SC, SC,"(Bar 1)TDiff vs. Coincident QDC");
    

    /*
    DeclareHistogram2D(DD_BETAWALK1, SC, SE,"Beta Walk Correction, QDC(bar 1) vs. TDIFF");
    DeclareHistogram2D(DD_BETAWALK2, SC, SE,"Beta Walk Correction, QDC(bar 2) vs. TDIFF");
    DeclareHistogram2D(DD_BETAWALK3, SC, SE,"Beta Walk Correction, QDC(bar 3) vs. TDIFF");
    DeclareHistogram2D(DD_BETAWALK4, SC, SE,"Beta Walk Correction, QDC(bar 4) vs. TDIFF");
    */
}

bool DoubleBetaProcessor::PreProcess(RawEvent &event) {
    if (!EventProcessor::PreProcess(event))
        return(false);
    
    static const vector<ChanEvent*> & events =
        event.GetSummary("beta:double")->GetList();

    BarBuilder builder(events);
    builder.BuildBars();

    map<unsigned int, pair<double,double> > lrtbars = builder.GetLrtBarMap();
    BarMap betas = builder.GetBarMap();

    double resolution = 2;
    double offset = 1500;

    for(map<unsigned int, pair<double,double> >::iterator it = lrtbars.begin();
	it != lrtbars.end(); it++) {
	stringstream place;
	place << "DoubleBeta" << (*it).first;
	EventData data((*it).second.first, (*it).second.second, (*it).first);
	TreeCorrelator::get()->place(place.str())->activate(data);
    }
    
    for(BarMap::const_iterator it = betas.begin(); it != betas.end(); it++) {
        unsigned int barNum = (*it).first.first;
        plot(DD_QDC, (*it).second.GetLeftSide().GetTraceQdc(), barNum * 2);
        plot(DD_QDC, (*it).second.GetRightSide().GetTraceQdc(), barNum * 2 + 1);
        plot(DD_TDIFF, (*it).second.GetTimeDifference()*resolution + offset, barNum);
        if(barNum == 4) {
            plot(DD_PP1, (*it).second.GetLeftSide().GetPhase()*resolution,
                        (*it).second.GetRightSide().GetPhase()*resolution);
            plot(DD_QDCTDIFF1, (*it).second.GetTimeDifference()*resolution+offset,
             (*it).second.GetLeftSide().GetTraceQdc());
	    //	    if((*it).second.GetLeftSide().GetTraceQdc() > 3000.0){
	    //plot(DD_BETAWALK1,  (*it).second.GetTimeDifference()*resolution+offset,
	    //	   (*it).second.GetRightSide().GetTraceQdc());
	    }
        
	
	if(barNum == 5) {
	     plot(DD_PP2, (*it).second.GetLeftSide().GetPhase()*resolution,
                        (*it).second.GetRightSide().GetPhase()*resolution);
             plot(DD_QDCTDIFF2, (*it).second.GetTimeDifference()*resolution+offset,
             (*it).second.GetLeftSide().GetTraceQdc());
	     //if((*it).second.GetLeftSide().GetTraceQdc() > 3000.0){
	     //plot(DD_BETAWALK2,  (*it).second.GetTimeDifference()*resolution+offset,
	     //	   (*it).second.GetRightSide().GetTraceQdc());
	    }
        
	if(barNum == 6) {
            plot(DD_PP3, (*it).second.GetLeftSide().GetPhase()*resolution,
                        (*it).second.GetRightSide().GetPhase()*resolution);
            plot(DD_QDCTDIFF3, (*it).second.GetTimeDifference()*resolution+offset,
             (*it).second.GetLeftSide().GetTraceQdc());
	    //if((*it).second.GetLeftSide().GetTraceQdc() > 3000.0){
	    //plot(DD_BETAWALK3,  (*it).second.GetTimeDifference()*resolution+offset,
	    //	   (*it).second.GetRightSide().GetTraceQdc());
	    }
        
	
	if(barNum == 7) {
	     plot(DD_PP4, (*it).second.GetLeftSide().GetPhase()*resolution,
                        (*it).second.GetRightSide().GetPhase()*resolution);
             plot(DD_QDCTDIFF4, (*it).second.GetTimeDifference()*resolution+offset,
             (*it).second.GetLeftSide().GetTraceQdc());
	     //if((*it).second.GetLeftSide().GetTraceQdc() > 3000.0){
	     //plot(DD_BETAWALK4,  (*it).second.GetTimeDifference()*resolution+offset,
	     //	   (*it).second.GetRightSide().GetTraceQdc());
	    }
        

        if(barNum == 1) {
            plot(DD_PP5, (*it).second.GetLeftSide().GetPhase()*resolution,
                        (*it).second.GetRightSide().GetPhase()*resolution);
            plot(DD_QDCTDIFF5, (*it).second.GetTimeDifference()*resolution+offset,
             (*it).second.GetLeftSide().GetTraceQdc());
	}
    }
    return(true);
}

bool DoubleBetaProcessor::Process(RawEvent &event) {
    if (!EventProcessor::Process(event))
        return(false);
    return(true);
}
