/** \file Anl1471Processor.hpp
 * \brief A class to process data from the LeRIBSS 2012 campaign.
 *
 *\author S. V. Paulauskas, update by S. Z. Taylor
 *\date September 19, 2015
 */
#ifndef __ANL1471PROCESSOR_HPP_
#define __ANL1471PROCESSOR_HPP_
#include <fstream>

#include "EventProcessor.hpp"


#ifdef useroot
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH1D.h>
#endif



/// Class to process ANL experiment related events
class Anl1471Processor : public EventProcessor {
public:
    /** Default Constructor */
    Anl1471Processor();
    /** Default Destructor */
    ~Anl1471Processor();
    /** Declare the plots used in the analysis */
    virtual void DeclarePlots(void);

    /** Constructor taking a list of detector types as an argument
     * \param [in] typeList : the list of bar types that are in the analysis
     * \param [in] res : The resolution of the DAMM histograms
     * \param [in] offset : The offset of the DAMM histograms */
     Anl1471Processor(const std::vector<std::string> &typeList,
                  const double &res, const double &offset,
                  const double &numStarts);

    /** PreProcess does nothing since this is solely dependent on results
	from other Processors*/
     virtual bool PreProcess(RawEvent &event);

    /** Process the event
    * \param [in] event : the event to process
    * \return Returns true if the processing was successful */
    virtual bool Process(RawEvent &event);
private:
    std::string fileName_;
    std::vector<std::string> fileNames_;

#ifdef useroot
    TFile *rootfile_;
    TTree *roottree_;
    TH2D *qdctof_;
    TH1D *Vsize;
    TH1D *Bsize;
    static double qdc, tof;
    static double ge_, GamEn;
    static double vandle_, VID, SNRVL, SNRVR, QDCVL, QDCVR;
    static double beta_, BID, SNRBL, SNRBR, QDCBL, QDCBR;
#endif
};
#endif
