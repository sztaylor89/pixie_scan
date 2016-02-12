/** \file ChanEvent.hpp
 * \brief A Class to define what a channel event is
 */
#ifndef __CHANEVENT_HPP
#define __CHANEVENT_HPP

#include <vector>

#include "ChannelEvent.hpp"

#include "DetectorLibrary.hpp"
#include "pixie16app_defs.h"
#include "Identifier.hpp"
#include "Globals.hpp"
#include "Trace.hpp"

/*! \brief A channel event
 *
 * All data is grouped together into channels.  For each pixie16 channel that
 * fires the energy, time (both trigger time and event time), and trace (if
 * applicable) are obtained.  Additional information includes the channels
 * identifier, calibrated energies, trace analysis information.
 * Note that this currently stores raw values internally through pixie word types
 *   but returns data values through native C types. This is potentially non-portable.
 */
class ChanEvent{
public:
    /** Default constructor */
    ChanEvent(){ }
    
    ChanEvent(ChannelEvent *event_){ 
    	event = event_; 
    	trace = Trace(event->trace); // Copy the trace from the ChannelEvent (messy, but needed for the Trace class).
    }
    
    ~ChanEvent(){ 
    	if(event != NULL){ delete event; } 
    }

    /** Set the raw energy in case we do not want to extract it from the trace
     * ourselves
     * \param [in] a : the energy to set */
    void SetEnergy(double a) {event->energy = a;}
    /** Set the calibrated energy
     * \param [in] a : the calibrated energy */
    void SetCalEnergy(double a) {event->calEnergy = a;}
    /** Set the time
     * \param [in] a : the time to set */
    void SetTime(double a) {event->time = a;}
    /** Set the Walk corrected time
     * \param [in] a : the walk corrected time */
    void SetCorrectedTime(double a) {event->correctedTime = a;}
    /** Set the high resolution time (Filter time + phase )
     * \param [in] a : the high resolution time */
    void SetHighResTime(double a) {event->hires_time =a;}

    bool GetCfdSourceBit() const {
		return(event->cfdTrigSource);
    }
    bool CfdForceTrig() const {
		return(event->cfdForceTrig);
    }
    double GetEnergy() const      {
        return event->energy;   /**< \return the raw energy */
    }
    double GetCalEnergy() const   {
        return event->calEnergy;   /**< \return the calibrated energy */
    }
    double GetTime() const        {
        return event->time;   /**< \return the raw time */
    }
    double GetCorrectedTime() const {
		return event->correctedTime;
    }
    double GetHighResTime() const {
        return event->hires_time;   /**< \return the high-resolution time */
    }
    double GetEventTime() const   {
        return event->eventTime;   /**< \return the event time */
    }
    const Trace& GetTrace() const {
        return trace;   /**< \return a reference to the trace */
    }
    Trace& GetTrace() {
        return trace;   /** \return a reference which can alter the trace */
    }
    unsigned long GetTrigTime() const {
        return event->trigTime;   /**< \return the channel trigger time */
    }
    unsigned long GetEventTimeLo() const {
        return event->eventTimeLo;   /**< \return the lower 32 bits of event time */
    }
    unsigned long GetEventTimeHi() const {
        return event->eventTimeHi;   /**< \return the upper 32 bits of event time */
    }
    bool IsPileup() const {
        return event->pileupBit;   //!< \return true if channel is pileup
    }
    bool IsSaturated() const { /**< \return whether the trace is saturated */
        return event->saturatedBit;
    }

	/// Save the pointer to the channel event associated with this ChanEvent.
	void SetChannelEvent(ChannelEvent *event_){ event = event_; }

    //! \return The identifier in the map for the channel event
    const Identifier& GetChanID() const;
    /** \return the channel id defined as pixie module # * 16 + channel number */
    int GetID() const;
    /** \return The Onboard QDC value at i
     * \param [in] i : the QDC number to obtain, possible values [0,7] */
    unsigned long GetQdcValue(int i) const;

    /** Channel event zeroing
     *
     * All numerical values are set to -1, and the trace,
     * and traceinfo vectors are cleared and the channel
     * identifier is zeroed using its identifier::zeroid method. */
    void ZeroVar();
private:
	ChannelEvent *event; /// The raw channel event.

    Trace trace; /**< Channel trace if present */

    bool virtualChannel; /**< Flagged if generated virtually in Pixie DSP */
    bool pileupBit;      /**< Pile-up flag from Pixie */
    bool saturatedBit;   /**< Saturation flag from Pixie */
    bool cfdForceTrig;   //!< CFD was forced to trigger
    bool cfdTrigSource;  //!< The ADC that the CFD/FPGA synched with
    
    void ZeroNums(void); /**< Zero members which do not have constructors associated with them */

    /** Make the front end responsible for reading the data able to set the
     * channel data directly from ReadBuffDataA - REVISION A */
    friend int ReadBuffDataA(pixie::word_t *, unsigned long *, std::vector<ChanEvent *> &);
    /** Make the front end responsible for reading the data able to set the
     * channel data directly from ReadBuffDataA - REVISION D */
    friend int ReadBuffDataD(pixie::word_t *, unsigned long *, std::vector<ChanEvent *> &);
    /** Make the front end responsible for reading the data able to set the
     * channel data directly from ReadBuffDataA - REVISION F */
    friend int ReadBuffDataF(pixie::word_t *, unsigned long *, std::vector<ChanEvent *> &);
};

/** Sort by increasing corrected time
 * \param [in] a : the left hand side for comparison
 * \param [in] b : the right hand side for comparison
 * \return True if LHS is less the RHS */
bool CompareCorrectedTime(const ChanEvent *a, const ChanEvent *b);

/** Sort by increasing raw time
 * \param [in] a : the left hand side for comparison
 * \param [in] b : the right hand side for comparison
 * \return True if LHS is less the RHS*/
bool CompareTime(const ChanEvent *a, const ChanEvent *b);

#endif
