/** \class HLTGetRaw
 *
 * See header file for documentation
 *
 *  $Date: 2008/02/19 13:58:31 $
 *  $Revision: 1.2.6.1 $
 *
 *  \author various
 *
 */

#include "HLTrigger/HLTanalyzers/interface/HLTGetRaw.h"

#include "DataFormats/Common/interface/Handle.h"

// system include files
#include <memory>
#include <vector>
#include <map>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DataFormats/Common/interface/DetSetVector.h"

#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/FEDRawData/interface/FEDRawData.h"
#include "DataFormats/FEDRawData/interface/FEDNumbering.h"

// using namespace edm;
// using namespace std;

//
// constructors and destructor
//
HLTGetRaw::HLTGetRaw(const edm::ParameterSet& ps)
{
  RawDataCollection_ = ps.getParameter<edm::InputTag>("RawDataCollection");
}

HLTGetRaw::~HLTGetRaw()
{ }

//
// member functions
//

// ------------ method called to produce the data  ------------
void
HLTGetRaw::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
//    using namespace edm;

    std::string errMsg("");
    edm::Handle<FEDRawDataCollection> RawDataHandle ; 
    //iEvent.getByLabel(RawDataCollection_, RawDataHandle );
    try {iEvent.getByLabel(RawDataCollection_, RawDataHandle);} catch (...) {errMsg=errMsg + " HLTGetRaw: No Raw Data Collection";}



    LogDebug("DigiInfo") << "Loaded Raw Data Collection: " << RawDataCollection_ ; 

    
}
