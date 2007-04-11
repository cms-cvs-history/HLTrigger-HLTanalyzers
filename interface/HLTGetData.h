#ifndef HLTGetData_h
#define HLTGetData_h

/** \class HLTGetData
 *
 *  
 *  This class is an EDAnalyzer implementing a "get data into RAM"
 *  functionality, to simulate online FF running/timimg.
 *
 *  $Date: 2007/03/28 12:58:54 $
 *  $Revision: 1.1 $
 *
 *  \author various
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class declaration
//

class HLTGetData : public edm::EDAnalyzer {

 public:
  explicit HLTGetData(const edm::ParameterSet&);
  ~HLTGetData();
  void analyze(const edm::Event&, const edm::EventSetup&);
  
 private:
  edm::InputTag EBdigiCollection_;
  edm::InputTag EEdigiCollection_;
  edm::InputTag ESdigiCollection_;      

};

#endif //HLTGetData_h
