#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string.h>

#include "HLTrigger/HLTanalyzers/interface/EventHeader.h"

EventHeader::EventHeader() {

  //set parameter defaults 
  _Debug=false;
}

EventHeader::~EventHeader() {

}

/*  Setup the analysis to put the branch-variables into the tree. */
void EventHeader::setup(TTree* HltTree) {

	fRun = -1;
	fEvent = -1;

  HltTree->Branch("run",&fRun,"run/I");
  HltTree->Branch("event",&fEvent,"event/I");

}

/* **Analyze the event** */
void EventHeader::analyze(edm::Event const& iEvent, TTree* HltTree) {
					
		fRun 		= iEvent.id().run();
		fEvent 	= iEvent.id().event();

    if (_Debug) {
		
			std::cout << "EventHeader -- run   = " << fRun << std::endl;
			std::cout << "EventHeader -- event = " << fEvent << std::endl;

		}

}
