#define OHltMenu_cxx
#include "OHltMenu.h"


void OHltMenu::AddHlt(TString trig, TString l1Bit, int prescale, TString threshold, TString desc)
{

	hlts.push_back(trig);
	hltL1Bit[trig] 					= l1Bit;
	hltThreshold[trig] 			= threshold;
	hltDescription[trig] 		= desc;
	hltPrescale[trig]				= prescale;

}
