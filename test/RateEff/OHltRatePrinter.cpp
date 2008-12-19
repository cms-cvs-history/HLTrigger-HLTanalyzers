#include <iostream>
#include <iomanip>
#include <fstream>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TString.h>
#include "OHltRatePrinter.h"
#include "OHltTree.h"

void OHltRatePrinter::SetupAll(vector<float> tRate,vector<float> tRateErr,vector<float> tspureRate,
			  vector<float> tspureRateErr,vector<float> tpureRate,
			  vector<float> tpureRateErr,vector< vector<float> >tcoMa) {
  Rate = tRate;
  RateErr = tRateErr;
  spureRate = tspureRate;
  spureRateErr = tspureRateErr;
  pureRate = tpureRate;
  pureRateErr = tpureRateErr;
  coMa = tcoMa;
}

/* ********************************************** */
// Print out rate as ascii
/* ********************************************** */
void OHltRatePrinter::printRatesASCII(OHltConfig *cfg, OHltMenu *menu) {

  cout.setf(ios::floatfield,ios::fixed);
  cout<<setprecision(3);

  cout << "\n";
  cout << "Trigger Rates [Hz] : " << "\n";
  cout << "         Name                       Prescale           Indiv.          Pure   Cumulative\n";
  cout << "----------------------------------------------------------------------------------------------\n";

  float cumulRate = 0.;
  float cumulRateErr = 0.;
  for (unsigned int i=0;i<menu->GetTriggerSize();i++) {
    cumulRate += spureRate[i];
    cumulRateErr += pow(spureRateErr[i],2.);
    cout<<setw(34)<<menu->GetTriggerName(i)<<" ("
	<<setw(8)<<menu->GetPrescale(i)<<")  "
	<<setw(8)<<Rate[i]<<" +- "
	<<setw(7)<<RateErr[i]<<"  "
	<<setw(8)<<spureRate[i]<<"  "
	<<setw(8)<<cumulRate
	<<endl;
  }

  cumulRateErr = sqrt(cumulRateErr);
  cout << "\n";
  cout << setw(60) << "TOTAL RATE : " << setw(5) << cumulRate << " +- " << cumulRateErr << " Hz" << "\n";
  cout << "----------------------------------------------------------------------------------------------\n";
  
}


/* ********************************************** */
// Fill histos
/* ********************************************** */
void OHltRatePrinter::writeHistos(OHltConfig *cfg, OHltMenu *menu) {
  TString tableFileName = GetFileName(cfg,menu);

  TFile *fr = new TFile(tableFileName+TString(".root"),"recreate");
  fr->cd();
  
  int nTrig = (int)menu->GetTriggerSize();
  TH1F *individual = new TH1F("individual","individual",nTrig,1,nTrig+1);
  TH1F *cumulative = new TH1F("cumulative","cumulative",nTrig,1,nTrig+1);
  TH1F *throughput = new TH1F("throughput","throughput",nTrig,1,nTrig+1);
  TH1F *eventsize = new TH1F("eventsize","eventsize",nTrig,1,nTrig+1);
  TH2F *overlap = new TH2F("overlap","overlap",nTrig,1,nTrig+1,nTrig,1,nTrig+1);


  float cumulRate = 0.;
  float cumulRateErr = 0.;
  float cuThru = 0.;
  float cuThruErr = 0.;
  for (unsigned int i=0;i<menu->GetTriggerSize();i++) {
    cumulRate += spureRate[i];
    cumulRateErr += pow(spureRateErr[i],2.);
    cuThru += spureRate[i] * menu->GetEventsize(i);
    cuThruErr += pow(spureRate[i]*menu->GetEventsize(i),2.);

    individual->SetBinContent(i+1,Rate[i]);
    individual->GetXaxis()->SetBinLabel(i+1,menu->GetTriggerName(i));
    cumulative->SetBinContent(i+1,cumulRate);
    cumulative->GetXaxis()->SetBinLabel(i+1,menu->GetTriggerName(i));

    throughput->SetBinContent(i+1,cuThru);
    throughput->GetXaxis()->SetBinLabel(i+1,menu->GetTriggerName(i)); 
    eventsize->SetBinContent(i+1,menu->GetEventsize(i)); 
    eventsize->GetXaxis()->SetBinLabel(i+1,menu->GetTriggerName(i));      
  }

  individual->SetStats(0); individual->SetYTitle("Rate (Hz)");
  individual->SetTitle("Individual trigger rate");
  cumulative->SetStats(0); cumulative->SetYTitle("Rate (Hz)");
  cumulative->SetTitle("Cumulative trigger rate"); 
  overlap->SetStats(0); overlap->SetTitle("Overlap");
  individual->Write();
  cumulative->Write();
  eventsize->Write();
  throughput->Write();
  overlap->Write();
  fr->Close();
}



/* ********************************************** */
// Generate basic file name
/* ********************************************** */
TString OHltRatePrinter::GetFileName(OHltConfig *cfg, OHltMenu *menu) {
  char sLumi[10],sEnergy[10]; 
  sprintf(sEnergy,"%1.0f",cfg->cmsEnergy); 
  sprintf(sLumi,"%1.1e",cfg->iLumi); 

  TString menuTag;
  if (menu->IsL1Menu()) menuTag = "l1menu_";
  else menuTag = "hltmenu_";
    
  TString tableFileName = menuTag  + sEnergy + TString("TeV_") + sLumi
    + TString("_") + cfg->alcaCondition + TString("_") + cfg->versionTag; 
  tableFileName.ReplaceAll("+","");

  return tableFileName;
}

/* ********************************************** */
// Print out corelation matrix
/* ********************************************** */
void OHltRatePrinter::printCorrelationASCII() {

  for (unsigned int i=0;i<Rate.size();i++) {
    for (unsigned int j=0;j<Rate.size();j++) {
      cout<<"("<<i<<","<<j<<") = "<<coMa[i][j]<<endl;
    }
  }
}

/* ********************************************** */
// Print out rates as tex
/* ********************************************** */
void OHltRatePrinter::printRatesTex(OHltConfig *cfg, OHltMenu *menu) {
  if (menu->IsL1Menu())
    printL1RatesTex(cfg,menu);
  else
    printHltRatesTex(cfg,menu);
    
}
/* ********************************************** */
// Print out L1 rates as tex
/* ********************************************** */
void OHltRatePrinter::printL1RatesTex(OHltConfig *cfg, OHltMenu *menu) {
  TString tableFileName = GetFileName(cfg,menu);
  
  char sLumi[10],sEnergy[10]; 
  sprintf(sEnergy,"%1.0f",cfg->cmsEnergy); 
  sprintf(sLumi,"%1.1e",cfg->iLumi); 

  TString texFile = tableFileName + TString(".tex");
  ofstream outFile(texFile.Data());
  if (!outFile){cout<<"Error opening output file"<< endl;}

  outFile <<setprecision(2);
  outFile.setf(ios::floatfield,ios::fixed);
  outFile << "\\documentclass[amsmath,amssymb]{revtex4}" << endl;
  outFile << "\\usepackage{longtable}" << endl;
  outFile << "\\usepackage{color}" << endl;
  outFile << "\\usepackage{lscape}" << endl;
  outFile << "\\begin{document}" << endl;
  outFile << "\\begin{landscape}" << endl;
  outFile << "\\newcommand{\\met}{\\ensuremath{E\\kern-0.6em\\lower-.1ex\\hbox{\\/}\\_T}}" << endl;

  outFile << "\\begin{footnotesize}" << endl;
  outFile << "\\begin{longtable}{|l|c|c|r|}" << endl;
  outFile << "\\caption[Cuts]{L1 bandwidth is 17 kHz. } \\label{CUTS} \\\\ " << endl;
  
  outFile << "\\hline \\multicolumn{4}{|c|}{\\bf \\boldmath L1 for L = "<< sLumi  << "}\\\\  \\hline" << endl;
  outFile << "{\\bf Path Name} & " << endl;
  outFile << "\\begin{tabular}{c} {\\bf L1} \\\\ {\\bf Prescale} \\end{tabular} & " << endl;
  outFile << "\\begin{tabular}{c} {\\bf L1 Rate} \\\\ {\\bf $[$Hz$]$} \\end{tabular} & " << endl;
  outFile << "\\begin{tabular}{c} {\\bf Total Rate} \\\\ {\\bf $[$Hz$]$} \\end{tabular} \\\\ \\hline" << endl;
  outFile << "\\endfirsthead " << endl;
  
  outFile << "\\multicolumn{4}{r}{\\bf \\bfseries --continued from previous page (L = " << sLumi << ")"  << "}\\\\ \\hline " << endl;
  outFile << "{\\bf Path Name} & " << endl;
  outFile << "\\begin{tabular}{c} {\\bf L1} \\\\ {\\bf Prescale} \\end{tabular} & " << endl;
  outFile << "\\begin{tabular}{c} {\\bf L1 Rate} \\\\ {\\bf $[$Hz$]$} \\end{tabular} & " << endl;
  outFile << "\\begin{tabular}{c} {\\bf Total Rate} \\\\ {\\bf $[$Hz$]$} \\end{tabular} \\\\ \\hline" << endl;
  outFile << "\\endhead " << endl;
  
  outFile << "\\hline \\multicolumn{4}{|r|}{{Continued on next page}} \\\\ \\hline " << endl;
  outFile << "\\endfoot " << endl;
  
  outFile << "\\hline " << endl;
  outFile << "\\endlastfoot " << endl;
  
  float cumulRate = 0.;
  float cumulRateErr = 0.;
  for (unsigned int i=0;i<menu->GetTriggerSize();i++) {
    cumulRate += spureRate[i];
    cumulRateErr += pow(spureRateErr[i],2.);
    
    TString tempTrigName = menu->GetTriggerName(i);
    tempTrigName.ReplaceAll("_","\\_");

    outFile << "\\color{blue}"  << tempTrigName
	    << " & " <<  menu->GetPrescale(i) 
	    << " & " << Rate[i] << " {$\\pm$ " << RateErr[i]
	    << "} & " << cumulRate << "\\\\" << endl;
  }
  
  cumulRateErr = sqrt(cumulRateErr);
  outFile << "\\hline \\multicolumn{2}{|l|}{\\bf \\boldmath Total L1 rate (Hz)} & \\multicolumn{2}{|r|} {\\bf "<<  cumulRate << " $\\pm$ " << cumulRateErr << "} \\\\  \\hline" << endl;  
  outFile << "\\end{longtable}" << endl;
    outFile << "\\end{footnotesize}" << endl;
    outFile << "\\clearpage" << endl;
    outFile << "\\end{landscape}" << endl;
    outFile << "\\end{document}";
    outFile.close();
}


/* ********************************************** */
// Print out Hlt rates as tex
/* ********************************************** */
void OHltRatePrinter::printHltRatesTex(OHltConfig *cfg, OHltMenu *menu) {
  TString tableFileName = GetFileName(cfg,menu);

  char sLumi[10],sEnergy[10]; 
  sprintf(sEnergy,"%1.0f",cfg->cmsEnergy); 
  sprintf(sLumi,"%1.1e",cfg->iLumi); 

  TString texFile = tableFileName + TString(".tex");
  ofstream outFile(texFile.Data());
  if (!outFile){cout<<"Error opening output file"<< endl;}

  outFile <<setprecision(2);
  outFile.setf(ios::floatfield,ios::fixed);
  outFile << "\\documentclass[amsmath,amssymb]{revtex4}" << endl;
  outFile << "\\usepackage{longtable}" << endl;
  outFile << "\\usepackage{color}" << endl;
  outFile << "\\usepackage{lscape}" << endl;
  outFile << "\\begin{document}" << endl;
  outFile << "\\begin{landscape}" << endl;
  outFile << "\\newcommand{\\met}{\\ensuremath{E\\kern-0.6em\\lower-.1ex\\hbox{\\/}\\_T}}" << endl;

  outFile << "\\begin{footnotesize}" << endl;
  outFile << "\\begin{longtable}{|l|l|c|c|c|r|c|r|}" << endl;
  outFile << "\\caption[Cuts]{Available HLT bandwith is 150 Hz = ((1 GB/s / 3) - 100 MB/s for AlCa triggers) / 1.5 MB/event. } \\label{CUTS} \\\\ " << endl;
  
  outFile << "\\hline \\multicolumn{8}{|c|}{\\bf \\boldmath HLT for L = "<< sLumi  << "}\\\\  \\hline" << endl;
  outFile << "{\\bf Path Name} & " << endl;
  outFile << "\\begin{tabular}{c} {\\bf L1} \\\\ {\\bf Condition} \\end{tabular} & " << endl;
  outFile << "\\begin{tabular}{c} {\\bf L1} \\\\ {\\bf Prescale} \\end{tabular} & " << endl;
  outFile << "\\begin{tabular}{c} {\\bf HLT} \\\\ {\\bf Prescale} \\end{tabular} & " << endl;
  outFile << "\\begin{tabular}{c} {\\bf HLT Rate} \\\\ {\\bf $[$Hz$]$} \\end{tabular} & " << endl;
  outFile << "\\begin{tabular}{c} {\\bf Total Rate} \\\\ {\\bf $[$Hz$]$} \\end{tabular} &" << endl;
  outFile << "\\begin{tabular}{c} {\\bf Avg. Size} \\\\ {\\bf $[$MB$]$} \\end{tabular} & " << endl;  
  outFile << "\\begin{tabular}{c} {\\bf Total} \\\\ {\\bf Throughput} \\\\ {\\bf $[$MB/s$]$} \\end{tabular} \\\\ \\hline" << endl;  
  outFile << "\\endfirsthead " << endl;
  
  outFile << "\\multicolumn{8}{r}{\\bf \\bfseries --continued from previous page (L = " << sLumi << ")"  << "}\\\\ \\hline " << endl;
  outFile << "{\\bf Path Name} & " << endl;
  outFile << "\\begin{tabular}{c} {\\bf L1} \\\\ {\\bf Condition} \\end{tabular} & " << endl;
  outFile << "\\begin{tabular}{c} {\\bf L1} \\\\ {\\bf Prescale} \\end{tabular} & " << endl;
  outFile << "\\begin{tabular}{c} {\\bf HLT} \\\\ {\\bf Prescale} \\end{tabular} & " << endl;
  outFile << "\\begin{tabular}{c} {\\bf HLT Rate} \\\\ {\\bf $[$Hz$]$} \\end{tabular} & " << endl;
  outFile << "\\begin{tabular}{c} {\\bf Total Rate} \\\\ {\\bf $[$Hz$]$} \\end{tabular} &" << endl;
  outFile << "\\begin{tabular}{c} {\\bf Avg. Size} \\\\ {\\bf $[$MB$]$} \\end{tabular} & " << endl;  
  outFile << "\\begin{tabular}{c} {\\bf Total} \\\\ {\\bf Throughput} \\\\ {\\bf $[$MB/s$]$} \\end{tabular} \\\\ \\hline" << endl;  
  outFile << "\\endhead " << endl;
  
  outFile << "\\hline \\multicolumn{8}{|r|}{{Continued on next page}} \\\\ \\hline " << endl;
  outFile << "\\endfoot " << endl;
  
  outFile << "\\hline " << endl;
  outFile << "\\endlastfoot " << endl;
  
  float cumulRate = 0.;
  float cumulRateErr = 0.;
  float cuThru = 0.;
  float cuThruErr = 0.;
  float physCutThru = 0.;
  float physCutThruErr = 0.;
  float cuPhysRate = 0.;
  float cuPhysRateErr = 0.;
  vector<TString> footTrigSeedPrescales;
  vector<TString> footTrigSeeds;
  vector<TString> footTrigNames;
  for (unsigned int i=0;i<menu->GetTriggerSize();i++) {
    cumulRate += spureRate[i];
    cumulRateErr += pow(spureRateErr[i],2.);
    cuThru += spureRate[i] * menu->GetEventsize(i);
    cuThruErr += pow(spureRateErr[i]*menu->GetEventsize(i),2.);

    if (!(menu->GetTriggerName(i).Contains("AlCa"))) {
      cuPhysRate += spureRate[i];
      cuPhysRateErr += pow(spureRateErr[i],2.);
      physCutThru += spureRate[i]*menu->GetEventsize(i);
      physCutThruErr += pow(spureRateErr[i]*menu->GetEventsize(i),2.);
    }
    
    TString tempTrigName = menu->GetTriggerName(i);
    tempTrigName.ReplaceAll("_","\\_");

    TString tempTrigSeedPrescales;
    TString tempTrigSeeds;
    std::map<TString, std::vector<TString> >
      mapL1seeds = menu->GetL1SeedsOfHLTPathMap(); // mapping to all seeds

    vector<TString> vtmp;
    vector<int> itmp;
    typedef map< TString, vector<TString> >  mymap;
    for(mymap::const_iterator it = mapL1seeds.begin();it != mapL1seeds.end(); ++it) {
      if (it->first.CompareTo(menu->GetTriggerName(i)) == 0) {
	vtmp = it->second;
	//cout<<it->first<<endl;
	for (unsigned int j=0;j<it->second.size();j++) {
	  itmp.push_back(menu->GetL1Prescale((it->second)[j]));
	  //cout<<"\t"<<(it->second)[j]<<endl;
	}
      }
    }
    // Faster, but crashes???:
    //vector<TString> vtmp = mapL1seeds.find(TString(menu->GetTriggerName(i)))->second;
    if (vtmp.size()>2) {
      for (unsigned int j=0;j<vtmp.size();j++) {
	tempTrigSeeds = tempTrigSeeds + vtmp[j];
	tempTrigSeedPrescales += itmp[j];
	if (j<(vtmp.size()-1)) {
	  tempTrigSeeds = tempTrigSeeds + ", ";
	  tempTrigSeedPrescales = tempTrigSeedPrescales + ", ";
	}
      }

      tempTrigSeeds.ReplaceAll("_","\\_");
      tempTrigSeedPrescales.ReplaceAll("_","\\_");
      footTrigSeedPrescales.push_back(tempTrigSeedPrescales);
      footTrigSeeds.push_back(tempTrigSeeds);
      TString tmpstr = menu->GetTriggerName(i);
      tmpstr.ReplaceAll("_","\\_");
      footTrigNames.push_back(tmpstr);

      
      tempTrigSeeds = "List Too Long";
      tempTrigSeedPrescales = "-";
    } else {
      for (unsigned int j=0;j<vtmp.size();j++) {
	tempTrigSeeds = tempTrigSeeds + vtmp[j];
	tempTrigSeedPrescales += itmp[j];
	if (j<(vtmp.size()-1)) {
	  tempTrigSeeds = tempTrigSeeds + ",";
	  tempTrigSeedPrescales = tempTrigSeedPrescales + ",";
	}
      }
    }
    tempTrigSeeds.ReplaceAll("_","\\_");
    tempTrigSeedPrescales.ReplaceAll("_","\\_");
    
    outFile << "\\color{blue}"  << tempTrigName
	    << " & " << tempTrigSeeds
	    << " & " << tempTrigSeedPrescales
	    << " & " <<  menu->GetPrescale(i) 
	    << " & " << Rate[i] << " {$\\pm$ " << RateErr[i]
	    << "} & " << cumulRate
      	    << " & " << menu->GetEventsize(i)
	    << " & " << cuThru  << "\\\\"
	    << endl;
  }
  
  cumulRateErr = sqrt(cumulRateErr);
  cuThruErr = sqrt(cuThruErr);
  physCutThruErr = sqrt(physCutThruErr);
  cuPhysRateErr = sqrt(cuPhysRateErr);
  
  outFile << "\\hline \\multicolumn{6}{|l|}{\\bf \\boldmath Total HLT Physics rate (Hz), AlCa triggers not included } &  \\multicolumn{2}{|r|} { \\bf "<<  cuPhysRate << " $\\pm$ " << cuPhysRateErr << "} \\\\  \\hline" << endl;
  outFile << "\\hline \\multicolumn{6}{|l|}{\\bf \\boldmath Total Physics HLT throughput (MB/s), AlCa triggers not included }  & \\multicolumn{2}{|r|} { \\bf   "<<  physCutThru<< " $\\pm$ " << physCutThruErr << "} \\\\  \\hline" << endl;
  outFile << "\\hline \\multicolumn{6}{|l|}{\\bf \\boldmath Total HLT rate (Hz) }  &  \\multicolumn{2}{|r|} { \\bf "<<  cumulRate << " $\\pm$ " << cumulRateErr << "} \\\\  \\hline" << endl;
  outFile << "\\hline \\multicolumn{6}{|l|}{\\bf \\boldmath Total HLT throughput (MB/s) } &  \\multicolumn{2}{|r|} { \\bf  "<<  cuThru << " $\\pm$ " << cuThruErr << "} \\\\  \\hline" << endl; 

  // Footer for remaining L1 seeds
  for (unsigned int i=0;i<footTrigNames.size();i++) {
    outFile << "\\hline  { \\begin{tabular}{p{1.6cm}} \\scriptsize "<<  footTrigNames[i]
	    << "\\end{tabular} } & \\multicolumn{4}{|l|}{ \\begin{tabular}{p{7cm}} \\scriptsize " << footTrigSeeds[i]
	    << " \\end{tabular} } & \\multicolumn{3}{|l|}{ \\begin{tabular}{p{2cm}} \\scriptsize " << footTrigSeedPrescales[i] << " \\end{tabular} } \\\\  \\hline" << endl; 
  } 
  
  outFile << "\\end{longtable}" << endl;
    outFile << "\\end{footnotesize}" << endl;
    outFile << "\\clearpage" << endl;
    outFile << "\\end{landscape}" << endl;
    outFile << "\\end{document}";
    outFile.close();
}