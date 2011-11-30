
cmsenv
source setup.sh

#for DS in  MET ; do
for DS in AlphaT DoubleMu EleHadEG12 HTMHT MET MuHad MultiJet PhotonDoubleEle PhotonDoublePhoton PhotonHad PhotonPhoton RMR SingleMu; do
		ls -l ~/Trigger/HighPU/cfgs/hltmenu_HighPU_r179828_${DS}_forHighPU_cfg.py
		ls -l ~/Trigger/HighPU/cfgs/hltmenu_3E33_r178479_${DS}_forHighPU_cfg.py
		nohup ./OHltRateEff ~/Trigger/HighPU/cfgs/hltmenu_HighPU_r179828_${DS}_forHighPU_cfg.py > ~/Trigger/HighPU/logs/hltmenu_HighPU_r179828_${DS}_forHighPU_cfg.log &
		nohup ./OHltRateEff ~/Trigger/HighPU/cfgs/hltmenu_3E33_r178479_${DS}_forHighPU_cfg.py > ~/Trigger/HighPU/logs/hltmenu_3E33_r178479_${DS}_forHighPU_cfg.log &
done


##for DS in AlphaT DoubleMu EleHadEG12 HTMHT MET MuHad MultiJet PhotonDoubleEle PhotonDoublePhoton PhotonHad PhotonPhoton RMR SingleMu; do 
##		echo "if (DS==\"${DS}\") {"; 
##		f1=`ls -1 *_r178479_3e33_${DS}*.root`; 
##		f2=`ls -1 *_r179828_HighPU_${DS}*.root`; 
##		echo "    f[fileCounter] = new TFile(\""HighPU/root/$f1"\"); vlumiSFperFile[fileCounter++]=1;";
##		echo "    f[fileCounter] = new TFile(\""HighPU/root/$f2"\"); vlumiSFperFile[fileCounter++]=131.8;"; 
##		echo }; 
##done
