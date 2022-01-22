#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_HistoProducer.h"
#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Variables.h"
#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_Cuts.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"

using namespace ana;
using namespace ICARUSNumuXsec;

void run_HistoProducer(int whichInputFile, TString outDirBase){

  //==== Inputfile

  string inputDef = "";
  bool isDataInput = false;
  bool isCosmicInput = false;
  double TargetPOT = 6.6e20;
  TString str_TargetPOT = "6.6e20 POT";
  //====   BNB, new central sample, IcarusProd2021B_BNB_Nu_Cosmics_v09_28_01_01_01_caf
  if(whichInputFile==0){
    cout << "[run_Efficiency] Running BNB, IcarusProd2021B_BNB_Nu_Cosmics_v09_28_01_01_01_caf" << endl;
    inputDef = "IcarusProd2021B_BNB_Nu_Cosmics_v09_28_01_01_01_caf";
  }
  //====   NuMI, new central sample, IcarusProd2021B_NuMI_Nu_Cosmics_v09_28_01_01_01_caf
  else if(whichInputFile==1){
    cout << "[run_Efficiency] Running NuMI, IcarusProd2021B_NuMI_Nu_Cosmics_v09_28_01_01_01_caf" << endl;
    inputDef = "IcarusProd2021B_NuMI_Nu_Cosmics_v09_28_01_01_01_caf";
    TargetPOT = 6e20;
    str_TargetPOT = "6e20 POT";
  }
  //====   Cosmic, IcarusProd2021B_Intime_Cosmic_v09_28_01_01_01_caf
  else if(whichInputFile==2){
    isCosmicInput = true;
    cout << "[run_Efficiency] Running In-time BNB cosmic" << endl;
    inputDef = "IcarusProd2021B_Intime_Cosmic_v09_28_01_01_01_caf";
  }
  //====   Cosmic Data, howard_howard_numi_cosmic_run6480_full_2_CAFMaker_20210913_224007_821563
  else if(whichInputFile==3){
    isCosmicInput = true;
    isDataInput = true;
    cout << "[run_Efficiency] Running Cosmic data" << endl;
    inputDef = "howard_howard_numi_cosmic_run6480_full_2_CAFMaker_20210913_224007_821563";
    TargetPOT = 6e20;
    str_TargetPOT = "6e20 POT";
  }
  //====   cosmics in-time (is a little old)
  else if(whichInputFile==4){
    isCosmicInput = true;
    cout << "[run_Efficiency] Running /pnfs/icarus/persistent/users/dmendez/Tests/cosmics_hadded0.flatcaf.root" << endl;
    inputDef = "/pnfs/icarus/persistent/users/dmendez/Tests/cosmics_hadded0.flatcaf.root";
  }
  //====   cosmics out-time (release v09_17_01):
  else if(whichInputFile==5){
    isCosmicInput = true;
    cout << "[run_Efficiency] Running /icarus/data/users/mmr/CAFMaker/v09_20_00/NuMIBeamWindow/icarus_numi_offaxis_cosmic_v09_17_01/icarus_numi_offaxis_cosmic_v09_17_01_CAFMaker_out_flat.root" << endl;
    inputDef = "/icarus/data/users/mmr/CAFMaker/v09_20_00/NuMIBeamWindow/icarus_numi_offaxis_cosmic_v09_17_01/icarus_numi_offaxis_cosmic_v09_17_01_CAFMaker_out_flat.root";
    TargetPOT = 6e20;
    str_TargetPOT = "6e20 POT";
  }
  //====   old NuMI sample
  else if(whichInputFile==6){
    cout << "[run_Efficiency] Running /icarus/data/users/mmr/CAFMaker/v09_20_00/NuMIBeamWindow/icarus_numi_offaxis_cosmic_v09_17_01/icarus_numi_offaxis_cosmic_v09_17_01_CAFMaker_out_flat.root" << endl;
    inputDef = "/icarus/data/users/mmr/CAFMaker/v09_20_00/NuMIBeamWindow/icarus_numi_offaxis_cosmic_v09_17_01/icarus_numi_offaxis_cosmic_v09_17_01_CAFMaker_out_flat.root";
    TargetPOT = 6e20;
    str_TargetPOT = "6e20 POT";
  }
  else if(whichInputFile==7){
    cout << "[run_Efficiency] Running ReCAF NuMI+Cosmic, jsmedley_jsmedley_numi_cafmaker_CAFMaker_20210913_224007_821563" << endl;
    inputDef = "jsmedley_jsmedley_numi_cafmaker_CAFMaker_20210913_224007_821563";
    inputDef = "/pnfs/icarus/scratch/users/jskim/mc/ICARUS_NuMI_Nu_Cosmics/FlatCAF/icaruscode/v09_28_01_01_01/flatcaf_*.root";
    TargetPOT = 6e20;
    str_TargetPOT = "6e20 POT";
  }

  cout << "[run_Efficiency] Input definition : " << inputDef << endl;
  SpectrumLoader loader(inputDef);

  //==== Define samples;
  vector<TString> baseSampleNames;
  vector<Cut> baseSampleCuts;
  vector<SpillCut> baseSampleSpillCuts;

  baseSampleNames.push_back("AllSamples");
  baseSampleCuts.push_back(kNoCut);
  baseSampleSpillCuts.push_back(kNoSpillCut);

  if(!isDataInput && !isCosmicInput){

/*
    baseSampleNames.push_back("NuMuCC");
    baseSampleCuts.push_back(kIsNuMuCC);
    baseSampleSpillCuts.push_back(kNoSpillCut);

    baseSampleNames.push_back("NuMuNC");
    baseSampleCuts.push_back(cutIsNuMuNC);
    baseSampleSpillCuts.push_back(kNoSpillCut);
*/

    baseSampleNames.push_back("QE");
    baseSampleCuts.push_back(!kIsCosmic && cutIsQE);
    baseSampleSpillCuts.push_back(kNoSpillCut);

    baseSampleNames.push_back("Res");
    baseSampleCuts.push_back(!kIsCosmic && cutIsRes);
    baseSampleSpillCuts.push_back(kNoSpillCut);

    baseSampleNames.push_back("DIS");
    baseSampleCuts.push_back(!kIsCosmic && cutIsDIS);
    baseSampleSpillCuts.push_back(kNoSpillCut);

    baseSampleNames.push_back("Coh");
    baseSampleCuts.push_back(!kIsCosmic && cutIsCoh);
    baseSampleSpillCuts.push_back(kNoSpillCut);

    baseSampleNames.push_back("CohElastic");
    baseSampleCuts.push_back(!kIsCosmic && cutIsCohElastic);
    baseSampleSpillCuts.push_back(kNoSpillCut);

    baseSampleNames.push_back("ElectronScattering");
    baseSampleCuts.push_back(!kIsCosmic && cutIsElectronScattering);
    baseSampleSpillCuts.push_back(kNoSpillCut);

    baseSampleNames.push_back("IMDAnnihilation");
    baseSampleCuts.push_back(!kIsCosmic && cutIsIMDAnnihilation);
    baseSampleSpillCuts.push_back(kNoSpillCut);

    baseSampleNames.push_back("InverseBetaDecay");
    baseSampleCuts.push_back(!kIsCosmic && cutIsInverseBetaDecay);
    baseSampleSpillCuts.push_back(kNoSpillCut);

    baseSampleNames.push_back("GlashowResonance");
    baseSampleCuts.push_back(!kIsCosmic && cutIsGlashowResonance);
    baseSampleSpillCuts.push_back(kNoSpillCut);

    baseSampleNames.push_back("AMNuGamma");
    baseSampleCuts.push_back(!kIsCosmic && cutIsAMNuGamma);
    baseSampleSpillCuts.push_back(kNoSpillCut);

    baseSampleNames.push_back("MEC");
    baseSampleCuts.push_back(!kIsCosmic && cutIsMEC);
    baseSampleSpillCuts.push_back(kNoSpillCut);

    baseSampleNames.push_back("Diffractive");
    baseSampleCuts.push_back(!kIsCosmic && cutIsDiffractive);
    baseSampleSpillCuts.push_back(kNoSpillCut);

    baseSampleNames.push_back("EM");
    baseSampleCuts.push_back(!kIsCosmic && cutIsEM);
    baseSampleSpillCuts.push_back(kNoSpillCut);

    baseSampleNames.push_back("WeakMix");
    baseSampleCuts.push_back(!kIsCosmic && cutIsWeakMix);
    baseSampleSpillCuts.push_back(kNoSpillCut);

    baseSampleNames.push_back("UnknownInteractionType1");
    baseSampleCuts.push_back(!kIsCosmic && cutIsUnknownInteractionType1);
    baseSampleSpillCuts.push_back(kNoSpillCut);

    baseSampleNames.push_back("UnknownInteractionType2");
    baseSampleCuts.push_back(!kIsCosmic && cutIsUnknownInteractionType2);
    baseSampleSpillCuts.push_back(kNoSpillCut);

    baseSampleNames.push_back("UnknownInteractionType3");
    baseSampleCuts.push_back(!kIsCosmic && cutIsUnknownInteractionType3);
    baseSampleSpillCuts.push_back(kNoSpillCut);

    baseSampleNames.push_back("Cosmic");
    baseSampleCuts.push_back(kIsCosmic);
    baseSampleSpillCuts.push_back(kNoSpillCut);

  }
  else{

    baseSampleNames.push_back("AllSamples");
    baseSampleCuts.push_back(kNoCut);
    baseSampleSpillCuts.push_back(kNoSpillCut);

  }

  //==== Define Cuts

  vector<Cut> cuts;
  vector<SpillCut> spillcuts;
  vector<TString> cutNames;

  //====   1) Numu-cc selection

  cuts.push_back( kNoCut );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "NoCut" );

  cuts.push_back( kNoCut && cutRFiducial );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "RecoFiducial" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "NotClearCosmic" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && kNuScore );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "NuScoreGT0p4" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && kNuScore && cutFMScore );
  spillcuts.push_back( kNoSpillCut );
  cutNames.push_back( "FMScoreLT6" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && kNuScore && cutFMScore );
  spillcuts.push_back( spillcutSideCRTHitVetoFD );
  cutNames.push_back( "SideCRTVeto" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && kNuScore && cutFMScore );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "CRTVeto" );

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && kNuScore && cutFMScore && cutHasMuon );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "HasMuon" );

  //====     Contained/Exiting

  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && kNuScore && cutFMScore && cutHasMuon && cutMuonContained );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "HasMuon_Contained" );
  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && kNuScore && cutFMScore && cutHasMuon && !cutMuonContained );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "HasMuon_Exiting" );

/*
  if(!isDataInput && !isCosmicInput){

    //====       Matched to muon/proton

    cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && cutMuonContained && cutMuonMatchedToMuon );
    spillcuts.push_back( kCRTHitVetoFD );
    cutNames.push_back( "HasMuon_Contained_MatchedToMuon" );

    cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && cutMuonContained && cutMuonMatchedToProton );
    spillcuts.push_back( kCRTHitVetoFD );
    cutNames.push_back( "HasMuon_Contained_MatchedToProton" );

    cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && !cutMuonContained && cutMuonMatchedToMuon );
    spillcuts.push_back( kCRTHitVetoFD );
    cutNames.push_back( "HasMuon_Exiting_MatchedToMuon" );

    cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && !cutMuonContained && cutMuonMatchedToProton );
    spillcuts.push_back( kCRTHitVetoFD );
    cutNames.push_back( "HasMuon_Exiting_MatchedToProton" );

    //====     Additional particles

    cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && !cutTruthNoNeutron );
    spillcuts.push_back( kCRTHitVetoFD );
    cutNames.push_back( "HasMuon_WithNeutron" );

    cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && !cutTruthNoChargedPion );
    spillcuts.push_back( kCRTHitVetoFD );
    cutNames.push_back( "HasMuon_WithChargedPion" );

    cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && !cutTruthNoPiZero );
    spillcuts.push_back( kCRTHitVetoFD );
    cutNames.push_back( "HasMuon_WithNeutralPion" );

    cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && cutFMScore && kNuScore && cutHasMuon && cutTruthNoNeutron && cutTruthNoChargedPion && cutTruthNoPiZero );
    spillcuts.push_back( kCRTHitVetoFD );
    cutNames.push_back( "HasMuon_WithoutNeutronOrPion" );

  }
*/
  //====     HasProton added
  cuts.push_back( kNoCut && cutRFiducial && kNotClearCosmic && kNuScore && cutFMScore && cutHasMuon && cutHasProton );
  spillcuts.push_back( kCRTHitVetoFD );
  cutNames.push_back( "HasMuon_HasProton" );

  //==== Declare HistoProducer

  HistoProducer m;
  m.outputDir = outDirBase;
  m.initialize();
  m.setSystematicWeights(); // TODO Run this to get systematic variations


  m.TargetPOT = TargetPOT;
  m.str_TargetPOT = str_TargetPOT;

  for(unsigned int is=0; is<baseSampleNames.size(); is++){

    TString baseSampleName = baseSampleNames.at(is);
    Cut baseSampleCut = baseSampleCuts.at(is);
    SpillCut baseSampleSpillCut = baseSampleSpillCuts.at(is);

    bool isCosmic = (baseSampleName=="Cosmic");

    for(unsigned int i=0; i<cuts.size(); i++){

      Cut myCut = baseSampleCut && cuts.at(i);
      SpillCut mySpillCut = baseSampleSpillCut && spillcuts.at(i);
      TString cutName = baseSampleName+"_"+cutNames.at(i);

      if( m.setCut(cutName) ){

        if(isCosmic){

          m.bookSlice(loader, myCut, mySpillCut);
          if(cutName.Contains("HasMuon")){
            m.bookRecoMuon(loader, myCut, mySpillCut);
          }

          //==== For cosmic, remove truth-related variables to prevent warning messages
          map< TString, vector<Spectrum *> > map_cutName_to_vec_Spectrums;
          map< TString, vector<EnsembleSpectrum *> > map_cutName_to_vec_EnsembleSpectrums;
          for(std::map< TString, vector<Spectrum *> >::iterator mapit=map_cutName_to_vec_Spectrums.begin(); mapit!=map_cutName_to_vec_Spectrums.end(); mapit++){
            bool isTruthVariable = false;
            if(mapit->first.Contains("True")) isTruthVariable = true;
            if(mapit->first.Contains("Residual")) isTruthVariable = true;
            if(mapit->first.Contains("Truth")) isTruthVariable = true;
            map_cutName_to_vec_Spectrums.erase(mapit->first);
          }
          for(std::map< TString, vector<EnsembleSpectrum *> >::iterator mapit=map_cutName_to_vec_EnsembleSpectrums.begin(); mapit!=map_cutName_to_vec_EnsembleSpectrums.end(); mapit++){
            bool isTruthVariable = false;
            if(mapit->first.Contains("True")) isTruthVariable = true;
            if(mapit->first.Contains("Residual")) isTruthVariable = true;
            if(mapit->first.Contains("Truth")) isTruthVariable = true;
            map_cutName_to_vec_EnsembleSpectrums.erase(mapit->first);
          }

        }
        else{

          m.bookSlice(loader, myCut, mySpillCut);
          m.bookTruth(loader, myCut, mySpillCut);
          if(cutName.Contains("HasMuon")){
            m.bookRecoMuon(loader, myCut, mySpillCut);
          }
          if(cutName.Contains("HasProton")){
            m.bookRecoProton(loader, myCut, mySpillCut);
            m.bookRecoNeutrino(loader, myCut, mySpillCut);
          }

        }

      } // END if setCut

    } // END loop cut

  } // END sample cut

  loader.Go();

  m.saveHistograms();

}
