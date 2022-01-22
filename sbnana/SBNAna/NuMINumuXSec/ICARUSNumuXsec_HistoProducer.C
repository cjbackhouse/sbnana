#include "sbnana/SBNAna/NuMINumuXSec/ICARUSNumuXsec_HistoProducer.h"

using namespace ICARUSNumuXsec;

HistoProducer::HistoProducer(){

  cout << "[HistoProducer::HistoProducer] called" << endl;
  TargetPOT = 6.6e20;
  str_TargetPOT = "6.6e20 POT";
  outputName = "output.root";
  currentCutName = "DefaultCutName";
  vec_cutNames.clear();
  map_cutName_to_vec_Spectrums.clear();
  map_cutName_to_vec_EnsembleSpectrums.clear();
  cout << "[HistoProducer::HistoProducer] Finished" << endl;

}

void HistoProducer::initialize(){

  cout << "[HistoProducer::initialize] outputDir = " << outputDir << endl;
  gSystem->mkdir(outputDir, kTRUE);
  TString outputName = "output.root";
  outputfile = new TFile(outputDir+"/"+outputName,"RECREATE");

}

bool HistoProducer::setCut(TString cutName){

  currentCutName = cutName;
  if( std::find(vec_cutNames.begin(), vec_cutNames.end(), cutName) != vec_cutNames.end() ){
    cout << "[HistoProducer::setCut] cutName = " << cutName << " already exist.." << endl;
    return false;
  }
  else{
    cout << "[HistoProducer::setCut] adding cutName = " << cutName << endl;
    vec_cutNames.push_back(cutName);
    return true;
  }

}

void HistoProducer::bookSlice(SpectrumLoader& loader, Cut cut, SpillCut spillCut){

  cout << "[HistoProducer::bookSlice] called" << endl;
  //==== some binnings
  double xFMScoreMin = 0.;
  double xFMScoreMax = 150.;
  double dxFMScore = 1.;
  const Binning binsFMScore = Binning::Simple( int( (xFMScoreMax-xFMScoreMin)/dxFMScore ), xFMScoreMin, xFMScoreMax );
  double xFMTimeMin = -50;
  double xFMTimeMax = 50.;
  double dxFMTime = 1.;
  const Binning binsFMTime = Binning::Simple( int( (xFMTimeMax-xFMTimeMin)/dxFMTime ), xFMTimeMin, xFMTimeMax );
  const Binning binsOne = Binning::Simple( 1, 0., 1.);
  const Binning binsIntCode = Binning::Simple( 17, -2, 15 );

  cout << "[HistoProducer::bookSlice] Delcaring HistAxis" << endl;

  //==== HistAxis
  const HistAxis axCountSlice("CountSlice", binsOne, varCountSlice);
  const HistAxis axFMScore("FMScore", binsFMScore, varFMScore);
  const HistAxis axFMTime("FMTime", binsFMTime, varFMTime);
  const HistAxis axGENIEIntCode("GENIEIntCode", binsIntCode, varGENIEIntCode);

  cout << "[HistoProducer::bookSlice] Delcaring Spectrum" << endl;

  //==== Spectrum
  Spectrum *sCountSlice = new Spectrum(loader, axCountSlice, spillCut, cut);
  Spectrum *sFMScore = new Spectrum(loader, axFMScore, spillCut, cut);
  Spectrum *sFMTime = new Spectrum(loader, axFMTime, spillCut, cut);
  //Spectrum *sFMScore_vs_FMTime = new Spectrum(loader, axFMScore, axFMTime, spillCut, cut);
  Spectrum *sGENIEIntCode = new Spectrum(loader, axGENIEIntCode, spillCut, cut);

  cout << "[HistoProducer::bookTruth] Saving spectrums" << endl;

  map_cutName_to_vec_Spectrums[currentCutName].push_back(sCountSlice);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sFMScore);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sFMTime);
  //map_cutName_to_vec_Spectrums[currentCutName].push_back(sFMScore_vs_FMTime);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sGENIEIntCode);

}

void HistoProducer::bookTruth(SpectrumLoader& loader, Cut cut, SpillCut spillCut){

  cout << "[HistoProducer::bookTruth] called" << endl;
  //==== some binnings
  double xEnergyMin = 0.; // GeV
  double xEnergyMax = 10.;
  double dxEnergy = 0.1;
  const Binning binsEnergy = Binning::Simple( int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax );
  double xCosineThetaMin = -1;
  double xCosineThetaMax = 1.;
  double dxCosineTheta = 0.05;
  const Binning binsCosineTheta = Binning::Simple( int( (xCosineThetaMax-xCosineThetaMin)/dxCosineTheta ), xCosineThetaMin, xCosineThetaMax );
  const Binning binsN = Binning::Simple(5, 0., 5.);
  const Binning binsLength = Binning::Simple(2000, 0., 2000.);
  const Binning binsChi2B = Binning::Simple(1500, 0., 150);
  const Binning binsReducedChi2B = Binning::Simple(1000, 0., 10);

  cout << "[HistoProducer::bookTruth] Delcaring HistAxis" << endl;

  //==== HistAxis
  const HistAxis axNeutrinoTruthE("NeutrinoTruthE", binsEnergy, varNeutrinoTruthE);
  const HistAxis axTruthQ2("TruthQ2", binsEnergy, varTruthQ2);
  const HistAxis axTruthq0_lab("Truthq0_lab", binsEnergy, varTruthq0_lab);
  const HistAxis axTruthmodq_lab("Truthmodq_lab", binsEnergy, varTruthmodq_lab);
  //==== Muon
  const HistAxis axMuonTruthP("MuonTruthP", binsEnergy, varMuonTruthP);
  const HistAxis axMuonTruthCosineTheta("MuonTruthCosineTheta", binsCosineTheta, varMuonTruthCosineTheta);
  //==== Proton
  const HistAxis axProtonTruthP("ProtonTruthP", binsEnergy, varProtonTruthP);
  const HistAxis axProtonTruthCosineTheta("ProtonTruthCosineTheta", binsCosineTheta, varProtonTruthCosineTheta);
  //==== Muon+Proton
  const HistAxis axTruthMuonProtonCosineTheta("TruthMuonProtonCosineTheta", binsCosineTheta, varTruthMuonProtonCosineTheta);
  //==== Truth match
  const HistAxis axTruthMuonMatchedTrackRangeP("TruthMuonMatchedTrackRangeP", binsEnergy, varTruthMuonMatchedTrackRangeP);
  const HistAxis axTruthMuonMatchedTrackMCSP("TruthMuonMatchedTrackMCSP", binsEnergy, varTruthMuonMatchedTrackMCSP);
  const HistAxis axTruthMuonMatchedTrackCombinedP("TruthMuonMatchedTrackCombinedP", binsEnergy, varTruthMuonMatchedTrackCombinedP);
  const HistAxis axTruthProtonMatchedTrackRangeP("TruthProtonMatchedTrackRangeP", binsEnergy, varTruthProtonMatchedTrackRangeP);
  const HistAxis axTruthProtonMatchedTrackMCSP("TruthProtonMatchedTrackMCSP", binsEnergy, varTruthProtonMatchedTrackMCSP);
  const HistAxis axTruthProtonMatchedTrackCombinedP("TruthProtonMatchedTrackCombinedP", binsEnergy, varTruthProtonMatchedTrackCombinedP);
  const HistAxis axTruthProtonMatchedTrackLength("TruthProtonMatchedTrackLength", binsLength, varTruthProtonMatchedTrackLength);
  const HistAxis axTruthProtonMatchedTrackChi2Proton("TruthProtonMatchedTrackChi2Proton", binsChi2B, varTruthProtonMatchedTrackChi2Proton);
  const HistAxis axTruthProtonMatchedTrackReducedChi2Proton("TruthProtonMatchedTrackReducedChi2Proton", binsReducedChi2B, varTruthProtonMatchedTrackReducedChi2Proton);
  //==== stub
  //const HistAxis axNStub("NStub", binsN, varNStub);

  cout << "[HistoProducer::bookTruth] Delcaring Spectrum" << endl;

  //==== Spectrum
  Spectrum *sNeutrinoTruthE = new Spectrum(loader, axNeutrinoTruthE, spillCut, cut);
  Spectrum *sTruthQ2 = new Spectrum(loader, axTruthQ2, spillCut, cut);
  Spectrum *sTruthq0_lab = new Spectrum(loader, axTruthq0_lab, spillCut, cut);
  Spectrum *sTruthmodq_lab = new Spectrum(loader, axTruthmodq_lab, spillCut, cut);
  Spectrum *sTruth_modq_lab_vs_q0_lab = new Spectrum(loader, axTruthmodq_lab, axTruthq0_lab, spillCut, cut);
  //==== Muon
  Spectrum *sMuonTruthP = new Spectrum(loader, axMuonTruthP, spillCut, cut);
  Spectrum *sMuonTruthCosineTheta = new Spectrum(loader, axMuonTruthCosineTheta, spillCut, cut);
  Spectrum *sMuonTruth_P_vs_CosineTheta = new Spectrum(loader, axMuonTruthP, axMuonTruthCosineTheta, spillCut, cut);
  //==== Proton
  Spectrum *sProtonTruthP = new Spectrum(loader, axProtonTruthP, spillCut, cut);
  Spectrum *sProtonTruthCosineTheta = new Spectrum(loader, axProtonTruthCosineTheta, spillCut, cut);
  //==== Muon+Proton
  Spectrum *sTruthMuonProtonCosineTheta = new Spectrum(loader, axTruthMuonProtonCosineTheta, spillCut, cut);
  //==== Truth match
  Spectrum *sTruthMuonMatchedTrackRangeP = new Spectrum(loader, axTruthMuonMatchedTrackRangeP, spillCut, cut);
  Spectrum *sTruthMuonMatchedTrackMCSP = new Spectrum(loader, axTruthMuonMatchedTrackMCSP, spillCut, cut);
  Spectrum *sTruthMuonMatchedTrackCombinedP = new Spectrum(loader, axTruthMuonMatchedTrackCombinedP, spillCut, cut);
  Spectrum *sTruthProtonMatchedTrackRangeP = new Spectrum(loader, axTruthProtonMatchedTrackRangeP, spillCut, cut);
  Spectrum *sTruthProtonMatchedTrackMCSP = new Spectrum(loader, axTruthProtonMatchedTrackMCSP, spillCut, cut);
  Spectrum *sTruthProtonMatchedTrackCombinedP = new Spectrum(loader, axTruthProtonMatchedTrackCombinedP, spillCut, cut);
  Spectrum *sTruthProtonMatchedTrackLength = new Spectrum(loader, axTruthProtonMatchedTrackLength, spillCut, cut);
  Spectrum *sTruthProtonMatchedTrackChi2Proton = new Spectrum(loader, axTruthProtonMatchedTrackChi2Proton, spillCut, cut);
  Spectrum *sTruthProtonMatchedTrackReducedChi2Proton = new Spectrum(loader, axTruthProtonMatchedTrackReducedChi2Proton, spillCut, cut);
  //Spectrum *sTruthProtonMatchedTrack_Length_vs_Chi2Proton = new Spectrum(loader, axTruthProtonMatchedTrackLength, axTruthProtonMatchedTrackChi2Proton, spillCut, cut);
  //Spectrum *sTruthProtonMatchedTrack_Length_vs_ReducedChi2Proton = new Spectrum(loader, axTruthProtonMatchedTrackLength, axTruthProtonMatchedTrackReducedChi2Proton, spillCut, cut);
  //Spectrum *sNStub = new Spectrum(loader, axNStub, spillCut, cut);

  cout << "[HistoProducer::bookTruth] Saving spectrums" << endl;

  map_cutName_to_vec_Spectrums[currentCutName].push_back(sNeutrinoTruthE);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthQ2);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthq0_lab);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthmodq_lab);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruth_modq_lab_vs_q0_lab);
  //==== Muon
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonTruthP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonTruthCosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonTruth_P_vs_CosineTheta);
  //==== Proton
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonTruthP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonTruthCosineTheta);
  //==== Muon+Proton
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthMuonProtonCosineTheta);
  //==== Truth match
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthMuonMatchedTrackRangeP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthMuonMatchedTrackMCSP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthMuonMatchedTrackCombinedP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthProtonMatchedTrackRangeP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthProtonMatchedTrackMCSP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthProtonMatchedTrackCombinedP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthProtonMatchedTrackLength);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthProtonMatchedTrackChi2Proton);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthProtonMatchedTrackReducedChi2Proton);
  //map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthProtonMatchedTrack_Length_vs_Chi2Proton);
  //map_cutName_to_vec_Spectrums[currentCutName].push_back(sTruthProtonMatchedTrack_Length_vs_ReducedChi2Proton);
  //map_cutName_to_vec_Spectrums[currentCutName].push_back(sNStub);

  if(ISysts.size()>0){

    EnsembleSpectrum *esNeutrinoTruthE = new EnsembleSpectrum(loader, axNeutrinoTruthE, spillCut, cut, systWs);
    EnsembleSpectrum *esTruthQ2 = new EnsembleSpectrum(loader, axTruthQ2, spillCut, cut, systWs);
    EnsembleSpectrum *esMuonTruthP = new EnsembleSpectrum(loader, axMuonTruthP, spillCut, cut, systWs);
    EnsembleSpectrum *esMuonTruthCosineTheta = new EnsembleSpectrum(loader, axMuonTruthCosineTheta, spillCut, cut, systWs);
    EnsembleSpectrum *esProtonTruthP = new EnsembleSpectrum(loader, axProtonTruthP, spillCut, cut, systWs);

    map_cutName_to_vec_EnsembleSpectrums[currentCutName].push_back(esNeutrinoTruthE);
    map_cutName_to_vec_EnsembleSpectrums[currentCutName].push_back(esTruthQ2);
    map_cutName_to_vec_EnsembleSpectrums[currentCutName].push_back(esMuonTruthP);
    map_cutName_to_vec_EnsembleSpectrums[currentCutName].push_back(esMuonTruthCosineTheta);
    map_cutName_to_vec_EnsembleSpectrums[currentCutName].push_back(esProtonTruthP);
  }

  cout << "[HistoProducer::bookTruth] Finished" << endl;

}

void HistoProducer::bookRecoMuon(SpectrumLoader& loader, Cut cut, SpillCut spillCut){

  cout << "[HistoProducer::bookRecoMuon] called" << endl;

  //==== some binnings
  double xEnergyMin = 0.; // GeV
  double xEnergyMax = 10.;
  double dxEnergy = 0.1;
  const Binning binsEnergy = Binning::Simple( int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax );
  double xCosineThetaMin = -1;
  double xCosineThetaMax = 1.;
  double dxCosineTheta = 0.05;
  const Binning binsCosineTheta = Binning::Simple( int( (xCosineThetaMax-xCosineThetaMin)/dxCosineTheta ), xCosineThetaMin, xCosineThetaMax );
  const Binning binsEResidual = Binning::Simple(600, -3., 3.);
  const Binning binsEResidualFraction = Binning::Simple(600, -3., 3.);
  const Binning binsPDG = Binning::Simple(100000, -5000., 5000.);
  const Binning binsLength = Binning::Simple(2000, 0., 2000.);
  const Binning binsChi2 = Binning::Simple(200, 0., 200);
  const Binning binsReducedChi2 = Binning::Simple(100, 0., 10);

  cout << "[HistoProducer::bookRecoMuon] Delcaring HistAxis" << endl;

  //==== HistAxis
  const HistAxis axMuonRecoP("MuonRecoP", binsEnergy, varMuonRecoP);
  const HistAxis axMuonCaloPlane0P("MuonCaloPlane0P", binsEnergy, varMuonCaloPlane0P);
  const HistAxis axMuonCaloPlane1P("MuonCaloPlane1P", binsEnergy, varMuonCaloPlane1P);
  const HistAxis axMuonCaloPlane2P("MuonCaloPlane2P", binsEnergy, varMuonCaloPlane2P);
  const HistAxis axMuonTrueP("MuonTrueP", binsEnergy, varMuonTrueP);
  const HistAxis axMuonTruePDG("MuonTruePDG", binsPDG, varMuonTruePDG);
  const HistAxis axMuonPResidual("MuonPResidual", binsEResidual, varMuonPResidual);
  const HistAxis axMuonPResidualFraction("MuonPResidualFraction", binsEResidualFraction, varMuonPResidualFraction);
  const HistAxis axMuonCaloPlane0PResidualFraction("MuonCaloPlane0PResidualFraction", binsEResidualFraction, varMuonCaloPlane0PResidualFraction);
  const HistAxis axMuonCaloPlane1PResidualFraction("MuonCaloPlane1PResidualFraction", binsEResidualFraction, varMuonCaloPlane1PResidualFraction);
  const HistAxis axMuonCaloPlane2PResidualFraction("MuonCaloPlane2PResidualFraction", binsEResidualFraction, varMuonCaloPlane2PResidualFraction);
  const HistAxis axMuonRecoCosineTheta("MuonRecoCosineTheta", binsCosineTheta, varMuonRecoCosineTheta);
  const HistAxis axMuonTrueCosineTheta("MuonTrueCosineTheta", binsCosineTheta, varMuonTrueCosineTheta);
  const HistAxis axMuonCosineThetaResidual("MuonCosineThetaResidual", binsEResidual, varMuonCosineThetaResidual);
  const HistAxis axMuonLength("MuonLength", binsLength, varMuonLength);
  const HistAxis axMuonChi2Muon("MuonChi2Muon", binsChi2, varMuonChi2Muon);
  const HistAxis axMuonReducedChi2Muon("MuonReducedChi2Muon", binsReducedChi2, varMuonReducedChi2Muon);

  cout << "[HistoProducer::bookRecoMuon] Delcaring Spectrum" << endl;

  //==== Spectrum
  Spectrum *sMuonRecoP = new Spectrum(loader, axMuonRecoP, spillCut, cut);
  Spectrum *sMuonCaloPlane0P = new Spectrum(loader, axMuonCaloPlane0P, spillCut, cut);
  Spectrum *sMuonCaloPlane1P = new Spectrum(loader, axMuonCaloPlane1P, spillCut, cut);
  Spectrum *sMuonCaloPlane2P = new Spectrum(loader, axMuonCaloPlane2P, spillCut, cut);
  Spectrum *sMuonTrueP = new Spectrum(loader, axMuonTrueP, spillCut, cut);
  Spectrum *sMuonTruePDG = new Spectrum(loader, axMuonTruePDG, spillCut, cut);
  Spectrum *sMuon_TrueP_vs_RecoP = new Spectrum(loader, axMuonTrueP, axMuonRecoP, spillCut, cut);
  Spectrum *sMuon_RecoP_vs_TrueP = new Spectrum(loader, axMuonRecoP, axMuonTrueP, spillCut, cut);
  Spectrum *sMuonPResidual = new Spectrum(loader, axMuonPResidual, spillCut, cut);
  Spectrum *sMuonPResidualFraction = new Spectrum(loader, axMuonPResidualFraction, spillCut, cut);
  Spectrum *sMuonCaloPlane0PResidualFraction = new Spectrum(loader, axMuonCaloPlane0PResidualFraction, spillCut, cut);
  Spectrum *sMuonCaloPlane1PResidualFraction = new Spectrum(loader, axMuonCaloPlane1PResidualFraction, spillCut, cut);
  Spectrum *sMuonCaloPlane2PResidualFraction = new Spectrum(loader, axMuonCaloPlane2PResidualFraction, spillCut, cut);
  Spectrum *sMuonRecoCosineTheta = new Spectrum(loader, axMuonRecoCosineTheta, spillCut, cut);
  Spectrum *sMuonTrueCosineTheta = new Spectrum(loader, axMuonTrueCosineTheta, spillCut, cut);
  Spectrum *sMuon_TrueCosineTheta_vs_RecoCosineTheta = new Spectrum(loader, axMuonTrueCosineTheta, axMuonRecoCosineTheta, spillCut, cut);
  Spectrum *sMuon_RecoCosineTheta_vs_TrueCosineTheta = new Spectrum(loader, axMuonRecoCosineTheta, axMuonTrueCosineTheta, spillCut, cut);
  Spectrum *sMuonCosineThetaResidual = new Spectrum(loader, axMuonCosineThetaResidual, spillCut, cut);
  Spectrum *sMuon_True_P_vs_CosineTheta = new Spectrum(loader, axMuonTrueP, axMuonTrueCosineTheta, spillCut, cut);
  Spectrum *sMuon_Reco_P_vs_CosineTheta = new Spectrum(loader, axMuonRecoP, axMuonRecoCosineTheta, spillCut, cut);
  //==== x : cosine theta, y : reco-true / true
  Spectrum *sMuon_TrueCosineTheta_vs_TrueCosineThetaResidual = new Spectrum(loader, axMuonTrueCosineTheta, axMuonCosineThetaResidual, spillCut, cut);
  Spectrum *sMuonLength = new Spectrum(loader, axMuonLength, spillCut, cut);
  Spectrum *sMuonChi2Muon = new Spectrum(loader, axMuonChi2Muon, spillCut, cut);
  Spectrum *sMuonReducedChi2Muon = new Spectrum(loader, axMuonReducedChi2Muon, spillCut, cut);
  //Spectrum *sMuon_Length_vs_Chi2Muon = new Spectrum(loader, axMuonLength, axMuonChi2Muon, spillCut, cut);
  //Spectrum *sMuon_Length_vs_ReducedChi2Muon = new Spectrum(loader, axMuonLength, axMuonReducedChi2Muon, spillCut, cut);

/*
  //==== Systematics
  for(unsigned int i=0; i<systWs.size(); i++){
    Spectrum *sMuonRecoPSyst = new Spectrum(loader, axMuonRecoP, spillCut, cut, kNoShift, systWs.at(i));
    std::string baseLabel = sMuonRecoPSyst->GetLabels()[0];
    std::string newLabel = baseLabel+"_"+ISysts.at(i)->ShortName();
    sMuonRecoPSyst->SetLabel(0, newLabel);
    map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonRecoPSyst);
  }
*/

  cout << "[HistoProducer::bookRecoMuon] Saving spectrums" << endl;

  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonRecoP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonCaloPlane0P);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonCaloPlane1P);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonCaloPlane2P);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonTrueP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonTruePDG);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuon_TrueP_vs_RecoP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuon_RecoP_vs_TrueP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonPResidual);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonPResidualFraction);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonCaloPlane0PResidualFraction);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonCaloPlane1PResidualFraction);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonCaloPlane2PResidualFraction);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonRecoCosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonTrueCosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuon_TrueCosineTheta_vs_RecoCosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuon_RecoCosineTheta_vs_TrueCosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonCosineThetaResidual);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuon_True_P_vs_CosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuon_Reco_P_vs_CosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuon_TrueCosineTheta_vs_TrueCosineThetaResidual);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonLength);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonChi2Muon);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuonReducedChi2Muon);
  //map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuon_Length_vs_Chi2Muon);
  //map_cutName_to_vec_Spectrums[currentCutName].push_back(sMuon_Length_vs_ReducedChi2Muon);

  if(ISysts.size()>0){

    EnsembleSpectrum *esMuonRecoP = new EnsembleSpectrum(loader, axMuonRecoP, spillCut, cut, systWs);
    EnsembleSpectrum *esMuonTrueP = new EnsembleSpectrum(loader, axMuonTrueP, spillCut, cut, systWs);
    EnsembleSpectrum *esMuonRecoCosineTheta = new EnsembleSpectrum(loader, axMuonRecoCosineTheta, spillCut, cut, systWs);
    EnsembleSpectrum *esMuonTrueCosineTheta = new EnsembleSpectrum(loader, axMuonTrueCosineTheta, spillCut, cut, systWs);

    map_cutName_to_vec_EnsembleSpectrums[currentCutName].push_back(esMuonRecoP);
    map_cutName_to_vec_EnsembleSpectrums[currentCutName].push_back(esMuonTrueP);
    map_cutName_to_vec_EnsembleSpectrums[currentCutName].push_back(esMuonRecoCosineTheta);
    map_cutName_to_vec_EnsembleSpectrums[currentCutName].push_back(esMuonTrueCosineTheta);
  }

  cout << "[HistoProducer::bookRecoMuon] Finished" << endl;

}

void HistoProducer::bookRecoProton(SpectrumLoader& loader, Cut cut, SpillCut spillCut){

  cout << "[HistoProducer::bookRecoProton] called" << endl;

  //==== some binnings
  double xEnergyMin = 0.; // GeV
  double xEnergyMax = 10.;
  double dxEnergy = 0.1;
  const Binning binsEnergy = Binning::Simple( int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax );
  double xCosineThetaMin = -1;
  double xCosineThetaMax = 1.;
  double dxCosineTheta = 0.05;
  const Binning binsCosineTheta = Binning::Simple( int( (xCosineThetaMax-xCosineThetaMin)/dxCosineTheta ), xCosineThetaMin, xCosineThetaMax );
  const Binning binsEResidual = Binning::Simple(600, -3., 3.);
  const Binning binsEResidualFraction = Binning::Simple(600, -3., 3.);
  const Binning binsPDG = Binning::Simple(100000, -5000., 5000.);
  const Binning binsLength = Binning::Simple(2000, 0., 200.);
  const Binning binsChi2 = Binning::Simple(1500, 0., 150);
  const Binning binsReducedChi2 = Binning::Simple(1000, 0., 10);

  cout << "[HistoProducer::bookRecoProton] Delcaring HistAxis" << endl;

  //==== HistAxis
  const HistAxis axProtonRecoP("ProtonRecoP", binsEnergy, varProtonRecoP);
  const HistAxis axProtonCaloP("ProtonCaloP", binsEnergy, varProtonCaloP);
  const HistAxis axProtonTrueP("ProtonTrueP", binsEnergy, varProtonTrueP);
  const HistAxis axProtonTruePDG("ProtonTruePDG", binsPDG, varProtonTruePDG);
  const HistAxis axProtonPResidual("ProtonPResidual", binsEResidual, varProtonPResidual);
  const HistAxis axProtonPResidualFraction("ProtonPResidualFraction", binsEResidualFraction, varProtonPResidualFraction);
  const HistAxis axProtonRecoCosineTheta("ProtonRecoCosineTheta", binsCosineTheta, varProtonRecoCosineTheta);
  const HistAxis axProtonTrueCosineTheta("ProtonTrueCosineTheta", binsCosineTheta, varProtonTrueCosineTheta);
  const HistAxis axProtonCosineThetaResidual("ProtonCosineThetaResidual", binsEResidual, varProtonCosineThetaResidual);
  const HistAxis axProtonLength("ProtonLength", binsLength, varProtonLength);
  const HistAxis axProtonChi2Proton("ProtonChi2Proton", binsChi2, varProtonChi2Proton);
  const HistAxis axProtonReducedChi2Proton("ProtonReducedChi2Proton", binsReducedChi2, varProtonReducedChi2Proton);

  cout << "[HistoProducer::bookRecoProton] Delcaring Spectrum" << endl;

  //==== Spectrum
  Spectrum *sProtonRecoP = new Spectrum(loader, axProtonRecoP, spillCut, cut);
  Spectrum *sProtonCaloP = new Spectrum(loader, axProtonCaloP, spillCut, cut);
  Spectrum *sProtonTrueP = new Spectrum(loader, axProtonTrueP, spillCut, cut);
  Spectrum *sProtonTruePDG = new Spectrum(loader, axProtonTruePDG, spillCut, cut);
  Spectrum *sProton_TrueP_vs_RecoP = new Spectrum(loader, axProtonTrueP, axProtonRecoP, spillCut, cut);
  Spectrum *sProton_RecoP_vs_TrueP = new Spectrum(loader, axProtonRecoP, axProtonTrueP, spillCut, cut);
  Spectrum *sProtonPResidual = new Spectrum(loader, axProtonPResidual, spillCut, cut);
  Spectrum *sProtonPResidualFraction = new Spectrum(loader, axProtonPResidualFraction, spillCut, cut);
  Spectrum *sProtonRecoCosineTheta = new Spectrum(loader, axProtonRecoCosineTheta, spillCut, cut);
  Spectrum *sProtonTrueCosineTheta = new Spectrum(loader, axProtonTrueCosineTheta, spillCut, cut);
  Spectrum *sProton_TrueCosineTheta_vs_RecoCosineTheta = new Spectrum(loader, axProtonTrueCosineTheta, axProtonRecoCosineTheta, spillCut, cut);
  Spectrum *sProton_RecoCosineTheta_vs_TrueCosineTheta = new Spectrum(loader, axProtonRecoCosineTheta, axProtonTrueCosineTheta, spillCut, cut);
  Spectrum *sProtonCosineThetaResidual = new Spectrum(loader, axProtonCosineThetaResidual, spillCut, cut);
  Spectrum *sProton_True_P_vs_CosineTheta = new Spectrum(loader, axProtonTrueP, axProtonTrueCosineTheta, spillCut, cut);
  Spectrum *sProton_Reco_P_vs_CosineTheta = new Spectrum(loader, axProtonRecoP, axProtonRecoCosineTheta, spillCut, cut);
  Spectrum *sProtonLength = new Spectrum(loader, axProtonLength, spillCut, cut);
  Spectrum *sProtonChi2Proton = new Spectrum(loader, axProtonChi2Proton, spillCut, cut);
  Spectrum *sProtonReducedChi2Proton = new Spectrum(loader, axProtonReducedChi2Proton, spillCut, cut);
  //Spectrum *sProton_Length_vs_Chi2Proton = new Spectrum(loader, axProtonLength, axProtonChi2Proton, spillCut, cut);
  //Spectrum *sProton_Length_vs_ReducedChi2Proton = new Spectrum(loader, axProtonLength, axProtonReducedChi2Proton, spillCut, cut);

  cout << "[HistoProducer::bookRecoProton] Saving spectrums" << endl;

  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonRecoP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonCaloP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonTrueP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonTruePDG);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProton_TrueP_vs_RecoP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProton_RecoP_vs_TrueP);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonPResidual);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonPResidualFraction);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonRecoCosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonTrueCosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProton_TrueCosineTheta_vs_RecoCosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProton_RecoCosineTheta_vs_TrueCosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonCosineThetaResidual);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProton_True_P_vs_CosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProton_Reco_P_vs_CosineTheta);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonLength);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonChi2Proton);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sProtonReducedChi2Proton);
  //map_cutName_to_vec_Spectrums[currentCutName].push_back(sProton_Length_vs_Chi2Proton);
  //map_cutName_to_vec_Spectrums[currentCutName].push_back(sProton_Length_vs_ReducedChi2Proton);

  if(ISysts.size()>0){

    EnsembleSpectrum *esProtonRecoP = new EnsembleSpectrum(loader, axProtonRecoP, spillCut, cut, systWs);
    EnsembleSpectrum *esProtonTrueP = new EnsembleSpectrum(loader, axProtonTrueP, spillCut, cut, systWs);
    EnsembleSpectrum *esProtonRecoCosineTheta = new EnsembleSpectrum(loader, axProtonRecoCosineTheta, spillCut, cut, systWs);
    EnsembleSpectrum *esProtonTrueCosineTheta = new EnsembleSpectrum(loader, axProtonTrueCosineTheta, spillCut, cut, systWs);

    map_cutName_to_vec_EnsembleSpectrums[currentCutName].push_back(esProtonRecoP);
    map_cutName_to_vec_EnsembleSpectrums[currentCutName].push_back(esProtonTrueP);
    map_cutName_to_vec_EnsembleSpectrums[currentCutName].push_back(esProtonRecoCosineTheta);
    map_cutName_to_vec_EnsembleSpectrums[currentCutName].push_back(esProtonTrueCosineTheta);
  }

  cout << "[HistoProducer::bookRecoProton] Finished" << endl;

}

void HistoProducer::bookRecoNeutrino(SpectrumLoader& loader, Cut cut, SpillCut spillCut){

  cout << "[HistoProducer::bookRecoNeutrino] called" << endl;

  //==== some binnings
  double xEnergyMin = 0.; // GeV
  double xEnergyMax = 10.;
  double dxEnergy = 0.1;
  const Binning binsEnergy = Binning::Simple( int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax );
  const Binning binsEResidual = Binning::Simple(600, -3., 3.);
  const Binning binsEResidualFraction = Binning::Simple(600, -3., 3.);

  cout << "[HistoProducer::bookRecoNeutrino] Delcaring HistAxis" << endl;

  //==== HistAxis
  const HistAxis axNeutrinoCombinedEnergy("NeutrinoCombinedEnergy", binsEnergy, varNeutrinoCombinedEnergy);
  const HistAxis axNeutrinoTestEnergy("NeutrinoTestEnergy", binsEnergy, varNeutrinoTestEnergy);
  const HistAxis axNeutrinoQE("NeutrinoQE", binsEnergy, varNeutrinoQE);
  const HistAxis axNeutrinoCombinedEnergyResidual("NeutrinoCombinedEnergyResidual", binsEResidual, varNeutrinoCombinedEnergyResidual);
  const HistAxis axNeutrinoQEResidual("NeutrinoQEResidual", binsEResidual, varNeutrinoQEResidual);
  const HistAxis axNeutrinoCombinedEnergyResidualFraction("NeutrinoCombinedEnergyResidualFraction", binsEResidualFraction, varNeutrinoCombinedEnergyResidualFraction);
  const HistAxis axNeutrinoQEResidualFraction("NeutrinoQEResidualFraction", binsEResidualFraction, varNeutrinoQEResidualFraction);

  cout << "[HistoProducer::bookRecoNeutrino] Delcaring Spectrum" << endl;

  //==== Spectrum
  Spectrum *sNeutrinoCombinedEnergy = new Spectrum(loader, axNeutrinoCombinedEnergy, spillCut, cut);
  Spectrum *sNeutrinoTestEnergy = new Spectrum(loader, axNeutrinoTestEnergy, spillCut, cut);
  Spectrum *sNeutrinoQE = new Spectrum(loader, axNeutrinoQE, spillCut, cut);
  Spectrum *sNeutrinoCombinedEnergyResidual = new Spectrum(loader, axNeutrinoCombinedEnergyResidual, spillCut, cut);
  Spectrum *sNeutrinoQEResidual = new Spectrum(loader, axNeutrinoQEResidual, spillCut, cut);
  Spectrum *sNeutrinoCombinedEnergyResidualFraction = new Spectrum(loader, axNeutrinoCombinedEnergyResidualFraction, spillCut, cut);
  Spectrum *sNeutrinoQEResidualFraction = new Spectrum(loader, axNeutrinoQEResidualFraction, spillCut, cut);

  cout << "[HistoProducer::bookRecoNeutrino] Saving spectrums" << endl;

  map_cutName_to_vec_Spectrums[currentCutName].push_back(sNeutrinoCombinedEnergy);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sNeutrinoTestEnergy);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sNeutrinoQE);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sNeutrinoCombinedEnergyResidual);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sNeutrinoQEResidual);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sNeutrinoCombinedEnergyResidualFraction);
  map_cutName_to_vec_Spectrums[currentCutName].push_back(sNeutrinoQEResidualFraction);

  if(ISysts.size()>0){

    EnsembleSpectrum *esNeutrinoCombinedEnergy = new EnsembleSpectrum(loader, axNeutrinoCombinedEnergy, spillCut, cut, systWs);
    EnsembleSpectrum *esNeutrinoQE = new EnsembleSpectrum(loader, axNeutrinoQE, spillCut, cut, systWs);

    map_cutName_to_vec_EnsembleSpectrums[currentCutName].push_back(esNeutrinoCombinedEnergy);
    map_cutName_to_vec_EnsembleSpectrums[currentCutName].push_back(esNeutrinoQE);
  }

  cout << "[HistoProducer::bookRecoNeutrino] Finished" << endl;

}

void HistoProducer::saveHistograms(){

  outputfile->cd();

  cout << "[HistoProducer::saveHistograms] Number of cuts = " << vec_cutNames.size() << endl;
  const unsigned int nCutName = vec_cutNames.size();
  if(nCutName>0){

    TH1D *hist_CutNames = new TH1D("hist_CutNames", "", nCutName, 0., 1.*nCutName);

    for(unsigned int ic=0; ic<nCutName; ic++){

      TString cutName = vec_cutNames.at(ic);
      //TString dirName = "cutName"+TString::Itoa(ic,10);
      TString dirName = cutName; // I'd rather make the cutName simpler and use it as dirName too
      TDirectory *dir = outputfile->GetDirectory(dirName);
      if(!dir){
        outputfile->mkdir(dirName);
        dir = outputfile->GetDirectory(dirName);
      }
      outputfile->cd(dirName);

      vector<Spectrum *> vec_Spectrums = map_cutName_to_vec_Spectrums[cutName];
      vector<EnsembleSpectrum *> vec_EnsembleSpectrums = map_cutName_to_vec_EnsembleSpectrums[cutName];

      //==== BinContetn = number of Spectrum
      //==== BinError = number of EnsembleSpectrum
      hist_CutNames->SetBinContent(ic+1, vec_Spectrums.size());
      hist_CutNames->SetBinError(ic+1, vec_EnsembleSpectrums.size());
      hist_CutNames->GetXaxis()->SetBinLabel(ic+1, cutName);

      cout << "[HistoProducer::saveHistograms] cutName = " << cutName << endl;
      cout << "[HistoProducer::saveHistograms]   Directory name = " << dirName << endl;

      //==== Spectrum
      cout << "[HistoProducer::saveHistograms]   Number of Spectrum = " << vec_Spectrums.size() << endl;
      for(unsigned int i=0; i<vec_Spectrums.size(); i++){

        cout << "[HistoProducer::saveHistograms]   " << i << "-th Spectrum.." << endl;

        if(i==0){
          cout << "[HistoProducer::saveHistograms]     POT = " << vec_Spectrums.at(i)->POT() << endl;
          cout << "[HistoProducer::saveHistograms]     Livetime = " << vec_Spectrums.at(i)->Livetime() << endl;
          TH1D *hist_POT = new TH1D("POT_"+dirName, "POT", 1, 0., 1.);
          hist_POT->SetBinContent(1, vec_Spectrums.at(i)->POT());
          TH1D *hist_Livetime = new TH1D("Livetime_"+dirName, "Livetime", 1, 0., 1.);
          hist_Livetime->SetBinContent(1, vec_Spectrums.at(i)->Livetime());
          hist_POT->Write();
          hist_Livetime->Write();
        }

        if(vec_Spectrums.at(i)->GetBinnings().size()==1){
          TH1 *h = vec_Spectrums.at(i)->ToTH1(TargetPOT);
          TString hName = vec_Spectrums.at(i)->GetLabels()[0];
          cout << "[HistoProducer::saveHistograms]     Writing TH1, \"" << hName << "\"" << endl;
          h->SetName(hName+"_"+dirName);
          h->Write();
          cout << "[HistoProducer::saveHistograms]     --> Done" << endl;
        }
        else if(vec_Spectrums.at(i)->GetBinnings().size()==2){
          TH2 *h = vec_Spectrums.at(i)->ToTH2(TargetPOT);
          TString hName = vec_Spectrums.at(i)->GetLabels()[0];
          cout << "[HistoProducer::saveHistograms]     Writing TH2, \"" << hName << "\"" << endl;
          h->SetName(hName+"_"+dirName);
          h->Write();
        }

      }

      //==== EnsembleSpectrum
      cout << "[HistoProducer::saveHistograms]   Number of EnsembleSpectrum = " << vec_EnsembleSpectrums.size() << endl;

      for(unsigned int i=0; i<vec_EnsembleSpectrums.size(); i++){

        cout << "[HistoProducer::saveHistograms]   " << i << "-th EnsembleSpectrum.." << endl;

        //==== Write nominal
        Spectrum sNominal = vec_EnsembleSpectrums.at(i)->Nominal();
        TString baseLabel = sNominal.GetLabels()[0];
        TString newLabel = baseLabel+"_TestES_Nominal";
        TH1 *h = sNominal.ToTH1(TargetPOT);
        h->SetName(newLabel+"_"+dirName);
        h->Write();
        cout << "[HistoProducer::saveHistograms]     Nominal histogram = " << baseLabel << endl;
        //==== Write Universes
        for(unsigned int iu=0; iu<vec_EnsembleSpectrums.at(i)->NUniverses(); ++iu){
          TH1 *hU = vec_EnsembleSpectrums.at(i)->Universe(iu).ToTH1(TargetPOT);
          TString systName = ISysts.at(iu)->ShortName();
          hU->SetName(baseLabel+"_TestES_"+systName+"_"+dirName);
          hU->Write();
        }
        TGraphAsymmErrors *grErrorBand = vec_EnsembleSpectrums.at(i)->ErrorBand(TargetPOT);
        grErrorBand->SetName(baseLabel+"_TestES_ErrorBand_"+dirName);
        grErrorBand->Write();

      }

      outputfile->cd();

    } // END loop cutname

    outputfile->cd();
    hist_CutNames->Write();

  } // END if nCutName>0

  outputfile->cd();
  TH1D *hist_TargetPOT = new TH1D("hist_TargetPOT", "", 1, 0., 1.);
  hist_TargetPOT->SetBinContent(1, TargetPOT);
  hist_TargetPOT->Write();

}

void HistoProducer::setSystematicWeights(){

  ISysts = GetSBNGenieWeightSysts();
  for(unsigned int i=0; i<ISysts.size(); ++i){
    systWs.push_back(GetUniverseWeight(ISysts, i));
  }

}

HistoProducer::~HistoProducer(){

  cout << "[HistoProducer::~HistoProducer] called" << endl;

  outputfile->Close();

  cout << "[HistoProducer::~HistoProducer] output file : " << outputDir+"/"+outputName << endl;

  cout << "[HistoProducer::~HistoProducer] Finished" << endl;

}

