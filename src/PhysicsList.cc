//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file radioactivedecay/rdecay02/src/PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "PhysListParticles.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "PhysListEmLowEnergy.hh"
#include "PhysListHadron.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4UnitsTable.hh"
#include "G4LossTableManager.hh"

#include "G4HadronPhysicsQGSP_BERT.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"

#include "G4EmExtraPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4DecayPhysics.hh"
#include "G4EmProcessOptions.hh"

#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PhysicsList::SetCutScale(G4double SF)
{
/**********************************************\
  defaultCutValue     = 1.0*micrometer*SF;
  fCutForGamma         = defaultCutValue;
  fCutForElectron      = 0.1*defaultCutValue;
  fCutForPositron      = defaultCutValue;
  fCutForProton        = defaultCutValue;
  fCutForAlpha         = 1.e-3*defaultCutValue;
  fCutForGenericIon    = 1.e-6*defaultCutValue;
\**********************************************/
  fSF = SF;
  defaultCutValue     = CLHEP::nanometer*fSF;
  fCutForGamma         = 1.e4*defaultCutValue;
  fCutForElectron      = 1.e3*defaultCutValue;
  fCutForPositron      = 1.e3*defaultCutValue;
  fCutForProton        = 1.e3*defaultCutValue;
  fCutForAlpha         = 10.*defaultCutValue;
  fCutForGenericIon    = defaultCutValue;
  fCutForOpticalPhoton = 1.0*CLHEP::mm;
  SetCuts();
}

#include "G4MaxTimeCuts.hh"
#include "G4MinEkineCuts.hh"

void PhysicsList::AddTransportation() {
  
auto theParticleIterator = GetParticleIterator();

  G4VUserPhysicsList::AddTransportation();
  
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    // time cuts for ONLY neutrons:
    if(particleName == "neutron" ) 
    pmanager->AddDiscreteProcess(new G4MaxTimeCuts());
    // Energy cuts to kill charged (embedded in method) particles:
    pmanager->AddDiscreteProcess(new G4MinEkineCuts());
  }		      
}

void PhysicsList::PhysListInitialize(bool fActiveNeutrons, bool fActiveOptical) 
{
  G4LossTableManager::Instance();
  SetCutScale(1.);
  fradDecay = true;

  fDetectorCuts = 0;
  fTargetCuts   = 0;

  ActiveNeutrons = fActiveNeutrons;
  ActiveOptical = fActiveOptical;

  fPMessenger = new PhysicsListMessenger(this);

  SetVerboseLevel(1);

  //default physics
  fParticleList = new G4DecayPhysics();

  //default physics
  fRaddecayList = new G4RadioactiveDecayPhysics();

  // EM physics
  G4EmLivermorePhysics *fEmLivermore;
  fEmPhysicsList = fEmLivermore = new G4EmLivermorePhysics(aIonEner,aIonStep);

  G4EmProcessOptions emOptions;
  emOptions.SetBuildCSDARange(true);
  emOptions.SetDEDXBinningForCSDARange(10*10);
  emOptions.SetFluo(true); // To activate deexcitation processes and fluorescence
  emOptions.SetAuger(true); // To activate Auger effect if deexcitation is activated
  emOptions.SetPIXE(true); // To activate Particle Induced X-Ray Emission (PIXE) 
  emOptions.SetPIXECrossSectionModel("ecpssr_interpolated"); // To activate Particle Induced X-Ray Emission (PIXE) 
  
  // Had physics 
  fHadPhysicsList = new G4HadronPhysicsQGSP_BERT_HP("std-hadron");
  //aBertPhysicsList = static_cast<G4HadronPhysicsQGSP_BERT_HP*> (fHadPhysicsList);//MLCC
  //aBertPhysicsList->ActivateNeutrons(ActiveNeutrons);

//  fHadPhysicsList = new PhysListHadron("std-hadron");
//  aHadPhysicsList = static_cast<PhysListHadron*> (fHadPhysicsList);
//  aHadPhysicsList->PhysListHadronSetNeutrons(ActiveNeutrons);

  AddExtraBuilders(false);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
  delete fPMessenger;
  if( fradDecay ) delete fRaddecayList;
  delete fEmPhysicsList;
  if (fHadPhysicsList) delete fHadPhysicsList;
  if (fNhadcomp > 0) {
    for(G4int i=0; i<fNhadcomp; i++) {
      delete fHadronPhys[i];
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  fParticleList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4OpticalPhysics.hh"
#include "G4OpticalProcessIndex.hh"

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  // em
  fEmPhysicsList->ConstructProcess();
  // decays
  fParticleList->ConstructProcess();
  if( fradDecay  ) fRaddecayList->ConstructProcess();
  // had
  if (fNhadcomp > 0) {
    for(G4int i=0; i<fNhadcomp; i++) {
      (fHadronPhys[i])->ConstructProcess();
    }
  }
  if (fHadPhysicsList) fHadPhysicsList->ConstructProcess();

  if( ActiveOptical ) ConstructOp();
  G4cout << "### PhysicsList::ConstructProcess is done" << fEmPhysicsList->GetPhysicsName() << G4endl;
 
  if( ActiveOptical ) G4cout << "+++++++++ OPTICAL ACTIVE " << G4endl;
}

#include "G4OpticalPhysics.hh"
#include "G4OpticalProcessIndex.hh"

void PhysicsList::ConstructOp()
{

auto theParticleIterator = GetParticleIterator();

  G4Cerenkov*          theCerenkovProcess;
  G4Scintillation*     theScintillationProcess;
  G4OpAbsorption*      theAbsorptionProcess;
  G4OpRayleigh*        theRayleighScatteringProcess;
  G4OpMieHG*           theMieHGScatteringProcess;
  G4OpBoundaryProcess* theBoundaryProcess;

  theCerenkovProcess           = new G4Cerenkov("Cerenkov");
  theScintillationProcess      = new G4Scintillation("Scintillation");
  theAbsorptionProcess         = new G4OpAbsorption();
  theRayleighScatteringProcess = new G4OpRayleigh();
  theMieHGScatteringProcess    = new G4OpMieHG();
  theBoundaryProcess           = new G4OpBoundaryProcess();

//  theCerenkovProcess->DumpPhysicsTable();
//  theScintillationProcess->DumpPhysicsTable();
//  theRayleighScatteringProcess->DumpPhysicsTable();
  
  theCerenkovProcess->SetMaxNumPhotonsPerStep(100);
  theCerenkovProcess->SetMaxBetaChangePerStep(10.0);
  theCerenkovProcess->SetTrackSecondariesFirst(true);
  
  theScintillationProcess->SetScintillationYieldFactor(1.);
  theScintillationProcess->SetTrackSecondariesFirst(true);

  // Use Birks Correction in the Scintillation process

  G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
  theScintillationProcess->AddSaturation(emSaturation);

  theScintillationProcess->SetScintillationByParticleType(true);
 
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    if (theCerenkovProcess->IsApplicable(*particle)) {
      pmanager->AddProcess(theCerenkovProcess);
      pmanager->SetProcessOrdering(theCerenkovProcess,idxPostStep);
    }
    if (theScintillationProcess->IsApplicable(*particle)) {
      pmanager->AddProcess(theScintillationProcess);
      pmanager->SetProcessOrderingToLast(theScintillationProcess, idxAtRest);
      pmanager->SetProcessOrderingToLast(theScintillationProcess, idxPostStep);
    }
    if (particleName == "opticalphoton") {
      G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl;
      pmanager->AddDiscreteProcess(theAbsorptionProcess);
      pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
      pmanager->AddDiscreteProcess(theMieHGScatteringProcess);
      pmanager->AddDiscreteProcess(theBoundaryProcess);
    }
  }
  G4cout << "### PhysicsList::ConstructProcess OpticalPhyscs added" << G4endl; 
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SelectPhysicsList(const G4String& name)
{
  G4cout << "PhysicsList::SelectPhysicsList: <" << name << ">" << G4endl;

  // default  Had physics
  if (name == "Hadron" && !fHadPhysicsList) {
	fHadPhysicsList = new PhysListHadron("hadron");
	aHadPhysicsList = static_cast<PhysListHadron*> (fHadPhysicsList);
	aHadPhysicsList->ActivateNeutrons(ActiveNeutrons);
  } else if (name == "QGSP_BERT") {
    AddExtraBuilders(false);
    fHadPhysicsList = new G4HadronPhysicsQGSP_BERT("std-hadron");
  } else if (name == "QGSP_BIC" && !fHadPhysicsList) {
    AddExtraBuilders(false);
    fHadPhysicsList = new G4HadronPhysicsQGSP_BIC("std-hadron");
  } else if (name == "QGSP_BERT_HP"  && !fHadPhysicsList) {
    AddExtraBuilders(true);
    fHadPhysicsList = new G4HadronPhysicsQGSP_BERT_HP("std-hadron");
  } else if (name == "QGSP_BIC_HP"  && !fHadPhysicsList) {
    AddExtraBuilders(true);
    fHadPhysicsList = new G4HadronPhysicsQGSP_BIC_HP("std-hadron");
  } else if (name == "LowEnergy_EM") {
      delete fEmPhysicsList;
      fEmPhysicsList = new PhysListEmLowEnergy("lowe-em");
  } else if (name == "Standard_EM") {
      delete fEmPhysicsList;
      fEmPhysicsList = new G4EmStandardPhysics();
  } else if (name == "Livermore_EM") {
      delete fEmPhysicsList;
      G4EmLivermorePhysics *fEmLivermore;
      fEmPhysicsList = fEmLivermore = new G4EmLivermorePhysics(aIonEner,aIonStep);
  } else if (name == "Penelope_EM") {
      delete fEmPhysicsList;
      fEmPhysicsList = new G4EmPenelopePhysics();
  } else {
      G4cout << "PhysicsList WARNING wrong or unkonwn <" 
	     << name << "> Physics " << G4endl;
	  return;
  }

  G4EmProcessOptions emOptions;
  emOptions.SetFluo(true); // To activate deexcitation processes and fluorescence
  emOptions.SetAuger(true); // To activate Auger effect if deexcitation is activated
  emOptions.SetPIXE(true); // To activate Particle Induced X-Ray Emission (PIXE) 
  G4cout << "Fluorescence options activated!" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddExtraBuilders(G4bool flagHP)
{
  fNhadcomp = 5;
  fHadronPhys.push_back( new G4EmExtraPhysics("extra EM"));
  fHadronPhys.push_back( new G4HadronElasticPhysics(verboseLevel));
  fHadronPhys.push_back( new G4StoppingPhysics("stopping",verboseLevel));
  fHadronPhys.push_back( new G4IonBinaryCascadePhysics("ionBIC"));
  fHadronPhys.push_back( new G4NeutronTrackingCut("Neutron tracking cut"));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
  SetCutValue(fCutForGamma, "gamma");
  SetCutValue(fCutForElectron, "e-");
  SetCutValue(fCutForPositron, "e+");
  G4cout << "world cuts are set" << G4endl;

  if( !fTargetCuts ) SetTargetCut(fCutForElectron);
//  G4Region* region = (G4RegionStore::GetInstance())->GetRegion("Target");
//  region->SetProductionCuts(fTargetCuts);
//  G4cout << "Target cuts are set" << G4endl;

  if( !fDetectorCuts ) SetDetectorCut(fCutForElectron);
//  region = (G4RegionStore::GetInstance())->GetRegion("Detector");
//  region->SetProductionCuts(fDetectorCuts);
//  G4cout << "Detector cuts are set" << G4endl;

  SetCutValue(fCutForProton, "proton");
  SetCutValue(fCutForAlpha, "alpha");
  SetCutValue(fCutForGenericIon, "GenericIon");

  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForGamma(G4double cut)
{
  fCutForGamma = cut;
  SetParticleCuts(fCutForGamma, G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForElectron(G4double cut)
{
  fCutForElectron = cut;
  SetParticleCuts(fCutForElectron, G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForPositron(G4double cut)
{
  fCutForPositron = cut;
  SetParticleCuts(fCutForPositron, G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForProton(G4double cut)
{
  fCutForProton = cut;
  SetParticleCuts(fCutForProton, G4Proton::Proton());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForAlpha(G4double cut)
{
  fCutForAlpha = cut;
  SetParticleCuts(fCutForAlpha, G4Alpha::Alpha());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCutForGenericIon(G4double cut)
{
  fCutForGenericIon = cut;
  SetParticleCuts(fCutForGenericIon, G4GenericIon::GenericIon());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetTargetCut(G4double cut)
{
  if( !fTargetCuts ) fTargetCuts = new G4ProductionCuts();

  fTargetCuts->SetProductionCut(cut, idxG4GammaCut);
  fTargetCuts->SetProductionCut(cut, idxG4ElectronCut);
  fTargetCuts->SetProductionCut(cut, idxG4PositronCut);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetDetectorCut(G4double cut)
{
  if( !fDetectorCuts ) fDetectorCuts = new G4ProductionCuts();

  fDetectorCuts->SetProductionCut(cut, idxG4GammaCut);
  fDetectorCuts->SetProductionCut(cut, idxG4ElectronCut);
  fDetectorCuts->SetProductionCut(cut, idxG4PositronCut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Optical Processes ////////////////////////////////////////////////////////
