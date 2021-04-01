#include "G4ParticleGun.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "TH1D.h"
#include "TF1.h"

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  PrimaryGeneratorAction(G4String source_name);
  virtual ~PrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event*);

  virtual void SetParticle(G4double &x,G4double &y,G4double &z,G4double &px,G4double &py,G4double &pz,G4double &e,TH1D* h);//MLCC
  virtual void SetParticle(G4double &x,G4double &y,G4double &z,G4double &px,G4double &py,G4double &pz,G4double &e,TF1* f);//MLCC
  const G4ParticleGun* GetParticleGun() const {return particleGun;}

private:
  G4ParticleGun* particleGun;
  G4String sname;
};
