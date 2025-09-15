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

//GEANT4 - Depth-of-Interaction enabled Positron emission tomography (PET) advanced example 

//Authors and contributors

// Author list to be updated, with names of co-authors and contributors from National Institute of Radiological Sciences (NIRS)

// Abdella M. Ahmed (1, 2), Andrew Chacon (1, 2), Harley Rutherford (1, 2),
// Hideaki Tashima (3), Go Akamatsu (3), Akram Mohammadi (3), Eiji Yoshida (3), Taiga Yamaya (3)
// Susanna Guatelli (2), and Mitra Safavi-Naeini (1, 2)

// (1) Australian Nuclear Science and Technology Organisation, Australia
// (2) University of Wollongong, Australia
// (3) National Institute of Radiological Sciences, Japan


#include "doiPETPrimaryGeneratorAction.hh"
#include "DicomFileMgr.hh"
#include "DicomFilePET.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ChargedGeantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4GeneralParticleSource.hh"



doiPETPrimaryGeneratorAction::doiPETPrimaryGeneratorAction()
	: G4VUserPrimaryGeneratorAction(),
	fParticleGun(0),
	fDicomMgr(nullptr),
	fUsePETActivity(false),
	fScanTime(1800.0), // Default 10 minutes scan
	fTotalActivity(0.0),
	fTotalEmissions(0)
{
	// Definition of the General particle Source
	fParticleGun  = new G4GeneralParticleSource();
	
	// Initialize PET source if activity data is available
	InitializePETSource();
}

doiPETPrimaryGeneratorAction::~doiPETPrimaryGeneratorAction()
{
	delete fParticleGun;
}


void doiPETPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{           
	// If PET activity data is available, sample position from it
	if (fUsePETActivity) {
		G4ThreeVector position = SampleFromPETActivity();
		
		G4cout << "DEBUG: Sampled position from PET activity: (" 
		       << position.x()/mm << ", " << position.y()/mm << ", " << position.z()/mm << ") mm" << G4endl;
		
		// Configure the GPS for F-18 positron emission at sampled position
		fParticleGun->GetCurrentSource()->GetPosDist()->SetPosDisType("Point");
		fParticleGun->GetCurrentSource()->GetPosDist()->SetCentreCoords(position);
	} else {
		G4cout << "DEBUG: Using default GPS configuration (no PET activity data)" << G4endl;
	}
	
	//create vertex
	fParticleGun->GeneratePrimaryVertex(anEvent);
}

void doiPETPrimaryGeneratorAction::InitializePETSource()
{
	// Get the DICOM file manager
	fDicomMgr = DicomFileMgr::GetInstance();
	
	G4cout << "=== DEBUG: Checking PET Activity Data ===" << G4endl;
	G4cout << "DicomFileMgr instance: " << fDicomMgr << G4endl;
	
	if (fDicomMgr) {
		G4cout << "PETFileAll pointer: " << fDicomMgr->GetPETFileAll() << G4endl;
		if (fDicomMgr->GetPETFileAll()) {
			G4cout << "PET file found!" << G4endl;
		} else {
			G4cout << "No PET file found!" << G4endl;
		}
	}
	
	if (fDicomMgr && fDicomMgr->GetPETFileAll()) {
		fUsePETActivity = true;
		
		G4cout << "=== PET Activity Source Configuration ===" << G4endl;
		G4cout << "Found PET activity data - enabling activity-guided particle generation" << G4endl;
		
		// Configure GPS for F-18 ions (for positron emission)
		G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
		G4IonTable* ionTable = particleTable->GetIonTable();
		G4ParticleDefinition* f18 = ionTable->GetIon(9, 18, 0.0); // Z=9, A=18, E=0
		
		if (f18) {
			fParticleGun->GetCurrentSource()->SetParticleDefinition(f18);
			fParticleGun->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
			fParticleGun->GetCurrentSource()->GetEneDist()->SetMonoEnergy(0.0 * keV);
			
			G4cout << "Configured for F-18 ion emission (will decay to positrons)" << G4endl;
		} else {
			G4cout << "Warning: Could not find F-18 ion definition, fallback to e+" << G4endl;
			G4ParticleDefinition* positron = particleTable->FindParticle("e+");
			if (positron) {
				fParticleGun->GetCurrentSource()->SetParticleDefinition(positron);
				fParticleGun->GetCurrentSource()->GetEneDist()->SetEnergyDisType("Mono");
				fParticleGun->GetCurrentSource()->GetEneDist()->SetMonoEnergy(0.0 * keV);
			}
		}
		
		G4cout << "PET phantom dimensions: " 
		       << fDicomMgr->GetPETFileAll()->GetNoVoxelsX() << " x "
		       << fDicomMgr->GetPETFileAll()->GetNoVoxelsY() << " x "
		       << fDicomMgr->GetPETFileAll()->GetNoVoxelsZ() << " voxels" << G4endl;
		       
		// Calculate total emissions for the scan
		CalculateTotalEmissions();
		       
		G4cout << "==========================================" << G4endl;
		
	} else {
		fUsePETActivity = false;
		G4cout << "No PET activity data found - using default GPS configuration" << G4endl;
	}
}

G4ThreeVector doiPETPrimaryGeneratorAction::SampleFromPETActivity()
{
	if (!fDicomMgr || !fDicomMgr->GetPETFileAll()) {
		G4cout << "Error: No PET data available for sampling!" << G4endl;
		return G4ThreeVector(0, 0, 0);
	}
	
	DicomFilePET* petFile = fDicomMgr->GetPETFileAll();
	
	// Get phantom dimensions
	G4int nX = petFile->GetNoVoxelsX();
	G4int nY = petFile->GetNoVoxelsY(); 
	G4int nZ = petFile->GetNoVoxelsZ();
	
	G4double minX = petFile->GetMinX();
	G4double maxX = petFile->GetMaxX();
	G4double minY = petFile->GetMinY();
	G4double maxY = petFile->GetMaxY();
	G4double minZ = petFile->GetMinZ();
	G4double maxZ = petFile->GetMaxZ();
	
	G4double voxelSizeX = (maxX - minX) / nX;
	G4double voxelSizeY = (maxY - minY) / nY;
	G4double voxelSizeZ = (maxZ - minZ) / nZ;
	
	// Get activity data  
	const std::vector<int>& activities = petFile->GetHounsfieldV();
	
	// Find maximum activity for normalization (calculate once and store)
	static G4double maxActivity = 0.0;
	static G4bool maxActivityCalculated = false;
	if (!maxActivityCalculated) {
		for (size_t i = 0; i < activities.size(); i++) {
			if (activities[i] > maxActivity) {
				maxActivity = activities[i];
			}
		}
		maxActivityCalculated = true;
		G4cout << "INFO: Maximum activity value in PET data: " << maxActivity << " (for normalization)" << G4endl;
	}
	
	// Sample voxel weighted by activity (rejection sampling)
	G4int maxTries = 10000;
	for (G4int iTry = 0; iTry < maxTries; iTry++) {
		// Random voxel indices
		G4int iX = G4int(G4UniformRand() * nX);
		G4int iY = G4int(G4UniformRand() * nY);
		G4int iZ = G4int(G4UniformRand() * nZ);
		
		// Get activity at this voxel
		G4int voxelIndex = iX + iY * nX + iZ * nX * nY;
		
		if (voxelIndex < activities.size()) {
			G4double activity = activities[voxelIndex];
			
			// Accept with probability proportional to activity (normalized)
			if (activity > 0 && G4UniformRand() < (activity / maxActivity)) {
				// Convert voxel indices to world coordinates
				G4double x = minX + (iX + 0.5) * voxelSizeX;
				G4double y = minY + (iY + 0.5) * voxelSizeY;
				G4double z = minZ + (iZ + 0.5) * voxelSizeZ;
				
				// Debug: Print voxel info occasionally  
				static G4int sampleCount = 0;
				sampleCount++;
				if (sampleCount <= 10 || sampleCount % 100 == 0) {
					G4cout << "SPATIAL DEBUG [Sample #" << sampleCount 
					       << "]: Voxel(" << iX << "," << iY << "," << iZ 
					       << ") -> Position(" << x/mm << "," << y/mm << "," << z/mm 
					       << ")mm, Activity=" << activity 
					       << " (norm=" << activity/maxActivity << ")" << G4endl;
				}
				
				return G4ThreeVector(x, y, z);
			}
		}
	}
	
	// Fallback to center if no activity found
	G4cout << "Warning: No activity found after " << maxTries << " tries, using phantom center" << G4endl;
	return G4ThreeVector((minX + maxX)/2, (minY + maxY)/2, (minZ + maxZ)/2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void doiPETPrimaryGeneratorAction::CalculateTotalEmissions()
{
	if (!fDicomMgr || !fDicomMgr->GetPETFileAll()) {
		return;
	}
	
	DicomFilePET* petFile = fDicomMgr->GetPETFileAll();
	const std::vector<int>& activities = petFile->GetHounsfieldV();
	
	// Calculate total activity from all voxels
	fTotalActivity = 0.0;
	for (size_t i = 0; i < activities.size(); i++) {
		fTotalActivity += activities[i]; // Sum all voxel activities
	}
	
	// Convert to Bq (assuming DICOM values are in appropriate units)
	// You may need to adjust this conversion factor based on your DICOM data scaling
	fTotalActivity *= becquerel; 
	
	// Calculate total expected emissions during scan time
	// For F-18: each decay produces one positron -> two 511 keV gamma rays
	// fTotalActivity is in Bq (decays/s), fScanTime in seconds
	fTotalEmissions = static_cast<G4long>(fTotalActivity * fScanTime / becquerel);
	
	G4cout << "=== SCAN TIME CALCULATIONS ===" << G4endl;
	G4cout << "Scan time: " << fScanTime << " seconds (" << fScanTime/60.0 << " minutes)" << G4endl;
	G4cout << "Total activity in phantom: " << fTotalActivity/(1.e+6*becquerel) << " MBq" << G4endl;
	G4cout << "Expected emissions during scan: " << fTotalEmissions << " decays" << G4endl;
	G4cout << "Expected annihilation photons: " << 2*fTotalEmissions << " gamma rays (511 keV)" << G4endl;
	
	if (fTotalEmissions > 1e9) {
		G4cout << "INFO: Large number of expected emissions. Consider using:" << G4endl;
		G4cout << "  - Shorter scan time (current: " << fScanTime << "s)" << G4endl;
		G4cout << "  - Statistical scaling in post-processing" << G4endl;
		G4cout << "  - Activity concentration reduction" << G4endl;
	}
	
	G4cout << "===============================" << G4endl;
}


