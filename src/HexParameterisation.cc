//
// Created by mogu on 3/31/21.
//
#include "HexParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HexParameterisation::HexParameterisation(
        G4int ncells,
        G4int ncolumns,          //  Z of center of first
        G4int nrow,       //  Z spacing of centers
        G4double first_position_X,
        G4double first_position_Y,
        G4double first_position_Z,
        G4double pore_pitch)
        : G4VPVParameterisation()
{
    fncells = ncells;
    fncolumns = ncolumns;
    fnrow = nrow;
    ffirst_position_X = first_position_X;
    ffirst_position_Y = first_position_Y;
    position_Z = first_position_Z;
    //fcomb_radius = pore_pitch / 2;
    x_spacing = pore_pitch;
    y_spacing = pore_pitch * sqrt(3) /2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HexParameterisation::~HexParameterisation()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HexParameterisation::ComputeTransformation
        (const G4int copyNo, G4VPhysicalVolume* physVol) const
{
    // Note: copyNo will start with zero!spacing_X
    /*G4int position_column, position_row;*/
    G4int i = 0;
    G4int j = 0;
    G4double xoffset = 0;
    i = G4int(copyNo / fncolumns);
    j = G4int(copyNo % fnrow);
    if(i % 2 == 1){
        xoffset = 0.5 * x_spacing;
    }
    else{
        xoffset = 0;
    }
    G4double position_X = ffirst_position_X + j * x_spacing + xoffset;
    G4double position_Y = ffirst_position_Y + i * y_spacing;
    /*if((copyNo % (2 * fnrow -1)) < fnrow)
    {
        position_column = copyNo / (2 * fnrow - 1) * 2;
        position_row = copyNo % (2 * fnrow -1) * 2;
    }
    else
    {
        position_column = copyNo / (2 * fnrow - 1) * 2 + 1;
        position_row = (copyNo % (2 * fnrow -1) - fnrow) * 2 + 1;
    }*/


    /*G4double position_X = ffirst_position_X - fcomb_radius * sqrt(3) * position_column;
    G4double position_Y = ffirst_position_Y - fcomb_radius * position_row;*/
    //G4cout << copyNo << "," << position_column << "," << position_X << "; " << position_row << "," << position_Y << G4endl;
    G4ThreeVector origin(position_X, position_Y, position_Z);
    physVol->SetTranslation(origin);
    physVol->SetRotation(0);
}
