//
// Created by mogu on 3/31/21.
//

#ifndef B1_HEXPARAMETERISATION_HH
#define B1_HEXPARAMETERISATION_HH

#include "G4VPVParameterisation.hh"

class G4VPhysicalVolume;
class G4Box;

// Dummy declarations to get rid of warnings ...
class G4Trd;
class G4Trap;
class G4Cons;
class G4Orb;
class G4Sphere;
class G4Ellipsoid;
class G4Torus;
class G4Para;
class G4Hype;
class G4Tubs;
class G4Polycone;
class G4Polyhedra;

///  A parameterisation that describes a series of boxes along Z.
///
///  The boxes have equal width, & their lengths are a linear equation.
///  They are spaced an equal distance apart, starting from given location.

class HexParameterisation : public G4VPVParameterisation
{
public:

    HexParameterisation(G4int ncells,
                        G4int ncolumns,          //  Z of center of first
                        G4int nrow,       //  Z spacing of centers
                        G4double first_position_X,
                        G4double first_position_Y,
                        G4double first_position_Z,
                        G4double pore_pitch);

    virtual ~HexParameterisation();

    void ComputeTransformation (const G4int copyNo,
                                G4VPhysicalVolume* physVol) const;

private:
    G4int fncells;
    G4int fncolumns;    //  Z of center of first
    G4int fnrow;       //  Z spacing of centers
    G4double ffirst_position_X;
    G4double ffirst_position_Y;
    G4double position_Z;
    //G4double fcomb_radius;
    G4double x_spacing, y_spacing;

};



#endif //B1_HEXPARAMETERISATION_HH
