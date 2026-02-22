/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    JouleThomsonInitialTemperature

Group
    grpHeatTransferSolvers

Description
    Calculates initial temperature for flow of He II in a pipe
    due to Joule-Thomshon effect.
    It is assumed that the mesh is 1D and in the positive x direction
    also the inlet surface should be in plane y-z so it starts at xInlet

    Uses equation:
    \f[
        \frac{dT}{dx} = -\frac{1}{\rho C_p} \frac{dP}{dx} = \frac{1}{\rho C_p} \frac{\Delta P}{L}
    \f]

    where:
TODO: 22.02.2026
        \f$ rho_{k} \f$ = the effective (driving) kinematic density
        beta = thermal expansion coefficient [1/K]
        T = temperature [K]
        \f$ T_{ref} \f$ = reference temperature [K]


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseHeliumTransportModel.H"
#include "turbulentTransportModel.H"
#include "radiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Calculates initial temperature for flow of He II in a pipe\n"
        " due to Joule-Thomshon effect."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting calculation of temperature due to"
        << " Joule-Thomson effect." << endl;

    // approximation based on Fig. 3.4 from Fuzier's PhD
    // v is velocity in m/s
    // function returns pressure drop in Pa
    auto calcDp = [](const scalar v)
    {
    	const scalar dp{(-0.08369 + 0.02645*v + 0.08533*v*v)*1000};
    	return dimensionedScalar{dimPressure, dp};
    };
    
    const dimensionedScalar U_in{"U_in", dimVelocity, laminarTransport};
    const dimensionedScalar T_in{"T_in", dimTemperature, laminarTransport};
    const dimensionedScalar rho_in{laminarTransport.heThermProp(T_in,"rho")};
    const dimensionedScalar L{"L", dimLength, laminarTransport};
    dimensionedScalar dP{"dP", dimPressure, calcDp(U_in.value()).value()};
    Info<< "Calculated experimental pressure drop: " << dP << endl;
    
    // it is assumed that the mesh is 1D and in the positive x direction
    // also the inlet surface should be in plane y-z so it starts at xInlet
    Info<< "Initial temperature is calculated." << endl;
    volScalarField x
    (
        IOobject
        (
            "x",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh.C().component(vector::X)  // or mesh.C().x()
    );

    word inletPatch = "INLET";

    label inletPatchID = mesh.boundaryMesh().findPatchID(inletPatch);

    if (inletPatchID == -1)
    {
        FatalErrorInFunction
            << "Patch " << inletPatch << " not found!"
            << exit(FatalError);
    }

    word outletPatch = "OUTLET";

    label outletPatchID = mesh.boundaryMesh().findPatchID(outletPatch);

    if (outletPatchID == -1)
    {
        FatalErrorInFunction
            << "Patch " << outletPatch << " not found!"
            << exit(FatalError);
    }

    // check if inlet lays in YZ plane
    const vectorField& Cf =
        mesh.Cf().boundaryField()[inletPatchID];  // face centres of that patch
    scalar xmin = gMin(Cf.component(vector::X));
    scalar xmax = gMax(Cf.component(vector::X));
    
    if (mag(xmax - xmin) > 1e-6)
    {
    	FatalErrorInFunction
    		<< "Probably inlet do not lie in y-z plane." << abort(FatalError);
    }
    
    scalar xInlet = gAverage(Cf.component(vector::X));
    //Info << "Inlet x-coordinate = " << xInlet << endl;

    const dimensionedScalar cp_in{laminarTransport.heThermProp(T_in,"cp")};
    //Info<< "\nCalculated specific heat capacity at constant pressure SVP" 
    //    << " at temperature of " << T_in.value() << " K is " << cp_in << endl;
    //Info<< "\nCalculated density at constant pressure SVP" 
    //    << " at temperature of " << T_in.value() << " K is " << rho_in << endl;
    
    forAll (T, celli)
    {
    	T[celli] = T_in.value() 
    		+ (dP/rho_in/cp_in/L).value()*(x[celli]-xInlet);
    }
    
    // analytical outlet temperature due to JT effect
    const scalar T_out{(T_in+dP/rho_in/cp_in).value()};
    Info<< "\nAnalytical Tout due to Joule-Thomson effect is: " << T_out << " K" << endl;

    // set temperature at inlet to T_in
    T.boundaryFieldRef()[inletPatchID] ==
    scalarField
    (
        T.boundaryField()[outletPatchID].size(),
        T_in.value()
    );

    // set temperature at outlet to T_out
    T.boundaryFieldRef()[outletPatchID] ==
    scalarField
    (
        T.boundaryField()[outletPatchID].size(),
        T_out
    );

    T.correctBoundaryConditions();

    T.write();

    Info<< "Initial temperature has been applied.\n" << endl;
    Info<< "T boundary conditions have been modified and wrote. \n" << endl;

    T.write();
    

    return 0;
}


// ************************************************************************* //
