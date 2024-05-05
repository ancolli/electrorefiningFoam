/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is not part of OpenFOAM.

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

\*---------------------------------------------------------------------------*/

#include "wallFlux.H"

#include "surfaceInterpolate.H"
#include "fvcSnGrad.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wallFlux, 0);
    addToRunTimeSelectionTable(functionObject, wallFlux, dictionary);
}
}

void Foam::functionObjects::wallFlux::writeFileHeader(Ostream& os) const
{
    // Add headers to output data
    writeHeader(os, "Wall flux");
    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "integral");
    os << endl;
}

void Foam::functionObjects::wallFlux::calcFlux 
(
    const volScalarField& C_red_, 
    volScalarField& wallFlux
)

{
    
    surfaceScalarField flux
    (
      nu_e_/(1-t_)*F_*Dmol_*fvc::snGrad(C_red_) 
    );

    volScalarField::Boundary& wallFluxBf =
        wallFlux.boundaryFieldRef();

    const surfaceScalarField::Boundary& fluxBf =
        flux.boundaryField();

    forAll(wallFluxBf, patchi)
    {
        wallFluxBf[patchi] = fluxBf[patchi];
    }
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::functionObjects::wallFlux::wallFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),

    writeFile(obr_, name, typeName, dict),
    patchSet_()
 
   /* C_red_
    (
       IOobject
       (
          "rho_iMean",
           mesh_.time().timeName(),
           mesh_,
           IOobject::MUST_READ,
           IOobject::NO_WRITE
       ),
       mesh_
      
    )*/
      
{
    volScalarField* wallFluxPtr
    (
	new volScalarField
        (
            IOobject
            (
                type(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("0", dimMoles/dimVolume/dimTime*dimLength, 0.0)
        )
    );
    
    mesh_.objectRegistry::store(wallFluxPtr);
    
    read(dict);
    
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::wallFlux::~wallFlux()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wallFlux::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);//writeLocalObjects::read(dict);
     
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            wordReList(dict.lookupOrDefault("patches", wordReList()))
        );

    Info<< type() << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            label patchi = iter.key();
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall flux on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }
   
    IOdictionary transportProperties
    (
    	IOobject
    	(
        	"transportProperties",
		mesh_.time().constant(),//mesh_.time().timeName(),
        	mesh_,
        	IOobject::MUST_READ,
        	IOobject::NO_WRITE
   	)
   );

    transportProperties.lookup("nu_e") >> nu_e_;
    transportProperties.lookup("F") >> F_;
    transportProperties.lookup("Dmol") >> Dmol_;
    transportProperties.lookup("t") >> t_;

    return true;
}


bool Foam::functionObjects::wallFlux::execute()
{
    volScalarField& wallFlux = lookupObjectRef<volScalarField>(type());
    
    volScalarField C_red_
    (
       IOobject
       (
          "rho_iMean",
           mesh_.time().timeName(),
           mesh_,
           IOobject::MUST_READ,
           IOobject::NO_WRITE
       ),
       mesh_
      
    );
    
    calcFlux(C_red_, wallFlux);

    return true;
}

bool Foam::functionObjects::wallFlux::write()
{

    const volScalarField& wallFlux = lookupObject<volScalarField>(type());
    
    Log << type() << " " << name() << " write:" << nl
    	<< "    writing field " << wallFlux.name() << endl;

    wallFlux.write();//writeLocalObjects::write();

    //logFiles::write();

    //const volScalarField& wallFlux =
       // obr_.lookupObject<volScalarField>(type());

    const fvPatchList& patches = mesh_.boundary();

    const surfaceScalarField::Boundary& magSf =
        mesh_.magSf().boundaryField();

    for (const label patchi : patchSet_) //forAllConstIter(labelHashSet, patchSet_, iter)
    {
        //label patchi = iter.key();
        const fvPatch& pp = patches[patchi];

        const scalarField& hfp = wallFlux.boundaryField()[patchi];

        const scalar minHfp = gMin(hfp);
        const scalar maxHfp = gMax(hfp);
        const scalar integralHfp = gSum(magSf[patchi]*hfp)/gSum(magSf[patchi]);
	const scalar current = gSum((magSf[patchi]*hfp)); //absolute value

        if (Pstream::master())
        {
            file()
                << mesh_.time().value()
                << tab << pp.name()
                << tab << minHfp
                << tab << maxHfp
                << tab << integralHfp
		<< tab << current
                << endl;
        }

        Log << "    min/max/average/current(" << pp.name() << ") = "
            << minHfp << ", " << maxHfp << ", " << integralHfp << ", " << current << endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
