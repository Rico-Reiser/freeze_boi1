/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "codedFixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = d4b09f53d9a61667a3e655c8d3f8c05d49c41f07
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void bottomH_d4b09f53d9a61667a3e655c8d3f8c05d49c41f07(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchScalarField,
    bottomHFixedValueFvPatchScalarField
);


const char* const bottomHFixedValueFvPatchScalarField::SHA1sum =
    "d4b09f53d9a61667a3e655c8d3f8c05d49c41f07";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

bottomHFixedValueFvPatchScalarField::
bottomHFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct bottomH sha1: d4b09f53d9a61667a3e655c8d3f8c05d49c41f07"
            " from patch/DimensionedField\n";
    }
}


bottomHFixedValueFvPatchScalarField::
bottomHFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct bottomH sha1: d4b09f53d9a61667a3e655c8d3f8c05d49c41f07"
            " from patch/dictionary\n";
    }
}


bottomHFixedValueFvPatchScalarField::
bottomHFixedValueFvPatchScalarField
(
    const bottomHFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct bottomH sha1: d4b09f53d9a61667a3e655c8d3f8c05d49c41f07"
            " from patch/DimensionedField/mapper\n";
    }
}


bottomHFixedValueFvPatchScalarField::
bottomHFixedValueFvPatchScalarField
(
    const bottomHFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct bottomH sha1: d4b09f53d9a61667a3e655c8d3f8c05d49c41f07 "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

bottomHFixedValueFvPatchScalarField::
~bottomHFixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy bottomH sha1: d4b09f53d9a61667a3e655c8d3f8c05d49c41f07\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void bottomHFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs bottomH sha1: d4b09f53d9a61667a3e655c8d3f8c05d49c41f07\n";
    }

//{{{ begin code
    #line 23 "/home/rico-reiser/cases_OF10/new_cases/freeze4/0/h.contactFluid/boundaryField/bottom"
const fvPatch& patch = this->patch();
            const fvMesh&  mesh  = patch.boundaryMesh().mesh();

            const scalar Twall = 200.0;     // [K]
            const scalar TrefK = 273.65;    // MUSS zu Solver passen

            scalarField& hPatch = *this;

            if (mesh.foundObject<volScalarField>("cpMix"))
            {
                const volScalarField& cpMix =
                    mesh.lookupObject<volScalarField>("cpMix");

                forAll(hPatch, faceI)
                {
                    const label cellI = patch.faceCells()[faceI];
                    const scalar cpI  = max(cpMix[cellI], VSMALL);
                    hPatch[faceI] = cpI*(Twall - TrefK);
                }
            }
            else
            {
                IOdictionary solidProps
                (
                    IOobject
                    (
                        "solidificationProperties",
                        "constant",
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );

                const dictionary& contDict =
                    solidProps.subDict("contactFluid");
                const scalar cpC = readScalar(contDict.lookup("cp"));

                forAll(hPatch, faceI)
                {
                    hPatch[faceI] = cpC*(Twall - TrefK);
                }
            }
//}}} end code

    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

