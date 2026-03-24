// ============================================================================
// createFields.H  (Option C: getrennte Enthalpien pro Phase)
// - ersetzt Mischfeld h durch: h.air, h.water, h.contactFluid
// - erzeugt dazu T.air, T.water, T.contactFluid
// - cpMix/kEff/kW_eff bleiben für Debug/Utilities ok (kannst du später optional entfernen)
// - fs (Eis) wird NUR aus T.water initialisiert
// ============================================================================

Info<< "Reading field p_rgh\n" << endl;

volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

// --- multiphase mixture
multiphaseMixture mixture(U, phi);

// --- phase volume fractions
const volScalarField& alphaW =
    mesh.lookupObject<volScalarField>("alpha.water");

const volScalarField& alphaC =
    mesh.lookupObject<volScalarField>("alpha.contactFluid");

const volScalarField& alphaA =
    mesh.lookupObject<volScalarField>("alpha.air");



// --- mixture density (still useful for UEqn/pEqn)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mixture.rho()
);
rho.oldTime();

// --- turbulence
autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, mixture)
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

// --- pressure p from p_rgh
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

pressureReference pressureReference(p, p_rgh, pimple.dict());

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pressureReference.refValue()
      - getRefCellValue(p, pressureReference.refCell())
    );

    p_rgh = p - rho*gh;
}

mesh.schemes().setFluxRequired(p_rgh.name());

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"

// ============================================================================
// Thermophysical + solidification properties
// ============================================================================

Info<< "Reading solidificationProperties\n" << endl;

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

dimensionedScalar Tsolid
(
    "Tsolid",
    dimTemperature,
    solidProps.lookup("Tsolid")
);

dimensionedScalar Tliquid
(
    "Tliquid",
    dimTemperature,
    solidProps.lookup("Tliquid")
);

dimensionedScalar L
(
    "L",
    dimEnergy/dimMass,
    solidProps.lookup("latentHeat")
);

const dictionary& waterDict = solidProps.subDict("water");
const dictionary& contDict  = solidProps.subDict("contactFluid");
const dictionary& airDict   = solidProps.subDict("air");

// densities
// densities water liquid - solid
dimensionedScalar rhoWLiquid("rhoWLiquid", dimMass/dimVolume, waterDict.lookup("rhoLiquid"));
dimensionedScalar rhoWSolid("rhoWSolid", dimMass/dimVolume, waterDict.lookup("rhoSolid"));


dimensionedScalar rhoC("rhoC", dimMass/dimVolume, contDict.lookup("rho"));
dimensionedScalar rhoA("rhoA", dimMass/dimVolume, airDict.lookup("rho"));
// heat capacaties water (liquid - solid)
dimensionedScalar cpWSolid
(
    "cpWSolid",
    dimEnergy/dimMass/dimTemperature,
    waterDict.lookup("cpSolid")
);
dimensionedScalar cpWLiquid
(
    "cpWLiquid",
    dimEnergy/dimMass/dimTemperature,
    waterDict.lookup("cpLiquid")
);
// heat capacity contactFluid
dimensionedScalar cpC("cpC", dimEnergy/dimMass/dimTemperature, contDict.lookup("cp"));

// heat capacity air
dimensionedScalar cpA("cpA", dimEnergy/dimMass/dimTemperature, airDict.lookup("cp"));

// conductivities
dimensionedScalar lambdaW("lambdaW", dimEnergy/dimTime/dimLength/dimTemperature, waterDict.lookup("lambda"));
dimensionedScalar lambdaC("lambdaC", dimEnergy/dimTime/dimLength/dimTemperature, contDict.lookup("lambda"));
dimensionedScalar lambdaA("lambdaA", dimEnergy/dimTime/dimLength/dimTemperature, airDict.lookup("lambda"));
dimensionedScalar lambdaSolid
(
    "lambdaSolid",
    dimEnergy/dimTime/dimLength/dimTemperature,
    waterDict.lookup("lambdaSolid")
);
dimensionedScalar lambdaLiquid
(
    "lambdaLiquid",
    dimEnergy/dimTime/dimLength/dimTemperature,
    waterDict.lookup("lambdaLiquid")
);

// Reference temperature for sensible enthalpy definition
dimensionedScalar Tref
(
    "Tref",
    dimTemperature,
    solidProps.lookupOrDefault<scalar>("Tref", Tliquid.value())
);
// Regularisierungstemperatur
dimensionedScalar Treg
(
    "Treg",
    dimTemperature,
    readScalar(solidProps.lookup("Treg"))
);

// ============================================================================
// Heat transfer coefficients
// ============================================================================
const dictionary& htDict =
    solidProps.found("heatTransfer")
  ? solidProps.subDict("heatTransfer")
  : solidProps;

// heat transfer: water - contactFluid
dimensionedScalar KwcDefault
(
    "KwcDefault",
    dimEnergy/(dimVolume*dimTemperature*dimTime),
    htDict.lookupOrDefault<scalar>("K_wc", 0.0)
);

// heat transfer: water - air
dimensionedScalar KwaDefault
(
    "KwaDefault",
    dimEnergy/(dimVolume*dimTemperature*dimTime),
    htDict.lookupOrDefault<scalar>("K_wa", 0.0)
);

// heat transfer: contactFluid - air
dimensionedScalar KcaDefault
(
    "KcaDefault",
    dimEnergy/(dimVolume*dimTemperature*dimTime),
    htDict.lookupOrDefault<scalar>("K_ca", 0.0)
);

// heat transfer field: water - contactFluid
volScalarField K_wc
(
    IOobject
    (
        "K_wc",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    KwcDefault
);

// heat transfer field: water - air
volScalarField K_wa
(
    IOobject
    (
        "K_wa",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    KwaDefault
);

// heat transfer field: air - contactFluid
volScalarField K_ca
(
    IOobject
    (
        "K_ca",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    KcaDefault
);

K_wc.correctBoundaryConditions();
K_wa.correctBoundaryConditions();
K_ca.correctBoundaryConditions();

// ============================================================================
// Phase enthalpies
// ============================================================================
// enthalpy air
volScalarField hA
(
    IOobject
    (
        "h.air",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("hA0", dimEnergy/dimMass, 0.0)
);
// enthalpy water
volScalarField hW
(
    IOobject
    (
        "h.water",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("hW0", dimEnergy/dimMass, 0.0)
);
// enthalpy contactFluid
volScalarField hC
(
    IOobject
    (
        "h.contactFluid",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("hC0", dimEnergy/dimMass, 0.0)
);

// ============================================================================
// Phase temperatures
// ============================================================================
// Phasentemperaturen (water - contactFluid - air)
volScalarField Tair
(
    IOobject
    (
        "T.air",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    Treg
);

volScalarField Twater
(
    IOobject
    (
        "T.water",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    Treg
);

volScalarField Tcontact
(
    IOobject
    (
        "T.contactFluid",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    Treg
);
// Optional: keep a "mixture T" for legacy debug/flux output (can be removed later)
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    Treg
);

// ============================================================================
// solid fraction fs declaration
// ============================================================================
volScalarField fs
(
    IOobject
    (
        "fs",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);
// Initialise fs from Twater only where water exists
if (runTime.timeIndex() == 0 && !IOobject("fs", runTime.timeName(), mesh).headerOk())
{
    const scalar alphaCut = 1e-3;

    forAll(fs, celli)
    {
        if (alphaW[celli] > alphaCut)
        {
            const scalar Ti = Twater[celli];
            scalar fi;

            if (Ti <= Tsolid.value())
            {
                fi = 1.0;
            }
            else if (Ti >= Tliquid.value())
            {
                fi = 0.0;
            }
            else
            {
                const scalar dT = (Tliquid.value() - Tsolid.value()) + VSMALL;
                const scalar xi = (Ti - Tsolid.value())/dT;
                const scalar x  = min(max(xi, 0.0), 1.0);

                // cosine smoothing in mushy region
                fi = 0.5*(1.0 + Foam::cos(constant::mathematical::pi*x));
            }

            fs[celli] = fi;
        }
        else
        {
            fs[celli] = 0.0;
        }
    }

    fs.correctBoundaryConditions();
}
// ============================================================================
// Effective density of WATER
// fs = 1  -> solid  -> rhoSolid
// fs = 0  -> liquid -> rhoLiquid
// ============================================================================
volScalarField rhoW
(
    IOobject
    (
        "rhoW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fs*rhoWSolid + (scalar(1) - fs)*rhoWLiquid
);

rhoW.correctBoundaryConditions();

rho =
      alphaW*rhoW
    + alphaC*rhoC
    + alphaA*rhoA;

rho.correctBoundaryConditions();

// ============================================================================
// Effective heat capacity of WATER
// fs = 1  -> solid  -> cpWSolid
// fs = 0  -> liquid -> cpWLiquid
// ============================================================================
volScalarField cpW_eff
(
    IOobject
    (
        "cpW_eff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fs*cpWSolid + (scalar(1) - fs)*cpWLiquid
);

cpW_eff.max
(
    dimensionedScalar
    (
        "cpWMin",
        cpW_eff.dimensions(),
        1.0
    )
);

cpW_eff.correctBoundaryConditions();

// ============================================================================
// Initialise phase enthalpies if not provided (time=0)
// Strategy:
// - If h.phase files exist -> keep them
// - else compute from T.phase if present
// - cpW_eff is now already available
// ============================================================================
const bool hAok = IOobject("h.air", runTime.timeName(), mesh).headerOk();
const bool hWok = IOobject("h.water", runTime.timeName(), mesh).headerOk();
const bool hCok = IOobject("h.contactFluid", runTime.timeName(), mesh).headerOk();

if (runTime.timeIndex() == 0)
{
    if (!hAok)
    {
        hA = cpA*(Tair - Tref);
        hA.correctBoundaryConditions();
    }

    if (!hWok)
    {
        hW = cpW_eff*(Twater - Tref);
        hW.correctBoundaryConditions();
    }

    if (!hCok)
    {
        hC = cpC*(Tcontact - Tref);
        hC.correctBoundaryConditions();
    }

    // Recompute temperatures consistently from h (single source of truth)
    Tair     = Tref + hA/cpA;
    Twater   = Tref + hW/cpW_eff;
    Tcontact = Tref + hC/cpC;

    Tair.correctBoundaryConditions();
    Twater.correctBoundaryConditions();
    Tcontact.correctBoundaryConditions();
}


// ============================================================================
// Optional mixture temperature for debug only
// ============================================================================
T =
    alphaA*Tair
  + alphaW*Twater
  + alphaC*Tcontact;
T.correctBoundaryConditions();

// ============================================================================
// Thermodynamic mixture properties
// ============================================================================
volScalarField cpMix
(
    IOobject
    (
        "cpMix",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alphaW*cpW_eff + alphaC*cpC + alphaA*cpA
);
cpMix.max

(
    dimensionedScalar
    (
        "cpMin", 
        cpMix.dimensions(), 
        500.0
    )
);
cpMix.correctBoundaryConditions();

volScalarField lambdaMix
(
    IOobject
    (
        "lambdaMix",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alphaW*lambdaW + alphaC*lambdaC + alphaA*lambdaA
);
lambdaMix.max(dimensionedScalar("lambdaMin", lambdaMix.dimensions(), 1e-12));
lambdaMix.correctBoundaryConditions();

// ============================================================================
// Effective conductivity of WATER in mushy region
// ============================================================================
volScalarField lambdaW_eff
(
    IOobject
    (
        "lambdaW_eff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fs*lambdaSolid + (scalar(1) - fs)*lambdaLiquid
);
lambdaW_eff.max
(
    dimensionedScalar
    (
        "lambdaMinW", 
        lambdaW_eff.dimensions(), 
        1e-12
    )
);
lambdaW_eff.correctBoundaryConditions();

// ============================================================================
// Effective total conductivity
// ============================================================================
volScalarField lambdaEff
(
    IOobject
    (
        "lambdaEff",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alphaW*lambdaW_eff + alphaC*lambdaC + alphaA*lambdaA
);
lambdaEff.max
(
    dimensionedScalar
    (
        "lambdaMinEff", 
        lambdaEff.dimensions(), 
        1e-12
    )
);
lambdaEff.correctBoundaryConditions();
