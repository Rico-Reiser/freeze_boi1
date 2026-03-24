MRF.correctBoundaryVelocity(U);

// --- 1) Sicherheit
fs.max(0.0);
fs.min(1.0);

// --- 2) mushy - zone damping
dimensionedScalar Amushy
(
    "Amushy",
    dimensionSet(0, 0, -1, 0, 0, 0, 0),  // 1/s
    1e3
);

volScalarField DarcyCoeff
(
    "DarcyCoeff",
    Amushy * sqr(fs) / (pow3(1.0 - fs) + scalar(1e-3))
);
   // DEBUG fs und Darcy-damping
   Info<< "fs[min,max]=(" << gMin(fs) << "," << gMax(fs)
    << ") Darcy[max]=" << gMax(DarcyCoeff) << endl;

// --- 3) Kopplung Dary Term an alpha.water
if (mesh.foundObject<volScalarField>("alpha.water"))
{
    const volScalarField& alphaW =
        mesh.lookupObject<volScalarField>("alpha.water");

    DarcyCoeff *= alphaW;                                      
}
   DarcyCoeff *= pos(fs - 0.05);
   DarcyCoeff.min(1e6);


// --- 5) Momentum equation
fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevTau(rho, U)
  + fvm::Sp(rho*DarcyCoeff, U)                                          
 ==
    fvModels.source(rho, U)
);

UEqn.relax();
fvConstraints.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
                mixture.surfaceTensionForce()
              - ghf*fvc::snGrad(rho)
              - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
    );

    fvConstraints.constrain(U);
}
