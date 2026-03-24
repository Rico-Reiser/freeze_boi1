// Solver Konfigurationsdatei: hEqn.H 
// Erzeugung eines Enthalpie Feldes für alle beiteiligten Phasen

// Inhaltsverzeichnis
//
// --- 1)  binary mask for water                                             Zeile:
// --- 2)  Temperatur konsistent aus h                                       Zeile: 
// --- 3)  fs_old explizit sichern                                           Zeile:
// --- 4)  Solid fraction aus T aktualisieren                                Zeile:
// --- 5)  kEff aus fs                                                       Zeile:
// --- 6)  Änderung des Phasenanteils über der Zeit d(fs)/dt berechnen       Zeile:
// --- 7)  Definition: Mushy mask                                            Zeile:
// --- 8)  Definition: Latent heat source (water only)                       Zeile:
// --- 9)  mass flux for mixture enthalpy transport (alias of rhoPhi)        Zeile:
// --- 10) Stützterm gegen degenerierte Zellen (alphaW=0)                    Zeile:
// --- 11) Mixture sensible enthalpy equation                                Zeile:
// --- 12) Debug des Heat-Flux                                               Zeile:

// --- 1)   binary mask for water
// --- 1.1) Defintion/ Klarstellung: Beginn des Water Interface
// ============================================================================

{
    //==========================================================================
    // Robust Option C hEqn.H (OF10-safe + konservative Phasen-Konvektion)
    //==========================================================================

    const scalar alphaCut    = 1e-6;  
    const scalar alphaCutWC  = 1e-6;   
    const scalar alphaCutE   = 5e-2;   

    //======================================================================
    // 1) Temperaturen aus Enthalpien
    //======================================================================
    Tair     = Tref + hA/cpA;
    Twater   = Tref + hW/cpW;
    Tcontact = Tref + hC/cpC;

    Tair.correctBoundaryConditions();
    Twater.correctBoundaryConditions();
    Tcontact.correctBoundaryConditions();

    T = alphaA*Tair + alphaW*Twater + alphaC*Tcontact;
    T.correctBoundaryConditions();

    //======================================================================
    // 2) fs_old
    //======================================================================
    volScalarField fsOld
    (
        IOobject
        (
            "fsOld", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        fs
    );

    //======================================================================
    // 3) fs-Update (nur Wasser!)
    //======================================================================
    #include "updateFs.H"
    fs.max(0.0);
    fs.min(1.0);
    fs.correctBoundaryConditions();

    //======================================================================
    // 4) effektive Leitfähigkeiten
    //======================================================================
    #include "updatelambdaEff.H"

    //======================================================================
    // 5) Eis-front dFs(t)/dt
    //======================================================================
    volScalarField dFsdt
    (
        IOobject
        (
            "dFsdt", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        (fs - fsOld)/runTime.deltaT()
    );
    dFsdt *= pos(alphaW - alphaCut);

    dimensionedScalar dFsdtMax
    (
        "dFsdtMax", 
        dimless/dimTime, 
        0.5
    );
    dFsdt.min(dFsdtMax);
    dFsdt.max(-dFsdtMax);

    //======================================================================
    // 6) Latent heat source
    //======================================================================
    volScalarField ShW
    (
        IOobject
        (
            "ShW", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        -rhoW * L * dFsdt
    );

    //======================================================================
    // 7) Energie-Masken 
    //======================================================================
    // 7.2)   Deklarierung Energie Masken
    // 7.1.1) Maske air
    volScalarField maskA
    (
        IOobject
        (
            "maskA", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        pos(alphaA - alphaCutE)
    );
    // 7.1.2) Maske Wasser
    volScalarField maskW
    (
        IOobject
        (
            "maskW", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        pos(alphaW - alphaCutE)
    );
    // 7.1.3) Maske contactFluid
    volScalarField maskC
    (
        IOobject
        (
            "maskC", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        pos(alphaC - alphaCutE)
    );

    // 7.2) Hilfsfelder für Masken 
    // 7.2.1) Hilfsfeld Luft
    volScalarField alphaAeff
    (
        IOobject
        ( 
            "alphaAeff", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );
    // 7.2.2) Hilfsfeld Wasser
    volScalarField alphaWeff
    (
        IOobject
        (
            "alphaWeff", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ, 
            IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );
    // 7.2.3) Hilfsfeld contactFluid
    volScalarField alphaCeff
    (
        IOobject
        ("alphaCeff", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    alphaAeff = max(alphaA, scalar(alphaCutE))*maskA;
    alphaWeff = max(alphaW, scalar(alphaCutE))*maskW;
    alphaCeff = max(alphaC, scalar(alphaCutE))*maskC;

    //======================================================================
    // 7) Regularisierung
    //======================================================================
    dimensionedScalar hReg
    (
        "hReg",
        dimensionSet(1,-3,-1,0,0,0,0), // ~ rho/dt
        1e-6
    );

    const scalar T0 = 290.0;
    volScalarField hARef
    (
        IOobject
        (
            "hARef",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("hARef0", hA.dimensions(), cpA.value()*(T0 - Tref.value()))
    );
    volScalarField hWRef
    (
        IOobject
        (
            "hWRef",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("hWRef0", hW.dimensions(), cpW.value()*(T0 - Tref.value()))
    );
    volScalarField hCRef
    (
        IOobject
        (
            "hCRef",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("hCRef0", hC.dimensions(), cpC.value()*(T0 - Tref.value()))
    );

    //======================================================================
    // 8) PHASEN-KONSERVATIVE FLÜSSE
    //    h ist J/kg -> Konvektion muss mit phasen-MASSENFLUSS laufen:
    //      div(rho_i * phi_i, h_i)
    //======================================================================
    surfaceScalarField phiA
    (
        IOobject
        (
            "phiA", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        fvc::interpolate(alphaAeff) * phi
    );
    surfaceScalarField phiW
    (
        IOobject
        (
            "phiW", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        fvc::interpolate(alphaWeff) * phi
    );
    surfaceScalarField phiC
    (
        IOobject
        (
            "phiC", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        fvc::interpolate(alphaCeff) * phi
    );

    surfaceScalarField rhoPhiA
    (
        IOobject
        (
            "rhoPhiA", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        rhoA * phiA
    );
    surfaceScalarField rhoPhiW
    (
        IOobject
        (
            "rhoPhiW", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        rhoW * phiW
    );
    surfaceScalarField rhoPhiC
    (
        IOobject
        (
            "rhoPhiC", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        rhoC * phiC
    );

    //======================================================================
    // 9) interphase heat-transfer
    //======================================================================

    // 9.1)   overlap indicators

    // 9.1.1) overlap indicators: water - contactFluid
    volScalarField maskWC
    (
        IOobject
        (
            "maskWC", 
            runTime.timeName(), 
            mesh
        ),
        pos(alphaW - alphaCutWC)*pos(alphaC - alphaCutWC));

    volScalarField Iwc
    (
        IOobject
        (
            "Iwc", 
            runTime.timeName(), 
            mesh
        ), 
        4.0*alphaW*alphaC
    );
        Iwc *= maskWC;
        Iwc.max(0.0);
        Iwc.min(1.0);


    // 9.1.2) overlap indicators: water - air 
    // 9.1.2.1)   Hilfsfelder
    // 9.1.2.1.1) aus Interface Kopplung volumetrische Kopllung erzielen
    volScalarField deltaCell
    (
       IOobject
       (
          "deltaCell",
          runTime.timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::NO_WRITE
       ),
    mesh,
    dimensionedScalar("zero", dimLength, 0.0)
    );

    // 9.1.2.1.2) Defintion primitiveField
    deltaCell.primitiveFieldRef() = cbrt(mesh.V());
    // 9.1.2.1.3) einfache, harmlose Randwerte
    deltaCell.correctBoundaryConditions();

    // 9.1.2.2) Interface Indikator
    volScalarField Iwa
    (
        IOobject
        (
            "Iwa", 
            runTime.timeName(), 
            mesh,
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        deltaCell * mag(fvc::grad(alphaW))
    );
    
    Iwa *= pos(alphaW - alphaCutWC)*pos(alphaA - alphaCutWC)*maskW*maskA;
    Iwa.max(0.0);
    Iwa.min(1.0);

    // 9.2)  effektive couplings
    // 9.2.1) heat transfer water-contactFluid
    volScalarField Kwc
    (
        IOobject
        (
            "Kwc", 
            runTime.timeName(), 
            mesh
        ), 
        K_wc * Iwc
    );
    // 9.2.2) heat transfer water - air 
    volScalarField Kwa
    (
        IOobject
        (
            "Kwa", 
            runTime.timeName(), 
            mesh
        ), 
        K_wa * Iwa
    );

    // 9.3) Limiter
    Kwc.max(dimensionedScalar("zeroK", Kwc.dimensions(), 0.0));
    Kwa.max(dimensionedScalar("zeroK", Kwa.dimensions(), 0.0));

    Kwc.min(dimensionedScalar("KwcMax", Kwc.dimensions(), 1e9));
    Kwa.min(dimensionedScalar("KwaMax", Kwa.dimensions(), 1e9));

    // 9.4) Phasengewichtete Kopplung (Kopplung nur dort, wo Phase wirklich Masse besitzt)

    // 9.4.1) wirksame Kopplung: wasserseitig
    volScalarField KwcW
    (
        IOobject
        (
            "KwcW",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphaWeff*Kwc
    );
    volScalarField KwaW
    (
        IOobject
        (
            "KwaW",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphaWeff*Kwa
    );

    // 9.4.2) wirksame Kopplung: luftseitig
    volScalarField KwaA
    (
        IOobject
        (
            "KwaA",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphaAeff*Kwa
    );
    // 9.4.2) wirksame Kopplung: contactFluid - seitig
    volScalarField KwcC
    (
        IOobject
        (
            "KwcC",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphaCeff*Kwc
    );

    // 9.5) Debug
    dimensionedScalar KwcSum = fvc::domainIntegrate(Kwc);
    dimensionedScalar TwWc   = fvc::domainIntegrate(Kwc*Twater);
    dimensionedScalar TcWc   = fvc::domainIntegrate(Kwc*Tcontact);

    Info<< "KwcSum = " << KwcSum.value() << " [W/K]\n";
    if (KwcSum.value() > SMALL)
    {
        Info<< "Twater_on_WC  = " << TwWc.value()/KwcSum.value() << " [K]\n";
        Info<< "Tcontact_on_WC= " << TcWc.value()/KwcSum.value() << " [K]\n";
        Info<< "dT_on_WC      = "
            << (TcWc.value()-TwWc.value())/KwcSum.value() << " [K]\n";
    }
    Info<< "KwaW[min,max] = (" << gMin(KwaW) << ", " << gMax(KwaW) << ")\n";
    Info<< "KwcW[min,max] = (" << gMin(KwcW) << ", " << gMax(KwcW) << ")\n";

    dimensionedScalar QwcDom =
        fvc::domainIntegrate(Kwc*(Tcontact - Twater));
    Info<< "Qwc(domain) = " << QwcDom.value()
        << " [W] (pos => heat to water)\n";
    //======================================================================
    // 10) Energiegleichungen
    //======================================================================

    // 10.1) Hilfsfelder h_old 
    // 10.1.1) h_old air
    const volScalarField hAold
    (
        IOobject
        (
            "hAold", 
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        hA
    );
    // 10.1.2) h_old water
    const volScalarField hWold
    (
        IOobject
        (
            "hWold", 
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        hW
    );
    // 10.1.3) h_old contactFluid
    const volScalarField hCold
    (
        IOobject
        (
            "hCold", 
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        hC
    );

    // 10.2) Energiegleichungen

    // 10.2.1) Energiegleichung LUFT (AIR)
    fvScalarMatrix hEqnAir
    (
        fvm::ddt(alphaAeff*rhoA, hA)
      + fvm::div(rhoPhiA, hA)
      - fvm::laplacian(alphaAeff*lambdaA/cpA, hA)
      + fvm::Sp(hReg*(scalar(1)-maskA), hA)
      + fvm::Sp(KwaA/cpA, hA)
     ==
        (hReg*(scalar(1)-maskA))*hARef
      + KwaA*(hWold/cpW)   
    );
    hEqnAir.relax();
    hEqnAir.solve();

    // 10.2.2) Energiegleichung CONTACTFLUID 
    fvScalarMatrix hEqnContact
    (
        fvm::ddt(alphaCeff*rhoC, hC)
      + fvm::div(rhoPhiC, hC)
      - fvm::laplacian(alphaCeff*lambdaC/cpC, hC)
      + fvm::Sp(hReg*(scalar(1)-maskC), hC)
      + fvm::Sp(KwcC/cpC, hC)   
     ==
        (hReg*(scalar(1)-maskC))*hCRef
      + KwcC*(hWold/cpW)
    );
    hEqnContact.relax();
    hEqnContact.solve();

    // 10.2.3) Energiegleichung WATER 
    fvScalarMatrix hEqnWater
    (
        fvm::ddt(alphaWeff*rhoW, hW)
      + fvm::div(rhoPhiW, hW)
      - fvm::laplacian(alphaWeff*lambdaW_eff/cpW, hW)
      + fvm::Sp(hReg*(scalar(1)-maskW), hW)
      + fvm::Sp((KwcW+KwaW)/cpW, hW) 
     ==
        ShW
      + (hReg*(scalar(1)-maskW))*hWRef
      + KwcW*(hCold/cpC) 
      + KwaW*(hAold/cpA) 
    );
    hEqnWater.relax();
    hEqnWater.solve();

    //======================================================================
    // 11) HARD PARKING
    //======================================================================
    hA = maskA*hA + (scalar(1)-maskA)*hARef;
    hW = maskW*hW + (scalar(1)-maskW)*hWRef;
    hC = maskC*hC + (scalar(1)-maskC)*hCRef;

    hA.correctBoundaryConditions();
    hW.correctBoundaryConditions();
    hC.correctBoundaryConditions();

    //======================================================================
    // 12) Temperaturen
    //======================================================================
    Tair     = Tref + hA/cpA;
    Twater   = Tref + hW/cpW;
    Tcontact = Tref + hC/cpC;

    Tair.correctBoundaryConditions();
    Twater.correctBoundaryConditions();
    Tcontact.correctBoundaryConditions();

    //======================================================================
    // 13) Debug
    //======================================================================
    Info<< "T.air[min,max]          = (" << gMin(Tair)     << ", " << gMax(Tair)     << ")\n";
    Info<< "T.water[min,max]        = (" << gMin(Twater)   << ", " << gMax(Twater)   << ")\n";
    Info<< "T.contactFluid[min,max] = (" << gMin(Tcontact) << ", " << gMax(Tcontact) << ")\n";
}
