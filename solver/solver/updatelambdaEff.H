// ============================================================================
// updateKEff.H  (OPTION C – getrennte Phasen)
// kW_eff : effektive Leitfähigkeit der Wasserphase (Eis + Flüssig)
// kEff   : OPTIONAL (nur Debug / Postprocessing)
// ============================================================================

{
    // ------------------------------------------------------------------------
    // 1) Effektive Wärmeleitfähigkeit der Wasserphase (fest + flüssig)
    //     fs = Festanteil des Wassers (0…1)
    // ------------------------------------------------------------------------
    lambdaW_eff =
          fs*lambdaSolid
        + (scalar(1) - fs)*lambdaLiquid;

    // Numerische Sicherheit
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

    // ------------------------------------------------------------------------
    // 2) OPTIONALES Misch-kEff (NICHT für Energiegleichung verwenden!)
    //    Nur für:
    //      - ParaView
    //      - Debug
    //      - evtl. Darcy / Diagnostik
    // ------------------------------------------------------------------------
    lambdaEff =
          alphaW*lambdaW_eff
        + alphaC*lambdaC
        + alphaA*lambdaA;

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

    // ------------------------------------------------------------------------
    // 3) Debug-Ausgabe
    // ------------------------------------------------------------------------
    Info<< "lambdaW_eff[min,max] = ("
        << gMin(lambdaW_eff) << ", " << gMax(lambdaW_eff) << ")\n";

    Info<< "lambdaEff[min,max]   = ("
        << gMin(lambdaEff) << ", " << gMax(lambdaEff) << ")\n";
}
