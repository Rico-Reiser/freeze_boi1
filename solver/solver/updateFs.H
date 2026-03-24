// Solver Konfigurationsdatei: updateFs.H
// Steuerung des Phasenanteils flüssig - fest in Wasser

// Inhaltsverzeichnis
//
// --- 1)       Deklarierung des festen und flüssigen Bereichs                                     Zeile:
// --- 1.1)     Solidos-Temperatur                                                                 Zeile: 
// --- 1.2)     iquidos-Temperatur                                                                 Zeile:
// --- 1.3)     mushy-zone                                                                         Zeile:
// --- 2. )     Deklarierung der Volumenanteile innerhalb der aktuellen Zelle                      Zeile:
// --- 2.1)     Verhindert Gefrieren von contactFluid oder Luft                                    Zeile:
// --- 2.2)     lokale Temperatur in der Zelle                                                     Zeile:
// --- 2.3)     temperaturbasierter Phasenindikator                                                Zeile:
// --- 2.4)     Grenzfälle der Temperaturen                                                        Zeile:
// --- 2.3)     temperaturbasierter Phasenindikator                                                Zeile:
// --- 2.4.1)   unterhalb Solidus-Temperatur: Wasser gefriert  (fi = 1 /vollständig gefroren)      Zeile:
// --- 2.4.2)   unterhalb Liquidos-Temperatur: Wasser gefriert (fi = 0 /vollständig flüssig)       Zeile:
// --- 2.4.3)   normierte Temperatur in der Mushy-Zone                                             Zeile:
// --- 2.4.3.1) Cosinus-Interpolation des Phasenübergangs                                          Zeile:
// --- 2.5)     nur Wasseranteil kann gefrieren                                                    Zeile:
// --- 3.)      numerische Sicherheit: nur Wasser kann gefrieren                                   Zeile:
// --- 3.1)     global bonds                                                                       Zeile:
// --- 3.3)     Relaxation  (0.2 - 0.3 for mushy models)                                           Zeile:

{
    // --- 1. ) Deklarierung des festen und flüssigen Bereichs
    // --- 1.1) Solidos-Temperatur 
    const scalar Tsol = Tsolid.value();
    const scalar Tliq = Tliquid.value();
    const scalar dT   = (Tliq - Tsol) + VSMALL;

// --- 2) Schleife über Zellen
forAll(fs, celli)
{
    const scalar aW = alphaW[celli];

    // --- 2.1) Kein Wasser → kein Eis
    if (aW <= alphaCut)
    {
        fs[celli] = 0.0;
        continue;
    }

    // --- 2.2) Wassertemperatur
    const scalar Ti = Twater[celli];

    scalar fi;

    // --- 2.3) Phasenlogik
    if (Ti <= Tsol)
    {
        // vollständig gefroren
        fi = 1.0;
    }
    else if (Ti >= Tliq)
    {
        // vollständig flüssig
        fi = 0.0;
    }
    else
    {
        // Mushy-Zone (Cosinus-Interpolation)
        scalar xi = (Ti - Tsol)/dT;
        xi = min(max(xi, 0.0), 1.0);

        fi = 0.5*(1.0 + Foam::cos(constant::mathematical::pi*xi));
    }

    // --- 2.4) fs = Festanteil DES WASSERS
    fs[celli] = fi;
}

// --- 3) Numerische Sicherheit
fs.max(0.0);
fs.min(1.0);

// --- 4) Relaxation
const scalar fsRelax = 0.25;
// 4.1.1) globale Herangehensweise
// fs = (scalar(1) - fsRelax)*fsOld + fsRelax*fs;
// 4.1.2) lokale Herangehensweise
fs = fsOld + fsRelax*(fs - fsOld)*pos(alphaW - alphaCutE);



fs.correctBoundaryConditions();

}
