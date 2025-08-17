{
    const scalarField& alpha1I = alpha1;
    const scalarField& TI = T;
    scalarField& kappaI = kappa;
    scalarField& cpI = cp;

    forAll(kappaI, celli)
    {
        const scalar Tloc = TI[celli];

        // Gas phase (unchanged): always polynomial
        const scalar kGas = polykappa_g.value(Tloc);
        const scalar cpGas = polycp_g.value(Tloc);

        // ---- Metal kappa(T) ----
        scalar kMetal = 0.0;
        if (kappaModel == "table")
        {
            if (!kappaTab_mPtr.valid())
            {
                FatalErrorInFunction
                    << "kappaModel=table but 'table_kappa' missing in physicalProperties.metal" << nl
                    << exit(FatalError);
            }
            kMetal = (*kappaTab_mPtr)(Tloc);
        }
        else if (kappaModel == "polynomial")
        {
            if (!polykappa_mPtr.valid())
            {
                FatalErrorInFunction
                    << "kappaModel=polynomial but 'poly_kappa' missing in physicalProperties.metal" << nl
                    << exit(FatalError);
            }
            kMetal = polykappa_mPtr->value(Tloc);
        }
        else if (kappaModel == "hybrid")
        {
            if (!kappaTab_mPtr.valid() || !polykappa_mPtr.valid())
            {
                FatalErrorInFunction
                    << "kappaModel=hybrid requires both 'table_kappa' and 'poly_kappa' in physicalProperties.metal" << nl
                    << exit(FatalError);
            }
            const scalar kTab = (*kappaTab_mPtr)(Tloc);
            const scalar kPol = polykappa_mPtr->value(Tloc);
            kMetal = wTableKappa * kTab + (1.0 - wTableKappa) * kPol;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown kappaModel: " << kappaModel
                << " (expected 'polynomial', 'table', or 'hybrid')." << nl
                << exit(FatalError);
        }

        // ---- Metal cp(T) ----
        scalar cpMetal = 0.0;
        if (cpModel == "table")
        {
            if (!cpTab_mPtr.valid())
            {
                FatalErrorInFunction
                    << "cpModel=table but 'table_cp' missing in physicalProperties.metal" << nl
                    << exit(FatalError);
            }
            cpMetal = (*cpTab_mPtr)(Tloc);
        }
        else if (cpModel == "polynomial")
        {
            if (!polycp_mPtr.valid())
            {
                FatalErrorInFunction
                    << "cpModel=polynomial but 'poly_cp' missing in physicalProperties.metal" << nl
                    << exit(FatalError);
            }
            cpMetal = polycp_mPtr->value(Tloc);
        }
        else if (cpModel == "hybrid")
        {
            if (!cpTab_mPtr.valid() || !polycp_mPtr.valid())
            {
                FatalErrorInFunction
                    << "cpModel=hybrid requires both 'table_cp' and 'poly_cp' in physicalProperties.metal" << nl
                    << exit(FatalError);
            }
            const scalar cpTab = (*cpTab_mPtr)(Tloc);
            const scalar cpPol = polycp_mPtr->value(Tloc);
            cpMetal = wTableCp * cpTab + (1.0 - wTableCp) * cpPol;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown cpModel: " << cpModel
                << " (expected 'polynomial', 'table', or 'hybrid')." << nl
                << exit(FatalError);
        }

        const scalar a = max(min(alpha1I[celli], 1.0), 0.0);
        kappaI[celli] = a*kMetal + (1.0 - a)*kGas;
        cpI[celli]    = a*cpMetal + (1.0 - a)*cpGas;
    }

    kappa.correctBoundaryConditions();
    cp.correctBoundaryConditions();
}
