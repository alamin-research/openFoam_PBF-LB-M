{
    const scalarField& alpha1I = alpha1;
    const scalarField& TI = T;
    scalarField& kappaI = kappa;
    scalarField& cpI = cp;

    forAll(kappaI, celli)
    {
        const scalar Tloc = TI[celli];

        // ---- Metal kappa(T): selectable model ----
        const scalar kMetal_poly  = polykappa_m.value(Tloc);

        // default to polynomial if no table present
        scalar kMetal_table = kMetal_poly;
        if (kappaTab_mPtr.valid())
        {
            kMetal_table = (*kappaTab_mPtr)(Tloc);
        }

        scalar kMetal = kMetal_poly;
        if (kappaModel == "table" && kappaTab_mPtr.valid())
        {
            kMetal = kMetal_table;
        }
        else if (kappaModel == "hybrid" && kappaTab_mPtr.valid())
        {
            // blend: wTableKappa * table + (1-w) * polynomial
            kMetal = wTableKappa * kMetal_table + (1.0 - wTableKappa) * kMetal_poly;
        }
        // else: "polynomial" (default)

        // ---- Metal cp(T): selectable model ----
        const scalar cpMetal_poly  = polycp_m.value(Tloc);

        // default to polynomial if no table present
        scalar cpMetal_table = cpMetal_poly;
        if (cpTab_mPtr.valid())
        {
            cpMetal_table = (*cpTab_mPtr)(Tloc);
        }

        scalar cpMetal = cpMetal_poly;
        if (cpModel == "table" && cpTab_mPtr.valid())
        {
            cpMetal = cpMetal_table;
        }
        else if (cpModel == "hybrid" && cpTab_mPtr.valid())
        {
            // blend: wTableCp * table + (1-w) * polynomial
            cpMetal = wTableCp * cpMetal_table + (1.0 - wTableCp) * cpMetal_poly;
        }
        // else: "polynomial" (default)

        // ---- Phase mixing (metal/gas) unchanged ----
        const scalar a = max(min(alpha1I[celli], 1.0), 0.0);

        kappaI[celli] =
              a       * kMetal
            + (1.0-a) * polykappa_g.value(Tloc);

        cpI[celli] =
              a       * cpMetal
            + (1.0-a) * polycp_g.value(Tloc);
    }

    kappa.correctBoundaryConditions();
    cp.correctBoundaryConditions();
}
