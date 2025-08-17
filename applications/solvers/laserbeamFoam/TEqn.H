{
    int iter = 0;
    scalarField residual(epsilon1.primitiveField().size(), 1.0);
    scalar maxResidual = 1.0;
    scalar meanResidual = 1.0;
    // scalar Tresidual = 1.0;

    TRHS *= 0.0;

    rhoCp = rho*cp;
    rhophicp = fvc::interpolate(cp)*rhoPhi;

    // const volTensorField gradU(fvc::grad(U));
    // volTensorField tau
    // (
    //     ((mixture.mu()*gradU) + (mixture.mu()*gradU.T()))
    //   - ((2.0/3.0)*mixture.mu()*(fvc::div(phi))*I)
    // );

    do
    {
        epsilon1.storePrevIter();
        T.storePrevIter();
        TRHS =
            LatentHeat*((fvc::ddt(rho, epsilon1) + fvc::div(rhoPhi, epsilon1)));

        Qv =
            0.82*LatentHeatVap*Mm*p0
           *Foam::exp(LatentHeatVap*Mm*((T - Tvap)/(R*T*Tvap)))
           /Foam::pow(2*M_PI*Mm*R*T,0.5);

        if (damperSwitch)
        {
            scalarField& thermalDamperI = thermalDamper;
            const scalarField& rhoCpI = rhoCp;
            const scalarField& TI = T;

            forAll(thermalDamperI, celli)
            {
                const scalar Tloc = TI[celli];

                // Select metal cp(T) according to model
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

                // Original formula, gas cp is polynomial
                thermalDamperI[celli] =
                    2.0 * rhoCpI[celli]
                /(
                    rho1.value()*cpMetal
                    + rho2.value()*polycp_g.value(Tloc)
                    );
            }

            thermalDamper.correctBoundaryConditions();
        }



        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCp,T)
          + fvm::div(rhophicp, T)
          - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhophicp), T)
          - fvm::laplacian(kappa, T)
          //- (tau && gradU) // Viscous dissipation
         ==
            fvModels.source(rhoCp,T)
          + laser.deposition()
          - Qv*mag(gradAlpha)*thermalDamper
          - TRHS
        );

        TEqn.relax();

        fvConstraints.constrain(TEqn);

        TEqn.solve();

        T.relax();

        Tcorr = (TLiquidus - TSolidus)*epsilon1 + TSolidus;

        epsilon1 =
            max
            (
                min
                (
                    epsilon1 + (epsilonRel*cp/LatentHeat)*(T - Tcorr),
                    scalar(1)
                ),
                scalar(0)
            );

        residual =
            mag
            (
                epsilon1.primitiveField()
              - epsilon1.prevIter().primitiveField()
            );

        maxResidual = gMax(residual);

        meanResidual = gSum(residual*mesh.V())/gSum(mesh.V().field());

        //ddte1 = fvc::ddt(epsilon1);

        Info<< "Correcting epsilon1, mean residual = " << meanResidual
            << ", max residual = " << maxResidual
            << endl;
    }
    while
    (
        (
            ++iter < minTCorr
         // || (meanResidual > epsilonTol.value())
         || (maxResidual > epsilonTol.value())
        )
     && iter <= maxTCorr
    );

    gradT = fvc::grad(T);
}





