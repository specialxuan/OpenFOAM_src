// Shared PIMPLE loop body included by the partitioned, integrated and explicit
// FSI branches of myRhoPimpleFoam.C; branch-specific logic stays in the main file.

                if (pimple.firstIter() || moveMeshOuterCorrectors)
                {
                    // Store momentum to set rhoUf for introduced faces.
                    autoPtr<volVectorField> rhoU;
                    if (rhoUf.valid())
                    {
                        rhoU.reset(new volVectorField("rhoU", rho*U));
                    }

                    // Do any mesh changes
                    mesh.controlledUpdate();

                    if (mesh.changing())
                    {
                        MRF.update();

                        if (correctPhi)
                        {
                            // Calculate absolute flux
                            // from the mapped surface velocity
                            phi = mesh.Sf() & rhoUf();

                            #include "correctPhi.H"

                            // Make the fluxes relative to the mesh-motion
                            fvc::makeRelative(phi, rho, U);
                        }

                        if (checkMeshCourantNo)
                        {
                            #include "meshCourantNo.H"
                        }
                    }
                }

                if (pimple.firstIter() && !pimple.SIMPLErho())
                {
                    #include "rhoEqn.H"
                }

                #include "UEqn.H"
                #include "EEqn.H"

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    if (pimple.consistent())
                    {
                        #include "pcEqn.H"
                    }
                    else
                    {
                        #include "pEqn.H"
                    }
                }

                if (pimple.turbCorr())
                {
                    turbulence->correct();
                }
