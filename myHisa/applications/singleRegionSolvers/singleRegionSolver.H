/*---------------------------------------------------------------------------*\

    HiSA: High Speed Aerodynamic solver

    Copyright (C) 2014 Oliver Oxtoby - CSIR, South Africa
    Copyright (C) 2020 Oliver Oxtoby
    Copyright (C) 2014 Johan Heyns - CSIR, South Africa
    Copyright (C) 2004-2011 OpenCFD Ltd.
    Copyright (C) 2016 OpenFOAM Foundation

-------------------------------------------------------------------------------
License
    This file is part of HiSA.

    HiSA is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HiSA is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with HiSA.  If not, see <http://www.gnu.org/licenses/>.

Application
    singleRegionSolver

Description
    Loads and calls the appropriate solver module for single region
    cases

See also
    multiRegionSolver

Authors
    Oliver Oxtoby
    Johan Heyns
    Council for Scientific and Industrial Research, South Africa

\*---------------------------------------------------------------------------*/

#include "pimpleControl.H"
#include "pseudotimeControl.H"
#include "solverModule.H"
#include "argList.H"
#include "timeSelector.H"
#include "fsiCoupling.H"

using namespace Foam; 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addBoolOption("version", "Print version information");

    argList::addBoolOption
    (
        argList::postProcessOptionName,
        "Execute functionObjects only"
    );

    if (argList::postProcess(argc, argv))
    {
        Foam::timeSelector::addOptions();
        #include "addRegionOption.H"
        #include "addFunctionObjectOptions.H"

        // Set functionObject post-processing mode
        functionObject::postProcess = true;

        #include "setRootCase.H"

        #if OPENFOAM >= 1906
        if (args.found("list"))
        #else
        if (args.optionFound("list"))
        #endif
        {
            #if FOUNDATION >= 11
            Info<< nl
                << "Available configured functionObjects:"
                << listAllConfigFiles
                   (
                       functionEntries::includeFuncEntry::functionObjectDictPath
                   )
                << endl;
            #else            
            functionObjectList::list();
            #endif
            return 0;
        }

        #include "createTime.H"
        Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

        // Externally stored dictionary for functionObjectList
        // if not constructed from runTime
        dictionary functionsDict;

        #if OPENFOAM >= 1906
        HashSet<wordRe> selectedFields;
        #else
        HashSet<word> selectedFields;
        #endif

        // Construct functionObjectList
        autoPtr<functionObjectList> functionsPtr;

        autoPtr<solverModule> solverPtr;

        forAll(timeDirs, timei)
        {
            runTime.setTime(timeDirs[timei], timei);

            Info<< "Time = " 
                #if FOUNDATION >= 12
                << runTime.name() 
                #else
                << runTime.timeName() 
                #endif
                << endl;

            functionsPtr.clear();

            dictionary dict;
            dict.add("solver", MODULE);
            solverPtr = solverModule::New(Foam::polyMesh::defaultRegion, runTime, dict);

            solverPtr->initialise();

            FatalIOError.throwExceptions();

            try
            {
                functionsPtr = functionObjectList::New
                (
                    args,
                    runTime
                    #if !(FOUNDATION >= 12)
                    , functionsDict
                    #endif
                    #if !(FOUNDATION >= 10)
                    , selectedFields
                    #endif
                );
                functionsPtr->execute();

                // Execute the functionObject 'end()' function for the last time
                if (timei == timeDirs.size()-1)
                {
                    functionsPtr->end();
                }
            }
            catch (IOerror& err)
            {
                Warning<< err << endl;
            }

            Info<< endl;
        }

        Info<< "End\n" << endl;

        return 0;
    }

    #include "setRootCase.H"

#if OPENFOAM >= 1906
    if (args.found("version"))
#else
    if (args.optionFound("version"))
#endif
    {
        Info<<
            #include "../../versionInfo.H"
            << endl;
        return 0;
    }

    #include "createTime.H"

    dictionary dict;
    dict.add("solver", MODULE);
    autoPtr<solverModule> solverPtr =
        solverModule::New(Foam::polyMesh::defaultRegion, runTime, dict);

    solverModule& solver = solverPtr();

    solver.initialise();

    fvMesh& mesh = solver.mesh();

    bool steadyState = solver.steadyState();

    if (!steadyState)
    {
        #include "createTimeControls.H"

        scalar scale = solver.timeStepScaling(maxCo);
        //Generate a pretend Courant number so we can use setInitialDeltaT.H unmodified
        scalar CoNum = maxCo/scale;
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    if (steadyState)
    {
        Info << "\nStarting pseudotime iteration loop\n" << endl;
    }
    else
    {
        Info<< "\nStarting time loop\n" << endl;
    }

    scalar prevCpuTime = runTime.elapsedCpuTime();
    while (runTime.run())
    {
        if (!steadyState)
        {
            #include "createTimeControls.H"

            // Scale deltaT
            scalar scale = solver.timeStepScaling(maxCo);
            if (adjustTimeStep)
            {
                scalar deltaT = scale*runTime.deltaTValue();
                deltaT = min
                (
                    min(deltaT, runTime.deltaTValue() + 0.1*deltaT),
                    1.2*runTime.deltaTValue()
                );
                runTime.setDeltaT(min(deltaT, maxDeltaT));

                Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
            }
        }

        solver.preTimeIncrement();

        runTime++;
        Info<< "Time = " 
            #if FOUNDATION >= 12
            << runTime.name() 
            #else
            << runTime.timeName() 
            #endif
            << nl << endl;

        solver.beginTimeStep();

        // --- Outer corrector loop
        bool notFinished = false;

        // FSI coupling controls from constant/dynamicMeshDict; shared
        // contract with the fastDynamicFvMesh modal dynamic-mesh library
        const fsiCouplingControls fsiControls =
            readFsiCouplingControls(runTime, mesh);

        if (isPartitionedFsiCoupling(fsiControls.coupling) && !steadyState)
        {
            // Partitioned FSI: the mesh stays fixed inside each pseudotime
            // sub-loop and moves once per FSI iteration (on the first
            // pseudotime iteration, requiring moveMeshOuterCorrectors 0)
            label fsiIter = 0;
            bool fsiConverged = false;
            pointField points0 = mesh.points();

            while (fsiIter < fsiControls.maxIter && !fsiConverged)
            {
                fsiIter++;

#if FOUNDATION >= 6
                while (true)
                {
                    if (isA<pimpleControl>(solver.solnControl()))
                    {
                        notFinished =
                            refCast<pimpleControl>
                                (solver.solnControl()).loop();
                    }
                    else
                    {
                        notFinished =
                            refCast<pseudotimeControl>
                                (solver.solnControl()).loop();
                    }
                    if (!notFinished)
                    {
                        break;
                    }
#else
                while ( (notFinished = solver.solnControl().loop()) )
                {
#endif
                    solver.outerIteration();
                }

                const scalar maxDiff = fsiMaxDisplacementChange
                (
                    mesh,
                    points0,
                    fsiIter,
                    "FSI Iteration"
                );
                const scalar fsiRes = readFsiResidual(mesh);
                fsiConverged = (fsiRes < fsiControls.tolerance);
                points0 = mesh.points();

                logFsiIteration("FSI Iteration", fsiIter, fsiRes, maxDiff);
            }

            if (!fsiConverged)
            {
                Info<< "Warning: FSI force residual did not converge within "
                    << fsiControls.maxIter << " FSI iterations." << endl;
            }
        }
        else
        {
            // Integrated/implicit FSI: the pseudotime control itself gates
            // convergence on the modal force residual; log it each iteration
            const bool integratedFsi =
                !steadyState
             && isIntegratedFsiCoupling(fsiControls.coupling);

            label fsiIter = 0;
            pointField points0;
            if (integratedFsi)
            {
                points0 = mesh.points();
            }

#if FOUNDATION >= 6
        while (true)
        {
            // loop() function was removed from the base class
            if (isA<pimpleControl>(solver.solnControl()))
            {
                notFinished = refCast<pimpleControl>(solver.solnControl()).loop();
            }
            else
            {
                notFinished = refCast<pseudotimeControl>(solver.solnControl()).loop();
            }
            if (!notFinished)
            {
                break;
            }
#else
        while ( (notFinished = solver.solnControl().loop()) )
        {
#endif
            solver.outerIteration();

            if (integratedFsi)
            {
                fsiIter++;
                const scalar maxDiff = fsiMaxDisplacementChange
                (
                    mesh,
                    points0,
                    fsiIter,
                    "FSI/Pseudo Iteration"
                );
                const scalar fsiRes = readFsiResidual(mesh);
                points0 = mesh.points();

                logFsiIteration
                (
                    "FSI/Pseudo Iteration",
                    fsiIter,
                    fsiRes,
                    maxDiff
                );
            }

            if (steadyState)
            {
                break;
            }
        }

            if
            (
                integratedFsi
             && readFsiResidual(mesh) >= fsiControls.tolerance
            )
            {
                Info<< "Warning: FSI force residual did not converge within "
                    << "pseudotime iterations." << endl;
            }
        }
        if (steadyState && !notFinished)
        {
            runTime.writeAndEnd();
        }
        else
        {
            runTime.write();
        }

        scalar cpuTime = runTime.elapsedCpuTime();
        Info<< "ExecutionTime = " << cpuTime << " s"
            << "  ExecutionStepTime = " << cpuTime-prevCpuTime << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        prevCpuTime = runTime.elapsedCpuTime();

    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
