/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 AUTHOR,AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    frozenRANS

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wallDist.H"
#include "turbulenceModel.H"
#include "turbulentTransportModel.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "bound.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volScalarField> blend(
    const volScalarField& F1,
    const dimensionedScalar& psi1,
    const dimensionedScalar& psi2
)
{
    return F1*(psi1 - psi2) + psi2;
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "createMesh.H"
    #include "createFields.H"

    volScalarField divU (fvc::div(phi));
    volSymmTensorField S ("S",symm(fvc::grad(U)));

    #pragma region lambdas

    auto CDkOmega_ = [&] (){
        return max (
            2*sigmaOmega2/omega*(fvc::grad(k) & fvc::grad(omega)),
            dimensionedScalar(dimless/sqr(dimTime), 1.0e-10)
    );
    };

    CDkOmega= CDkOmega_();

    auto F1_ = [&] () {
        tmp<volScalarField> arg1 = min
        (
            min
            (
                max
                (
                    (scalar(1)/betaStar)*sqrt(k)/(omega*y),
                    scalar(500)*(nu)/(sqr(y)*omega)
                ),
                (4*sigmaOmega2)*k/(CDkOmega*sqr(y))
            ),
            scalar(10)
        );
    return tanh(pow4(arg1));
    };

    F1 = F1_();

    auto updateBlended =  [&] () mutable{
        sigmaK = blend(F1,sigmaK1,sigmaK2);
        sigmaOmega = blend(F1,sigmaOmega1,sigmaOmega2);
        gamma = blend(F1,gamma1,gamma2);
        beta = blend(F1,beta1,beta2);
    };

    updateBlended();

    auto F2_ = [&] () {
        tmp<volScalarField> arg2 = min
        (
            max
            (
                (scalar(2)/betaStar)*sqrt(k)/(omega*y),
                scalar(500)*(nu)/(sqr(y)*omega)
            ),
            scalar(100)
        );

        return tanh(sqr(arg2));
    };

    F2 = F2_();

    auto nut_ = [&] () {
        return a1*k/max(a1*omega,F2*mag(S));
    };
    
    nut = nut_();

    auto Pk_ = [&] () {
        return min(2.0*nut*magSqr(S),10*betaStar*omega*k);
    };

    Pk = Pk_();
    #pragma endregion
    
    while(runTime.loop()){
        Info<<"\nIteration "<<runTime.timeName()<<nl<<nl;
        CDkOmega = CDkOmega_();
        F1 = F1_();
        updateBlended();
        F2 = F2_();
        nut = nut_();
        Pk = Pk_();

        tmp<fvScalarMatrix> kEqn
        (
            fvm::div(phi, k)
            - fvm::laplacian(sigmaK*nut+nu, k)
            ==
            Pk
            - fvm::SuSp((2.0/3.0)*divU, k)
            - fvm::Sp(betaStar*omega, k)
            //TODO: add omegaInf and kInf (0 in the base case)
            // +betaStar_*omegaInf_*kInf_
            //+ kSource()                   //should be also 0 (?)
            //+ fvOptions(alpha, rho, k_)   // this too (?)
        );

        R.primitiveFieldRef() = kEqn().residual();

    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::div(phi, omega)
        - fvm::laplacian(sigmaOmega*nut+nu, omega)
        ==
        gamma*(Pk+R)/nut
        - fvm::SuSp((2.0/3.0)*gamma*divU, omega)
        - fvm::Sp(beta*omega, omega)
        - fvm::SuSp
        (
            (F1 - scalar(1))*CDkOmega/omega,
            omega
        )
        //+ beta_*sqr(omegaInf_)
        //+ Qsas(2*magSqr(S_), gamma_, beta_) // should be 0
        //+ omegaSource()
        //+ fvOptions(alpha, rho, omega_)
    );
    omegaEqn.ref().relax();
    omegaEqn.ref().boundaryManipulate(omega.boundaryFieldRef());
    solve(omegaEqn);
    bound(omega,dimensionedScalar("omegaMin",dimless/dimTime,1e-10)); 
    nut = nut_();  
    runTime.write();
    }
    Info<< nl;
    runTime.printExecutionTime(Info);
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
