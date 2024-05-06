singlePhaseTransportModel laminarTransport(U, phi);

// Solutal expansion coefficient [1/K]
dimensionedScalar beta_T
(
    "beta_T",
    dimless/dimTemperature,
    laminarTransport
);

// Solutal expansion coefficient [m3/mol]
dimensionedScalar beta_c
(
    "beta_c",
    dimVolume/dimMoles,
    laminarTransport
);

// Reference Temperature [K]
dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);

// Reference concentration [mol/m3]
dimensionedScalar rho_iRef("rho_iRef", dimMoles/dimVolume, laminarTransport);

// Laminar Pr number
dimensionedScalar Pr("Pr", dimless, laminarTransport);

// 1/rho/cp
dimensionedScalar rho_cp("rho_cp", dimMass/dimLength/dimTemperature/dimTime/dimTime, laminarTransport);

// Diffusion Coefficient
dimensionedScalar Dmol("Dmol", dimViscosity, laminarTransport);

// Turbulent Schmidt number
dimensionedScalar Sct("Sct", dimless, laminarTransport);

// Turbulent Pr number
dimensionedScalar Prt("Prt", dimless, laminarTransport);

// Electrical condictivity of the fluid
dimensionedScalar kf("kf", dimensionSet(-1,-3,3,0,0,2,0), laminarTransport); //S/m == A/V/m == A^2s^3/kg/m^3



