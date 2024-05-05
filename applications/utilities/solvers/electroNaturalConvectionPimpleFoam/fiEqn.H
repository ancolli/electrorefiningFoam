{

bool repeat(false);
// check if wall_anode exist
string patchNameNeeded = "wall_anode";
scalar patchOK = 0;

    forAll(mesh.boundaryMesh(), p_ith)
    {     
      string patchName = mesh.boundaryMesh()[p_ith].name();	
      if (patchNameNeeded == patchName)  // checks if "wall_anode" is in the current patchName
      {
        patchOK = 1;
      }
    }  

if (patchOK==1)  // checks if "wall_anode" is in the current patchName
      {
// if wall_anode exists then: Calculates flux on the anode
	
	// reading and writing probeLocation (reference electrode)

        label cellNo = mesh.findCell(vector(controlProperties.lookup("refLocation")));
	
        //Results.set("Vref", Vref);

// calculate cell potential
        scalar mult = 0.0;
        scalar Vcell_old = 0.0;
        scalar var = 0.0;
         
        int control = readScalar(controlProperties.lookup("control"));  
        scalar I = readScalar(controlProperties.lookup("I"));
        scalar maxPot = readScalar(controlProperties.lookup("maxPot"));
        scalar alpha = readScalar(controlProperties.lookup("alpha"));
        scalar tol_I = readScalar(controlProperties.lookup("tol_I"));

        scalar ba = readScalar(transportProperties.lookup("ba"));

// potentiostatic mode
	scalar ref = readScalar(controlProperties.lookup("ref"));
        scalar maxvar = readScalar(controlProperties.lookup("maxvar"));
        scalar tol_ref = readScalar(controlProperties.lookup("tol_ref"));       

	scalar Vcell = readScalar(Results.lookup("Vcell"));

	do {  // permite cumplir con modo galvanostatico o potenciostatico luego de cada timeStep
		
// Calculates flux on the anode and cathode
        label patchi_A = mesh.boundaryMesh().findPatchID("wall_anode"); 
	label patchi_C = mesh.boundaryMesh().findPatchID("wall_cathode");  
        Flux.boundaryFieldRef()[patchi_A] = keff.boundaryFieldRef()[patchi_A]*fi.boundaryFieldRef()[patchi_A].snGrad(); //gradients are defined inwards
	Flux.boundaryFieldRef()[patchi_C] = keff.boundaryFieldRef()[patchi_C]*fi.boundaryFieldRef()[patchi_C].snGrad(); //gradients are defined inwards
// Calculates average Flux
	scalar sumField = 0.0;
	scalar sumField_C = 0.0;
        {
        	sumField = gSum
                (
                	mesh.magSf().boundaryField()[patchi_A]
                      * Flux.boundaryField()[patchi_A]
                );
		sumField_C = gSum
                (
                	mesh.magSf().boundaryField()[patchi_C]
                      * Flux.boundaryField()[patchi_C]
                );
        }

        Info<< "Total current anode = " << sumField << " A "<< endl;
	Info<< "Total current cathode = " << sumField_C << " A "<< endl;	
	
// potentiostatic mode
	scalar Vref = fi[cellNo];
	var=(ref < 0.0) ? -ref-Vref : ref-(Vcell-Vref); // <0 cathodic overpotential, >0 anodic overpotential
        var=(var > 0.0) ? min(var,maxvar) : max(var,-maxvar);
        var=(mag(var/ref) > tol_ref) ? var : 0.0; // limito ref a x%

// galvanostatic
    
        mult = (mag(sumField-I)/I > tol_I) ? (sumField-I)/(I/ba) : 0.0; // limito corriente a x%
               
// overrelaxation       
        Vcell_old = Vcell;
	Vcell = (control==1) ? Vcell_old+var : ((control==2) ? maxPot : Vcell_old-mult);
        Vcell = Vcell_old+alpha*(Vcell-Vcell_old);
        //Info<< "var:" << var << endl;
	Vcell = max(min(Vcell, maxPot),mag(ref)); // limited potentials

// writing Vcell, I and Vref       
        Results.set("Vcell", Vcell);
        Results.set("I", sumField);
	Results.set("Vref", (ref < 0.0) ? -Vref : (Vcell-Vref));//
          
// potential field
        fi.storePrevIter();
    	solve
    	(
       		- fvm::laplacian(keff, fi)
    	);
    	fi.relax();
    	
    	fi.correctBoundaryConditions();
    	
    	Info<< "Vcell:" << Vcell << endl;

	Info<< "Min/max fi:" << min(fi).value() << ' '
    		<< max(fi).value() << endl;
    	   	
    	switch(control) {
        	case 1: repeat = (var != 0.0) ? true : false;   break; // potentiostatic mode
        	case 2: repeat = ((mag(sumField+sumField_C)/sumField) > tol_I) ? true : false; break; //fixed cell potential
        	case 3: repeat = (mult != 0.0) ? true : false;  break; // galvanostatic mode
        	default: Info<< "control = 1-potentiostatic mode" 
        			<< "control = 2-fixed cell potential"
        			<< "control = 3-galvanostatic mode"
        			<< endl
            			<< exit(FatalError); 
    	}
   		
    	} while (repeat);	
   }
   else
   {
	FatalErrorIn
        (
            " Foam::patchProbes::findElements(const fvMesh&)"
        )   << " No patch name: wall_anode found "
            << " the patch name: wall_anode should exist "
            << endl
            << exit(FatalError);     
   }
}

// ************************************************************************* //
