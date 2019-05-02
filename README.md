# YellowRust_Wheat

We present a model that simulates epidemics of a pathogen on a crop, where the life-cycle parameters of the pathogen are driven by weather.

The model is parameterised for yellow rust (Puccinia striiformis) on wheat, as published in Kirtphaiboon et al., 2019 (if all goes to plan), but can be generalised to other foliar crop pathogens.

The model uses a simple wheat canopy model coupled with an SLIR model of a pathogen epidemic. The life-cycle parameters of the epidemic are driven by weather, specifically temperature, leaf wetness and light.
Using this model we were able to demonstrate that the temperature and leaf wetness are both limiting factors at different times in the epidemic life-cycle.

The code in this repository is written in the R programming language.

## Files

* ModelEquations.r - contains all the functions needed for running the model. Key among these functions are:
	- runModel(), which is the top-level function that runs the model over multiple years
	- modelEquations2(), which contains the derivatives of the state variables
* Parameters.r - contains the parameters used in the Kirtphaiboon paper
* WData.csv - an example weather dataset downloaded from the Enviromental Change Network in 2001 and 2002 at Rothamsted.
* Script.r - a script to run the model for Rothamsted 2016 and plot graphs
