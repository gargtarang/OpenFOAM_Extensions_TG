# OpenFOAM_Extensions_TG
Extensions used by me in OpenFOAM. I use OpenFOAM foundation version.

All the extensions given below are working with OpenFOAM-6. Although I have not
tested, but I believe they will work with OpenFOAM-7 also, as there are no
major changes between them.

All working function objects are merged into a single library. 
Those which are under testing will made a separate library

Turbulence model is still work in progress for code cleanup and references.
Functionally it is working.

External Libraries to be imported in controlDict file:

1. ebeamTurbulentHeatFlux --> "libebeamCondition.so"
2. marangoniCondition --> "libmarangoniGradientCondition.so"
3. All Function Objects --> "libmyFunctionObjects.so"
4. kEpsilonMelting Turbulence Model --> "libmyTurbulenceModels.so"

