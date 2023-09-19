# The Geometric Mental Line
Matlab and Mathematica scripts used to perform the simulations and the analyses show in 

Di Antonio G, Raglio S and Mattia M. 2023.
_Ranking and serial thinking: A geometric solution_.
**bioRxiv**, 2023.08.03.551859: 1-21. doi:[10.1101/2023.08.03.551859v1](https://www.biorxiv.org/content/10.1101/2023.08.03.551859v1)

## Guide to reproduce the results
The following codes are used to reproduce the results in Figure 5-6-7.

The "MatlabLibrary" folder must be added to the path in order to run all other scripts. It contains codes related to network settings and dynamics.

The "Single_Run" folder contains codes to train and test a single network on the Transitive inference task and check it inside.

The "Multiple_Runs" folder contains codes related to results averaged over multiple task realizations and network structures.

For Figure 3 and 4 the code used to produce the presented results in the following folders:

"Fig3_GMLlearning", simulations and related plots (PDF files) of the Perceptron-like shallow networks learning via classical conditioning (pseudo-inverse approach equivalent to _delta rule_ learning) a transitive inference task.
Each subfolder is associated to a different sensory noise.
The scripts to run are **computeAndTestGMLOnAverage.m**.

"Fig4_IntegratorOnGML", simulations and related plots (PDF files) of the stochastic dynamics associated to the projection on the geometric mental line of the neural state of the network receiving as input the symbolic distance of the presented pairs of items.
The script to run is **plotIntegratorOnGML.m**.
In the subfolder "Mathematica" can be found the script **FPTofWienerProcWithTwoAbsorbBarr.nb** used to computed and plot the theoretical derivations related to the average reaction times and accuracy.
