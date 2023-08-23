# MA_thesis_2023
Example code for numerical continuation of a physical model of a trumpet-player system, used in Margaret Hopkins' 2023 MA thesis for comparisons between trumpets

This code repository provides example code for the numerical continuation used in Margaret Hopkins' MA thesis, "Comparisons between trumpet models using numerical continuation methods" (Hopkins, 2023), to assist with future research. Some instructions are provided below to assist in using the example code in this repository. For additional information and example files, see the user guides and example files provided with the Manlab package download on the [Manlab homepage](url). System B, described in Equations 3.17 and 3.18 in (Hopkins, 2023), is used for all numerical continuation computations in the thesis. The code used in this example was created using Manlab version 4.0, by modifying the example files included with the download. The following sections give instructions for using the files provided to locate a bifurcatuion point and perform continuation of the periodic solution. 

## Location of a Hopf bifurcation point

To locate a Hopf bifurcation point, the Manlab "SystAQ" class is used. First, create a directory within your Manlab directory with files called "equations.m," "point_display.m," "global_display.m," and a launch file with any title. The examples and notes below describe the content of the "Launch" and "Equations" files. The "point_display.m" and "global_display.m" can be used without modification.

The system used for this computation is System B from (Hopkins, 2023), with all time derivatives set to zero so that the system is at equilibrium.

### Launch File
The launch file is used to run the Manlab code, specify options, declare parameters, and initialize the system. To locate the bifurcation, run the launch file, and then, once the Manlab graphical user interface (GUI) has opened, use the "forward" button to compute the solution in steps as the $p_0$ value is varied. Continue to press the "forward" button until the curve changes from solid to dotted, at the point of the Hopf bifurcation, marked with a star. Next, click "Section" in the "Export" section of the GUI, then click on the bifurcation point in the diagram. Right click on the "Section" variable in the Matlab workspace, and save the section as a ".mat" file, to be used later to initialize the oscillating solution. 

### Equations file
The "Equations" file is used to define the equations of the system. For the location of the Hopf bifurcation point, we use the system with all terms with time derivatives set to zero, as the system is at equilibrium. In the file provided, the vector "u" includes the main variables in the system, in the order in which they are defined in the launch file, "lambda" is the continuation parameter, and the "Ua" vector includes the auxiliary variables in the order in which they are defined. 

## Continuation of the oscillating solution
In this section, the section representing the Hopf bifurcation point found in the previous section is used to initialize the periodic solution, using the "SystODE" class. All time derivatives are included, since the system is no longer at equilibrium. The directory is set up in the same way as for the location of the Hopf bifurcation point, with a launch file, an "equations.m" file, a "point_display" file, and a "global_display" file. Once again, the "point_display" and "global_display" files do not need to be modified.

### Launch file
As before, the launch file is run to perform the computations. The "Section.mat" file used should correspond to the one saved when locating the Hopf bifurcation. To create the full bifurcation diagram, a few different steps are executed after the launch file is run. 

First, click "reverse tangent," in order to follow the solution that continues towards the fold, and not the one that continues towards negative amplitude. After this, click the "forward" button a few times to compute a few steps. 

Next, change the value in the "correction" box to 1, so that correction is turned on and the results will be corrected to follow the solution. The reason that correction is not enabled previously is so that the solution is not close enough to equilibrium to correct towards that solution instead. After enabling correction, click "forward" one more time to compute another step. 

Next, click "cancel all" to remove all previously computed steps so that the diagram can be generated again, starting with the correct solution path. Now, continue the solution to the desired maximum pressure (keeping in mind that the dimensionless pressure is scaled from the actual blowing pressure by a factor of *P_M*, then click "reverse tangent" to compute the solution backwards until the y-axis reaches zero. A complete bifurcation diagram has now been computed. The "stability computation" box can be checked during these computations in order to plot the solutions with the stability indicated, but the computation will be much slower and the parameters associated with stability computation may need to be adjusted. 

### Equations file

The "equations.m" file once again contains the definitions of all main and auxiliary equations, except that the time derivatives are now included. 
