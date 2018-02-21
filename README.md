# Local-rainfall-model
This code allows for model calibration and stochastic simulation using a space-time stochastic rainfall model specially designed for local scale applications.

Written by: Lionel Benoit - University of Lausanne - Institute of Earth Surface Dynamics

How to use it?
(1) Example datasets are available in the Data folder.
(2) For first use, simply run (in Matlab) the main function called script_total. This script is split into sections that can be run  sequentially. 
(3) To change the data used in the script, change the variable input_file. To use your own dataset you have to create new data files as well as a new metedata file. Just mimic the examples and format your files following the same formats.
(4) Default values should work for most of the cases. However, if you wish, you can tune the parameters of the algorithm. Tunable parameters are grouped in the sections called 'parameters to set'.
(5) Advanced users can also modify the parameters of the priors used in the Metropolis Hastings sampler. To do it, you have to edit the function called set_Metropolis_hastings.

Enjoy the rain! :-)
Lionel
