# Purpose

Perform lightweight fragment analysis on [Amsterdam Modeling Suite](https://www.scm.com/amsterdam-modeling-suite/) xyz movies, independent of the engine (adf, band, etc...).
This script was developed specifically to work alongside with slurm on a hpc cluster, but can also be run locally. 

Reports the following EDA terms over the span of a reaction (or movement):
 * Strain of first fragment 
 * Strain of second fragment
 * Strain of system
 * Interaction of system

 * Pauli repulsion
 * Dispersion
 * Elecrostatic interction
 * Orbital interaction
   
### Installation
1. Download/clone the repository
2. Navigate to the project

`  cd LightFrag/`

3. Run the installation script
   
`  bash installation.sh`

4. Test your installation by running the sample input (LightFrag/sample_input/sample_bash)
   
`  sbatch sample_bash` on slurm or `bash sample_bash` (not recommended)

### Usage
1. Copy the sample_bash input file to the desired location for your calculations, as long as LightFrag is located in the root folder (cd ~).
2. Change the slurm, ams and engine settings in your input file to your desired settings. 
3. Add your .ams.amv (xyz movie) to the same directory as your input file
4. Configure the LightFrag settings in your input file
5. Run your input file with `sbatch whatever_your_input_file_is`
6. Analyse results. Raw results are displayed in the .out and a .csv file. Graphs are printed as a .ipynb in the same directory

### Script workflow
To perform the fragment analysis, the input_file is split into `fragment_1`, `fragment_2` and `fragment_full` for every frame, using the .ams.amv coordinates and frames, and using the fragment_mapping 
defined in the input file.

After running every fragment for every frame, the script compiles all results, using the strain_fragment_1 and strain_fragment_2 energies as basis, plotted over the bond_length defined in the input file. 
Running a LightFrag job also produces the following:
* `scripts/`: Folder containing the scripts that were run for every frame
* `frame_n`: For every frame, contains the ams output for the fragments. By default, these are pruned to only contain the ams.rkf and the ams.log for storage efficiency. 
