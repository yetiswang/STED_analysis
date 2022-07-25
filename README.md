# Matlab code accompanying ACS paper "Multicolor Super-resolution microscopy of protein corona on single nanoparticles"

These codes were developed by Yuyang Wang to analyze the data presented in the ACS paper Multicolor Super-resolution microscopy of protein corona on single nanoparticles.
The following codes involve full image analysis of STED and confocal microscopy data collected from Abberior Expertline STED microscope and Imspector software. 
A simple outline of the script is as below:
1. Image import with Bioformat Toolbox for Matlab, extraction of useful metadata for automated plotting and calculation
2. Multicolor image plotting 
3. Single particle recognition based on circle detection algorithm and near-neighbour distance filtering. A visualization handle is used for ease of evaluation
4. Cross-talk correction based on experimentally determined parameters and generation of corrected results.


For any questions, please contact Yuyang at <y.wang8@tue.nl>. 
