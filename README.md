# General information

In this repository, the analysis of smFISH data, used in *Selective dendritic localization of mRNA in Drosophila mushroom body output neurons* (https://elifesciences.org/articles/62770), is committed. The detection and spot analysis is based on the literature:

- Smith, C., Joseph, N., Rieger, B. et al. Fast, single-molecule localization that achieves theoretically minimum uncertainty. Nat Methods 7, 373–375 (2010). https://doi.org/10.1038/nmeth.1449
- Smith CS, Stallinga S, Lidke KA, Rieger B, Grunwald D. Probability-based particle detection that enables threshold-free and robust in vivo single-molecule tracking. Mol Biol Cell. 2015;26(22):4057-4062. doi:10.1091/mbc.E15-06-0448

The analysis consists of three parts, the spot analysis of detections inside the dendrites, calyx and the analysis of the transcription focus. 
All code is working and tested on MATLAB R2019b with CUDA 10.2, Visual Studio 2019 and the toolbox DIPimage. Those can be found at:

https://developer.nvidia.com/cuda-10.2-download-archive

https://visualstudio.microsoft.com/vs/

http://www.diplib.org/download

For help on running the code, please read the manual. Added comments with/in the different function m-files will also help a user to find his/her way in the code. The example data can be downloaded from https://figshare.com/articles/dataset/Example_data/13568438.
