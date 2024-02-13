# Omega-K Radar Image Reconstruction Algorithm

## Introduction

The Omega-K algorithm is a radar image reconstruction algorithm that consists of several key steps: Residual Video Phase compensation, Fourier Transform, Referent Function Multiplication, Stolt Interpolation, and Inverse Fourier Transform. This Python code implements the Omega-K algorithm for Ground Based Syntetic Aperture Radar (GBSAR) image reconstruction. Algorithm is based on [1], while its implementation is described in [2]. 


## Dependencies

Make sure you have the following Python libraries installed before running the code:

- matplotlib
- numpy
- seaborn
- scipy

## Usage

Set the filename variable to the path of the radar data file. 'Pastic_Glass_Metal_object_Example.txt' is an example data obtained using GBSAR-Pi [2]. More similar files are given in [3,4]. 
GBSAR-Pi works in stop-and-go mode. GBSAR sensor is FMCW (Frequency Modulated Continuous Wave) radar with central frequency of 24 GHz. 
Adjust the radar parameters such as fc (central frequency), bw (chirp bandwidth), step_size (GBSAR step size), num_of_steps (number of GBSAR steps), gbsar_direction (direction of sensor movement), tsweep (chirp sweep time), etc. In given examples central frequency is set to 24 GHz, chirp bandwidth to 1.3 GHz, step size to 0.5 cm, number of steps to 120 (meaning that the total aperture is 60 cm), and tsweep to 155 ms.





## References

[1] S. Guo and X. Dong, "Modified Omega-K algorithm for ground-based FMCW SAR imaging," 2016 IEEE 13th International Conference on Signal Processing (ICSP), Chengdu, China, 2016, pp. 1647-1650, doi: 10.1109/ICSP.2016.7878107.

[2] Kačan, M.; Turčinović, F.; Bojanjac, D.; Bosiljevac, M. Deep Learning Approach for Object Classification on Raw and Reconstructed GBSAR Data. Remote Sens. 2022, 14, 5673. https://doi.org/10.3390/rs14225673 

[3] Turčinović, Filip (2024), “Ground Based SAR Data Obtained With Different Polarizations”, Mendeley Data, V1, doi: 10.17632/nbc9xpwv96.1

[4] Turčinović, Filip (2023), “Ground Based SAR Data for Classification - 9 Objects in the Near Distance”, Mendeley Data, V1, doi: 10.17632/p2yhyx7335.1

