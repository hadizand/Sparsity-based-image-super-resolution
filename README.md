# Sparsity-based-image-super-resolution
We have proposed a novel method to resize input images through compressed sensing recovery. We have used super resolution in framework of compressed sensing, that is, conversion of a low resolution image to a high resolution image using compressed recovery technique.



The main code proposed by the authors is in Interpolation_CS_based.m file. Just give the local path of the images used in your work.

The rest of the .m files, (SL0, estimate_SNR, and sparseSigGen4plusNoise) are only for fast sparse recovery. Other sparse representation methods like OMP, L1-magic and etc can be used instead of SL0. 
