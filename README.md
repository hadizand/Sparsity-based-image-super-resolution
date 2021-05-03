# Sparsity-based-image-super-resolution
we have proposed a novel method to resize input images through compressed sensing recovery. We have used super resolution in framework of compressed sensing, that is, conversion of a low resolution image to a high resolution image using compressed recovery technique.



The main code proposed by the author is in Interpolation_CS_based.m file. Just give the local path of image.

The rest of .m files, (SL0, estimate_SNR, and sparseSigGen4plusNoise) are only for fast sparse recovery. Other sparse representation methods like OMP, L1-magic and etc can be used instead of SL0. 
