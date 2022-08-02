This is the code to reproduce the fusion algorithms reported in 

[1] Mitianoudis N., Stathaki T., " Pixel-based and Region-based Image Fusion schemes using ICA bases", Information Fusion 8 (2),  pp. 131-142, April 2007.
[2] Mitianoudis N., Stathaki T., "Optimal Contrast Correction for  ICA-based Fusion of Multimodal Images" , IEEE Sensors Journal, Vol. 8, No. 12, pp. 2016 - 2026, Dec. 2008. 
[3] Mitianoudis N., Antonopoulos S.A, Stathaki T., "Region-based ICA Image Fusion using Textural Information ", 18th Int. Conf on DSP (DSP2013), July 1-3, 2013, Santorini, Greece. 

Remarks
-------

1. Please refer to [1], if you are using the basic ICA fusion schemes and
refer to [2], if you are using the optimal contrast scheme.

2. Please refer to [3], if you are using the regional ICA fusion scheme.

3. If you need to train new ICA bases, please put the images you wish to train in the folder 'training_images' and run ICA_bases_train.m. You can use already trained bases in 'icabases.mat'

4. The main fusion function is 'ica_fusion_prog_mn.m'. You can learn its functionality through the test file 'ICA_fusion_test.m', where you can find several fusion examples.

If you have queries, concerning the code, don't hesitate to write to nmitiano@ee.duth.gr



Nikolaos Mitianoudis, 2016.