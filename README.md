# Multi-mic-AFC
This is the newest virtual test environment for PEM-AFC2 in MATLAB
# AFC Simulation

A virtual feedback environment in MATLAB R2014a or later, for PEM-AFC2 algorithm test based on HASQI sound quality evaluation

- PEM-AFC2 reference: https://ieeexplore.ieee.org/document/8270355
- HASQI source code: https://github.com/dhimnsen/OpenSpeechPlatform-UCSD
- HASQI evaluation of AFC algorithm, reference:  https://doi.org/10.1121/1.5007278

The .wav files in directory <clean_signal> are 7 different groups of signal recorded in anechoic test setup by mic1 and mic2 to simulate the clean audio without feedback. You can define the feedback path to generate echo and use AFC algorithm to cancel it.


### Run simulation

- get AFC score under constant input audio for different mu - run <Run_AFC_mutest.m>
- get AFC score under constant mu for different input audio - run <Run_AFC_audiotest.m>
- You can write the audio with echo, cancelled or not cancelled by AFC, to wav files and listen 


### Change AFC version

You can change AFC version by defining "% version" in function <AFC_processing.m> or <AFC_processing_VGSS.m>

##### Available AFC version: 

- PEM_AFC-2mic
- PEM_AFC
- AFC-2mic
- AFC

##### Available update method: 

- fix
- norm
- VGSS, ref: https://ieeexplore.ieee.org/document/7252068
- norm+VGSS


### Change test condition

You can change test condition by defining "% test parameters" and "% AFC parameters" in <Run_AFC_mutest.m> or <Run_AFC_audiotest.m>. Mention that

- For 2mic AFC test, input signal of 2 microphones should both be given and have the same length.
- Feedback path F1 and F2 must be given as two 1xN vectors.







# FIR filter Calibration

Calibration process is used for weight length estimation. Before testing in a new setup, the proper weight length should be choose with calibration. 

### Filter estimation method

Given the relationship between signal x(n) and y(n) is 

		y(n)=a0x(n)+a1x(n-1)+a2x(n-2)+...akx(n-k)
		or in z domain, y(z)=F(z)x(z)

The function <compute_filter.m> can estimate the coefficients a0-ak based on given x(n),y(n). The estimated a0-ak are based on Least Mean Square solution of the algebraic equation:

		y(n)=[x(n) x(n-1) ... x(n-k)]*[a0 a1 ... ak]^T=X(n)F(z),     n=k,k+1,....N

so that ||y_predicted - y|| is minimized under this F(z). The method to get Least Mean Square solution of an algebraic equation is "generalized inverse matrix".


### Run calibration

Get calibration result under different weight length - run <Run_calibration.m>. The .wav file in directory <signal_F> is the signals received by mic0,mic1,mic2 under hearing-aids speaking. The ones in directory <signal_H> is the signals received by mic1,mic2 under external speaker.

- The filter is estimated every signal frame to get a more general result, with user-defined frame size. Recommended frame number is 10 to 20.

- You can set the weight length as 1xN vector to find the best weight length
- You must use different audio file for F(z) and H(z) 

