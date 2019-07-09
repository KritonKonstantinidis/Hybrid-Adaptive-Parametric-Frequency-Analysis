# Hybrid-Autoregressive-Toolbox
Toolbox for parametric analysis of non-stationary univariate time series data
A description of how to use the routines is given below.
For the mathematical models, see pdf files in master: My MSc. thesis and the paper published in IEEE BIBE 2019 related to this work.

# Normal Hybrid Kalman Filtering
Call HybridAlgorithm. In order to generate the spectrogram of a time-series, the user has to just input the time series and run this routine. The output is the time-varying spectral estimation of the time-series. Output also inlcudes metrics including AIC scores for the different model orders and EM algorithm convergence. Finally, the temporal evolution of the autoregressive coefficients is plotted.

The user has the liberty to modify the following parameters:

Exponential smoothing of AR coefficients (default:Off) Sampling rate (default: 250Hz)
Frequency resolution (default: 4 index/Hz)
Maximum of frequency range (default: 50 Hz)
Minimum of frequency range (default: 0.0001 Hz)
Maximum number of EM algorithm iterations (default:50)
Batch size to run EM on (default: 10 seconds)
Tolerance for EM convergence (default: 10−3)
Initial observation noise variance for EM iterations (default: 0.2)
Starting coefficient of diagonal covariance matrix for EM iterations (default: 1) Minimum AR order to be considered in model selection (default: 2)
Maximum AR order to be considered in model selection (default: 20) Information criterion (default: AIC)
Remove outliers or not (default: On)
Run online EM or not (default: Off)
Buffer size in case of online EM (default 5 min)

# Cauchy Hybrid Sequential Monte Carlo Filtering
Run SO_SMC_Hybrid_Spectrogram. Generates a spectrogram of the input data based on sequential importance resampling with adaptive estimation of the scale parameters.
Input Argumetns 
y: Univariate Time-series EEG data
p: Autoregressive model order
N: Number of particles for the SMC
R: Variance of Gaussian observation noise
c: Scale parameter of the Cauchy distribution of the AR parameters for SMC_Hybrid_Spectrogram a,b: interval to uniformly initialize the scale parameters
delta: Number of discretization steps between the observations
B: Number of Bootstrap realizations

In case of MCMC routine, the following 2 arguments are added.
crandom: Variance of the local random walk
L: Number of steps of MCMC at each time point

# EMD + HAR(2+2) for propofol
Run EMD_Hybrid_Spectrogram. This function takes as input the raw EEG time series and outputs the calculated spectrogram.
