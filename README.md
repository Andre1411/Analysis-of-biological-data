# Analysis of biological data
This folder contains the code and results of two individual projects required for the course "Analysis of biomedical data" at the University of Padova:
- Bin smoothing.
- Local Field Potential.

More details are provided in the respective sections.

# HOMEWORK 1 - Bin smoothing

## Text

Consider the test signal A. Apply smoothing to it on a virtual uniform grid at 1-minute intervals (1:1:ts(end)) by invoking a function implementing the bin-smoother discussed in class. The function should take the data vector ys, the sampling time vector ts, the virtual grid tv, the number of bins (assumed to have constant support) as input, and return the smoothed vector uhat (relative to tv) and the residuals vector res. Illustrate the results (fit and residuals) for three exemplary bin numbers representing typical situations: undersmoothing, oversmoothing, and reasonable smoothing. Experimentally find the bin number p that approximately satisfies the discrepancy criterion and display the corresponding profile graph (fit and residuals). Calculate the norm of the absolute estimation error for various p values and check if the minimum error occurs for the recommended smoothing criterion.

## Solution

### Introduction and Methodology

To implement the bin-smoother, the time axis needs to be divided into intervals equal to the number of bins. For each interval $R_k$, the smoother value is the average of the signal samples collected at times ti falling within that interval:

$$ u_k = mean \[ y(t_i) \] $$

To find the bin number satisfying the discrepancy criterion, an iterative procedure using a bisection method is employed. It starts with $p_{min}=1$ and $p_{max}=200$, calculates the optimal value as the integer closest to the average of the two and updates the bounds based on the value of the sum of squared absolute residuals.

### Code Presentation

Two code files are provided:

- `bin_smoother.m`: A function implementing the bin-smoother with input parameters as described in the prompt. It calculates the smoother value in a bin by averaging the values of samples between the current bin's upper bound and the upper bound of the previous bin. It computes the smoother value at the time instants on the virtual grid where the signal is sampled and uses them to calculate the residuals.

- `h1_31_main.m`: Initializes the virtual grid, calculates three examples of smoothing (over, under, and reasonable) using the above function, and plots them. Calculates the optimal bin number according to the discrepancy criterion: starting from $p_{min}=1$ and $p_{max}=200$, it computes p as

$$ p = \text{round}\left(10^{\frac{\log_{10}(p_{min}) + \log_{10}(p_{max})}{2}}\right) $$

and uses it to build the smoother and residuals, iterating until the ratio

$$ \bigg|\frac{\text{ARSS} - \text{ns} \times \text{sd}^2}{\text{ns} \times \text{sd}^2}\bigg| $$

becomes less than a tolerance value. The $p_{max}$ bound is updated to the calculated value of $p$ if $ARSS < ns \times \text{sd}^2$ otherwise, the bound $p_{min}$ is updated. At the end of the loop, the optimal bin number is calculated using a different method (to ensure the two methods match) that explores all bin values between 1 and 200. For each value, it calculates residuals and ARSS, compares each ARSS value to the target value $ns \times \sigma^2$, and selects the optimal nbin as the one deviating least from that value. Plots the fit and residuals for the optimal value and the ARSS trend as the number of bins varies.

### Results and Discussion

The graphs show:
- **Figure 1:** The value of the sum of squared absolute residuals (ARSS) as the number of bins varies.<br />
  ![h1_arss_neg](https://github.com/Andre1411/Analysis-of-biological-data/assets/107708093/6e6e6c4f-ce7c-47c7-a326-eb912ba24f89)
- **Figure 2:** Oversmoothing case: too few bins do not capture the data trend, and two distinct peaks are not recognized. Residuals exhibit clear trends.<br />
  ![h1_oversmoothing_neg](https://github.com/Andre1411/Analysis-of-biological-data/assets/107708093/a1bf313f-ea62-4a91-b67b-fd32c7af580f)
- **Figure 3:** Undersmoothing case: too many bins, and the smoother follows the noise. Residuals are uncorrelated but with reduced variability.<br />
  ![h1_undersmoothing_neg](https://github.com/Andre1411/Analysis-of-biological-data/assets/107708093/9af15ac1-3c0b-410f-bd1c-5d2fbc7c198d)
- **Figure 4:** Reasonable smoothing: the smoother follows the noise less (approximately constant at the ends) but can recognize two distinct peaks. Residuals are uncorrelated and have more variable amplitude than the previous case.<br />
  ![h1_oksmoothing_neg](https://github.com/Andre1411/Analysis-of-biological-data/assets/107708093/8808f79a-f4df-4d3d-a685-414cc936ded5)
- **Figure 5:** Discrepancy: optimal smoothing, even though it slightly follows the noise at the ends.<br />
  ![h1_bestsmoothing_neg](https://github.com/Andre1411/Analysis-of-biological-data/assets/107708093/da280988-8276-4591-964d-74dcbf9752a3)

The minimum value of the norm of the absolute estimation error is not achieved for the bin number found using the discrepancy criterion. In such cases, the minimum is reached using 182 bins or more, as seen from the graph. In these instances, only one sample falls in each interval, so the smoother value corresponds to the sample value itself, leading to overfitting, undersmoothing, and a zero ARSS.

# HOMEWORK 2 - Local Field Potential

## Text

The file data_LFP.mat contains a noisy recording of a Local Field Potential measured in the primary sensory area of a rat during an experiment, after a whisker stimulus, to study the mechanisms of processing tactile stimuli in the telencephalic cortex. The data vector contains samples in mV, and the time vector contains the corresponding sampling times in ms. It is known that the noise standard deviation can be assumed to be constant, σ=0.005. From electrophysiology, it is expected that, with time 0 corresponding to the stimulus, the noise-free LFP profile shows a very pronounced peak, inflection point, and trough in the time interval [5 25] ms.

Create a script that, using a single regularized deconvolution procedure implemented with appropriate numerical algorithms, estimates and shows:

- Regularized second derivative of the LFP
- Regularized first derivative of the LFP
- Regularized LFP
- Residues

Using the two regularized derivatives obtained, write code that then identifies the coordinates (time, amplitude) of the first absolute peak, the next absolute trough, and the inflection point between them.

**Hints:** Set up the deconvolution problem by observing that a causal signal can be interpreted as the double integration of its second derivative, and that, in modeling terms, a double integration of a signal can be interpreted as the output of an LTI system with its input being the second derivative of the signal and its impulse response being the ramp function. The regularized LFP and latencies/amplitudes should resemble what is reported in the figure. [...]

## Solution

### Introduction and Methodology

Estimates of the regularized LFP and its regularized derivatives were developed according to the suggestion. An approach of deconvolution of a signal viewed as the output of a system performing double integration is exploited: the input is, therefore, the second derivative of the output signal. The transfer matrix G is obtained as a Toeplitz matrix and is lower triangular with the same elements along the diagonals: discretizations of the integral of the impulse response, which in this case is a ramp function. The energy of the second derivative is minimized through the F matrix. To obtain the first derivative, a second transfer matrix G1 is constructed, not used in the regularized estimation, which implements the single integration of the second derivative. To optimize computational cost, a change of coordinates was implemented by diagonalizing the problem through single value decomposition, allowing more efficient search for gamma before returning to the original coordinates.

### Code Presentation

There are two code files:

- `consistency.m`: implements the procedure of regularized deconvolution via a consistency criterion, optimizing the algorithm using diagonalization through a change of coordinates via svd. It has the following input arguments:

  - sigma2: scaling factor of matrix B
  - B: weight matrix (sigma2*B = covariance matrix of measurement error)
  - F: regularization matrix
  - g_max and g_min: search bounds for the regularization parameter gamma
  - ys: measurements (system output affected by noise)
  - ns and nv: length of sampling and virtual grids
  - crit: number of the consistency criterion to be implemented

  Defines the matrix $B^{-1/2}$ used later for the estimation calculation. Decomposes the matrix $B^{-1/2}⋅G⋅F^{-1}$ into matrices $U$, $D$, and $V$ through single value decomposition. Converts the measurements $ys$ into new coordinates (psi vector) Implements the three consistency criteria, trying to achieve equality:

  - criterion 1: WESS(γ) = σ^2⋅q(γ)/γ
  - criterion 2: WRSS(γ) = σ^2⋅(n-q(γ))
  - criterion 3: γ⋅WESS(γ)/q(γ) = WRSS(γ)/(n-q(γ))

  calculating the gamma value by averaging the values of the two search bounds that are updated based on the difference from the target value. Stops when this difference is less than a tolerance value. Returns the input vector to the original coordinates by implementing the product $F^{-1}⋅V⋅$ ni_hat as filtering the $V⋅$ ni_hat vector through a filter that implements the second derivative (second difference in discrete). Computes the regularized output and the residuals with respect to the initial measurements ys.

- `quinto_h2_11_main.m`: loads the data, defines the virtual grid, the degree of the derivative to minimize energy (m=2), and matrix $B$. Constructs matrix $F$ for the implementation of this minimization. Builds the Toeplitz matrix $G$ as a discretization of the integral of the impulse response of the system, which is a ramp since the system implements double integration.

Defines a second matrix $G1$ obtained from the discretization of the integral of the impulse response of a system implementing single integration. This matrix will be used for the calculation of the first derivative.

Uses the function `consistency.m` to compute the regularized output and input and the residuals. Calculates the regularized output on the virtual grid y_reg, the first derivative dy_dt as the product of the regularized output and the matrix $Gv$, and the second derivative d2y_dt2 as the reconstructed input.

Searches for the minimum and maximum values as points where the first derivative nullifies in the interval [5 25], knowing that the maximum should occur first and then the minimum. Uses the maximum and minimum values as extremes within which to find the inflection point, where the second derivative nullifies.


## Results and Discussion

The graphs, obtained using the consistency criterion 2, show:

- **Figure 1:** First derivative of the regularized LFP, search extremes, zero level, and nullification points within the extremes, corresponding to the maximum and minimum points of the LFP.<br />
![h2_fig_1_neg](https://github.com/Andre1411/Analysis-of-biological-data/assets/107708093/a0775161-33db-4166-86cc-e8663f87e56d)
- **Figure 2:** Second derivative of the regularized LFP, search extremes constituted by the maximum and minimum points found with the first derivative, zero level, and nullification point corresponding to the inflection point of the LFP.<br />
![h2_fig_2_neg](https://github.com/Andre1411/Analysis-of-biological-data/assets/107708093/f2af646e-59b7-4c75-94ed-73c03c30f7d1)
- **Figure 3:** Regularized LFP with respective points of maximum, inflection, and minimum.<br />
![h2_fig_3_neg](https://github.com/Andre1411/Analysis-of-biological-data/assets/107708093/c8cdb805-5fda-40c4-8667-bc91720938d8)
- **Figure 4:** Comparison between measured and regularized LFP at the top and respective absolute residuals at the bottom.<br />
![h2_fig_4_neg](https://github.com/Andre1411/Analysis-of-biological-data/assets/107708093/897151d3-8916-4967-8863-e88d135ea2e1)

The residuals exhibit trends attributable to the high frequency at which the signal is sampled. This phenomenon is particularly noticeable in the initial part of the signal, where the most abrupt change occurs (negative peak of the first derivative). Despite this, the estimate appears not to follow the noise in the data, and the first and second derivatives allow identifying maximum, minimum, and inflection points, reflecting physiological expectations.


