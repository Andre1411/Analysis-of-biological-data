# Analysis of biological data
This folder contains the code and results of two individual projects required for the course "Analysis of biomedical data" at the University of Padova:
- Bin smoothing.
- Local Field Potential.

More details are provided in the respective sections.

# HOMEWORK 1 - Bin smoothing

## Exercise H1 31

Consider the test signal A. Apply smoothing to it on a virtual uniform grid at 1-minute intervals (1:1:ts(end)) by invoking a function implementing the bin-smoother discussed in class. The function should take the data vector ys, the sampling time vector ts, the virtual grid tv, the number of bins (assumed to have constant support) as input, and return the smoothed vector uhat (relative to tv) and the residuals vector res. Illustrate the results (fit and residuals) for three exemplary bin numbers representing typical situations: undersmoothing, oversmoothing, and reasonable smoothing. Experimentally find the bin number p that approximately satisfies the discrepancy criterion and display the corresponding profile graph (fit and residuals). Calculate the norm of the absolute estimation error for various p values and check if the minimum error occurs for the recommended smoothing criterion.

## Solution

### Introduction and Methodology

To implement the bin-smoother, the time axis needs to be divided into intervals equal to the number of bins. For each interval R$_k$, the smoother value is the average of the signal samples collected at times ti falling within that interval:

$$ u_k = mean \[ y(t_i) \] $$

To find the bin number satisfying the discrepancy criterion, an iterative procedure using a bisection method is employed. It starts with $p_{min}=1$ and $p_{max}=200$, calculates the optimal value as the integer closest to the average of the two and updates the bounds based on the value of the sum of squared absolute residuals.

### Code Presentation

Two code files are provided:

- **bin_smoother.m:** A function implementing the bin-smoother with input parameters as described in the prompt. It calculates the smoother value in a bin by averaging the values of samples between the current bin's upper bound and the upper bound of the previous bin. It computes the smoother value at the time instants on the virtual grid where the signal is sampled and uses them to calculate the residuals.

- **h1_31_main.m:** Initializes the virtual grid, calculates three examples of smoothing (over, under, and reasonable) using the above function, and plots them. Calculates the optimal bin number according to the discrepancy criterion: starting from $p_{min}=1$ and $p_{max}=200$, it computes p as

$$ p = \text{round}\left(10^{\frac{\log_{10}(p_{min}) + \log_{10}(p_{max})}{2}}\right) $$

and uses it to build the smoother and residuals, iterating until the ratio

$$ \bigg|\frac{\text{ARSS} - \text{ns} \times \text{sd}^2}{\text{ns} \times \text{sd}^2}\bigg| $$

becomes less than a tolerance value. The $p_{max}$ bound is updated to the calculated value of $p$ if $ARSS < ns \times \text{sd}^2$ otherwise, the bound $p_{min}$ is updated. At the end of the loop, the optimal bin number is calculated using a different method (to ensure the two methods match) that explores all bin values between 1 and 200. For each value, it calculates residuals and ARSS, compares each ARSS value to the target value $ns \times \sigma^2$, and selects the optimal nbin as the one deviating least from that value. Plots the fit and residuals for the optimal value and the ARSS trend as the number of bins varies.

### Results and Discussion

The graphs show:
- **Figure 1:** The value of the sum of squared absolute residuals (ARSS) as the number of bins varies.![h1_arss_neg](https://github.com/Andre1411/Analysis-of-biological-data/assets/107708093/6e6e6c4f-ce7c-47c7-a326-eb912ba24f89)
- **Figure 2:** Oversmoothing case: too few bins do not capture the data trend, and two distinct peaks are not recognized. Residuals exhibit clear trends.![h1_oversmoothing_neg](https://github.com/Andre1411/Analysis-of-biological-data/assets/107708093/a1bf313f-ea62-4a91-b67b-fd32c7af580f)
- **Figure 3:** Undersmoothing case: too many bins, and the smoother follows the noise. Residuals are uncorrelated but with reduced variability.![h1_undersmoothing_neg](https://github.com/Andre1411/Analysis-of-biological-data/assets/107708093/9af15ac1-3c0b-410f-bd1c-5d2fbc7c198d)
- **Figure 4:** Reasonable smoothing: the smoother follows the noise less (approximately constant at the ends) but can recognize two distinct peaks. Residuals are uncorrelated and have more variable amplitude than the previous case.![h1_oksmoothing_neg](https://github.com/Andre1411/Analysis-of-biological-data/assets/107708093/8808f79a-f4df-4d3d-a685-414cc936ded5)
- **Figure 5:** Discrepancy: optimal smoothing, even though it slightly follows the noise at the ends.![h1_bestsmoothing_neg](https://github.com/Andre1411/Analysis-of-biological-data/assets/107708093/da280988-8276-4591-964d-74dcbf9752a3)

The minimum value of the norm of the absolute estimation error is not achieved for the bin number found using the discrepancy criterion. In such cases, the minimum is reached using 182 bins or more, as seen from the graph. In these instances, only one sample falls in each interval, so the smoother value corresponds to the sample value itself, leading to overfitting, undersmoothing, and a zero ARSS.
