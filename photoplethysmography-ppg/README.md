Event detection (time domain + SpO₂)
Implements time-domain event detection on reflective PPG data from a MAX30102 sensor to estimate heart rate and SpO₂. Detects systolic peaks using derivative thresholding and refractory periods, computes interbeat intervals and HRV, and calculates SpO₂ via AC/DC ratio-of-ratios.

Frequency-domain analysis (FFT)
Estimates heart rate from reflective PPG signals using frequency-domain analysis. Removes DC offset, applies windowing and FFT, and identifies the dominant spectral peak within a physiological heart-rate band. Demonstrates robustness to waveform variability and noise compared to time-domain methods.

Autocorrelation-based heart rate
Computes heart rate from reflective PPG signals using autocorrelation. Removes DC offset and identifies the dominant lag corresponding to the cardiac period within physiologically valid bounds. This approach emphasizes waveform periodicity and demonstrates improved robustness to amplitude variation and noise.
