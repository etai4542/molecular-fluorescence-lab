# Quantitative Molecular Fluorescence Lab Analysis

This repository contains a comprehensive MATLAB-based data processing and regression analysis suite designed to evaluate the operational boundaries of the Beer-Lambert law in molecular fluorescence. The toolkit processes experimental broadband emission spectra and spatial digital imaging data across various concentration series for three organic fluorophores: **Fluorescein**, **Rhodamine B**, and **Rhodamine 6G**.

The project is structured into two independent working directories representing the two main analytical tracks of the investigation.

---

## Directory Workspace Structure

To ensure the automated execution of the scripts, organize the data paths, raw spectra, and uploaded beam photographs according to the following layout matching your MATLAB Drive environment:

```text
molecular-fluorescence-lab/
├── Part A/                         # Isotropic Spectral Analysis Track
│   ├── beer_lambert_analysis_A.m   # Main processing and fitting script for Part A
│   ├── F/                          # Fluorescein raw spectral data (.csv files)
│   ├── RB/                         # Rhodamine B raw spectral data (.csv files)
│   ├── R6G/                        # Rhodamine 6G raw spectral data (.csv files)
│   └── Figures/                    # Automatically generated output folder for Part A plots
│
└── Part B/                         # Spatial Beam Attenuation Track
    ├── beer_lambert_analysis_B.m   # Main decay profiling and imaging script for Part B
    ├── Fluorescein/                # Upload path for Fluorescein beam photos (e.g., F-0_1.jpg)
    ├── Rhodamine B/                # Upload path for Rhodamine B beam photos (e.g., RB-0_05.jpg)
    ├── Rhodamine 6G/               # Upload path for Rhodamine 6G beam photos (e.g., R6G-0_01.jpg)
    └── Figures_PartB/              # Automatically generated output folder for Part B plots
```

---

## Technical Specifications

### Part A: Spectral Analysis Framework (`/Part A/beer_lambert_analysis_A.m`)

#### Inputs & Core Arrays
The script dynamically processes structured streams of matrix data from `.csv` files organized by fluorophore type subdirectories. For each material, the workspace constructs three corresponding vector components:
* **Concentration Vector ($c$):** The complete series of sample concentrations prepared for observation.
* **Integration Time Vector ($t_{\text{int}}$):** Tracks the manually adjusted exposure durations utilized during detector collection to prevent photon clipping or underexposure.
* **Spectral Truncation Guard ($wl_{\text{min}}$):** Defines a strict lower-bound cutoff at $460\text{ nm}$ to filter out scattering residue from the light source and isolate the true emission profiles.

#### Computational Flow & Loops
The execution logic processes data sequentially under a main double-run configuration controlled via a boolean framework:
1. **Raw Intensity Loop:** Analyzes data streams directly using raw, uncorrected photon count readings.
2. **Temporally Normalized Intensity Loop:** Scales the integrated intensity entries by their active exposure parameter ($t_{\text{int}}$), yielding consistent units of $\text{counts}\cdot\text{ms}^{-1}\cdot\text{nm}$.
3. **Numerical Integration:** Applies a trapezoidal integration rule (`trapz`) across the active emission boundaries to evaluate the cumulative experimental photon yield parameter, denoted as $S$.
4. **Data-Pair Mapping:** Compiles the computed parameters into unified experimental $(c, S)$ coordinate pairs.

#### Regression & Fitting Architecture
The compiled dataset is fitted across three distinct analytical tracks:
* **Global Linear Fit:** An unconstrained proportional regression model ($S = a \cdot c$) forced through the origin, evaluated across the entire concentration range ($0$ to $0.1\text{ mM}$).
* **Dilute Linear Fit:** A constrained proportional model that filters out high-concentration boundaries to isolate the ideal linear bounds ($c \le 0.01\text{ mM}$ for Fluorescein, $c \le 0.025\text{ mM}$ for the Rhodamine series) where the signal follows a monotonic upward scaling trend.
* **3rd-Degree Polynomial Fit:** An empirical framework matching up to degree 3, designed to track non-linear trends, sub-linear curvature, and high-concentration saturation profiles.

The script outputs comprehensive goodness-of-fit benchmarks ($R^2$, $\text{RMSE}$) for all approaches and explicitly extracts the calibration slope ($a$) from the dilute range framework.

---

### Part B: Spatial Beam Attenuation Framework (`/Part B/beer_lambert_analysis_B.m`)

#### Inputs & Visual Interface
The spatial attenuation pipeline utilizes digital photography tracks captured across an extended $10\text{-cm}$ physical cell container. The imaging data is read using native matrix parsing routines:
```matlab
A = imread('Rhodamine B/RB-0_05.jpg');
figure; imagesc(im2double(A(:,:,1))); colorbar; axis image;
```
*Note: Users must upload their raw beam attenuation photographs into the respective material folders utilizing the localized format notation (e.g., `F-[conc].jpg`, `RB-[conc].jpg`, `R6G-[conc].jpg`).*

To optimize signal-to-noise separation, the framework isolates targeted R/G/B intensity matrices depending on the respective fluorescence color profile of the target dye (e.g., extracting channel 1 for the red-emitting Rhodamine compounds or channel 2 for Fluorescein).

#### Coordinate Profiling & Path Definition
Using interactive data visualization panels, the coordinate limits are defined to isolate physical boundaries:
* $x_{\text{start}}$ and $x_{\text{end}}$ coordinates register the internal entrance and exit boundaries of the liquid container.
* The $y$ coordinate isolates the specific pixel row alignment corresponding to the center of the propagating excitation beam path.

#### Spatial Decay & Molar Absorption Extraction
The pixel indices are converted into a continuous metric path ($x$ in $\text{cm}$) to model intensity profile trends:
1. **Logarithmic Evaluation:** The spatial intensity track is transformed to its natural logarithm ($\ln I$) to linearize the continuous spatial decay function.
2. **Local Attenuation Extraction:** The script solves for the continuous natural spatial attenuation coefficient ($\alpha$, in $\text{cm}^{-1}$) across the clean propagation profile using the downward sloping relationship:
   $$\ln(I) = -\alpha \cdot x$$
3. **Molar Absorption Mapping:** The extracted discrete $\alpha$ values are plotted against the sample concentration vector ($c$). The framework evaluates the data using both a forced-zero linear regression based on ideal conventions:
   $$\alpha = \ln(10) \cdot \varepsilon_{\lambda_{\text{ex}}} \cdot c \approx 2.303 \cdot \varepsilon_{\lambda_{\text{ex}}} \cdot c$$
   to directly isolate the base-10 molar absorption coefficient $\varepsilon_{\lambda_{\text{ex}}}$, and a comparative 3rd-degree polynomial model to capture high-concentration sub-linear saturation trends.

#### Hardware Performance Boundaries
The script explicitly identifies and accounts for two systematic hardware artifacts that distort raw spatial tracking:
* **Pixel Clipping (Saturation):** Occurs at highly concentrated entry zones ($x = 0$ to $1.5\text{ cm}$), where the intense emission hits the digital collection ceiling of the camera sensor.
* **Sensor Noise Floor:** Occurs at downstream locations when the beam undergoes near-complete attenuation, flattening the logarithmic decay profile into a baseline floor dominated by pixel granularity.

---

## Runtime Execution & Controls (Part A)

The computational execution of the spectral track is driven by boolean flags located at the initialization block of `beer_lambert_analysis_A.m`:

```matlab
norm_flag        = true;   % Enables/disables exposure time normalization (counts/ms)
GENERATE_LINEAR  = true;   % Exports independent global and dilute range linear plots (Fig1a, Fig1b)
GENERATE_POLY    = true;   % Generates comparative 3rd-degree polynomial curves (Fig2)
GENERATE_SPECTRA = false;  % Saves standalone raw and normalized emission spectra panels
```

### Generated File Outputs

#### Part A Outputs (Saved automatically into `/Part A/Figures/`)
* `Fig1a_GlobalLinearFits_norm.png` — Visual breakdown of full-range linear evaluation for total emission.
* `Fig1b_DiluteLinearFits_norm.png` — Proportional linear verification inside the isolated, optically dilute bounds.
* `Fig2_PolyFits_norm.png` — High-precision empirical tracking of non-linear integrated fluorescence curves and peak suppression.

#### Part B Outputs (Saved automatically into `/Part B/Figures_PartB/`)
* Analytical spatial decay graphs ($\ln I$ versus distance) for individual concentrations.
* Comprehensive $\alpha$ versus concentration summary regression panels (Linear and Polynomial models).
```
