%% Beer-Lambert Analysis - Fluorescein, Rhodamine B, Rhodamine 6G
% Compares linear (Beer-Lambert) vs poly4 (IFE) fits
% Normalizes each spectrum by its integration time before integration
clear; clc; close all;

%% ============================================================
%  1. DATA DEFINITION - Files, concentrations, integration times
% =============================================================

% --- Fluorescein ---
FL.name   = 'Fluorescein';
FL.color  = [0.10, 0.60, 0.20];   % green
FL.files  = {'F-0_1mM.csv', 'F-0_05.csv', 'F-0_025.csv', 'F-0_01.csv', ...
             'F-0_005.csv', 'F-0_0025.csv', 'F-0_001.csv', 'F-0_0008.csv'};
FL.conc   = [0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001, 8e-4];  % mM
FL.tint   = [4000, 4000, 4000, 5000, 8000, 10000, 12000, 15000];    % ms
FL.wl_min = 460;   % nm - all datasets start at 460nm

% --- Rhodamine B ---
RB.name   = 'Rhodamine B';
RB.color  = [0.85, 0.10, 0.40];   % magenta-red
RB.files  = {'RB-0_1mM.csv', 'RB-0_05.csv', 'RB-0_025.csv', ...
             'RB-0_01.csv', 'RB-0_005.csv'};
RB.conc   = [0.1, 0.05, 0.025, 0.01, 0.005];   % mM
RB.tint   = [6000, 8000, 8000, 10000, 15000];   % ms
RB.wl_min = 460;   % nm - all datasets start at 460nm

% --- Rhodamine 6G ---
R6G.name  = 'Rhodamine 6G';
R6G.color = [1.00, 0.50, 0.00];   % orange
R6G.files = {'R6G-0_1mM.csv', 'R6G-0_05.csv', 'R6G-0_025.csv', 'R6G-0_01.csv', ...
             'R6G-0_005.csv', 'R6G-0_0025.csv', 'R6G-0_001.csv'};
R6G.conc  = [0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001];   % mM
R6G.tint  = [4000, 4000, 4000, 6000, 10000, 15000, 15000];     % ms
R6G.wl_min = 460;  % nm - all datasets start at 460nm

% Collect all materials into a cell array for looping
materials = {FL, RB, R6G};
n_mat     = length(materials);

%% ============================================================
%  2. LOAD, NORMALIZE, INTEGRATE - loop over all materials
% =============================================================

for m = 1:n_mat
    mat       = materials{m};
    n_conc    = length(mat.conc);
    S         = zeros(1, n_conc);
    wl_all    = cell(1, n_conc);   % store wavelengths for spectra plot
    int_all   = cell(1, n_conc);   % store normalized intensities for spectra plot

    for i = 1:n_conc
        data       = readmatrix(mat.files{i});
        wavelength = data(:, 1);
        intensity  = data(:, 2);

        % --- Normalize by integration time ---
        intensity_norm = intensity / mat.tint(i);   % counts / ms

        % --- Filter wavelengths above threshold ---
        idx    = wavelength >= mat.wl_min;
        w_fit  = wavelength(idx);
        i_fit  = intensity_norm(idx);

        % --- Numerical integration (trapezoid rule) ---
        S(i) = trapz(w_fit, i_fit);

        % Store full normalized spectrum for plotting
        wl_all{i}  = wavelength;
        int_all{i} = intensity_norm;
    end

    % Save back into struct
    materials{m}.S       = S;
    materials{m}.wl_all  = wl_all;
    materials{m}.int_all = int_all;
end

%% ============================================================
%  3. FITS - linear (Beer-Lambert) and poly4 (IFE) per material
% =============================================================

linear_model = fittype('a*x');   % forced zero intercept

for m = 1:n_mat
    c_col = materials{m}.conc(:);
    S_col = materials{m}.S(:);

    % Linear fit
    [lf, lg]              = fit(c_col, S_col, linear_model);
    materials{m}.lin_fit  = lf;
    materials{m}.lin_gof  = lg;

    % Poly4 fit (degree chosen = min(4, n_points-1) to avoid overfitting)
    deg = min(4, length(c_col) - 1);
    poly_type = sprintf('poly%d', deg);
    [pf, pg]              = fit(c_col, S_col, poly_type);
    materials{m}.poly_fit = pf;
    materials{m}.poly_gof = pg;
    materials{m}.poly_deg = deg;
end

%% ============================================================
%  4. PRINT GOODNESS-OF-FIT STATISTICS
% =============================================================

fprintf('=======================================================\n');
fprintf('          GOODNESS-OF-FIT SUMMARY\n');
fprintf('=======================================================\n');
for m = 1:n_mat
    mat = materials{m};
    fprintf('\n--- %s ---\n', mat.name);
    fprintf('  Linear fit  |  R^2 = %.6f  |  RMSE = %.4e\n', ...
        mat.lin_gof.rsquare, mat.lin_gof.rmse);
    fprintf('  Poly%d fit   |  R^2 = %.6f  |  RMSE = %.4e\n', ...
        mat.poly_deg, mat.poly_gof.rsquare, mat.poly_gof.rmse);
end
fprintf('\n');

%% ============================================================
%  5. FIGURE 1 - LINEAR FITS (one subplot per material)
% =============================================================

fig1 = figure('Name', 'Beer-Lambert Linear Fits', 'NumberTitle', 'off', ...
              'Position', [100, 100, 1200, 900]);

for m = 1:n_mat
    mat   = materials{m};
    c_col = mat.conc(:);
    S_col = mat.S(:);
    c_smooth = linspace(0, max(c_col)*1.15, 200)';

    subplot(3, 1, m);
    scatter(c_col, S_col, 90, 'filled', 'MarkerFaceColor', mat.color, ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.8); hold on;
    plot(c_smooth, mat.lin_fit(c_smooth), '--k', 'LineWidth', 2.0);

    xlabel('Concentration [mM]', 'FontSize', 11);
    ylabel('S [counts ms^{-1} nm^{-1}]', 'FontSize', 11);
    title(sprintf('%s — Linear Fit (Beer-Lambert)   R^2 = %.4f', ...
          mat.name, mat.lin_gof.rsquare), 'FontSize', 12);
    legend('Experimental Data', 'Linear Fit (S = a \cdot c)', 'Location', 'best');
    grid on; box on;

    % Annotate slope
    cf = coeffvalues(mat.lin_fit);
    text(0.05, 0.85, sprintf('Slope a = %.3e', cf(1)), ...
         'Units', 'normalized', 'FontSize', 10, 'Color', 'k');
end
sgtitle('Linear (Beer-Lambert) Fits — All Materials', 'FontSize', 14, 'FontWeight', 'bold');

%% ============================================================
%  6. FIGURE 2 - POLY4 FITS (one subplot per material)
% =============================================================

fig2 = figure('Name', 'Poly4 IFE Fits', 'NumberTitle', 'off', ...
              'Position', [150, 150, 1200, 900]);

for m = 1:n_mat
    mat   = materials{m};
    c_col = mat.conc(:);
    S_col = mat.S(:);
    c_smooth = linspace(0, max(c_col)*1.15, 200)';

    subplot(3, 1, m);
    scatter(c_col, S_col, 90, 'filled', 'MarkerFaceColor', mat.color, ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.8); hold on;
    plot(c_smooth, mat.poly_fit(c_smooth), '-', 'Color', mat.color*0.7, 'LineWidth', 2.0);
    plot(c_smooth, mat.lin_fit(c_smooth), '--k', 'LineWidth', 1.2);   % keep linear for reference

    xlabel('Concentration [mM]', 'FontSize', 11);
    ylabel('S [counts ms^{-1} nm^{-1}]', 'FontSize', 11);
    title(sprintf('%s — Poly%d Fit (IFE curve)   R^2 = %.4f', ...
          mat.name, mat.poly_deg, mat.poly_gof.rsquare), 'FontSize', 12);
    legend('Experimental Data', sprintf('Poly%d Fit (IFE)', mat.poly_deg), ...
           'Linear Fit (reference)', 'Location', 'best');
    grid on; box on;
end
sgtitle('Polynomial (IFE) Fits — All Materials', 'FontSize', 14, 'FontWeight', 'bold');

%% ============================================================
%  7. FIGURES 3-5 - EMISSION SPECTRA (one figure per material)
% =============================================================

for m = 1:n_mat
    mat    = materials{m};
    n_conc = length(mat.conc);

    fig_s = figure('Name', sprintf('Emission Spectra — %s', mat.name), ...
                   'NumberTitle', 'off', 'Position', [200+m*30, 200+m*30, 900, 550]);

    % Create a colormap going from light to dark in the material's hue
    cmap = zeros(n_conc, 3);
    for i = 1:n_conc
        frac      = (i - 1) / max(n_conc - 1, 1);   % 0 = lightest, 1 = darkest
        cmap(i,:) = mat.color * (0.35 + 0.65*(1-frac)) + [1,1,1]*(0.65*frac*0);
        % Simpler: interpolate from light tint to full color
        cmap(i,:) = mat.color * frac + (mat.color*0.25 + [0.75,0.75,0.75]) * (1-frac);
        cmap(i,:) = min(max(cmap(i,:), 0), 1);   % clamp to [0,1]
    end

    hold on;
    legend_entries = cell(1, n_conc);
    for i = 1:n_conc
        plot(mat.wl_all{i}, mat.int_all{i}, 'Color', cmap(i,:), 'LineWidth', 1.6);
        legend_entries{i} = sprintf('%.4g mM', mat.conc(i));
    end

    % Mark the integration cutoff wavelength
    xline(mat.wl_min, '--k', 'LineWidth', 1.2, ...
          'Label', sprintf('Integration start (%d nm)', mat.wl_min), ...
          'LabelVerticalAlignment', 'bottom');

    xlabel('Wavelength [nm]', 'FontSize', 12);
    ylabel('Normalized Intensity [counts ms^{-1}]', 'FontSize', 12);
    title(sprintf('%s — Emission Spectra (normalized by t_{int})', mat.name), ...
          'FontSize', 13, 'FontWeight', 'bold');
    legend(legend_entries, 'Location', 'best', 'FontSize', 9);
    grid on; box on;
    xlim([min(mat.wl_all{1}), max(mat.wl_all{1})]);
end

fprintf('Done. Figures generated:\n');
fprintf('  Fig 1 — Linear fits (3 subplots)\n');
fprintf('  Fig 2 — Poly4 fits  (3 subplots)\n');
fprintf('  Fig 3 — Fluorescein emission spectra\n');
fprintf('  Fig 4 — Rhodamine B emission spectra\n');
fprintf('  Fig 5 — Rhodamine 6G emission spectra\n');
