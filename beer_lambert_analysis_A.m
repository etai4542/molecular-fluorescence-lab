%% Beer-Lambert Analysis - Fluorescein, Rhodamine B, Rhodamine 6G
% Compares global linear, dilute linear, and poly3 fits (Fixed Range Cutting)
clear; clc; close all;
set(0, 'DefaultFigureWindowStyle', 'normal');

%% ============================================================
%  0. OUTPUT FOLDER FOR SAVED FIGURES
% =============================================================
output_dir = 'Figures';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% ============================================================
%  1. DATA DEFINITION (With Explicit Dilute Limits)
% =============================================================
% --- Fluorescein ---
FL.name       = 'Fluorescein';
FL.color      = [0.10, 0.60, 0.20];
FL.files      = {'F/F-0_1mM.csv', 'F/F-0_05mM.csv', 'F/F-0_025mM.csv', 'F/F-0_01mM.csv', ...
                 'F/F-0_005mM.csv', 'F/F-0_0025mM.csv', 'F/F-0_001mM.csv', 'F/F-0_0008mM.csv'};
FL.conc       = [0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001, 8e-4];
FL.tint       = [4000, 4000, 4000, 5000, 8000, 10000, 12000, 15000];
FL.wl_min     = 460;
FL.max_linear = 0.01; % Truncating the linear range
% --- Rhodamine B ---
RB.name       = 'Rhodamine B';
RB.color      = [0.85, 0.10, 0.40];
RB.files      = {'RB/RB-0_1mM.csv', 'RB/RB-0_05mM.csv', 'RB/RB-0_025mM.csv', ...
                 'RB/RB-0_01mM.csv', 'RB/RB-0_005mM.csv'};
RB.conc       = [0.1, 0.05, 0.025, 0.01, 0.005];
RB.tint       = [6000, 8000, 8000, 10000, 15000];
RB.wl_min     = 460;
RB.max_linear = 0.025; % Truncating the linear range

% --- Rhodamine 6G ---
R6G.name       = 'Rhodamine 6G';
R6G.color      = [1.00, 0.50, 0.00];
R6G.files      = {'R6G/R6G-0_1mM.csv', 'R6G/R6G-0_05mM.csv', 'R6G/R6G-0_025mM.csv', 'R6G/R6G-0_01mM.csv', ...
                  'R6G/R6G-0_005mM.csv', 'R6G/R6G-0_0025mM.csv', 'R6G/R6G-0_001mM.csv'};
R6G.conc       = [0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001];
R6G.tint       = [4000, 4000, 4000, 6000, 10000, 15000, 15000];
R6G.wl_min     = 460;
R6G.max_linear = 0.025; % Truncating the linear range

materials = {FL, RB, R6G};
n_mat     = length(materials);

%% ============================================================
%  2. MAIN LOOP - runs twice: once normalized, once raw
% =============================================================
GENERATE_LINEAR  = true;   
GENERATE_POLY    = false;  
GENERATE_SPECTRA = false;  

spectral_colors = [
    0.894  0.102  0.110;
    0.216  0.494  0.722;
    0.302  0.686  0.290;
    1.000  0.498  0.000;
    0.596  0.306  0.639;
    1.000  1.000  0.200;
    0.651  0.337  0.157;
    0.969  0.506  0.749;
];

for norm_flag = [true, false]
    if norm_flag
        norm_label  = 'Normalized Intensity [counts ms^{-1}]';
        S_label     = 'S [counts ms^{-1} nm]';
        norm_suffix = '(normalized by t_{int})';
        file_tag    = 'norm';
    else
        norm_label  = 'Raw Intensity [counts]';
        S_label     = 'S [counts \cdot nm]';
        norm_suffix = '(raw, NOT normalized)';
        file_tag    = 'raw';
    end
    
    fprintf('\n=======================================================\n');
    fprintf('  PROCESSING: %s\n', upper(file_tag));
    fprintf('=======================================================\n');
    
    %% --- Load, process, integrate ---
    for m = 1:n_mat
        mat     = materials{m};
        n_conc  = length(mat.conc);
        S       = zeros(1, n_conc);
        wl_all  = cell(1, n_conc);
        int_all = cell(1, n_conc);
        
        for i = 1:n_conc
            data       = readmatrix(mat.files{i});
            wavelength = data(:, 1);
            intensity  = data(:, 2);
            
            if norm_flag
                intensity_proc = intensity / mat.tint(i);
            else
                intensity_proc = intensity;
            end
            
            idx  = wavelength >= mat.wl_min;
            S(i) = trapz(wavelength(idx), intensity_proc(idx));
            wl_all{i}  = wavelength;
            int_all{i} = intensity_proc;
        end
        materials{m}.S       = S;
        materials{m}.wl_all  = wl_all;
        materials{m}.int_all = int_all;
    end
    
    %% --- Fits ---
    linear_model = fittype('a*x');
    
    for m = 1:n_mat
        c_col = materials{m}.conc(:);
        S_col = materials{m}.S(:);
        
        % 1. Global Linear Fit (All Points)
        [lf, lg]              = fit(c_col, S_col, linear_model, 'StartPoint', 1);
        materials{m}.lin_fit  = lf;
        materials{m}.lin_gof  = lg;
        
        % 2. Dilute Linear Fit (Strictly Constrained Linear Range Only)
        dilute_idx = c_col <= materials{m}.max_linear;
        c_dilute   = c_col(dilute_idx);
        S_dilute   = S_col(dilute_idx);
        
        if length(c_dilute) >= 2
            [dlf, dlg] = fit(c_dilute, S_dilute, linear_model, 'StartPoint', 1);
            materials{m}.dilute_ok       = true;
            materials{m}.dilute_lin_fit  = dlf;
            materials{m}.dilute_lin_gof  = dlg;
            materials{m}.max_linear_conc = materials{m}.max_linear;
        else
            materials{m}.dilute_ok = false;
        end
        
        % 3. Polynomial Fit (Degree 3)
        deg = min(3, length(c_col) - 1);
        [pf, pg] = fit(c_col, S_col, sprintf('poly%d', deg));
        materials{m}.poly_fit = pf;
        materials{m}.poly_gof = pg;
        materials{m}.poly_deg = deg;
    end
    
    %% --- SEPARATE LINEAR PLOTS GENERATION ---
    if GENERATE_LINEAR
        
        % ---------------------------------------------------------
        % GRAPH 1: Global Linear Fits (All Prepared Concentrations)
        % ---------------------------------------------------------
        figGlobal = figure('Visible', 'off');
        for m = 1:n_mat
            mat      = materials{m};
            c_col    = mat.conc(:);
            S_col    = mat.S(:);
            c_smooth = linspace(0, max(c_col)*1.15, 200)';
            
            subplot(3, 1, m);
            scatter(c_col, S_col, 90, 'filled', 'MarkerFaceColor', mat.color, ...
                    'MarkerEdgeColor', 'k', 'LineWidth', 0.8); hold on;
            
            plot(c_smooth, mat.lin_fit(c_smooth), '-', 'Color', mat.color*0.6, 'LineWidth', 2.0);
            
            xlabel('Concentration [mM]', 'FontSize', 11);
            ylabel(S_label, 'FontSize', 11);
            title(sprintf('%s — Unconstrained Global Linear Fit', mat.name), 'FontSize', 12);
            legend('Experimental Data', sprintf('Global Linear Fit (R^2 = %.4f)', mat.lin_gof.rsquare), ...
                   'Location', 'best');
            grid on; box on;
        end
        sgtitle(sprintf('Global Beer-Lambert Linear Fits (All Points)  %s', norm_suffix), ...
                'FontSize', 14, 'FontWeight', 'bold');
        saveFigure(figGlobal, output_dir, sprintf('Fig1a_GlobalLinearFits_%s', file_tag));
        
        % ---------------------------------------------------------
        % GRAPH 2: Dilute Linear Fits (Fixed Monotonic Low Range Only)
        % ---------------------------------------------------------
        figDilute = figure('Visible', 'off');
        for m = 1:n_mat
            mat      = materials{m};
            c_col    = mat.conc(:);
            S_col    = mat.S(:);
            c_smooth = linspace(0, mat.max_linear*1.15, 200)'; % Graphing up to the dilute range limit        
            
            subplot(3, 1, m);
            % Presenting the datapoints inside the dilute linear range
            dilute_mask = c_col <= mat.max_linear;
            scatter(c_col(dilute_mask), S_col(dilute_mask), 90, 'filled', 'MarkerFaceColor', mat.color, ...
                    'MarkerEdgeColor', 'k', 'LineWidth', 0.8); hold on;
            
            if mat.dilute_ok
                plot(c_smooth, mat.dilute_lin_fit(c_smooth), '-', 'Color', mat.color, 'LineWidth', 2.2);
                legend_str = sprintf('Dilute Linear Range \\leq %g mM (R^2 = %.4f)', ...
                                     mat.max_linear_conc, mat.dilute_lin_gof.rsquare);
            else
                legend_str = 'No Dilute Range Found';
            end
            
            xlabel('Concentration [mM]', 'FontSize', 11);
            ylabel(S_label, 'FontSize', 11);
            title(sprintf('%s — Constrained Dilute Range Linear Fit', mat.name), 'FontSize', 12);
            legend('Experimental Data', legend_str, 'Location', 'best');
            grid on; box on;
            xlim([0, mat.max_linear*1.2]); % focusing the axis onto the lower linear range
        end
        sgtitle(sprintf('Dilute Range Beer-Lambert Linear Fits (Ideal Steps)  %s', norm_suffix), ...
                'FontSize', 14, 'FontWeight', 'bold');
        saveFigure(figDilute, output_dir, sprintf('Fig1b_DiluteLinearFits_%s', file_tag));
        
    end
    
    %% --- Figure 2: Polynomial fits ---
    if GENERATE_POLY
        fig = figure('Visible', 'off');
        for m = 1:n_mat
            mat      = materials{m};
            c_col    = mat.conc(:);
            S_col    = mat.S(:);
            c_smooth = linspace(0, max(c_col)*1.15, 200)';
            
            subplot(3, 1, m);
            scatter(c_col, S_col, 90, 'filled', 'MarkerFaceColor', mat.color, ...
                    'MarkerEdgeColor', 'k', 'LineWidth', 0.8); hold on;
            plot(c_smooth, mat.poly_fit(c_smooth), '-',  'Color', mat.color*0.7, 'LineWidth', 2.0);
            plot(c_smooth, mat.lin_fit(c_smooth),  '--k', 'LineWidth', 1.2);
            
            xlabel('Concentration [mM]', 'FontSize', 11);
            ylabel(S_label, 'FontSize', 11);
            title(sprintf('%s — Poly%d Fit (IFE)   R^2 = %.4f', ...
                  mat.name, mat.poly_deg, mat.poly_gof.rsquare), 'FontSize', 12);
            legend('Experimental Data', sprintf('Poly%d Fit (IFE)', mat.poly_deg), ...
                   'Global Linear Reference', 'Location', 'best');
            grid on; box on;
        end
        sgtitle(sprintf('Polynomial (IFE) Fits — All Materials  %s', norm_suffix), ...
                'FontSize', 14, 'FontWeight', 'bold');
        saveFigure(fig, output_dir, sprintf('Fig2_PolyFits_%s', file_tag));
    end
    
    %% --- Print stats ---
    for m = 1:n_mat
        mat = materials{m};
        fprintf('\n--- %s ---\n', mat.name);
        fprintf('  Global Linear Fit |  R^2 = %.6f  |  RMSE = %.4e\n', ...
            mat.lin_gof.rsquare, mat.lin_gof.rmse);
        if mat.dilute_ok
            cf_dilute = coeffvalues(mat.dilute_lin_fit);
            fprintf('  Dilute Linear Fit |  R^2 = %.6f  |  RMSE = %.4e  |  Max Cutoff = %g mM  |  Slope = %.3e\n', ...
                mat.dilute_lin_gof.rsquare, mat.dilute_lin_gof.rmse, mat.max_linear_conc, cf_dilute(1));
        end
        fprintf('  Poly%d fit         |  R^2 = %.6f  |  RMSE = %.4e\n', ...
            mat.poly_deg, mat.poly_gof.rsquare, mat.poly_gof.rmse);
    end
end
fprintf('\nAll done. Separate plots fixed to ideal dilute ranges.\n');

function saveFigure(fig, folder, filename)
    fig.Position = [100, 100, 1200, 900];
    print(fig, fullfile(folder, filename), '-dpng', '-r150');
    close(fig);
end
