%% Part B - Spatial Beer-Lambert Analysis
% Extracts attenuation coefficient alpha from beam decay photos
% Alpha = -slope of ln(I) vs x
% Then fits alpha vs concentration with linear and polynomial models (cf. Part A)
clear; clc; close all;
set(0, 'DefaultFigureWindowStyle', 'normal');

output_dir = 'Figures_PartB';
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

%% ============================================================
%  1. DATA DEFINITION - coordinates found via ginput
% ============================================================

% --- Fluorescein (channel 2 = green) ---
FL.name    = 'Fluorescein';
FL.color   = [0.10, 0.60, 0.20];
FL.channel = 2;
FL.files   = {'Fluorescein/F-0_1.jpg',   'Fluorescein/F-0_05.jpg', ...
              'Fluorescein/F-0_025.jpg',  'Fluorescein/F-0_01.jpg', ...
              'Fluorescein/F-0_005.jpg',  'Fluorescein/F-0_0025.jpg', ...
              'Fluorescein/F-0_001.jpg',  'Fluorescein/F-0_0008.jpg', ...
              'Fluorescein/F-0_0001.jpg'};
FL.conc    = [0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001, 0.0008, 0.0001];
FL.row     = [1889, 1656, 1758, 1828, 1470, 2391, 2028, 2096, 1962];
FL.x_start = [ 506,  530,  544,  339,  246,  418,  569,  472,  497];
FL.x_end   = [1955, 2011, 1922, 1900, 1853, 2470, 2079, 2024, 1850];

% --- Rhodamine 6G (channel 1 = red/orange) ---
R6G.name    = 'Rhodamine 6G';
R6G.color   = [1.00, 0.50, 0.00];
R6G.channel = 1;
R6G.files   = {'Rhodamine 6G/R6G-0_1.jpg',   'Rhodamine 6G/R6G-0_05.jpg', ...
               'Rhodamine 6G/R6G-0_01.jpg',   'Rhodamine 6G/R6G-0_005.jpg', ...
               'Rhodamine 6G/R6G-0_0025.jpg'};
R6G.conc    = [0.1, 0.05, 0.01, 0.005, 0.0025];
R6G.row     = [1850, 2055, 2001, 2080, 2053];
R6G.x_start = [ 380,  522,  498,  693,  467];
R6G.x_end   = [2533, 2298, 2246, 2348, 2280];

% --- Rhodamine B (channel 1 = red) ---
RB.name    = 'Rhodamine B';
RB.color   = [0.85, 0.10, 0.40];
RB.channel = 1;
RB.files   = {'Rhodamine B/RB-0_1.jpg',   'Rhodamine B/RB-0_05.jpg', ...
              'Rhodamine B/RB-0_025.jpg',  'Rhodamine B/RB-0_01.jpg', ...
              'Rhodamine B/RB-0_005.jpg'};
RB.conc    = [0.1, 0.05, 0.025, 0.01, 0.005];
RB.row     = [1831, 1969, 1894, 1998, 2072];
RB.x_start = [ 494,  378,  338,  448,  328];
RB.x_end   = [2622, 2189, 2258, 2132, 1964];

materials = {FL, R6G, RB};
n_mat     = length(materials);

%% ============================================================
%  2. EXTRACT ALPHA FROM EACH IMAGE + SAVE ln(I) vs x PLOTS
% ============================================================

for m = 1:n_mat
    mat   = materials{m};
    n     = length(mat.files);
    alpha = zeros(1, n);

    fprintf('\n--- %s ---\n', mat.name);

    for i = 1:n
        % Load correct channel
        A  = imread(mat.files{i});
        A1 = im2double(A);
        Av = double(A1(mat.row(i), mat.x_start(i):mat.x_end(i), mat.channel));

        % x axis in cm
        x = linspace(0, 10, length(Av));

        % Smooth and take log
        Av_smooth  = movmean(Av, 50);
        Avl_smooth = log(Av_smooth);

        % Linear fit to log profile (exclude non-finite values)
        valid = isfinite(Avl_smooth) & Av_smooth > 0;
        if sum(valid) > 2
            p        = polyfit(x(valid), Avl_smooth(valid), 1);
            alpha(i) = -p(1);
        else
            alpha(i) = NaN;
        end

        fprintf('  %s  ->  alpha = %.4f cm^-1\n', mat.files{i}, alpha(i));

        % --- ln(I) vs x figure (instructor output) ---
        fig = figure('Visible', 'off');
        scatter(x, log(Av), 5, 'filled', ...
                'MarkerFaceColor', mat.color, 'MarkerFaceAlpha', 0.3); hold on;
        if sum(valid) > 2
            plot(x, polyval(p, x), '-k', 'LineWidth', 2);
        end
        xlabel('x [cm]', 'FontSize', 12);
        ylabel('ln(Intensity) [AU]', 'FontSize', 12);
        title(sprintf('%s  %g mM — ln(I) vs x   \\alpha = %.4f cm^{-1}', ...
              mat.name, mat.conc(i), alpha(i)), 'FontSize', 12);
        legend('Data points', sprintf('Linear fit  \\alpha = %.4f cm^{-1}', alpha(i)), ...
               'Location', 'best');
        grid minor; box on;
        fname = strrep(mat.files{i}, '/', '_');
        fname = strrep(fname, '.jpg', '');
        saveFigure(fig, output_dir, sprintf('logprofile_%s', fname));
    end

    materials{m}.alpha = alpha;
end

%% ============================================================
%  3. PRINT GOODNESS-OF-FIT SUMMARY
% ============================================================
linear_model = fittype('a*x');   % forced zero intercept, same as Part A
fprintf('\n=======================================================\n');
fprintf('       ALPHA vs CONCENTRATION — FIT SUMMARY\n');
fprintf('=======================================================\n');
for m = 1:n_mat
    mat   = materials{m};
    valid = isfinite(mat.alpha);
    c_col = mat.conc(valid)';
    a_col = mat.alpha(valid)';
    [lf, lg] = fit(c_col, a_col, linear_model, 'StartPoint', 1);
    deg      = min(3, length(c_col) - 1);
    [pf, pg] = fit(c_col, a_col, sprintf('poly%d', deg));
    
    % --- שלב תיקון הבסיס למציאת ה-epsilon של הרקע התיאורטי ---
    cf = coeffvalues(lf);
    empirical_slope = cf(1); 
    decadic_epsilon = empirical_slope / 2.302585; % diving by 2.303 to transition to a base-10 logarithm
    
    fprintf('\n--- %s ---\n', mat.name);
    fprintf('  Linear fit  |  epsilon (base 10) = %.4f cm^-1 mM^-1  |  R^2 = %.6f  |  RMSE = %.4e\n', ...
            decadic_epsilon, lg.rsquare, lg.rmse);
    fprintf('  Poly%d fit   |  R^2 = %.6f  |  RMSE = %.4e\n', deg, pg.rsquare, pg.rmse);
    
    % Store fits for plotting
    materials{m}.lin_fit  = lf;       % Keeps the original fit in the natural range
    materials{m}.lin_gof  = lg;
    materials{m}.poly_fit = pf;
    materials{m}.poly_gof = pg;
    materials{m}.poly_deg = deg;
    materials{m}.epsilon  = decadic_epsilon; % Stores the epsolon values in base-10 logarithm
end

%% ============================================================
%  4. FIGURE — alpha vs concentration, linear fit
% ============================================================

fig = figure('Visible', 'off');
hold on;
legend_entries = {};

for m = 1:n_mat
    mat      = materials{m};
    c_smooth = linspace(0, max(mat.conc)*1.1, 200)';

    scatter(mat.conc, mat.alpha, 100, 'filled', ...
            'MarkerFaceColor', mat.color, 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
    plot(c_smooth, mat.lin_fit(c_smooth), '--', 'Color', mat.color, 'LineWidth', 1.8);

    legend_entries{end+1} = sprintf('%s data', mat.name);
    legend_entries{end+1} = sprintf('%s linear (\\epsilon=%.3f cm^{-1}mM^{-1})', ...
                                    mat.name, mat.epsilon);
end

xlabel('Concentration [mM]', 'FontSize', 12);
ylabel('\alpha [cm^{-1}]', 'FontSize', 12);
title('\alpha vs Concentration — Linear Fit (Beer-Lambert)', ...
      'FontSize', 13, 'FontWeight', 'bold');
legend(legend_entries, 'Location', 'best', 'FontSize', 9);
grid on; box on;
saveFigure(fig, output_dir, 'Summary_alpha_LinearFit');

%% ============================================================
%  5. FIGURE — alpha vs concentration, polynomial fit
% ============================================================

fig = figure('Visible', 'off');
hold on;
legend_entries = {};

for m = 1:n_mat
    mat      = materials{m};
    c_smooth = linspace(0, max(mat.conc)*1.1, 200)';

    scatter(mat.conc, mat.alpha, 100, 'filled', ...
            'MarkerFaceColor', mat.color, 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
    plot(c_smooth, mat.poly_fit(c_smooth), '-', 'Color', mat.color*0.6, 'LineWidth', 1.8);
    plot(c_smooth, mat.lin_fit(c_smooth),  '--', 'Color', mat.color,    'LineWidth', 1.2);

    legend_entries{end+1} = sprintf('%s data', mat.name);
    legend_entries{end+1} = sprintf('%s poly%d  R^2=%.4f', mat.name, mat.poly_deg, mat.poly_gof.rsquare);
    legend_entries{end+1} = sprintf('%s linear  R^2=%.4f', mat.name, mat.lin_gof.rsquare);
end

xlabel('Concentration [mM]', 'FontSize', 12);
ylabel('\alpha [cm^{-1}]', 'FontSize', 12);
title('\alpha vs Concentration — Polynomial vs Linear Fit', ...
      'FontSize', 13, 'FontWeight', 'bold');
legend(legend_entries, 'Location', 'best', 'FontSize', 9);
grid on; box on;
saveFigure(fig, output_dir, 'Summary_alpha_PolyFit');

fprintf('\nDone. Figures saved to %s/\n', output_dir);
fprintf('  Per-concentration: ln(I) vs x with alpha\n');
fprintf('  Summary Fig 1: alpha vs c, linear fit\n');
fprintf('  Summary Fig 2: alpha vs c, poly vs linear\n');

%% ============================================================
function saveFigure(fig, folder, filename)
    fig.Position = [100, 100, 900, 500];
    print(fig, fullfile(folder, filename), '-dpng', '-r150');
    close(fig);
end
