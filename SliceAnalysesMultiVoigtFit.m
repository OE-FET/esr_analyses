function [out_struct, coefs_con] = SliceAnalysesMultiVoigtFit(varargin)
%SLICEANALYSESMULTIVOIGTFIT performs normalization and spin-counting of an ESR
%signal by fitting it to a voigt function.
%
%   Performs a 1D fit of a cw-EPR spectrum to determine spin cohrence time
%   T2, the magnetic susceptibility and the number of spins in the sample.
%   This program assumes a Voigtian resonance line which is the convolution
%   of multiple Lorentzian spin-ensembles with identifical lifetimes T1 and
%   T2 with a Gaussian distribution of resonance fields.
%
%   Field modulation and the 1st harmonic detection of the cw-EPR signal,
%   together with possible distortions from overmodulation, are explicitly
%   accounted for. The MW field distribution in the cavity is taken from
%   the Bruker DSC file and is used when calculating the MW field amplitude
%   over the sample volume.
%
%   INPUT(S):
%   SLICEANALYSESMULTIVOIGTFIT() - prompts user for spectrum file
%   ...FIT('Path')               - path to file with ESR data
%   ...FIT('PathSIG', 'PathBG')  - path to signal, path to background
%   ...FIT(dset)                 - data set from BrukerRead
%
%   OUTPUT(S):
%	out_struct  - structure containing the measurement data and fit results
%   out_table   - fit results in table format
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   Updated by Remington Carey to fit multiple Voigtians in a curve

%% Initialize =============================================================

%
import esr_analyses.*
import esr_analyses.utils.*

close all

% Load data
dset = load_spectrum_dialog(varargin);
[x,y,pars] = dset_to_tuple(dset);

% Number of Voigtians to fit
N = 2;

% Calculate MW field
Bmw = get_mw_fields(pars);

%% Get starting points for fit ============================================

% Fit to pseudo-Voigtian
fit = pseudo_voigt_fit(x, y, 'deriv', 1);

% Set initial parameters from pseudo-Voigtian fit
FWHM_lorentz  = fit.FWHM_lorentz;       % in Gauss
FWHM_gauss    = fit.FWHM_gauss;         % in Gauss
A0   = fit.a/(pars.B0MA*1e4 * 1e4/8 * Bmw);
B0   = fit.x0;                          % in Gauss
T1   = 1e-07;                           % in sec, assume not saturated
T2   = 1/(gmratio * FWHM_lorentz*1E-4); % in sec

% Assign parameters to array
var0 = ones(N,1)*[A0/2 B0 T2 FWHM_gauss];


%% Perform voigt fit ======================================================

% Create single fit function
func_single = @(var, x) abs(var(1))*esr_voigt_simulation(x, abs(var(2)), T1, ...
    abs(var(3)), abs(var(4)), Bmw, pars.B0MA*1e4, 1);

% Expand to multiple peaks
multi_fit_func = @(v, x) to_multi(func_single, v, x);

% Create fit options
opt = optimset('TolFun', 1e-9, 'TolX', 1e-9, 'MaxFunEvals', 1e10,...
    'MaxIter', 1e10);

% Fit model to data with Nelder Mead algorithm (unconstrained fit)
fitres   = nelder_mead_fit(multi_fit_func, x, y, var0, opt);

% Estimate confidence intervals
conf_int = standarderror(fitres);

% Convert fit results to EPR parameters
A     = abs(fitres.coef(:,1));
B0    = abs(fitres.coef(:,2));
T2    = abs(fitres.coef(:,3));
Brms  = abs(fitres.coef(:,4));

dA    = full(abs(conf_int(1:N)));
dB0   = full(abs(conf_int(N+1:2*N)));
dT2   = full(abs(conf_int(3*N+1:3*N)));
dBrms = full(abs(conf_int(4*N+1:4*N)));

% Plot
p = plot(fitres);
xlabel('Magnetic field [G]')
ylabel('ESR signal [a.u.]')
legend_texts = ['Data', 'Fit' cell(1,N)];
for i=1:N
    c = fitres.coef(i,:);
    plot(x, func_single(c, x));
    legend_texts{i+2} = ['Peak ' num2str(i)];
end
legend(legend_texts);

%% Susceptibility calculation =============================================

% Get g-factor and scaling factor for pseudo-modulation
pars.GFactor    = b2g(mean(B0)*1e-4, pars.MWFQ);
modScaling      = pars.B0MA*1e4 * 1e4/8;

% Get susceptibility & num. spins
Chi = zeros(size(A));
dChi = zeros(size(A));
NSpin = zeros(size(A));
dNSpin = zeros(size(A));

for i=1:length(A) % calculate for each peak
    areaDI = modScaling * Bmw .* A(i);
    areaDIerror = modScaling * Bmw .* dA(i);

    % get 'maximum' value, even though all values are equal...
    [Chi(i), dChi(i)]   = max(susceptibility_calc(areaDI, pars, 'dA', areaDIerror));
    [NSpin(i), dNSpin(i)] = max(spincounting(areaDI, pars, 'dA', areaDIerror));
end

%% Create output structure ================================================
out_struct = struct(...
    'x', x, 'y', y, 'pars', pars, 'fitres', fitres, ...
    'A', A, 'dA', dA, 'B0', B0, 'dB0', dB0, 'T2', T2, 'dT2', dT2, ...
    'Brms', Brms, 'dBrms', dBrms, 'Chi', Chi(1), 'dChi', dChi(1), ...
    'NSpin', NSpin(1), 'dNSpin', dNSpin(1));


%% Perform constrained fit=================================================

% Fit model w/ parameter constraints
coefs_con = constrained_fit(x, y, fitres);

% Plot constrained fit
figure('Name','Constrained fit')
hold on
plot(x,y,'ko','MarkerSize',1)
plot(x,fitres.fitfunc(coefs_con,x),'-r');
xlabel('Magnetic field [G]')
ylabel('ESR signal [a.u.]')
set(gca,'xlim',[min(x) max(x)],'ylim',[min(y) max(y)]);
legend_texts = ['Data', 'Fit' cell(1,N)];
for i=1:N
    c = coefs_con(i,:);
    plot(x, func_single(c, x));
    legend_texts{i+2} = ['Peak ' num2str(i)];
end
legend(legend_texts);

end

%% Multi fit function =====================================================
function y = to_multi(func_single, variables, x)

    [N, ~] = size(variables);

    result1 = func_single(variables(1,:), x);
    shape = num2cell(size(result1));

    results = zeros(N, shape{:});
    results(1,:,:) = result1;

    for i=2:N
        results(i,:,:) = func_single(variables(i,:), x);
    end

    y = sum(results, 1).';

end

%% Constrained fit function ===============================================
function coefs_con = constrained_fit(x, y, fitres)

    % Constraints hard-programmed for N = 2 (Voigtians to fit)

    coefs = fitres.coef;

    % Inequality constraints (none)
    con_problem.Aineq     = [];
    con_problem.bineq     = [];

    % Equality constraints
    con_problem.Aeq       = [];
    con_problem.Beq       = [];
%     con_problem.Aeq       = [1 -1 0 0 0 0 0 0]; % equal areas
%     con_problem.Beq       = 0;                  % equal areas

    % Lower bounds
    con_problem.lb        = zeros(size(coefs));
    con_problem.lb(:,2)   = min(x); % No resonance centers outside x-axis
    con_problem.lb(:,3)   = 10^-9;  % min T2

    % Upper bounds
    con_problem.ub        = inf(size(coefs));
    con_problem.ub(:,2)   = max(x); % No resonance centers outside x-axis
    con_problem.ub(:,3)   = 10^-6;  % max T2
    con_problem.ub(:,4)   = 10;     % Restrict Brms to be 10 G at max

    % Solver, objective function, initial parameter estimates, and fit
    % options
    con_problem.solver    = 'fmincon';
    con_problem.objective = @(coefs) sum((y-fitres.fitfunc(coefs,x)).^2);
    con_problem.x0        = coefs;
    con_problem.options   = optimoptions('fmincon', 'ConstraintTolerance',...
                            1e-14,'MaxFunctionEvaluations',1e14, ...
                            'MaxIterations',1e14,'OptimalityTolerance',...
                            1e-18, 'StepTolerance',1e-18,'Display','final-detailed');

    % Perform the fit
    coefs_con = fmincon(con_problem);

end