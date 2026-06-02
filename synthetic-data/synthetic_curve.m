%% PLANT-LIKE GROWTH (REALISTIC SYNTHETIC DATA GENERATOR)
clear; close all; clc;

%% Growth parameters
Lmax = 10;
k    = 1;
t0   = 5;
L = @(t) Lmax ./ (1 + exp(-k*(t - t0)));

t_samples = [1 3 5 7 9];

%% Discretization (ground truth curve)
N = 300;                      % dense underlying curve
M = 200;                      % observed points (dense like images)

s = linspace(0,1,N);

%% Randomness controls
rng(1)

noise_window = 15;
noise_sigmaA = 1.0;
noise_sigmaP = 1.0;

drift_sigma_alpha = 0.01;
drift_sigma_phi   = 0.01;

%% Geometry parameters
alpha0 = 0.10;
bend_strength  = 10.0;
twist_strength = 2.0;

alpha_clip = 1.30;

%% Envelope (growth region)
a_env = 1.0;
b_env = 0.3;
w = (s.^a_env) .* ((1 - s).^b_env);
w = w / max(w);

%% CSV output
fid = fopen('synthetic_data_independent_views.csv','w');
fprintf(fid,'gen,time,view,point,coord1,coord2\n');

%% Number of generations
nGen = 10;

for g = 1:nGen

    %% --- Smooth random fields
    etaA = smoothdata(randn(size(s)), 'gaussian', noise_window);
    etaA = etaA - mean(etaA);
    etaA = etaA / std(etaA);

    etaP = smoothdata(randn(size(s)), 'gaussian', noise_window);
    etaP = etaP - mean(etaP);
    etaP = etaP / std(etaP);

    f_alpha = (1 + 0.15*randn) * (w .* etaA);
    f_phi   = (1 + 0.15*randn) * (w .* etaP);

    phi0 = (2*rand - 1)*pi;

    for i = 1:length(t_samples)

        Lt  = L(t_samples(i));

        % --- Slightly perturb arc-length (break perfect parametrisation)
        ell = s * Lt + 0.02 * Lt * sin(2*pi*s);

        %% curvature drivers
        kappa_alpha = (bend_strength / Lt) * f_alpha;
        kappa_phi   = (twist_strength / Lt) * f_phi;

        %% integrate angles
        alpha = alpha0 + cumtrapz(ell, kappa_alpha);
        phi   = phi0   + cumtrapz(ell, kappa_phi);

        if drift_sigma_alpha > 0
            alpha = alpha + drift_sigma_alpha * cumsum(randn(size(s)));
        end
        if drift_sigma_phi > 0
            phi   = phi   + drift_sigma_phi   * cumsum(randn(size(s)));
        end

        alpha = min(max(alpha, -alpha_clip), alpha_clip);

        %% tangent
        Tx = sin(alpha).*cos(phi);
        Ty = sin(alpha).*sin(phi);
        Tz = cos(alpha);

        %% integrate 3D curve
        x = cumtrapz(ell, Tx);
        y = cumtrapz(ell, Ty);
        z = cumtrapz(ell, Tz);

        %% --- ADD STRUCTURED NOISE (image-like)
        noise_level = 0.01;

        z = z + noise_level * smoothdata(randn(size(z)), 'gaussian', 10);
        x = x + noise_level * randn(size(x));
        y = y + noise_level * randn(size(y));

        %% ===== NON-UNIFORM SAMPLING =====

        base = linspace(1, N, M);

        % small local irregularity
        perturb_xz = 5 * randn(size(base));
        perturb_yz = 5 * randn(size(base));

        idx_xz = round(base + perturb_xz);
        idx_yz = round(base + perturb_yz);

        idx_xz = max(min(idx_xz, N), 1);
        idx_yz = max(min(idx_yz, N), 1);

        % keep ordering (important)
        idx_xz = sort(idx_xz);
        idx_yz = sort(idx_yz);

        %% --- projection-specific distortion (small)
        x_proj = x + 0.01 * x.^2;
        y_proj = y + 0.01 * y.^2;

        %% ===== EXPORT XZ =====
        for j = 1:M
            jj = idx_xz(j);
            fprintf(fid, '%d,%.3f,xz,%d,%.6f,%.6f\n', ...
                g, t_samples(i), j, ...
                x_proj(jj), z(jj));
        end

        %% ===== EXPORT YZ =====
        for j = 1:M
            jj = idx_yz(j);
            fprintf(fid, '%d,%.3f,yz,%d,%.6f,%.6f\n', ...
                g, t_samples(i), j, ...
                y_proj(jj), z(jj));
        end

    end
end

fclose(fid);

disp('synthetic dataset generated');