%% SYNTHETIC PLANT IMAGE GENERATOR (CENTERLINE → MASK)

clear; close all; clc;

%% ==========================================
%% Ground-truth data storage
%% ==========================================

trueCurveData = table();
viewCurveData = table();

%% ===============================
%% PHYSICAL SCALE
%% ===============================

Lmax = 7.0;      % cm
k    = 1;
t0   = 2;

L = @(t) Lmax ./ (1 + exp(-k*(t - t0)));

t_samples = [1 3 5 7 9];

%% Image resolution
imgH = 512;
imgW = 512;

%% ===============================
%% DISCRETISATION
%% ===============================

N = 400;
s = linspace(0,1,N);

rng(1)

%% ===============================
%% NOISE / GEOMETRY
%% ===============================

noise_window = 15;

alpha0 = 0.10;
bend_strength  = 10.0;
twist_strength = 2.0;

alpha_clip = 1.3;

%% Envelope
a_env = 1.0;
b_env = 0.3;
w = (s.^a_env) .* ((1 - s).^b_env);
w = w / max(w);

%% ===============================
%% OUTPUT DIRECTORY
%% ===============================

outputDir = 'synthetic_curves';

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% ===============================
%% LOOP
%% ===============================

nGen = 5;

for g = 1:nGen

    %% Smooth random fields
    etaA = smoothdata(randn(size(s)), 'gaussian', noise_window);
    etaA = etaA / std(etaA);

    etaP = smoothdata(randn(size(s)), 'gaussian', noise_window);
    etaP = etaP / std(etaP);

    f_alpha = (1 + 0.15*randn) * (w .* etaA);
    f_phi   = (1 + 0.15*randn) * (w .* etaP);

    phi0 = (2*rand - 1)*pi;

    %% ======================================
    %% FIXED CAMERA FROM FINAL TIME STEP
    %% ======================================

    Lt_final = L(t_samples(end));

    ell = s * Lt_final + 0.02 * Lt_final * sin(2*pi*s);

    kappa_alpha = (bend_strength / Lt_final) * f_alpha;
    kappa_phi   = (twist_strength / Lt_final) * f_phi;

    alpha = alpha0 + cumtrapz(ell, kappa_alpha);
    phi   = phi0   + cumtrapz(ell, kappa_phi);

    alpha = min(max(alpha, -alpha_clip), alpha_clip);

    Tx = sin(alpha).*cos(phi);
    Ty = sin(alpha).*sin(phi);
    Tz = cos(alpha);

    x_final = cumtrapz(ell, Tx);
    y_final = cumtrapz(ell, Ty);
    z_final = cumtrapz(ell, Tz);

    % enforce upward growth
    z_final = z_final - min(z_final);
    if z_final(end) < z_final(1)
        z_final = max(z_final) - z_final;
    end

    % thickness (worst-case envelope)
    radius_final = 0.02 * (1 + 0.2 * smoothdata(randn(size(s)),'gaussian',10));
    radius_final = max(radius_final, 0.02);

    %% --- include BOTH projections in bounding box
    horizontal_extent = max( ...
        max(x_final)-min(x_final), ...
        max(y_final)-min(y_final));

    vertical_extent = max(z_final)-min(z_final);

    margin = 0.85;

    scale = margin * min( ...
        imgW/horizontal_extent, ...
        imgH/vertical_extent);

    x_center = 0.5*(max(x_final)+min(x_final));
    y_center = 0.5*(max(y_final)+min(y_final));
    z_center = 0.5*(max(z_final)+min(z_final));

    x_offset = imgW/2 - scale*x_center;
    y_offset = imgH/2 - scale*y_center;
    z_offset = imgH/2 - scale*z_center;

    %% ===============================
    %% TIME LOOP
    %% ===============================

    for i = 1:length(t_samples)

        Lt = L(t_samples(i));

        ell = s * Lt + 0.02 * Lt * sin(2*pi*s);

        kappa_alpha = (bend_strength / Lt) * f_alpha;
        kappa_phi   = (twist_strength / Lt) * f_phi;

        alpha = alpha0 + cumtrapz(ell, kappa_alpha);
        phi   = phi0   + cumtrapz(ell, kappa_phi);

        alpha = min(max(alpha, -alpha_clip), alpha_clip);

        Tx = sin(alpha).*cos(phi);
        Ty = sin(alpha).*sin(phi);
        Tz = cos(alpha);

        x = cumtrapz(ell, Tx);
        y = cumtrapz(ell, Ty);
        z = cumtrapz(ell, Tz);

        % enforce bottom → top growth
        z = z - min(z);
        if z(end) < z(1)
            z = max(z) - z;
        end

        %% ==========================================
        %% TRUE 3D CURVE
        %% ==========================================

        nPts = length(x);

        trueBlock = table( ...
            repmat(g,nPts,1), ...
            repmat(t_samples(i),nPts,1), ...
            (1:nPts)', ...
            x(:), ...
            y(:), ...
            z(:), ...
            'VariableNames', ...
            {'gen','time','point','x','y','z'});

        trueCurveData = [trueCurveData; trueBlock];

        %% ==========================================
        %% XZ VIEW
        %% ==========================================

        xzBlock = table( ...
            repmat(g,nPts,1), ...
            repmat(t_samples(i),nPts,1), ...
            repmat("xz",nPts,1), ...
            (1:nPts)', ...
            x(:), ...
            z(:), ...
            'VariableNames', ...
            {'gen','time','view','point','coord1','coord2'});

        %% ==========================================
        %% YZ VIEW
        %% ==========================================

        yzBlock = table( ...
            repmat(g,nPts,1), ...
            repmat(t_samples(i),nPts,1), ...
            repmat("yz",nPts,1), ...
            (1:nPts)', ...
            y(:), ...
            z(:), ...
            'VariableNames', ...
            {'gen','time','view','point','coord1','coord2'});

        %% ==========================================
        %% XY VIEW
        %% ==========================================
        xyBlock = table( ...
            repmat(g,nPts,1), ...
            repmat(t_samples(i),nPts,1), ...
            repmat("xy",nPts,1), ...
            (1:nPts)', ...
            x(:), ...
            y(:), ...
            'VariableNames', ...
            {'gen','time','view','point','coord1','coord2'});

        viewCurveData = [viewCurveData; xzBlock; yzBlock; xyBlock];

        %% thickness
        radius_cm = 0.02 * (1 + 0.2 * smoothdata(randn(size(s)),'gaussian',10));
        radius_cm = max(radius_cm, 0.02);

        %% render with FIXED camera
        mask_xz = render_mask_fixed(x, z, radius_cm, scale, x_offset, z_offset, imgH, imgW, true);
        mask_yz = render_mask_fixed(y, z, radius_cm, scale, y_offset, z_offset, imgH, imgW, true);
        mask_xy = render_mask_fixed(x, y, radius_cm, scale, x_offset, y_offset, imgH, imgW, false);
        
        %% save
        genDir = fullfile(outputDir, sprintf('gen_%02d', g));
        if ~exist(genDir, 'dir')
            mkdir(genDir);
        end

        cameraTable = table( ...
            g,...
            scale,...
            x_offset,...
            y_offset,...
            z_offset,...
            imgH,...
            imgW,...
            'VariableNames', ...
            {'gen','scale','x_offset','y_offset','z_offset','imgH','imgW'});

        writetable( ...
            cameraTable, ...
            fullfile(genDir,'camera_params.csv'));

        fname_xz = sprintf('plant_gen%02d_t%02d_xz.png', g, i);
        fname_yz = sprintf('plant_gen%02d_t%02d_yz.png', g, i);
        fname_xy = sprintf('plant_gen%02d_t%02d_xy.png', g, i);

        imwrite(mask_xz, fullfile(genDir, fname_xz));
        imwrite(mask_yz, fullfile(genDir, fname_yz));
        imwrite(mask_xy, fullfile(genDir, fname_xy));

        fprintf('Saved %s, %s, and %s\n', fname_xz, fname_yz, fname_xy);

    end
end

%% ==========================================
%% SAVE CSV FILES
%% ==========================================

writetable( ...
    viewCurveData, ...
    fullfile(outputDir,'view_curves.csv'));

writetable( ...
    trueCurveData, ...
    fullfile(outputDir,'true_curve.csv'));

fprintf('Saved view_curves.csv\n');
fprintf('Saved true_curve.csv\n');

disp('Done.');

%% ============================================================
%% FUNCTION: CAMERA RENDERING
%% ============================================================

function mask = render_mask_fixed(x1_raw, x2_raw, radius_cm, scale, x1_offset,x2_offset, imgH, imgW, flip_second_axis)

    x1 = x1_raw*scale + x1_offset;

    if flip_second_axis
        x2 = imgH - (x2_raw*scale + x2_offset);
    else
        x2 = x2_raw*scale + x2_offset;
    end

    radius_px = radius_cm*scale;

    [X1,X2] = meshgrid(1:imgW,1:imgH);

    mask = false(imgH,imgW);

    N = length(x1);

    for k = 1:(N-1)

        segLength = hypot( ...
            x1(k+1)-x1(k), ...
            x2(k+1)-x2(k));

        nInterp = max(2,ceil(segLength));

        sInterp = linspace(0,1,nInterp);

        xs = x1(k) + sInterp*(x1(k+1)-x1(k));
        zs = x2(k) + sInterp*(x2(k+1)-x2(k));

        rs = radius_px(k) ...
            + sInterp*(radius_px(k+1)-radius_px(k));

        for j = 1:nInterp

            dx1 = X1 - xs(j);
            dx2 = X2 - zs(j);

            mask = mask | ...
                (dx1.^2 + dx2.^2 <= rs(j)^2);

        end
    end

    mask = imgaussfilt(double(mask),0.6);
    mask = mask > 0.10;
    mask = imclose(mask,strel('disk',2));
    mask = bwareaopen(mask,5);
end