%% SYNTHETIC PLANT IMAGE GENERATOR (CENTERLINE → MASK)

clear; close all; clc;

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

outputDir = 'images';

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% ===============================
%% LOOP
%% ===============================

nGen = 5;

for g = 1:nGen + 2

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
    if g  <= nGen
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
    elseif g == nGen + 1
        % helix (for testing)
        ell = s * Lt;

        R = 0.4;             % helix radius
        n_turns = 2.5;       % controls complexity

        % helix
        x = R * cos(2*pi*n_turns*s);
        y = R * sin(2*pi*n_turns*s);
        z = Lt * s;

        % optional: slight vertical modulation (more realistic)
        z = z + 0.1 * R * sin(2*pi*n_turns*s);
    else
        % loop de loop (for testing)
        ell = s * Lt + 0.02 * Lt * sin(2*pi*s);   % keep your original variation

        center = 0.4;
        width  = 0.2;

        bump = exp(-((s - center).^2)/width^2);

        % --- LOOP CONTROL ---
        % full rotation localized in s
        alpha_loop = 2*pi * bump;

        % add upward bias so curve keeps going up
        alpha_bias = 0.3;     % small tilt upward

        alpha = alpha_loop + alpha_bias;

        % --- Y-direction motion
        phi = 0.2 * bump;

        % --- compute tangent ---
        Tx = sin(alpha).*cos(phi);
        Ty = sin(alpha).*sin(phi) +  0.3;
        Tz = cos(alpha);

        % integrate
        x = cumtrapz(ell, Tx);
        y = cumtrapz(ell, Ty);
        z = cumtrapz(ell, Tz);

        % normalize z
        z = z - min(z);
        if z(end) < z(1)
            z = max(z) - z;
        end
    end
    

    % thickness (worst-case envelope)
    radius_final = 0.02 * (1 + 0.2 * smoothdata(randn(size(s)),'gaussian',10));
    radius_final = max(radius_final, 0.02);

    %% --- include BOTH projections in bounding box
    xmin = min([x_final - radius_final, y_final - radius_final], [], 'all');
    xmax = max([x_final + radius_final, y_final + radius_final], [], 'all');

    zmin = min(z_final - radius_final);
    zmax = max(z_final + radius_final);

    width_cm  = xmax - xmin;
    height_cm = zmax - zmin;

    % LARGE margin to prevent clipping
    margin = 0.5;

    width_cm  = width_cm  * (1 + margin);
    height_cm = height_cm * (1 + margin);

    sx = imgW / width_cm;
    sz = imgH / height_cm;

    scale = min(sx, sz);

    % enforce x=0 center
    x_offset = imgW / 2;

    % bottom alignment
    margin_px = 20;   % small buffer in pixels

    z_offset = margin_px - zmin * scale;

    % ensure top fits
    top_pixel = zmax * scale + z_offset;

    if top_pixel > imgH - margin_px
        scale = (imgH - 2*margin_px) / (zmax - zmin);
        z_offset = margin_px - zmin * scale;
    end

    %% ===============================
    %% SAVE CAMERA PARAMETERS
    %% ===============================

    genDir = fullfile(outputDir, sprintf('gen_%02d', g));
    if ~exist(genDir, 'dir')
        mkdir(genDir);
    end

    camera_file = fullfile(genDir, 'camera_params.csv');

    % column names
    headers = {'scale','x_offset','z_offset','imgH','imgW'};
    % values (1 row)
    values = [scale, x_offset, z_offset, imgH, imgW];
    % write CSV
    fid = fopen(camera_file, 'w');
    % header row
    fprintf(fid, '%s,%s,%s,%s,%s\n', headers{:});
    % data row
    fprintf(fid, '%.10f,%.10f,%.10f,%d,%d\n', values);
    fclose(fid);

    %% ===============================
    %% TIME LOOP
    %% ===============================

    for i = 1:length(t_samples)
        if g <= nGen
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
        elseif g == nGen + 1
            % helix (for testing)
            Lt = L(t_samples(i));
            ell = s * Lt;

            R = 0.4;             % helix radius
            n_turns = 2.5;       % controls complexity

            % helix
            x = R * cos(2*pi*n_turns*s);
            y = R * sin(2*pi*n_turns*s);
            z = Lt * s;

            % optional: slight vertical modulation (more realistic)
            z = z + 0.1 * R * sin(2*pi*n_turns*s);
        else
            % loop de loop (for testing) 
            Lt = L(t_samples(i));
            ell = s * Lt + 0.02 * Lt * sin(2*pi*s);   % keep your original variation

            center = 0.4;
            width  = 0.2;

            bump = exp(-((s - center).^2)/width^2);

            % --- LOOP CONTROL ---
            % full rotation localized in s
            alpha_loop = 2*pi * bump;

            % add upward bias so curve keeps going up
            alpha_bias = 0.3;     % small tilt upward

            alpha = alpha_loop + alpha_bias;

            % --- Y-direction motion (small!!)
            phi = 0.2 * bump;

            % --- compute tangent ---
            Tx = sin(alpha).*cos(phi);
            Ty = sin(alpha).*sin(phi) + 0.3; 
            Tz = cos(alpha);

            % integrate
            x = cumtrapz(ell, Tx);
            y = cumtrapz(ell, Ty);
            z = cumtrapz(ell, Tz);

            % normalize z
            z = z - min(z);
            if z(end) < z(1)
                z = max(z) - z;
            end
        end 
        

        %% thickness
        radius_cm = 0.02 * (1 + 0.2 * smoothdata(randn(size(s)),'gaussian',10));
        radius_cm = max(radius_cm, 0.02);

        %% save curve data
        curve_table = table( ...
            s(:), ...
            x(:), ...
            y(:), ...
            z(:), ...
            radius_cm(:), ...
            'VariableNames', ...
            {'s','x','y','z','radius'} ...
        );

        curve_file = fullfile( ...
            genDir, ...
            sprintf('curve_gen%02d_t%02d.csv',g,i) ...
        );

        ds = sqrt(diff(x).^2 + diff(y).^2 + diff(z).^2);
        arc = [0 cumsum(ds)];
        arc = arc / arc(end);

        arc_uniform = linspace(0,1,length(s));

        x_uniform = interp1(arc, x, arc_uniform, 'pchip');
        y_uniform = interp1(arc, y, arc_uniform, 'pchip');
        z_uniform = interp1(arc, z, arc_uniform, 'pchip');

        T = table( ...
            arc_uniform(:), ...
            x_uniform(:), ...
            y_uniform(:), ...
            z_uniform(:), ...
            'VariableNames', ...
            {'s','x','y','z'} );

        writetable( ...
            T, ...
            fullfile(genDir, ...
            sprintf('U_gen%02d_t%02d.csv',g,i)) ...
        );


        %% render with FIXED camera
        mask_xz = render_mask_fixed(x, z, radius_cm, scale, x_offset, z_offset, imgH, imgW);
        mask_yz = render_mask_fixed(y, z, radius_cm, scale, x_offset, z_offset, imgH, imgW);

        %% save
        genDir = fullfile(outputDir, sprintf('gen_%02d', g));
        if ~exist(genDir, 'dir')
            mkdir(genDir);
        end

        fname_xz = sprintf('plant_gen%02d_t%02d_xz.png', g, i);
        fname_yz = sprintf('plant_gen%02d_t%02d_yz.png', g, i);

        imwrite(mask_xz, fullfile(genDir, fname_xz));
        imwrite(mask_yz, fullfile(genDir, fname_yz));

        fprintf('Saved %s and %s\n', fname_xz, fname_yz);

    end
end

disp('Done.');

%% ============================================================
%% FUNCTION: CAMERA RENDERING
%% ============================================================

function mask = render_mask_fixed(xp_raw, zp_raw, radius_cm, scale, x_offset, z_offset, imgH, imgW)

    xp = xp_raw * scale + x_offset;
    zp = imgH - (zp_raw * scale + z_offset);

    radius_px = radius_cm * scale;

    [X, Z] = meshgrid(1:imgW, 1:imgH);
    mask = false(imgH, imgW);

    N = length(xp);

    for k = 1:N
        dx = X - xp(k);
        dz = Z - zp(k);

        mask = mask | (dx.^2 + dz.^2 <= radius_px(k)^2);
    end

    % smoothing
    % mask = imgaussfilt(double(mask), 1.2);
    % mask = mask > 0.35;

    % speckle noise
    % noise = rand(imgH, imgW) < 0.002;
    % mask = mask | noise;

    % clean
    % mask = bwareaopen(mask, 30);
end