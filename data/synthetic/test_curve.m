function test_curve()

clear; clc;

imgW = 512;
imgH = 512;

%% ==========================================
%% OUTPUT DIRECTORY
%% ==========================================

outputDir = 'test_curves';

if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

%% ==========================================
%% STORAGE TABLES
%% ==========================================

viewCurveData = table();
trueCurveData = table();

%% ==========================================
%% CURVE 1 : TRANSITIONAL SPIRAL
%% ==========================================

curve_id = 1;

t = linspace(0,1,200)';

x = t .* cos(2*pi*t.^2);
y = t .* sin(2*pi*t.^2);
z = 4*t;

camera = auto_camera(x, y, z, imgW, imgH);
save_curve(camera);


%% ==========================================
%% CURVE 2 : S-Bend
%% ==========================================

curve_id = 2;

t = linspace(0,0.75,300)';

x = 2*sin(pi*t);
y = 4*t;
z = sin(2*pi*t);

camera = auto_camera(x, y, z, imgW, imgH);
save_curve(camera);


%% ==========================================
%% CURVE 3 : CUBIC BEZIER
%% ==========================================

curve_id = 3;

t = linspace(0,1,1000)';

P0 = [ 0 0 0];
P1 = [ 2 5 -3];
P2 = [-2 5 3];
P3 = [ 0 10 0];

B0 = (1-t).^3;
B1 = 3*(1-t).^2.*t;
B2 = 3*(1-t).*t.^2;
B3 = t.^3;

x = B0*P0(1) + B1*P1(1) + B2*P2(1) + B3*P3(1);
y = B0*P0(2) + B1*P1(2) + B2*P2(2) + B3*P3(2);
z = B0*P0(3) + B1*P1(3) + B2*P2(3) + B3*P3(3);

camera = auto_camera(x, y, z, imgW, imgH);
save_curve(camera);


%% ==========================================
%% SAVE AGGREGATED CSV FILES
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
%% NESTED FUNCTION : SAVE CURVE
%% ============================================================

    function save_curve(camera)

        n = length(x);

        %% ------------------------------------------
        %% TRUE CURVE CSV
        %% ------------------------------------------

        Ttrue = table( ...
            repmat(curve_id,n,1), ...
            ones(n,1), ...
            (1:n)', ...
            x, y, z, ...
            'VariableNames', ...
            {'gen','time','point','x','y','z'});

        trueCurveData = [trueCurveData; Ttrue];

        %% ------------------------------------------
        %% VIEW CSV
        %% ------------------------------------------

        Txz = table( ...
            repmat(curve_id,n,1), ...
            ones(n,1), ...
            repmat("xz",n,1), ...
            (1:n)', ...
            x, ...
            z, ...
            'VariableNames', ...
            {'gen','time','view','point','coord1','coord2'});

        Tyz = table( ...
            repmat(curve_id,n,1), ...
            ones(n,1), ...
            repmat("yz",n,1), ...
            (1:n)', ...
            y, ...
            z, ...
            'VariableNames', ...
            {'gen','time','view','point','coord1','coord2'});

        viewCurveData = [viewCurveData; Txz; Tyz];

        %% ------------------------------------------
        %% THICKNESS
        %% ------------------------------------------

        s = linspace(0,1,n)';

        radius_cm = ...
            0.02 * ...
            (1 + 0.2 * smoothdata(randn(size(s)),'gaussian',10));

        radius_cm = max(radius_cm,0.02);

        %% ------------------------------------------
        %% RENDER
        %% ------------------------------------------

        mask_xz = render_mask_fixed( ...
            x, z, radius_cm, ...
            camera.scale, camera.x_offset, camera.z_offset, ...
            imgH, imgW);

        mask_yz = render_mask_fixed( ...
            y, z, radius_cm, ...
            camera.scale, camera.x_offset, camera.z_offset, ...
            imgH, imgW);

        %% ------------------------------------------
        %% SAVE PNGS
        %% ------------------------------------------

        genDir = fullfile( ...
            outputDir, ...
            sprintf('gen_%02d',curve_id));

        if ~exist(genDir,'dir')
            mkdir(genDir);
        end

        fname_xz = sprintf( ...
            'plant_gen%02d_t01_xz.png', ...
            curve_id);

        fname_yz = sprintf( ...
            'plant_gen%02d_t01_yz.png', ...
            curve_id);

        imwrite(mask_xz, ...
            fullfile(genDir,fname_xz));

        imwrite(mask_yz, ...
            fullfile(genDir,fname_yz));

        fprintf('Saved %s\n',fname_xz);
        fprintf('Saved %s\n',fname_yz);

        cameraTable = table( ...
            curve_id,...
            camera.scale,...
            camera.x_offset,...
            camera.z_offset,...
            imgH,...
            imgW,...
            'VariableNames',...
            {'gen','scale','x_offset','z_offset','imgH','imgW'});

        writetable( ...
            cameraTable, ...
            fullfile(genDir,'camera_params.csv'));

        fprintf('Saved camera_params.csv\n');

    end

end

%% ============================================================
%% FUNCTION : CAMERA RENDERING
%% ============================================================

function mask = render_mask_fixed( ...
    xp_raw, ...
    zp_raw, ...
    radius_cm, ...
    scale, ...
    x_offset, ...
    z_offset, ...
    imgH, ...
    imgW)

    xp = xp_raw * scale + x_offset;

    zp = imgH - ...
         (zp_raw * scale + z_offset);

    fprintf("xp range = [%f, %f]\n", min(xp), max(xp));
    fprintf("zp range = [%f, %f]\n", min(zp), max(zp));

    radius_px = radius_cm * scale;

    [X,Z] = meshgrid(1:imgW,1:imgH);

    mask = false(imgH,imgW);

    N = length(xp);

    for k = 1:(N-1)

        segLength = hypot( ...
            xp(k+1)-xp(k), ...
            zp(k+1)-zp(k));

        nInterp = max(2,ceil(segLength));

        s = linspace(0,1,nInterp);

        xs = xp(k) + s*(xp(k+1)-xp(k));
        zs = zp(k) + s*(zp(k+1)-zp(k));

        rs = radius_px(k) + ...
            s*(radius_px(k+1)-radius_px(k));

        for j = 1:nInterp

            dx = X - xs(j);
            dz = Z - zs(j);

            mask = mask | ...
                (dx.^2 + dz.^2 <= rs(j)^2);

        end
    end

    mask = imgaussfilt(double(mask),0.6);
    mask = mask > 0.10;

    mask = imclose(mask, strel('disk',2));

    mask = bwareaopen(mask,5);

end

function camera = auto_camera(x,y,z,imgW,imgH)

    margin = 0.85;

    horizontal_extent = max( ...
        max(x)-min(x), ...
        max(y)-min(y));

    vertical_extent = max(z)-min(z);

    camera.scale = margin * min( ...
        imgW / horizontal_extent, ...
        imgH / vertical_extent);

    % centre both x and y around image centre
    xy_centre = 0.5 * ...
        (max([x;y]) + min([x;y]));

    z_centre = 0.5 * ...
        (max(z) + min(z));

    camera.x_offset = ...
        imgW/2 - camera.scale*xy_centre;

    camera.z_offset = ...
        imgH/2 - camera.scale*z_centre;

    camera.imgW = imgW;
    camera.imgH = imgH;

end