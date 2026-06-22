

function run_torsion()
    clc; close all;

    E_PAPER      = 35000;
    nu_PAPER     = 0.20;
    ft_PAPER     = 3.0;
    GF_PAPER     = 0.08;
    kappa0_PAPER = 6.0e-5;
    k_PAPER      = 10.0;

    base_prefix = 'Job-1';
    out_dir     = 'out_torsion_LIVE_ONLY_OLIVER';
    if ~exist(out_dir,'dir'), mkdir(out_dir); end

    load_mode      = 'torsion_brokenshire';
    theta_snapshot = 2.8e-3;
    a_arm          = 100.0;

    AB_points = [ 200,  0, +25;
                  200,  0, -25];

    opts = struct( ...
        'n_increments',        140, ...
        'max_iter',            18,  ...
        'tol_damage',          5.0e-4, ...
        'tol_residual',        1.0e-4, ...
        'min_stiffness_ratio', 1.0e-4, ...
        'relax_damage',        0.35, ...
        'load_path',           load_mode, ...
        'E',                   E_PAPER, ...
        'nu',                  nu_PAPER, ...
        'GF',                  GF_PAPER, ...
        'ft',                  ft_PAPER, ...
        'k',                   k_PAPER, ...
        'kappa0',              kappa0_PAPER, ...
        'Uy_end',              0.30, ...
        'theta_end',           3.0e-3, ...
        'eeq_model',           'mod_vm', ...
        'use_crack_band',      true, ...
        'crack_band_method',   'oliver_directional', ...
        'principal_dir_iters', 10, ...
        'defscale_3d',         10.0, ...
        'crack_threshold',     0.70, ...
        'do_visualization',    true,  ...
        'viz_every',           1,  ...
        'save_video',          true,  ...
        'video_rotate_camera', true, ...
        'video_elevation',     18, ...
        'video_revolutions',   1.25, ...
        'video_framerate',     20, ...
        'camera_view_angles',  [  0 20;  90 20; 180 20; 270 20; ...
                                  45 30; 135 30; 225 30; 315 30; ...
                                   0 90], ...
        'camera_view_names',   {{'front','right','back','left','iso1','iso2','iso3','iso4','top'}}, ...
        'save_multiview_snapshots', false, ...
        'n_snapshots',         5, ...
        'save_peak_picture',    true, ...
        'damage_save_levels',   [], ...
        'save_peak_multiview',  false, ...
        'save_live_snapshots', true, ...
        'n_live_snaps',        8, ...
        'live_snap_view',      [45 30], ...
        'live_snap_capture_fig', false, ...
        'cmod_AB',             AB_points, ...
        'lever_arm',           a_arm, ...
        'Fmax',                2000, ...
        'slice_trigger',       inf, ...
        'check_residual_every',0, ...
        'save_full_results',   false, ...
        'save_images',         false, ...
        'save_summary_plots',  true, ...
        'save_tables',         true);

    run_prefix = fullfile(out_dir, sprintf('%s_StaticFast_%s_OLIVER', base_prefix, opts.eeq_model));
    fprintf('\n--- Running FAST QUASI-STATIC Torsion Damage Solver: %s + full Oliver bandwidth ---\n', opts.eeq_model);
    damage_static_vectorized(run_prefix, base_prefix, opts);
    fprintf('Done. CSV + purple snapshots + video saved in: %s\n', out_dir);
end

function damage_static_vectorized(prefix_out, prefix_mesh, opts)
    get = @(f,def) local_get(opts,f,def);

    do_viz      = get('do_visualization', true);
    save_video  = get('save_video', true);
    E           = get('E',35000.0);
    nu          = get('nu',0.20);
    GF          = get('GF',0.060);
    ft          = get('ft',E*get('kappa0',6.0e-5));
    k_tc        = get('k',10.0);
    kappa0      = get('kappa0',6.0e-5);
    Uy_end      = get('Uy_end',0.10);
    path        = get('load_path','torsion_brokenshire');
    eeq_model   = get('eeq_model','mod_vm');
    use_cb      = get('use_crack_band',true);
    cb_method   = char(lower(string(get('crack_band_method','oliver_directional'))));
    principal_dir_iters = get('principal_dir_iters',10);
    theta_opt   = get('theta_end',[]);
    nsnap       = get('n_snapshots',6);
    save_peak_picture = get('save_peak_picture',true);
    damage_save_levels = get('damage_save_levels',[]);
    damage_save_levels = sort(damage_save_levels(:));
    damage_level_saved = false(numel(damage_save_levels),1);
    save_peak_multiview = get('save_peak_multiview',true);
    AB          = get('cmod_AB',[]);
    a_arm       = get('lever_arm',100.0);
    Fmax        = get('Fmax',2000);
    slice_trig  = get('slice_trigger', inf);
    ninc        = get('n_increments',240);
    max_iter    = get('max_iter',20);
    tolD        = get('tol_damage',2e-4);
    tolR        = get('tol_residual',1e-5);
    check_res_every = get('check_residual_every',0);
    kmin_ratio  = get('min_stiffness_ratio',1e-6);
    relaxD      = get('relax_damage',1.0);
    viz_every   = get('viz_every',5);
    rotate_cam  = get('video_rotate_camera',true);
    cam_elev    = get('video_elevation',18);
    cam_revs    = get('video_revolutions',1.0);
    vid_fps     = get('video_framerate',20);
    cam_views   = get('camera_view_angles',[0 20;90 20;180 20;270 20;45 30;135 30;0 90]);
    cam_names   = get('camera_view_names',{'front','right','back','left','iso1','iso2','top'});
    save_mv     = get('save_multiview_snapshots',true);
    save_images = get('save_images',false);
    save_tables = get('save_tables',true);
    save_summary = get('save_summary_plots',false);
    save_live   = get('save_live_snapshots',true);
    n_live      = get('n_live_snaps',8);
    live_view   = get('live_snap_view',[45 30]);
    live_cap_fig= get('live_snap_capture_fig',false);

    nd = readmatrix([prefix_mesh '_nodes.txt']);
    ids = nd(:,1);
    p = nd(:,2:4);
    np = size(p,1);

    ed = readmatrix([prefix_mesh '_elements.txt']);
    [ok,T] = ismember(ed(:,2:5),ids);
    if any(~ok(:)), error('Some element connectivity node IDs were not found in node file.'); end
    ne = size(T,1);
    T = uint32(T);

    Nleft  = uint32(readmatrix([prefix_mesh '_left_nodes.txt']));
    Nright = uint32(readmatrix([prefix_mesh '_right_nodes.txt']));

    if ~isempty(AB)
        XA = AB(1,:); XB = AB(2,:);
        idA = nearest_node(p, XA);
        idB = nearest_node(p, XB);
        rAB0 = p(idB,:) - p(idA,:);
        nAB  = rAB0 / max(norm(rAB0),eps);
    else
        idA = []; idB = []; nAB = [0 0 1];
    end

    fprintf('Precomputing TET4 B matrices, gradients, and element stiffness...\n');
    [B3,V3,gradN3,h_min] = precompute_TET4_vector_fast(p,T);
    [Dunit,~,~] = iso_3D_D(1.0,nu);
    [D0,~,~]    = iso_3D_D(E,nu);

    Ke_unit3 = pagemtimes(permute(B3,[2 1 3]), pagemtimes(repmat(Dunit,[1 1 ne]),B3));
    Ke_unit3 = bsxfun(@times,Ke_unit3,reshape(V3,1,1,ne));

    ndof = 3*np;
    Tdouble = double(T).';
    EDOF = uint32(reshape([3*Tdouble(:)-2, 3*Tdouble(:)-1, 3*Tdouble(:)].',12,ne));
    [Iglob,Jglob] = make_sparse_indices(EDOF);
    idx_glob = double(EDOF(:));

    eps0  = ft/E;
    j0    = eps0;
    kappa = max(kappa0,j0)*ones(ne,1);
    De    = zeros(ne,1);
    h_band = ones(ne,1);
    ef_param = GF./(h_band*ft) + eps0/2;
    ef_param = max(ef_param, eps0 + 1e-12);

    fprintf('Mesh: nodes=%d | elements=%d | dof=%d | min edge=%.4g mm\n', np, ne, ndof, h_min);
    fprintf('Static setup: increments=%d | max_iter=%d | eps0=%.4e | model=%s | crack_band=%d | method=%s\n', ...
        ninc, max_iter, eps0, eeq_model, use_cb, cb_method);
    fprintf('Oliver bandwidth formula: h_n = 2 / sum_a |grad(N_a) dot n_crack|\n');

    dof = @(nid,comp) 3*(double(nid)-1)+comp;

    switch lower(path)
        case 'torsion_brokenshire'
            fix_left    = [dof(Nleft,1); dof(Nleft,2); dof(Nleft,3)];
            presc_right = [dof(Nright,1); dof(Nright,2); dof(Nright,3)];
            bc_dofs     = unique([fix_left; presc_right]);
            free_dofs   = setdiff((1:ndof).', bc_dofs);

            yr = p(double(Nright),2);
            zr = p(double(Nright),3);
            yc = mean(yr); zc = mean(zr);
            Ravg = mean(sqrt((yr-yc).^2 + (zr-zc).^2));
            if isempty(theta_opt)
                theta_end = Uy_end/max(Ravg,eps);
            else
                theta_end = theta_opt;
            end

        case 'torsion_forces'
            fix_left  = [dof(Nleft,1); dof(Nleft,2); dof(Nleft,3)];
            bc_dofs   = unique(fix_left);
            free_dofs = setdiff((1:ndof).', bc_dofs);

            yr = p(double(Nright),2); yc = mean(yr);
            Nright_pos = Nright(yr >= yc);
            Nright_neg = Nright(yr <  yc);
            theta_end = 0;

        otherwise
            error('Unknown load_path: %s', path);
    end

    theta_hist  = zeros(ninc,1);
    torque_hist = zeros(ninc,1);
    cmod_hist   = zeros(ninc,1);
    fapp_hist   = zeros(ninc,1);
    maxD_hist   = zeros(ninc,1);
    iter_hist   = zeros(ninc,1);
    res_hist    = zeros(ninc,1);
    hmean_hist  = zeros(ninc,1);
    hmin_hist   = zeros(ninc,1);
    hmax_hist   = zeros(ninc,1);

    [extFaces, extOwner] = external_hull(double(T));

    fig = []; ax1 = []; ax2 = []; hHull = []; hCrack = []; hCurve = []; v = [];
    if do_viz
        try
            fig = figure('Color','w','Name','LIVE Abaqus-style Torsion Damage Viewer - Oliver bandwidth','Position',[60 60 1250 620]);
            tl = tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');

            ax1 = nexttile(tl,1); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
            axis(ax1,'equal'); axis(ax1,'vis3d'); daspect(ax1,[1 1 1]);
            camproj(ax1,'perspective'); camup(ax1,[0 0 1]);
            title(ax1,'Quasi-static torsion damage - Oliver bandwidth', 'Interpreter','none');
            colormap(ax1,parula(256)); caxis(ax1,[0 1]); colorbar(ax1);
            mins = min(p,[],1); maxs = max(p,[],1); pad = 0.05*(maxs-mins+eps);
            xlim(ax1,[mins(1)-pad(1), maxs(1)+pad(1)]);
            ylim(ax1,[mins(2)-pad(2), maxs(2)+pad(2)]);
            zlim(ax1,[mins(3)-pad(3), maxs(3)+pad(3)]);
            view(ax1,cam_views(1,:)); camva(ax1,8); camlight(ax1,'headlight'); lighting(ax1,'gouraud');
            [hHull,hCrack] = init_live_3d(ax1, double(p), extFaces);

            ax2 = nexttile(tl,2); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
            xlabel(ax2,'Twist angle theta [rad]'); ylabel(ax2,'Reaction torque [N-mm]');
            hCurve = plot(ax2,0,0,'k-','LineWidth',1.8);
            title(ax2,'Torque-twist curve');

            if save_video
                try
                    v = VideoWriter([prefix_out '_animation.mp4'],'MPEG-4');
                    v.FrameRate = vid_fps; v.Quality = 95; open(v);
                catch MEv
                    warning('Video disabled: %s', MEv.message);
                    v = [];
                end
            end
        catch ME
            warning('Visualization disabled: %s', ME.message);
            fig = []; ax1 = []; ax2 = []; v = [];
        end
    end

    slice_saved = false;
    if save_images && nsnap > 0
        save_steps = unique(round(linspace(1,ninc,nsnap)));
    else
        save_steps = [];
    end

    if save_live && ~isempty(fig) && n_live > 0
        live_snap_steps = unique(round(linspace(1,ninc,n_live)));
    else
        live_snap_steps = [];
    end

    U = zeros(ndof,1);
    perm = [];

    peak = struct('absTorque',-inf,'inc',0,'theta',0,'torque',0,'fapp',0,'maxD',0,'U',[],'De',[],'h_band',[]);

    fprintf('Starting static increments...\n');
    tic;
    for inc = 1:ninc
        lambda_load = inc/ninc;
        current_theta = 0;

        Ubc = zeros(ndof,1);
        Fext = zeros(ndof,1);

        switch lower(path)
            case 'torsion_brokenshire'
                current_theta = lambda_load*theta_end;
                Uy_rot = -current_theta*(zr - zc);
                Uz_rot = +current_theta*(yr - yc);

                Ubc(fix_left) = 0;
                Ubc(dof(Nright,1)) = 0;
                Ubc(dof(Nright,2)) = Uy_rot;
                Ubc(dof(Nright,3)) = Uz_rot;

            case 'torsion_forces'
                Fnow = Fmax*lambda_load;
                if ~isempty(Nright_pos), Fext(dof(Nright_pos,2)) = +Fnow/numel(Nright_pos); end
                if ~isempty(Nright_neg), Fext(dof(Nright_neg,2)) = -Fnow/numel(Nright_neg); end
                Ubc(fix_left) = 0;
        end

        U(bc_dofs) = Ubc(bc_dofs);

        rel_res = inf;
        dDmax = inf;
        for it = 1:max_iter
            De_start_iter = De;

            K = assemble_global_K(Iglob,Jglob,Ke_unit3,E,De,kmin_ratio,ndof);

            Kfb = K(free_dofs,bc_dofs);
            rhs = Fext(free_dofs) - Kfb*Ubc(bc_dofs);
            Kff = K(free_dofs,free_dofs);

            if isempty(perm)
                perm = amd(Kff);
            end
            Kp = Kff(perm,perm);
            rp = rhs(perm);
            [R,flag] = chol(Kp);
            if flag == 0
                xp = R \ (R.' \ rp);
                uf = zeros(numel(free_dofs),1);
                uf(perm) = xp;
                U(free_dofs) = uf;
            else
                U(free_dofs) = Kff \ rhs;
            end
            U(bc_dofs) = Ubc(bc_dofs);

            ue12 = reshape(U(EDOF),12,1,ne);
            epsv = pagemtimes(B3,ue12);
            eeq = compute_equivalent_strain(eeq_model,epsv,Dunit,D0,E,nu,k_tc);
            eeq = reshape(eeq,[],1);

            if use_cb
                switch cb_method
                    case {'oliver','oliver_directional','directional'}
                        [h_band,~,~] = oliver_bandwidth_TET4_from_strain(epsv,gradN3,principal_dir_iters);
                    otherwise
                        error('Unknown crack_band_method: %s', cb_method);
                end
                ef_param = GF./(h_band(:)*ft) + eps0/2;
            else
                h_band = ones(ne,1);
                ef_param = GF./(1.0*ft) + eps0/2;
            end
            ef_param = max(ef_param, eps0 + 1e-12);

            kappa_trial = max(kappa,eeq);
            De_trial = damage_exp_stress_strain(kappa_trial,eps0,ef_param);

            De_trial = max(De,De_trial);
            De = (1-relaxD)*De + relaxD*De_trial;
            De = min(max(De,0),0.999999);
            kappa = max(kappa,kappa_trial);

            dDmax = max(abs(De - De_start_iter));

            if check_res_every > 0 && (mod(it,check_res_every)==0 || dDmax < tolD)
                s_e = E*max(1-De,kmin_ratio);
                fe12 = pagemtimes(Ke_unit3,ue12);
                fe12 = bsxfun(@times,fe12,reshape(s_e,1,1,ne));
                Fint_iter = assemble_accum(idx_glob,fe12(:),ndof);
                Rfree = Fext(free_dofs) - Fint_iter(free_dofs);
                denom = max(norm(Fext(free_dofs)), max(norm(Fint_iter(free_dofs)),1.0));
                rel_res = norm(Rfree)/denom;
            end

            if dDmax < tolD || (~isnan(rel_res) && rel_res < tolR)
                break;
            end
        end

        iter_hist(inc) = it;
        res_hist(inc) = rel_res;
        hmean_hist(inc) = mean(h_band);
        hmin_hist(inc)  = min(h_band);
        hmax_hist(inc)  = max(h_band);

        s_e = E*max(1-De,kmin_ratio);
        ue12 = reshape(U(EDOF),12,1,ne);
        fe12 = pagemtimes(Ke_unit3,ue12);
        fe12 = bsxfun(@times,fe12,reshape(s_e,1,1,ne));
        Fint = assemble_accum(idx_glob,fe12(:),ndof);

        switch lower(path)
            case 'torsion_brokenshire'
                right_dofs = [dof(Nright,1); dof(Nright,2); dof(Nright,3)];
                [Tx,~] = torque_about_x(p,right_dofs,Fint-Fext);
            case 'torsion_forces'
                right_dofs = [dof(Nright,1); dof(Nright,2); dof(Nright,3)];
                [Tx,~] = torque_about_x(p,right_dofs,Fext);
                current_theta = estimate_theta_from_right_end(p,U,Nright);
        end

        theta_hist(inc) = current_theta;
        torque_hist(inc) = Tx;
        maxD_hist(inc) = max(De);
        fapp_hist(inc) = Tx/a_arm;

        if ~isempty(idA)
            uA = [U(dof(idA,1)), U(dof(idA,2)), U(dof(idA,3))];
            uB = [U(dof(idB,1)), U(dof(idB,2)), U(dof(idB,3))];
            cmod_hist(inc) = dot((uB-uA),nAB);
        end

        if abs(Tx) > peak.absTorque
            peak.absTorque = abs(Tx);
            peak.inc = inc;
            peak.theta = current_theta;
            peak.torque = Tx;
            peak.fapp = fapp_hist(inc);
            peak.maxD = maxD_hist(inc);
            peak.U = U;
            peak.De = De;
            peak.h_band = h_band;
        end

        for ilev = 1:numel(damage_save_levels)
            if save_images && ~damage_level_saved(ilev) && maxD_hist(inc) >= damage_save_levels(ilev)
                lvl = damage_save_levels(ilev);
                save_damage_level_image(prefix_out,lvl,inc,current_theta,p,double(T),U,De,extFaces,extOwner, ...
                    get('defscale_3d',1.0),cam_views,cam_names,save_mv);
                damage_level_saved(ilev) = true;
            end
        end

        if save_images && ~slice_saved && current_theta >= slice_trig
            save_slice_figure(prefix_out,current_theta,p,double(T),U,De,get('crack_threshold',0.9),get('defscale_3d',1.0));
            slice_saved = true;
        end

        if save_images && ismember(inc,save_steps)
            save_damage_image(prefix_out,inc,current_theta,p,double(T),U,De,extFaces,extOwner, ...
                get('defscale_3d',1.0),cam_views,cam_names,save_mv);
        end

        is_snap = ismember(inc,live_snap_steps);
        if ~isempty(fig) && isvalid(fig) && (mod(inc,viz_every)==0 || inc==1 || inc==ninc || is_snap)
            set(hCurve,'XData',theta_hist(1:inc),'YData',torque_hist(1:inc));
            update_live_3d(ax1,hHull,hCrack,double(p),double(U),get('defscale_3d',1.0),double(T),double(De), ...
                get('crack_threshold',0.9),extFaces,extOwner);
            title(ax1,sprintf('inc=%d/%d | theta=%.3e | maxD=%.3f | it=%d | mean h_O=%.3g mm', ...
                inc,ninc,current_theta,maxD_hist(inc),it,hmean_hist(inc)));

            if rotate_cam
                az = 360*cam_revs*(inc-1)/max(ninc-1,1);
                view(ax1,[az cam_elev]);
            end
            drawnow limitrate nocallbacks;
            if ~isempty(v)
                try, writeVideo(v,getframe(fig)); catch, end
            end

            if is_snap
                if live_cap_fig
                    save_live_snapshot(fig,prefix_out,inc,current_theta,maxD_hist(inc),[]);
                else
                    save_live_snapshot(ax1,prefix_out,inc,current_theta,maxD_hist(inc),live_view);
                end
            end
        end

        if mod(inc,max(1,round(ninc/20)))==0 || inc==1 || inc==ninc
            fprintf('inc %4d/%4d | theta=% .4e | Tx=% .4e | maxD=%.4f | it=%2d | dD=%.2e | hOmean=%.3g | res=%.2e\n', ...
                inc,ninc,current_theta,Tx,maxD_hist(inc),it,dDmax,hmean_hist(inc),rel_res);
        end

        if maxD_hist(inc) > 0.995 && inc > 10
            fprintf('Severe damage at inc %d: maxD = %.4f, continuing to full theta_end...\n', inc,maxD_hist(inc));
        end
    end
    elapsed = toc;
    fprintf('Static solve time: %.2f s\n',elapsed);

    if ~isempty(v)
        try, close(v); fprintf('Saved video: %s_animation.mp4\n',prefix_out); catch, end
    end

    Tout = table(theta_hist(:),torque_hist(:),cmod_hist(:),fapp_hist(:),maxD_hist(:),iter_hist(:),res_hist(:), ...
        hmean_hist(:),hmin_hist(:),hmax_hist(:), ...
        'VariableNames',{'theta_rad','Torque_Nmm','CMOD_mm','Fapp_N','max_damage','iterations','relative_residual', ...
                         'mean_h_oliver_mm','min_h_oliver_mm','max_h_oliver_mm'});
    if save_tables
        writetable(Tout,[prefix_out '_T_theta_CMOD_static.csv']);
    end

    if save_summary
        save_summary_plots(prefix_out,theta_hist,torque_hist,cmod_hist,fapp_hist,maxD_hist,iter_hist,res_hist);
    end

    if save_peak_picture && ~isempty(peak.De)
        fprintf('Peak load at inc=%d: theta=%.4e rad, torque=%.4e N-mm, P=%.4e N, maxD=%.4f\n', ...
            peak.inc,peak.theta,peak.torque,peak.fapp,peak.maxD);
        save_peak_damage_image(prefix_out,peak,p,double(T),extFaces,extOwner, ...
            get('defscale_3d',1.0),cam_views,cam_names,save_peak_multiview);
        save_peak_curve_marker(prefix_out,theta_hist,torque_hist,fapp_hist,peak);
        [~,name_only] = fileparts(prefix_out);
        writematrix(peak.De,fullfile(fileparts(prefix_out),sprintf('%s_PEAK_DAMAGE_VECTOR.csv',name_only)));
        if ~isempty(peak.h_band)
            writematrix(peak.h_band,fullfile(fileparts(prefix_out),sprintf('%s_PEAK_OLIVER_BANDWIDTH_VECTOR.csv',name_only)));
        end
    end

    if save_images
        save_damage_image(prefix_out,numel(theta_hist),theta_hist(end),p,double(T),U,De,extFaces,extOwner, ...
            get('defscale_3d',1.0),cam_views,cam_names,save_mv);
    end

    if get('save_full_results',false)
        save([prefix_out '_final_state.mat'],'p','T','U','De','kappa','theta_hist','torque_hist','cmod_hist','maxD_hist', ...
            'iter_hist','res_hist','hmean_hist','hmin_hist','hmax_hist','h_band','opts','-v7.3');
    end
end

function [h_band,ncrack,lambda_max] = oliver_bandwidth_TET4_from_strain(epsv,gradN3,n_power_iter)

    if nargin < 3 || isempty(n_power_iter)
        n_power_iter = 10;
    end

    [ncrack,lambda_max] = max_principal_strain_direction_power_vec(epsv,n_power_iter);
    nx = ncrack(:,1); ny = ncrack(:,2); nz = ncrack(:,3);

    g1x = reshape(gradN3(1,1,:),[],1); g1y = reshape(gradN3(1,2,:),[],1); g1z = reshape(gradN3(1,3,:),[],1);
    g2x = reshape(gradN3(2,1,:),[],1); g2y = reshape(gradN3(2,2,:),[],1); g2z = reshape(gradN3(2,3,:),[],1);
    g3x = reshape(gradN3(3,1,:),[],1); g3y = reshape(gradN3(3,2,:),[],1); g3z = reshape(gradN3(3,3,:),[],1);
    g4x = reshape(gradN3(4,1,:),[],1); g4y = reshape(gradN3(4,2,:),[],1); g4z = reshape(gradN3(4,3,:),[],1);

    den = abs(g1x.*nx + g1y.*ny + g1z.*nz) + ...
          abs(g2x.*nx + g2y.*ny + g2z.*nz) + ...
          abs(g3x.*nx + g3y.*ny + g3z.*nz) + ...
          abs(g4x.*nx + g4y.*ny + g4z.*nz);

    h_band = 2.0 ./ max(den,1e-14);
    bad = ~isfinite(h_band) | h_band <= 0;
    if any(bad)
        h_band(bad) = median(h_band(~bad));
        if any(~isfinite(h_band) | h_band <= 0)
            h_band(~isfinite(h_band) | h_band <= 0) = 1.0;
        end
    end
end

function [nvec,lambda_max] = max_principal_strain_direction_power_vec(epsv,niter)

    exx = reshape(epsv(1,1,:),[],1);
    eyy = reshape(epsv(2,1,:),[],1);
    ezz = reshape(epsv(3,1,:),[],1);
    exy = 0.5*reshape(epsv(4,1,:),[],1);
    eyz = 0.5*reshape(epsv(5,1,:),[],1);
    ezx = 0.5*reshape(epsv(6,1,:),[],1);
    ne = numel(exx);

    frob = sqrt(exx.^2 + eyy.^2 + ezz.^2 + 2*(exy.^2 + eyz.^2 + ezx.^2));
    shift = frob + 1e-30;

    nx = ones(ne,1)/sqrt(3);
    ny = ones(ne,1)/sqrt(3);
    nz = ones(ne,1)/sqrt(3);

    for k = 1:max(3,niter)
        wx = (exx+shift).*nx + exy.*ny + ezx.*nz;
        wy = exy.*nx + (eyy+shift).*ny + eyz.*nz;
        wz = ezx.*nx + eyz.*ny + (ezz+shift).*nz;
        wn = sqrt(wx.^2 + wy.^2 + wz.^2);
        wn = max(wn,1e-300);
        nx = wx./wn; ny = wy./wn; nz = wz./wn;
    end

    Ax = exx.*nx + exy.*ny + ezx.*nz;
    Ay = exy.*nx + eyy.*ny + eyz.*nz;
    Az = ezx.*nx + eyz.*ny + ezz.*nz;
    lambda_max = nx.*Ax + ny.*Ay + nz.*Az;

    bad = ~isfinite(nx) | ~isfinite(ny) | ~isfinite(nz);
    if any(bad)
        nx(bad) = 1; ny(bad) = 0; nz(bad) = 0;
        lambda_max(bad) = 0;
    end
    nvec = [nx,ny,nz];
end

function De = damage_exp_stress_strain(kappa,eps0,ef_param)
    De = zeros(size(kappa));
    act = kappa > eps0;
    if any(act)
        tau_decay = ef_param(act) - eps0;
        tau_decay = max(tau_decay,1e-12);
        De(act) = 1 - (eps0./kappa(act)).*exp(-(kappa(act)-eps0)./tau_decay);
    end
    De = min(max(De,0),0.999999);
end

function K = assemble_global_K(Iglob,Jglob,Ke_unit3,E,De,kmin_ratio,ndof)
    s_e = E*max(1-De,kmin_ratio);
    Ke = bsxfun(@times,Ke_unit3,reshape(s_e,1,1,[]));
    K = sparse(Iglob,Jglob,Ke(:),ndof,ndof);
    K = 0.5*(K+K.');
end

function [Iglob,Jglob] = make_sparse_indices(EDOF)
    ed = double(EDOF);
    Iglob = repmat(ed,12,1);
    Jglob = repelem(ed,12,1);
    Iglob = Iglob(:);
    Jglob = Jglob(:);
end

function Fglob = assemble_accum(idx_glob,vals_col,ndof)
    Fglob = accumarray(idx_glob,vals_col,[ndof,1],@sum,0);
end

function eeq = compute_equivalent_strain(model,epsv,Dunit,D0,E,nu,k_tc)
    switch lower(model)
        case 'mazars'
            eeq = eqv_strain_mazars_vec(epsv);
        case 'mod_vm'
            eeq = eqv_strain_mod_vm_vec(epsv,nu,k_tc);
        case 'energy_norm'
            eeq = eqv_strain_energy_norm_vec(epsv,Dunit);
        case 'rankine'
            eeq = eqv_strain_rankine_vec(epsv,D0,E);
        case 'smooth_rankine'
            eeq = eqv_strain_smooth_rankine_vec(epsv,D0,E);
        otherwise
            error('Unknown equivalent strain model: %s',model);
    end
end

function eeq = eqv_strain_mazars_vec(epsv)
    exx = squeeze(epsv(1,1,:)); eyy = squeeze(epsv(2,1,:)); ezz = squeeze(epsv(3,1,:));
    exy = 0.5*squeeze(epsv(4,1,:)); eyz = 0.5*squeeze(epsv(5,1,:)); ezx = 0.5*squeeze(epsv(6,1,:));
    [l1,l2,l3] = eig3x3_sym_vec(exx,eyy,ezz,exy,eyz,ezx);
    p1 = max(0,l1); p2 = max(0,l2); p3 = max(0,l3);
    eeq = reshape(sqrt(p1.^2 + p2.^2 + p3.^2),1,1,[]);
end

function eeq = eqv_strain_mod_vm_vec(epsv,nu,k)
    exx = squeeze(epsv(1,1,:)); eyy = squeeze(epsv(2,1,:)); ezz = squeeze(epsv(3,1,:));
    gxy = squeeze(epsv(4,1,:)); gyz = squeeze(epsv(5,1,:)); gzx = squeeze(epsv(6,1,:));
    exy = 0.5*gxy; eyz = 0.5*gyz; ezx = 0.5*gzx;

    I1 = exx + eyy + ezz;
    sxx = exx - I1/3; syy = eyy - I1/3; szz = ezz - I1/3;
    J2 = 0.5*(sxx.^2 + syy.^2 + szz.^2 + 2*(exy.^2 + eyz.^2 + ezx.^2));
    a1 = (k-1)./(1-2*nu);
    a2 = 12*k./(1+nu).^2;
    term1 = (k-1)./(2*k*(1-2*nu)).*I1;
    term2 = (1./(2*k)).*sqrt((a1.^2).*I1.^2 + a2.*J2);
    eeq = reshape(max(0,term1+term2),1,1,[]);
end

function eeq = eqv_strain_energy_norm_vec(epsv,Dunit)
    term1 = pagemtimes(epsv,'transpose',Dunit,'none');
    enorm = pagemtimes(term1,'none',epsv,'none');
    eeq = sqrt(max(enorm,0));
end

function eeq = eqv_strain_rankine_vec(epsv,D0,E)
    E6 = reshape(epsv,6,[]);
    SV = D0*E6;
    [l1,~,~] = eig3x3_sym_vec(SV(1,:).',SV(2,:).',SV(3,:).',SV(4,:).',SV(5,:).',SV(6,:).');
    eeq = reshape(max(0,l1)/E,1,1,[]);
end

function eeq = eqv_strain_smooth_rankine_vec(epsv,D0,E)
    E6 = reshape(epsv,6,[]);
    SV = D0*E6;
    [l1,l2,l3] = eig3x3_sym_vec(SV(1,:).',SV(2,:).',SV(3,:).',SV(4,:).',SV(5,:).',SV(6,:).');
    p1 = max(0,l1); p2 = max(0,l2); p3 = max(0,l3);
    eeq = reshape(sqrt(p1.^2 + p2.^2 + p3.^2)/E,1,1,[]);
end

function [l1,l2,l3] = eig3x3_sym_vec(a11,a22,a33,a12,a23,a13)
    a11=a11(:); a22=a22(:); a33=a33(:); a12=a12(:); a23=a23(:); a13=a13(:);
    p1 = a12.^2 + a13.^2 + a23.^2;
    q  = (a11 + a22 + a33)/3;
    p2 = (a11-q).^2 + (a22-q).^2 + (a33-q).^2 + 2*p1;
    p  = sqrt(max(p2,0)/6);
    pe = max(p,1e-300);
    b11=(a11-q)./pe; b22=(a22-q)./pe; b33=(a33-q)./pe;
    b12=a12./pe;     b13=a13./pe;     b23=a23./pe;
    detB = b11.*(b22.*b33 - b23.^2) - b12.*(b12.*b33 - b23.*b13) + b13.*(b12.*b23 - b22.*b13);
    r = min(max(detB/2,-1),1);
    phi = acos(r)/3;
    l1 = q + 2*p.*cos(phi);
    l3 = q + 2*p.*cos(phi + 2*pi/3);
    l2 = 3*q - l1 - l3;
    diagmask = p1 < 1e-300;
    if any(diagmask)
        ld = sort([a11(diagmask), a22(diagmask), a33(diagmask)],2,'descend');
        l1(diagmask) = ld(:,1); l2(diagmask) = ld(:,2); l3(diagmask) = ld(:,3);
    end
end

function [B3,V3,gradN3,hmin] = precompute_TET4_vector_fast(nodes,tets)
    ne = size(tets,1);
    n1 = double(tets(:,1)); n2 = double(tets(:,2)); n3 = double(tets(:,3)); n4 = double(tets(:,4));
    x1 = nodes(n1,:); x2 = nodes(n2,:); x3 = nodes(n3,:); x4 = nodes(n4,:);

    a = x2 - x1; b = x3 - x1; c = x4 - x1;
    bc = cross(b,c,2); ca = cross(c,a,2); ab = cross(a,b,2);

    detJ = sum(a.*bc,2);
    V = abs(detJ)/6;
    if any(V <= 1e-16)
        bad = find(V <= 1e-16,1,'first');
        error('Non-positive or extremely small volume at element %d.',bad);
    end

    g2 = bc ./ detJ; g3 = ca ./ detJ; g4 = ab ./ detJ;
    g1 = -(g2 + g3 + g4);

    B3 = zeros(6,12,ne);

    B3(1,1,:)=reshape(g1(:,1),1,1,ne); B3(2,2,:)=reshape(g1(:,2),1,1,ne); B3(3,3,:)=reshape(g1(:,3),1,1,ne);
    B3(4,1,:)=reshape(g1(:,2),1,1,ne); B3(4,2,:)=reshape(g1(:,1),1,1,ne);
    B3(5,2,:)=reshape(g1(:,3),1,1,ne); B3(5,3,:)=reshape(g1(:,2),1,1,ne);
    B3(6,1,:)=reshape(g1(:,3),1,1,ne); B3(6,3,:)=reshape(g1(:,1),1,1,ne);

    B3(1,4,:)=reshape(g2(:,1),1,1,ne); B3(2,5,:)=reshape(g2(:,2),1,1,ne); B3(3,6,:)=reshape(g2(:,3),1,1,ne);
    B3(4,4,:)=reshape(g2(:,2),1,1,ne); B3(4,5,:)=reshape(g2(:,1),1,1,ne);
    B3(5,5,:)=reshape(g2(:,3),1,1,ne); B3(5,6,:)=reshape(g2(:,2),1,1,ne);
    B3(6,4,:)=reshape(g2(:,3),1,1,ne); B3(6,6,:)=reshape(g2(:,1),1,1,ne);

    B3(1,7,:)=reshape(g3(:,1),1,1,ne); B3(2,8,:)=reshape(g3(:,2),1,1,ne); B3(3,9,:)=reshape(g3(:,3),1,1,ne);
    B3(4,7,:)=reshape(g3(:,2),1,1,ne); B3(4,8,:)=reshape(g3(:,1),1,1,ne);
    B3(5,8,:)=reshape(g3(:,3),1,1,ne); B3(5,9,:)=reshape(g3(:,2),1,1,ne);
    B3(6,7,:)=reshape(g3(:,3),1,1,ne); B3(6,9,:)=reshape(g3(:,1),1,1,ne);

    B3(1,10,:)=reshape(g4(:,1),1,1,ne); B3(2,11,:)=reshape(g4(:,2),1,1,ne); B3(3,12,:)=reshape(g4(:,3),1,1,ne);
    B3(4,10,:)=reshape(g4(:,2),1,1,ne); B3(4,11,:)=reshape(g4(:,1),1,1,ne);
    B3(5,11,:)=reshape(g4(:,3),1,1,ne); B3(5,12,:)=reshape(g4(:,2),1,1,ne);
    B3(6,10,:)=reshape(g4(:,3),1,1,ne); B3(6,12,:)=reshape(g4(:,1),1,1,ne);

    V3  = reshape(V,1,1,ne);

    gradN3 = zeros(4,3,ne);
    gradN3(1,1,:) = reshape(g1(:,1),1,1,ne); gradN3(1,2,:) = reshape(g1(:,2),1,1,ne); gradN3(1,3,:) = reshape(g1(:,3),1,1,ne);
    gradN3(2,1,:) = reshape(g2(:,1),1,1,ne); gradN3(2,2,:) = reshape(g2(:,2),1,1,ne); gradN3(2,3,:) = reshape(g2(:,3),1,1,ne);
    gradN3(3,1,:) = reshape(g3(:,1),1,1,ne); gradN3(3,2,:) = reshape(g3(:,2),1,1,ne); gradN3(3,3,:) = reshape(g3(:,3),1,1,ne);
    gradN3(4,1,:) = reshape(g4(:,1),1,1,ne); gradN3(4,2,:) = reshape(g4(:,2),1,1,ne); gradN3(4,3,:) = reshape(g4(:,3),1,1,ne);

    e12 = vecnorm(x1-x2,2,2); e13 = vecnorm(x1-x3,2,2); e14 = vecnorm(x1-x4,2,2);
    e23 = vecnorm(x2-x3,2,2); e24 = vecnorm(x2-x4,2,2); e34 = vecnorm(x3-x4,2,2);
    hmin = min([e12; e13; e14; e23; e24; e34]);
end

function [D,lambda,G] = iso_3D_D(E,nu)
    G = E/(2*(1+nu));
    lambda = E*nu/((1+nu)*(1-2*nu));
    D = [lambda+2*G, lambda,     lambda,     0, 0, 0;
         lambda,     lambda+2*G, lambda,     0, 0, 0;
         lambda,     lambda,     lambda+2*G, 0, 0, 0;
         0,          0,          0,          G, 0, 0;
         0,          0,          0,          0, G, 0;
         0,          0,          0,          0, 0, G];
end

function [Tx,Rvec] = torque_about_x(p,dof_idx,F)
    nid = unique(ceil(double(dof_idx)/3));
    Fy = F(3*(nid-1)+2);
    Fz = F(3*(nid-1)+3);
    r = p(nid,:);
    y = r(:,2); z = r(:,3);
    Tx = sum(y.*Fz - z.*Fy);
    if nargout > 1, Rvec = [y,Fz,z,Fy]; end
end

function theta = estimate_theta_from_right_end(p,U,Nright)
    nid = double(Nright(:));
    y = p(nid,2); z = p(nid,3);
    yc = mean(y); zc = mean(z);
    uy = U(3*(nid-1)+2); uz = U(3*(nid-1)+3);
    r2 = (y-yc).^2 + (z-zc).^2;
    theta_i = ((y-yc).*uz - (z-zc).*uy)./max(r2,eps);
    theta = mean(theta_i(isfinite(theta_i)));
end

function id = nearest_node(P,Xq)
    [~,id] = min(vecnorm(P-Xq,2,2));
end

function [extFaces,extOwner] = external_hull(T)
    allF = [T(:,[1 2 3]); T(:,[1 2 4]); T(:,[1 3 4]); T(:,[2 3 4])];
    own = repelem((1:size(T,1)).',4,1);
    sF = sort(allF,2);
    [uf,~,ic] = unique(sF,'rows');
    cnt = accumarray(ic,1);
    owner_per_uf = accumarray(ic,own,[],@(x)x(1));
    mask = cnt == 1;
    extFaces = uf(mask,:);
    extOwner = owner_per_uf(mask);
end

function val = local_get(s,field,def)
    if ~isfield(s,field) || isempty(s.(field))
        val = def;
    else
        val = s.(field);
    end
end

function [hHull,hCrack] = init_live_3d(ax,p,extFaces)
    hHull = patch(ax,'Faces',extFaces,'Vertices',p,'FaceColor','flat', ...
        'FaceVertexCData',zeros(size(extFaces,1),1), ...
        'EdgeColor','none','FaceAlpha',0.30,'Visible','on');
    hCrack = patch(ax,'Faces',zeros(0,3),'Vertices',p,'FaceColor','flat', ...
        'FaceVertexCData',zeros(0,1), ...
        'EdgeColor','none','FaceAlpha',0.90,'Visible','off');
    caxis(ax,[0 1]);
end

function update_live_3d(ax,hHull,hCrack,p,U,defscale,T,De,thr,extFaces,extOwner)
    Ux = U(1:3:end); Uy = U(2:3:end); Uz = U(3:3:end);
    pdef = p + defscale*[Ux,Uy,Uz];
    set(hHull,'Vertices',pdef,'FaceVertexCData',De(extOwner),'Visible','on');

    hi = find(De >= thr);
    if isempty(hi), hi = find(De >= max(0.70,thr-0.20)); end
    if isempty(hi)
        set(hCrack,'Visible','off');
    else
        faces_hd = [T(hi,[1 2 3]); T(hi,[1 2 4]); T(hi,[1 3 4]); T(hi,[2 3 4])];
        own = repelem(hi,4,1);
        sF = sort(faces_hd,2);
        [uf,~,ic] = unique(sF,'rows');
        cnt = accumarray(ic,1);
        owner_per_uf = accumarray(ic,own,[],@(x)x(1));
        fb = uf(cnt==1,:);
        cface = De(owner_per_uf(cnt==1));
        set(hCrack,'Faces',fb,'Vertices',pdef,'FaceVertexCData',cface,'Visible','on');
    end
    caxis(ax,[0 1]);
end

function save_live_snapshot(handle,prefix,inc,theta,maxD,snapview)
    try
        is_ax = isa(handle,'matlab.graphics.axis.Axes');
        v_old = [];
        if is_ax && nargin >= 6 && ~isempty(snapview)
            v_old = get(handle,'View');
            view(handle,snapview);
        end
        drawnow;
        fn = sprintf('%s_LIVE_snap_inc_%04d_theta_%.3e.png',prefix,inc,theta);
        exportgraphics(handle,fn,'Resolution',200,'BackgroundColor','white');
        fprintf('Saved purple snapshot: %s (maxD=%.3f)\n',fn,maxD);
        if is_ax && ~isempty(v_old), view(handle,v_old); end
    catch ME
        warning('Live snapshot failed at inc %d: %s',inc,ME.message);
    end
end

function save_slice_figure(prefix,theta,p,T,U,De,thr,defscale)
    try
        hFig = figure('Color','w','Visible','off','Position',[100 100 1100 850]);
        ax = axes(hFig); hold(ax,'on'); axis(ax,'equal'); axis(ax,'vis3d'); grid(ax,'on'); box(ax,'on');
        colormap(ax,parula(256)); caxis(ax,[0 1]);
        Ux = U(1:3:end); Uy = U(2:3:end); Uz = U(3:3:end);
        Pdef = double(p) + defscale*[Ux,Uy,Uz];
        ne = size(T,1);
        Xe = reshape(double(p(T',:)),4,ne,3);
        centers = squeeze(mean(Xe,1));
        if size(centers,1) ~= ne, centers = centers.'; end
        keep_idx = centers(:,2) < 0;
        T_half = T(keep_idx,:); De_half = De(keep_idx);
        if isempty(T_half)
            [faces,owner] = external_hull(T); cdata = De(owner);
        else
            [faces,owner] = external_hull(T_half); cdata = De_half(owner);
        end
        patch(ax,'Faces',faces,'Vertices',Pdef,'FaceColor','flat','FaceVertexCData',cdata, ...
            'EdgeColor','none','FaceAlpha',1.0);
        if ~isempty(T_half)
            hi_local = find(De_half >= max(0.70,thr-0.20));
            if ~isempty(hi_local)
                faces_hd = [T_half(hi_local,[1 2 3]); T_half(hi_local,[1 2 4]); T_half(hi_local,[1 3 4]); T_half(hi_local,[2 3 4])];
                patch(ax,'Faces',faces_hd,'Vertices',Pdef,'FaceColor','none','EdgeColor',[0.1 0.1 0.1],'LineWidth',0.25);
            end
        end
        view(ax,[90 15]); camlight(ax,'headlight'); lighting(ax,'gouraud');
        title(ax,sprintf('Slice Y < 0 at theta = %.3e rad',theta),'Interpreter','none'); colorbar(ax);
        saveas(hFig,sprintf('%s_Slice_Theta_%.3e.png',prefix,theta));
        close(hFig);
    catch ME
        warning('Failed to save slice snapshot: %s',ME.message);
    end
end

function save_damage_image(prefix,inc,theta,p,T,U,De,extFaces,extOwner,defscale,cam_views,cam_names,save_multiview)
    try
        if nargin < 11 || isempty(cam_views), cam_views = [0 20;90 20;180 20;270 20;45 30;135 30;0 90]; end
        if nargin < 12 || isempty(cam_names), cam_names = {'front','right','back','left','iso1','iso2','top'}; end
        if nargin < 13 || isempty(save_multiview), save_multiview = true; end
        if numel(cam_names) < size(cam_views,1)
            for k = numel(cam_names)+1:size(cam_views,1), cam_names{k} = sprintf('view_%02d',k); end
        end
        hFig = figure('Color','w','Visible','off','Position',[100 100 1200 900]);
        ax = axes(hFig); hold(ax,'on'); axis(ax,'equal'); axis(ax,'vis3d'); grid(ax,'on'); box(ax,'on');
        colormap(ax,parula(256)); caxis(ax,[0 1]);
        Ux = U(1:3:end); Uy = U(2:3:end); Uz = U(3:3:end);
        Pdef = double(p) + defscale*[Ux,Uy,Uz];
        patch(ax,'Faces',extFaces,'Vertices',Pdef,'FaceColor','flat','FaceVertexCData',De(extOwner), ...
            'EdgeColor','none','FaceAlpha',0.95);
        camlight(ax,'headlight'); lighting(ax,'gouraud');
        title(ax,sprintf('Damage at increment %d, theta=%.3e rad, maxD=%.3f',inc,theta,max(De)),'Interpreter','none');
        colorbar(ax);
        view(ax,[45 30]); camva(ax,8); drawnow;
        saveas(hFig,sprintf('%s_damage_inc_%04d_iso.png',prefix,inc));
        if save_multiview
            for iv = 1:size(cam_views,1)
                view(ax,cam_views(iv,:)); camva(ax,8); drawnow;
                safe_name = regexprep(cam_names{iv},'[^A-Za-z0-9_\-]','_');
                saveas(hFig,sprintf('%s_damage_inc_%04d_%s.png',prefix,inc,safe_name));
            end
        end
        close(hFig);
    catch ME
        warning('Could not save damage image at inc %d: %s',inc,ME.message);
    end
end

function save_damage_level_image(prefix,lvl,inc,theta,p,T,U,De,extFaces,extOwner,defscale,cam_views,cam_names,save_multiview)
    try
        prefix_level = sprintf('%s_D%03d',prefix,round(100*lvl));
        save_damage_image(prefix_level,inc,theta,p,T,U,De,extFaces,extOwner,defscale,cam_views,cam_names,save_multiview);
    catch ME
        warning('Could not save damage-level image %.2f at inc %d: %s',lvl,inc,ME.message);
    end
end

function save_peak_damage_image(prefix,peak,p,T,extFaces,extOwner,defscale,cam_views,cam_names,save_multiview)
    try
        U = peak.U; De = peak.De;
        Ux = U(1:3:end); Uy = U(2:3:end); Uz = U(3:3:end);
        Pdef = double(p) + defscale*[Ux,Uy,Uz];
        hFig = figure('Color','w','Visible','off','Position',[100 100 1300 950]);
        ax = axes(hFig); hold(ax,'on'); axis(ax,'equal'); axis(ax,'vis3d'); grid(ax,'on'); box(ax,'on');
        colormap(ax,parula(256)); caxis(ax,[0 1]);
        patch(ax,'Faces',extFaces,'Vertices',Pdef,'FaceColor','flat','FaceVertexCData',De(extOwner), ...
            'EdgeColor','none','FaceAlpha',0.96);
        thr = 0.70; hi = find(De >= thr);
        if isempty(hi), hi = find(De >= max(0.50,max(De)-0.05)); end
        if ~isempty(hi)
            faces_hd = [T(hi,[1 2 3]); T(hi,[1 2 4]); T(hi,[1 3 4]); T(hi,[2 3 4])];
            own = repelem(hi,4,1);
            sF = sort(faces_hd,2);
            [uf,~,ic] = unique(sF,'rows');
            cnt = accumarray(ic,1);
            owner_per_uf = accumarray(ic,own,[],@(x)x(1));
            fb = uf(cnt==1,:); cface = De(owner_per_uf(cnt==1));
            patch(ax,'Faces',fb,'Vertices',Pdef,'FaceColor','flat','FaceVertexCData',cface, ...
                'EdgeColor',[0.05 0.05 0.05],'LineWidth',0.20,'FaceAlpha',1.0);
        end
        title(ax,sprintf('PEAK LOAD DAMAGE | inc=%d | theta=%.3e rad | T=%.3e N-mm | P=%.3e N | maxD=%.3f', ...
            peak.inc,peak.theta,peak.torque,peak.fapp,peak.maxD),'Interpreter','none');
        xlabel(ax,'X [mm]'); ylabel(ax,'Y [mm]'); zlabel(ax,'Z [mm]');
        colorbar(ax); camlight(ax,'headlight'); lighting(ax,'gouraud'); camva(ax,8);
        view(ax,[45 30]); drawnow;
        saveas(hFig,sprintf('%s_PEAK_LOAD_DAMAGE_iso.png',prefix));
        if save_multiview
            if numel(cam_names) < size(cam_views,1)
                for k = numel(cam_names)+1:size(cam_views,1), cam_names{k} = sprintf('view_%02d',k); end
            end
            for iv = 1:size(cam_views,1)
                view(ax,cam_views(iv,:)); camva(ax,8); drawnow;
                safe_name = regexprep(cam_names{iv},'[^A-Za-z0-9_\-]','_');
                saveas(hFig,sprintf('%s_PEAK_LOAD_DAMAGE_%s.png',prefix,safe_name));
            end
        end
        close(hFig);
    catch ME
        warning('Could not save peak-load damage image: %s',ME.message);
    end
end

function save_peak_curve_marker(prefix,theta,torque,fapp,peak)
    try
        h1 = figure('Color','w','Visible','off');
        plot(theta,torque,'k-','LineWidth',1.8); hold on; grid on; box on;
        plot(peak.theta,peak.torque,'ro','MarkerSize',8,'LineWidth',2.0);
        xlabel('Twist angle theta [rad]'); ylabel('Torque [N-mm]');
        title('Torque-angle response with peak-load marker');
        legend({'Torque-angle','Peak load'},'Location','best');
        saveas(h1,[prefix '_torque_theta_WITH_PEAK.png']); close(h1);

        h2 = figure('Color','w','Visible','off');
        plot(theta,fapp,'k-','LineWidth',1.8); hold on; grid on; box on;
        plot(peak.theta,peak.fapp,'ro','MarkerSize',8,'LineWidth',2.0);
        xlabel('Twist angle theta [rad]'); ylabel('Equivalent load P = T / a [N]');
        title('Load-angle response with peak-load marker');
        legend({'Load-angle','Peak load'},'Location','best');
        saveas(h2,[prefix '_load_theta_WITH_PEAK.png']); close(h2);
    catch ME
        warning('Could not save peak curve marker plots: %s',ME.message);
    end
end

function save_summary_plots(prefix,theta,torque,cmod,fapp,maxD,iters,resid)
    try
        h1 = figure('Color','w','Visible','off');
        plot(theta,torque,'k-','LineWidth',1.8); grid on; box on;
        xlabel('Twist angle theta [rad]'); ylabel('Torque [N-mm]'); title('Torque-twist response');
        saveas(h1,[prefix '_torque_theta.png']); close(h1);

        hload = figure('Color','w','Visible','off');
        plot(theta,fapp,'k-','LineWidth',1.8); grid on; box on;
        xlabel('Twist angle theta [rad]'); ylabel('Equivalent load P = T/a [N]'); title('Load-angle response');
        saveas(hload,[prefix '_load_theta.png']); close(hload);

        h2 = figure('Color','w','Visible','off');
        plot(cmod,torque,'k-','LineWidth',1.8); grid on; box on;
        xlabel('CMOD [mm]'); ylabel('Torque [N-mm]'); title('Torque-CMOD response');
        saveas(h2,[prefix '_torque_CMOD.png']); close(h2);

        h3 = figure('Color','w','Visible','off');
        plot(theta,maxD,'k-','LineWidth',1.8); grid on; box on;
        xlabel('Twist angle theta [rad]'); ylabel('Maximum damage'); title('Maximum damage evolution');
        ylim([0 1]);
        saveas(h3,[prefix '_max_damage.png']); close(h3);

        h4 = figure('Color','w','Visible','off');
        yyaxis left; plot(theta,iters,'-','LineWidth',1.5); ylabel('Iterations');
        yyaxis right; semilogy(theta,max(resid,eps),'-','LineWidth',1.5); ylabel('Relative residual');
        grid on; box on; xlabel('Twist angle theta [rad]'); title('Solver convergence history');
        saveas(h4,[prefix '_convergence.png']); close(h4);
    catch ME
        warning('Could not save summary plots: %s',ME.message);
    end
end
