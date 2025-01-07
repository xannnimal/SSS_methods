%% implementation testing leadfield calculated usin J. Sarvas cirrent dipole
% modify function given by Samu to take in a grid of dipole positions
% assume each grid point has unit mag in (x,y,z) direction, then need to do
% some regularization
clear all

%% load sensors from fif file
coordsys='device';
rawfile = 'sample_audvis_raw.fif';
[R,EX,EY,EZ] = fiff_getpos(rawfile, coordsys);
k=1;
for i=(1:size(EZ,2))
    if mod(i,3)==0 %every third is a magnetometer
        ch_types(i)=1;
        mags(k)=i;
        k=k+1;
    end
end

%field positions
x_sens = R(1,:)';
y_sens = R(2,:)';
z_sens = R(3,:)';

%% generate spherical grid with N_sph layers from rsph_min to rsph_max
spacing = 25;
N_sph = 7;
rsph_min = 0.01;
rsph_max = 0.07;

rsph = [rsph_min:(rsph_max-rsph_min)/(N_sph-1):rsph_max];

for nsph = 1:N_sph
    pos{nsph} = gensph2(rsph(nsph),spacing);
    r_grid(:,nsph) = pos{nsph}(:,3);
    th_grid(:,nsph) = pos{nsph}(:,2);
    phi_grid(:,nsph) = pos{nsph}(:,1);
end

%% convert both grids into matlab theta, phi convention, get cartesian coordinates

x_grid = zeros(size(r_grid));
y_grid = zeros(size(r_grid));
z_grid = zeros(size(r_grid));

th_temp = - th_grid + pi/2; 
r_temp = r_grid;
phi_temp = phi_grid;
for u = 1:size(phi_grid,2)
    tf = phi_grid(:,u) > pi & phi_grid(:,u) < 2*pi;
    phi_temp(tf,u) = - 2*pi + phi_grid(tf,u);

    [x1,y1,z1] = sph2cart(phi_temp(:,u),th_temp(:,u),r_temp(:,u));
    x_grid(:,u) = x1;
    y_grid(:,u) = y1;
    z_grid(:,u) = z1;
end

%% set dipole position (one for now)

dip_layer = 2; %ceil(N_sph/2);
dip_point = ceil(size(x_grid,1)/3);

r_dip = pos{dip_layer}(dip_point,3);
th_dip = pos{dip_layer}(dip_point,2);
phi_dip = pos{dip_layer}(dip_point,1);


x_dip = zeros(size(r_dip));
y_dip = zeros(size(r_dip));
z_dip = zeros(size(r_dip));

th_temp = - th_dip + pi/2; 
r_temp = r_dip;
phi_temp = phi_dip;
for u = 1:length(phi_dip)
    tf = phi_dip(:,u) > pi & phi_dip(:,u) < 2*pi;
    phi_temp(tf,u) = - 2*pi + phi_dip(tf,u);

    [x1,y1,z1] = sph2cart(phi_temp(:,u),th_temp(:,u),r_temp(:,u));
    x_dip(:,u) = x1;
    y_dip(:,u) = y1;
    z_dip(:,u) = z1;
end

% dipole strengths in theta, phi directions
dipmom_th = 1; 
dipmom_phi = 1;

%% calculate leadfield
% rs: center of the conducting sphere
% r0: location of the dipole in the coordinate system of the conductor
rs = [0,0,0];
r0 = {x_grid,y_grid,z_grid};
lf = dipole_field_sarvas_lf(rs',r0,R,EX,EY,EZ,mags);
[lf_svd,~,~]=svd(lf,'econ');

%% compare to SSS methods
% spheroidal harmonics
% [semi_major,semi_minor,origin]=find_ellipse_axis(R');
% [Sin_spm_p,Sout_spm_p] = spheroidIN_spheroidOUT(R',EZ',origin,semi_major,semi_minor,8,3);
% for j = 1:size(Sin_spm_p,2)
%   SNin_spm(:,j) = Sin_spm_p(:,j)/norm(Sin_spm_p(:,j));
% end
% for j = 1:size(Sout_spm_p,2)
%   SNout_spm(:,j) = Sout_spm_p(:,j)/norm(Sout_spm_p(:,j));
% end
% %calculate multi-vsh in and single-vsh out
% thresh = 0.005;
% center1= [-0.00350699, 0.01138051, 0.05947857] - [0,0,0.05]; 
% center2= [-0.00433911, 0.04081329, 0.05194245] - [0,0,0.05];
% [SNin_tot, SNout] = multi_sss(center1,center2,R,EX,EY,EZ,ch_types,8,3, thresh);
% %single SSS
% [Sin,SNin] = Sin_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,8);
% [Sout,SNout] = Sout_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,3);
% 
% %create SSS method matricies
% sVSH_sVSH=SNin; %[SNin SNout];
% mVSH_sVSH=SNin_tot; %[SNin_tot SNout];
% oid_oid=SNin_spm; %[SNin_spm SNout_spm];
% oid_sVSH=[SNin_spm SNout];
% for i=1:size(SNin,2)
%     angle_SSS(i) = subspace(sVSH_sVSH(:,i),lf_svd)*180/pi;
%     angle_mSSS(i) = subspace(mVSH_sVSH(:,i),lf_svd)*180/pi;
%     angle_oid(i) = subspace(oid_oid(:,i),lf_svd)*180/pi;
%     angle_oSSS(i) = subspace(oid_sVSH(:,i),lf_svd)*180/pi;
% end
% min_SSS = min(angle_SSS);
% max_SSS = max(angle_SSS);
% mean_SSS = mean(angle_SSS);
% 
% min_mSSS = min(angle_mSSS);
% max_mSSS = max(angle_mSSS);
% mean_mSSS = mean(angle_mSSS);
% 
% min_oid = min(angle_oid);
% max_oid = max(angle_oid);
% mean_oid = mean(angle_oid);
% 
% min_oSSS = min(angle_oSSS);
% max_oSSS = max(angle_oSSS);
% mean_oSSS = mean(angle_oSSS);