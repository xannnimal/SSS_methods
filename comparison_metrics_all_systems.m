%% quantify differences between SSS methods using subspace angles
%QUID, Sandia OPM (two sensing directions), Kernel OPM
clear
%% constant variables 
coordsys = 'device';
Lin = 8; % Truncation order of the internal VSH basis
Lout = 3; % Truncation order of the external VSH basis
dim_in = (Lin+1)^2 - 1; % Dimension of the internal SSS basis, should be 80

%% Sensor systems
%306cnah SQUID system 
rawfile = 'sample_audvis_raw.fif';
[R,EX,EY,~] = fiff_getpos(rawfile, coordsys);
EZn = load('normal_vectors.mat');
EZ = EZn.EZp;
point_mags = 1; %0: account for grads and mags, %1: only point mags
j=1;
for i=(1:size(R,2))
    if mod(i,3)==0 %every third is a magnetometer
        ch_types_mags(:,j)=1;
        R_mags(:,j)=R(:,i);
        EX_mags(:,j)=EX(:,i);
        EY_mags(:,j)=EY(:,i);
        EZ_mags(:,j)=EZ(:,i);
        j=j+1;
    end
end

R=R_mags;
EX=EX_mags;
EY=EY_mags;
EZ=EZ_mags;
ch_types=ch_types_mags;

% % % opm geometry from Peter at SANDIA
% filename="headwithsensors1.mat";
% [R,R_hat,theta,phi,~] = gen_opm_geometry(filename);
% % "EZ" is the sensing direction, change to EY or EX for "phi"/"theta" hat
% R=R'; EZ=theta'; EY=phi'; EX=R_hat';
% point_mags = 1; %1: only point mags

% Kernel opm data: AKCLEE_110 updated April 2024 correct sensor positions
% filename= 'C:/Users/xanmc/RESEARCH/audio_ERF_notebook_portal/audio_ERF_portal_raw.fif';
% info = fiff_read_meas_info(filename);
% nchan=info.nchan;
% for i=1:nchan
%     R(:,i)=info.chs(i).loc(1:3,:);
%     EX(:,i)=info.chs(i).loc(4:6,:);
%     EY(:,i)=info.chs(i).loc(7:9,:);
%     EZ(:,i)=info.chs(i).loc(10:12,:);
% end
% point_mags = 1;

%% specify point mags or no
k=1;
if point_mags ==0 
    for i=(1:size(EZ,2))
        if mod(i,3)==0 %every third is a magnetometer
            ch_types(i)=1;
            mags(k)=i;
            k=k+1;
        else
            ch_types(i)=0;
            k=k;
        end
    end
else
    for i=(1:size(EZ,2))
        ch_types(i)=1;
        mags(i)=i;
    end
end

%% SSS methods
%VSH SSS
[Sin,SNin] = Sin_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,Lin);
[Sout,SNout] = Sout_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,Lout);

% spheroidal harmonics
[semi_major,semi_minor,origin]=find_ellipse_axis(R');
[Sin_spm_p,Sout_spm_p] = spheroidIN_spheroidOUT(R',EZ',origin,semi_major,semi_minor,Lin,Lout);
for j = 1:size(Sin_spm_p,2)
  SNin_spm(:,j) = Sin_spm_p(:,j)/norm(Sin_spm_p(:,j));
end
for j = 1:size(Sout_spm_p,2)
  SNout_spm(:,j) = Sout_spm_p(:,j)/norm(Sout_spm_p(:,j));
end

%calculate multi-vsh in and single-vsh out
thresh = 0.005;
center1= [-0.00350699, 0.01138051, 0.05947857] - [0,0,0.05]; 
center2= [-0.00433911, 0.04081329, 0.05194245] - [0,0,0.05];
[SNin_tot, SNout] = multi_sss(center1,center2,R,EX,EY,EZ,ch_types,Lin, Lout, thresh);


%% calculate subspace angles
sVSH_sVSH=[SNin SNout];
mVSH_sVSH=[SNin_tot SNout];
oid_oid=[SNin_spm SNout_spm];
oid_sVSH=[SNin_spm SNout];
%calculations independent of the type of raw data
for i=(1:80)
    %mVSH In and Spheroid In
    angles_mVSH_oid(i)=subspace(SNin_tot(:,i),SNin_spm)*180/pi;
    angles_sVSH_mVSH(i)=subspace(SNin_tot(:,i),SNin)*180/pi;
    angles_sVSH_oid(i)=subspace(SNin(:,i),SNin_spm)*180/pi;
%sVSH/sVSH and mVSH/sVSH 
    angles_sVSHsVSH_mVSHsVSH(i)=subspace(sVSH_sVSH(:,i),mVSH_sVSH)*180/pi;
%sVSH/sVSH and Spheroid/Spheroid 
    angles_sVSHsVSH_oidoid(i)=subspace(sVSH_sVSH(:,i),oid_oid)*180/pi;
%sVSH/sVSH and Spheroid/sVSH 
    angles_sVSHsVSH_oidsVSH(i)=subspace(sVSH_sVSH(:,i),oid_sVSH)*180/pi;
%mVSH/sVSH and Spheroid/Spheroid  
    angles_mVSHsVSH_oidoid(i)=subspace(mVSH_sVSH(:,i),oid_oid)*180/pi;
%mVSH/sVSH and Spheroid/sVSH
    angles_mVSHsVSH_oidsVSH(i)=subspace(mVSH_sVSH(:,i),oid_sVSH)*180/pi;
%Spheroid/Spheroid and Spheroid/sVSH
    angles_oidoid_oidsVSH(i)=subspace(oid_oid(:,i),oid_sVSH)*180/pi;
end

%find min/max/average
%mVSH In and Spheroid In
max_mVSH_oid_t=max(angles_mVSH_oid);
min_mVSH_oid_t=min(angles_mVSH_oid);
av_mVSH_oid_t=mean(angles_mVSH_oid);

max_sVSH_mVSH_t=max(angles_sVSH_mVSH);
min_sVSH_mVSH_t=min(angles_sVSH_mVSH);
av_sVSH_mVSH_t=mean(angles_sVSH_mVSH);

max_sVSH_oid_t=max(angles_sVSH_oid);
min_sVSH_oid_t=min(angles_sVSH_oid);
av_sVSH_oid_t=mean(angles_sVSH_oid);

%sVSH/sVSH and mVSH/sVSH 
max_sVSHsVSH_mVSHsVSH_t = max(angles_sVSHsVSH_mVSHsVSH);
min_sVSHsVSH_mVSHsVSH_t = min(angles_sVSHsVSH_mVSHsVSH);
av_sVSHsVSH_mVSHsVSH_t = mean(angles_sVSHsVSH_mVSHsVSH);
%sVSH/sVSH and Spheroid/Spheroid 
max_sVSHsVSH_oidoid_t = max(angles_sVSHsVSH_oidoid);
min_sVSHsVSH_oidoid_t = min(angles_sVSHsVSH_oidoid);
av_sVSHsVSH_oidoid_t = mean(angles_sVSHsVSH_oidoid);
%sVSH/sVSH and Spheroid/sVSH 
max_sVSHsVSH_oidsVSH_t = max(angles_sVSHsVSH_oidsVSH);
min_sVSHsVSH_oidsVSH_t = min(angles_sVSHsVSH_oidsVSH);
av_sVSHsVSH_oidsVSH_t = mean(angles_sVSHsVSH_oidsVSH);
%mVSH/sVSH and Spheroid/Spheroid  
max_mVSHsVSH_oidoid_t = max(angles_mVSHsVSH_oidoid);
min_mVSHsVSH_oidoid_t = min(angles_mVSHsVSH_oidoid);
av_mVSHsVSH_oidoid_t = mean(angles_mVSHsVSH_oidoid);

%mVSH/sVSH and Spheroid/sVSH
max_mVSHsVSH_oidsVSH_t = max(angles_mVSHsVSH_oidsVSH);
min_mVSHsVSH_oidsVSH_t = min(angles_mVSHsVSH_oidsVSH);
av_mVSHsVSH_oidsVSH_t = mean(angles_mVSHsVSH_oidsVSH);

%Spheroid/Spheroid and Spheroid/sVSH
max_oidoid_oidsVSH_t = max(angles_oidoid_oidsVSH);
min_oidoid_oidsVSH_t = min(angles_oidoid_oidsVSH);
av_oidoid_oidsVSH_t = mean(angles_oidoid_oidsVSH);