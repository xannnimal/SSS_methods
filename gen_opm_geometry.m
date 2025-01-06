function [opm_matrix,R_hat,theta_hat,phi_hat,ch_types] = gen_opm_geometry(filename)
%% given a matrix (x,y,z) of OPM sensors, find r-hat, phi-hat, and theta-hat
% Xan McPherson, 2023
%INPUT filename of SANDIA helmet geometry, ".mat" file
%OUTPUT r_hat,theta_hat,phi_hat three nchanx3 matricies normal coil ori

geom = load(filename);
opm_geometry = geom.Headwithsensors1ScaledHeadPoints; %table
xsens=opm_geometry.x; 
ysens=opm_geometry.y;
zsens=opm_geometry.z;
opm_matrix=[xsens,ysens,zsens]; %matrix

% find r-hat vector
for i=(1:size(xsens,1))
    RX(i)=xsens(i)/sqrt(xsens(i)^2+ ysens(i)^2+ zsens(i)^2);
    RY(i)=ysens(i)/sqrt(xsens(i)^2+ ysens(i)^2+ zsens(i)^2);
    RZ(i)=zsens(i)/sqrt(xsens(i)^2+ ysens(i)^2+ zsens(i)^2);
end
R_hat=[RX',RY',RZ']; 

%find theta and phi hat
for i=(1:size(xsens,1))
    theta(i)=acos(zsens(i)/sqrt(xsens(i)^2+ ysens(i)^2+ zsens(i)^2));
    phi(i)=sign(ysens(i))*acos(xsens(i)/sqrt(xsens(i)^2+ ysens(i)^2));
end
phi_z=zeros(1,size(xsens,1));
for i=(1:size(xsens,1))
    theta_x(i)=cos(theta(i))*cos(phi(i));
    theta_y(i)=cos(theta(i))*sin(phi(i));
    theta_z(i)=-sin(theta(i));

    phi_x(i)= -sin(phi(i));
    phi_y(i)=cos(phi(i));
end
theta_hat=[theta_x',theta_y',theta_z'];
phi_hat=[phi_x',phi_y',phi_z'];

ch_types = ones(size(xsens,1),1); % 1 = magnetometers, 0 = gradiometers for S_in code

% check angles and plot sensor coil orientations
% opm_matrix=[xsens,ysens,zsens]
% figure(1)
% hold on;
% grid on;
% quiver3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3),R_hat(:,1),R_hat(:,2),R_hat(:,3))
% quiver3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3),theta_hat(:,1),theta_hat(:,2),theta_hat(:,3))
% quiver3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3),phi_hat(:,1),phi_hat(:,2),phi_hat(:,3))
% scatter3(opm_matrix(:,1),opm_matrix(:,2),opm_matrix(:,3),'r')
% scatter3(0,0,0,'g*')
% %scatter3(center2(1),center2(2),center2(3),'b*')
% title('OPM Helmet')
% xlabel('x axis (m)')
% ylabel('y axis (m)')
% zlabel('z axis (m)')
% hold off;
% rotate3d
% view(135, 20);


end