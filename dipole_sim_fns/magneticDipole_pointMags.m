function magneticFlux = magneticDipole_pointMags(chanpos,coilori,dip_pos, dip_mom)
% Iman and Xan 2024
%Calculating the magnetic field from a dipole using Griffiths EM eq. 5.89
% Assumes point magnetometers as the channel type ONLY- for OPM
% dip_mom: magnetic moment 1x3
% dip_pos: position vector, from the center of the magnetic dipole to where the magnetic field is measured
% chanpos:  3xnchan channel position matrix 
% coilori: the orientation/measuring direction of the sensors
scale = 1.0e-7; 

%calculate B for each channel position
nchan=size(chanpos,2);
for i=(1:nchan)
    rvec= chanpos(:,i)-dip_pos;
    r_hat=norm(rvec);
    Bfield = (3*dot(dip_mom,rvec)/r_hat^5)*rvec;
    Bfield2= dip_mom/r_hat^3;
    B(:,i)= scale*(Bfield-Bfield2);
    %projuct B into sensing direction- no integral needed
    magneticFlux(i)=dot(B(:,i), coilori(:,i));
end

end