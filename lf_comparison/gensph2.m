function sphcor = gensph2(r,spacing)
%generates (spacing*(spacing-1)) uniform points on a sphere.

phi = [0:2*pi/spacing:2*pi-2*pi/spacing]';
th = [pi/spacing:pi/spacing:pi-pi/spacing]';

sphcor = zeros(length(th)*length(phi),3);
j = 1;
for k = 1:length(phi)
    for l = 1:length(th)
        sphcor(j,:) = [phi(k),th(l),r];
        j = j + 1;
    end
end

%[x,y,z] = sph2cart(sphcor(:,1),sphcor(:,2),sphcor(:,3));

end