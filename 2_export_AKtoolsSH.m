addpath("~/projects_EXT/AKtools/2_Tools/SphericalHarmonics");

nmax = 5;

filename = ["gaussian", num2str(nmax), ".txt"];

cart = load(filename);
x = cart(:,1);
y = cart(:,2);
z = cart(:,3);

az = mod(atan2(y, x) * 180 / pi, 360);
el =  atan2(sqrt(abs(x).^2 + abs(y).^2), z) * 180 / pi;

Yreal = zeros(length(az), (nmax+1)^2);
Ycomplex = zeros(length(az), (nmax+1)^2);

i1 = 0;
for n = 0:nmax
    for m = -n:n
        i1 += 1;
        for ind = 1:length(az)
            Yreal(ind, i1) = AKsh(n, m, az(ind), el(ind), mode='real');
            Ycomplex(ind, i1) = AKsh(n, m, az(ind), el(ind), mode='complex');
        end
    end    
end

save(["AK_SHgaussian" num2str(nmax) ".mat"], "-v7", "Yreal", "Ycomplex");
