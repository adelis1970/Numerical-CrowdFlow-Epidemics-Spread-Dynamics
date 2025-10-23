% speed-density relationship, can be changed for a new simulation
function f = UV(rho)
[Ny,Nx] = size(rho);
um=1.4;
rm=6;

f= um.*exp(-7.5*(rho./rm).^2);

end
