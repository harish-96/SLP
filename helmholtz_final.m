clear all
load data.mat
x = linspace(0, grid_diff*Nx, Nx);
y = linspace(0, grid_diff*Ny, Ny);
z = linspace(0, grid_diff*Nz, Nz);
v_div = divergence(x, y, z, u, v, w);
freqs = fftfreq(v_div)*2*pi/grid_diff;

div_v_spec = fftn(v_div);
psi_f = zeros(Nx,Ny,Nz);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            if freqs(i)^2 + freqs(j)^2+freqs(k)^2 ~=0
                psi_f(i,j,k) = -div_v_spec(i,j,k) / (freqs(i)^2 +freqs(j)^2 + freqs(k)^2);
            else
                psi_f(i,j,k) = -div_v_spec(i,j,k) / 10^-8;
            end
        end
    end
end

psi = ifftn(psi_f);
[vx,vy,vz] = gradient(psi, grid_diff, grid_diff, grid_diff);
u_incomp = real(vx);
v_incomp = real(vy);
w_incomp = real(vz);
u_comp = u - u_incomp;
v_comp = v - v_incomp;
w_comp = w - w_incomp;


