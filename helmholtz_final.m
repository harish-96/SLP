clear all
load data.mat
x = linspace(0, grid_diff*Nx, Nx);
y = linspace(0, grid_diff*Ny, Ny);
z = linspace(0, grid_diff*Nz, Nz);
velocity = zeros(3,192,192,192);

for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            velocity(:,i,j,k) = [u(i,j,k), v(i,j,k), w(i,j,k)];
        end
    end
end


freqs = fftfreq(velocity)*2*pi/grid_diff;
v_spec = fftn(velocity);
div_v_spec = zeros(192,192,192);
psi_spec = zeros(Nx,Ny,Nz);

for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            div_v_spec(i,j,k) = dot(v_spec(:,i,j,k), [freqs(i), freqs(j), freqs(k)]);
        end
    end
end
vx_spec = zeros(192,192,192);
vy_spec = zeros(192,192,192);
vz_spec = zeros(192,192,192);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            if freqs(i)^2 + freqs(j)^2+freqs(k)^2 ~=0
                psi_spec(i,j,k) = -div_v_spec(i,j,k) / (freqs(i)^2 +freqs(j)^2 + freqs(k)^2);
            else
                psi_spec(i,j,k) = -div_v_spec(i,j,k) / 10^-10;
            end
            temp = psi_spec(i,j,k) * [freqs(i) freqs(j) freqs(k)];
            vx_spec(i,j,k) = temp(1);
            vy_spec(i,j,k) = temp(2);
            vz_spec(i,j,k) = temp(3);
        end
    end
end

% psi = ifftn(psi_spec);
% [vx,vy,vz] = gradient(psi, grid_diff, grid_diff, grid_diff);
vx = ifftn(vx_spec);
vy = ifftn(vx_spec);
vz = ifftn(vx_spec);

u_comp = abs(vx);
v_comp = abs(vy);
w_comp = abs(vz);
u_incomp = u - u_comp;
v_incomp = v - v_comp;
w_incomp = w - w_comp;

[curlx,curly,curlz,cav] = curl(x,y,z,u_comp,v_comp,w_comp);
div_incomp = divergence(x,y,z,u_comp,v_comp,w_comp);

u_incomp_00 = zeros(192, 192);
u_comp_00 = zeros(192, 192);
for i=1:Ny
    for j=1:Nz
        u_incomp_00(i,j) = u_incomp(1,i,j);
        u_comp_00(i,j) = u_comp(1,i,j);
    end
end
%---------------------------------------------------------------%
figure;
surf(x, y, u_incomp_00)
title('x-velocity variation with y and z for x = 0')
xlabel('y')
ylabel('z')
zlabel('u-Incompressible')
saveas(gcf,'incomp_x_0','epsc')
%---------------------------------------------------------------%
figure;
surf(x, y, u_comp_00)
title('x-velocity variation with y and z for x = 0')
xlabel('y')
ylabel('z')
zlabel('u-Compressible')
saveas(gcf,'comp_x_0','epsc')
%---------------------------------------------------------------%
f0 = figure;
plot(x, reshape(u_incomp(1,1,:), 192, 1))
hold on
plot(x, reshape(u(1,1,:), 192, 1))
title('x-velocity variation with x for y = z = 0')
legend('Incompressible velocity', 'Total velocity')
xlabel('x')
ylabel('Velocity')
saveas(gcf,'incomp_x_:00','epsc')
%---------------------------------------------------------------%
f0 = figure;
plot(x, reshape(u_incomp(1,1,:), 192, 1))
hold on
plot(x, reshape(u(1,1,:), 192, 1))
title('x-velocity variation with x for y = z = 0')
legend('Incompressible velocity', 'Total velocity')
xlabel('x')
ylabel('Velocity')
saveas(gcf,'incomp_x_:00','epsc')
%---------------------------------------------------------------%
f = figure;
plot(y, reshape(u_comp(1,:,1), 192, 1))
hold on
plot(y, reshape(u(1,:,1), 192, 1))
title('x-velocity variation with y for x = z = 0')
legend('Compressible velocity', 'Total velocity')
xlabel('y')
ylabel('Velocity')
saveas(gcf,'comp_y_0:0','epsc')
%---------------------------------------------------------------%

f2 = figure;
plot(x, reshape(v_incomp(1,1,:), 192, 1))
hold on
plot(x, reshape(u(1,1,:), 192, 1))
title('y-velocity variation with x for y = z = 0')
legend('Incompressible velocity', 'Total velocity')
xlabel('x')
ylabel('Velocity')
saveas(gcf,'incomp_y_:00','epsc')
%---------------------------------------------------------------%
f = figure;
plot(y, reshape(u_incomp(1,:,1), 192, 1))
hold on
plot(y, reshape(u(1,:,1), 192, 1))
title('x-velocity variation with y for x = z = 0')
legend('Incompressible velocity', 'Total velocity')
xlabel('y')
ylabel('Velocity')
saveas(gcf,'incomp_y_0:0','epsc')
%---------------------------------------------------------------%
