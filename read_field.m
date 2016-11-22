fid = fopen('iso192.prim','r','n') ;
prec = fread(fid,1,'integer*4') ;
version = fread(fid,1,'integer*4') ;
step  = fread(fid,1,'integer*4') ;
Nx = fread(fid,1,'integer*4') ;
Ny = fread(fid,1,'integer*4') ;
Nz = fread(fid,1,'integer*4') ;
nvar = fread(fid,1,'integer*4') ;
vars = fread(fid,nvar,'integer*4') ;
time = fread(fid,1,'double') ;
x = fread(fid,Nx,'double') ;
y = fread(fid,Ny,'double') ;
z = fread(fid,Nz,'double') ;
rho = fread(fid,Nx*Ny*Nz,'float') ;
u = fread(fid,Nx*Ny*Nz,'float') ;
v = fread(fid,Nx*Ny*Nz,'float') ;
w = fread(fid,Nx*Ny*Nz,'float') ;
p = fread(fid,Nx*Ny*Nz,'float') ;
Nx
diff(x)
%save filename.mat u v w Nx Ny Nz
