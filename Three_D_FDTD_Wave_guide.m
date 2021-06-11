
%3D FDTD for the air filled rectangular cavity.

%Dimension of the wave guide is given by
%x_direction=10 cm.
%y_direction=5 cm.
%z_direction=30 cm. (wave progation takes places in z direction).

% PEC boundary conditions:     
% ex(i,j,k)=0 on the j=1, j=jb, k=1, and k=kb planes
% ey(i,j,k)=0 on the i=1, i=ib, k=1, and k=kb planes
% ez(i,j,k)=0 on the i=1, i=ib, j=1, and j=jb planes
% These are the PEC boundaries for the cavity.


%clear all; %Clears the workspace.
clc; %Clears the command Window.
tic
%This are the fundamental constants.
cc=2.99792458e8;            %speed of light in free space
muz=4.0*pi*1.0e-7;          %permeability of free space
epsz=1.0/(cc*cc*muz);       %permittivity of free space

% waveguide dimension parameter.
Lx=0.1; %length of the waveguide in x direction (Metres).
Ly=0.05; %length of the waveguide in y direction (Metres).
Lz=0.3; %length of the waveguide in z direction (Metres).

dx=0.002;   %space increment of waveguide(m).
dt=dx/(2.0*cc); %time increment of waveguide(sec).

ie=ceil(Lx/dx);  %number of grid cells  along x-direction
je=ceil(Ly/dx);  %number of grid cells  along y-direction
ke=ceil(Lz/dx);  %number of grid cells  along z-direction

ib=ie+1;    %grid boundary coordinates
jb=je+1;    %grid boundary coordinates
kb=ke+1;    %grid boundary coordinates

is=ceil(ie/2);  %location of source.
ks=2;       

jobs=ceil(je/2); %location of observation.

nmax=3000;        %Number of time steps


%parameters of the differentiated Gaussian pluse.
rtau=50.0e-12;
tau=rtau/dt;
ndelay=3*tau;
J0=-1.0;

%Here we have considered air filled cavity with sigma=0.0.
eps=1.0;
ca=1.0;
cb=(dt/epsz/eps/dx);
da=1.0;
db=dt/muz/dx;

%intialising the matrix of electric and magnetic components.
ex=zeros(ie,jb,kb);
ey=zeros(ib,je,kb);
ez=zeros(ib,jb,ke);
hx=zeros(ib,je,ke);
hy=zeros(ie,jb,ke);
hz=zeros(ie,je,kb);

%runing the loop for the nmax time steps.
for n=1:nmax

%updating electric fields

ex(1:ie,2:je,2:ke)=ca*ex(1:ie,2:je,2:ke)+...
                   cb*(hz(1:ie,2:je,2:ke)-hz(1:ie,1:je-1,2:ke)+...
                       hy(1:ie,2:je,1:ke-1)-hy(1:ie,2:je,2:ke));

ey(2:ie,1:je,2:ke)=ca*ey(2:ie,1:je,2:ke)+...
                   cb*(hx(2:ie,1:je,2:ke)-hx(2:ie,1:je,1:ke-1)+...
                       hz(1:ie-1,1:je,2:ke)-hz(2:ie,1:je,2:ke));
                    
ez(2:ie,2:je,1:ke)=ca*ez(2:ie,2:je,1:ke)+...
                   cb*(hx(2:ie,1:je-1,1:ke)-hx(2:ie,2:je,1:ke)+...
                       hy(2:ie,2:je,1:ke)-hy(1:ie-1,2:je,1:ke));

%input pulse.                   
ez(is,1:je,ks)=ez(is,1:je,ks)+...
               10*J0*(n-ndelay)*exp(-((n-ndelay)^2/tau^2));

%updating magnetic fields.
hx(2:ie,1:je,1:ke)=hx(2:ie,1:je,1:ke)+...
                   db*(ey(2:ie,1:je,2:kb)-ey(2:ie,1:je,1:ke)+...
                       ez(2:ie,1:je,1:ke)-ez(2:ie,2:jb,1:ke));
                
hy(1:ie,2:je,1:ke)=hy(1:ie,2:je,1:ke)+...
                   db*(ex(1:ie,2:je,1:ke)-ex(1:ie,2:je,2:kb)+...
                       ez(2:ib,2:je,1:ke)-ez(1:ie,2:je,1:ke));
                
hz(1:ie,1:je,2:ke)=hz(1:ie,1:je,2:ke)+...
                   db*(ex(1:ie,2:jb,2:ke)-ex(1:ie,1:je,2:ke)+...
                       ey(1:ie,1:je,2:ke)-ey(2:ib,1:je,2:ke));

hz(is,1:je,ks)=hz(is,1:je,ks)+...
             10*J0*(n-ndelay)*exp(-((n-ndelay)^2/tau^2));

%plots
 timestep=int2str(n);
 tview=squeeze(ez(:,jobs,:));
 sview=squeeze(ez(is,:,:));
 
 figure(1)
 imagesc(tview')
 shading flat;
 caxis([-1.0 1.0]); 
 colorbar;
 axis image; axis xy; 
 title(['Ez(i=:,j=13,k=:), time step = ',timestep]);
 xlabel('x coordinate'); ylabel('z coordinate');
 
 figure(2)
 imagesc(sview');
 shading flat;
 caxis([-1.0 1.0]); 
 colorbar;
 axis image; axis xy;
 title(['Ez(i=25,j=:,k=:), time step = ',timestep]);
 xlabel('y coordinate'); ylabel('z coordinate');
 
pause(0.0000005)

end

toc
