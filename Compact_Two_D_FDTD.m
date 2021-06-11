
%2D compact for 3D FDTD for the air filled rectangular cavity.

%Dimension of the wave guide is given by
%x_direction=10 cm.
%y_direction=5 cm.
%z_direction=30 cm. (wave progation takes places in z direction).

% PEC boundary conditions:     
% ex(i,k)=0 on the k=1, and k=kb planes
% ey(i,k)=0 on the i=1, i=ib, k=1, and k=kb planes
% ez(i,k)=0 on the i=1, i=ib planes
% These are the PEC boundaries for the cavity.

clear all; %Clears the workspace.
clc; %Clears the command Window.
tic
%This are the fundamental constants.
cc=2.99792458e8;            %speed of light in free space
muz=4.0*pi*1.0e-7;          %permeability of free space
eps=1.0/(cc*cc*muz);       %permittivity of free space
epsz=1;

% waveguide dimension parameter.
Lx=0.1; %length of the waveguide in x direction (Metres).
Lz=0.3;  %length of the waveguide in z direction (Metres).

dx=0.002; %space increment of waveguide(m).
dt=dx/(2.0*cc); %time increment of waveguide(sec).

ie=ceil(Lx/dx);    %number of grid cells along x-direction
ke=ceil(Lz/dx);    %number of grid cells  along z-direction

is=ceil(ie/2);       %location of current source
ks=2;

ib=ie+1;    %grid boundary coordinates
kb=ke+1;

nmax=3000;          %number of time steps

%parameters of the differentiated Gaussian pluse.
 rtau=50.0e-12;
 tau=rtau/dt;
 ndelay=3*tau;
 J0=-1.0;

%N standing wave parameter 
l=1; %mode number l along y direction.
Ly=0.05; %length of the waveguide in z direction (Metres).

%coefficients and Here we have considered air filled cavity with sigma=0.
ca=(dt*l*pi)/(epsz*eps*Ly);
cb=dt/(epsz*eps*dx);
da=(dt*l*pi)/(muz*Ly);
db=dt/(muz*dx);

%intialising the matrix of electric and magnetic components.
ex=zeros(ie,kb);
ey=zeros(ib,kb);
ez=zeros(ib,ke);
hx=zeros(ib,ke);
hy=zeros(ie,ke);
hz=zeros(ie,kb);

%runing the loop for the nmax time steps.
for n=1:nmax
    
    %updating electric fields
    ex(1:ie,2:ke)=ex(1:ie,2:ke)-ca*(hz(1:ie,2:ke))+...
                       cb*(hy(1:ie,1:ke-1)-hy(1:ie,2:ke));
                   
    ey(2:ie,2:ke)=ey(2:ie,2:ke)+ cb*(hx(2:ie,2:ke)-hx(2:ie,1:ke-1)+...
                       hz(1:ie-1,2:ke)-hz(2:ie,2:ke));
                   
    ez(2:ie,1:ke)=ez(2:ie,1:ke)+ ca*(hx(2:ie,1:ke))+...
                       cb*(hy(2:ie,1:ke)-hy(1:ie-1,1:ke));
                   
    %input pulse.                
    ez(is,ks)=ez(is,ks)+...
               10*J0*(n-ndelay)*exp(-((n-ndelay)^2/tau^2));
   
   
          
    %updating Magnetic Field.
    
    hx(2:ie,1:ke)=hx(2:ie,1:ke)+db*(ey(2:ie,2:kb)-ey(2:ie,1:ke))-...
                       da*ez(2:ie,1:ke);
                   
    hy(1:ie,1:ke)=hy(1:ie,1:ke)+db*(ex(1:ie,1:ke)-ex(1:ie,2:kb))+...
                         db*(ez(2:ib,1:ke)-ez(1:ie,1:ke));
                   
    hz(1:ie,2:ke)=hz(1:ie,2:ke)+ da*ex(1:ie,2:ke)+...
                       db*(ey(1:ie,2:ke)-ey(2:ib,2:ke));

      %input pulse.                
    hz(is,ks)=hz(is,ks)+...
               10*J0*(n-ndelay)*exp(-((n-ndelay)^2/tau^2));
 %plots
    timestep=int2str(n);
    timestep
%   figure(1)
%   imagesc(ez);
%   shading flat;
%   caxis([-1.0 1.0]); 
%   colorbar;
%   axis image; axis xy; 
%   title(['ez, time step = ',timestep]);
%   xlabel('x coordinate'); ylabel('z coordinate');
  pause(0.0000005)
   
end

toc
