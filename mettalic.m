%Metallic Waveguide
%Solution of Maxwell equation for TE wave on XZ plane
%using 2D FDTD (XY)in a charge and current free region
%Source we can placed manually
%%
close all; clc; clear all;
tic
%%
%some constants
epsilon0=8.854187817e-12;  %free space permitivity
mu0=4*pi*10^-7;             %free space perimiability
c=1/sqrt(mu0*epsilon0); %velosity of light
%%
%dimension of region in micron
xdim=0.01;%total x length in m
ydim=0.05; %total y length in m

%%
%discritation of space and time
delta=0.0001; %space 
deltax=delta; %x step
deltay=delta; %y step
deltat=delta/(c*sqrt(2)); %COURANT STABILITY CRITERIA (time step)2.35ps
%%
imax=round(xdim/deltax);
jmax=round(ydim/deltay);
nmax=500;
time_total=nmax*deltat; %total time in sec
%initialization
Hz=zeros(imax,jmax,nmax);
Ex=zeros(imax,jmax,nmax);
Ey=zeros(imax,jmax,nmax);
%%
%material description
epsilon=epsilon0*ones(imax,jmax,nmax);
mu=mu0*ones(imax,jmax,nmax);
%%
%figure();
%Source of light
xsource=[imax/4:3*imax/4];
ysource=1;
%calculation
for n=1:nmax
   
    
        Ex(1:imax-1,1:jmax-1,n+1)=Ex(1:imax-1,1:jmax-1,n) ...
            +(deltat./(epsilon(1:imax-1,1:jmax-1,n).*delta)).*(Hz(1:imax-1,2:jmax,n)-Hz(1:imax-1,1:jmax-1,n));
        Ey(1:imax-1,1:jmax-1,n+1)=Ey(1:imax-1,1:jmax-1,n) ...
            -(deltat./(epsilon(1:imax-1,1:jmax-1,n).*delta)).*(Hz(2:imax,1:jmax-1,n)-Hz(1:imax-1,1:jmax-1,n));
       
    
     
        Hz(2:imax,2:jmax,n+1)=Hz(2:imax,2:jmax,n) ...
            +(deltat./(mu(2:imax,2:jmax,n).*delta)) ...
            .*(Ex(2:imax,2:jmax,n+1)-Ex(2:imax,1:jmax-1,n+1)-Ey(2:imax,2:jmax,n+1)+Ey(1:imax-1,2:jmax,n+1));
    
    
    %Boundary condition (PML)
    Hz(imax/4,:,:)=0; Hz(3*imax/4,:,:)=0; %along zaxis
    Hz(:,1,:)=0; Hz(:,jmax,:)=0; %along yaxis
    Hz(xsource,ysource,:)=10*sin(n/10);
    if mod(n,10)==0
    imagesc(Hz(:,:,n),[-8,8]);axis('xy'); colormap('jet'); colorbar;
    title(num2str(n))
    end
    pause(0.0001);
end
figure();
imagesc(Hz(:,:,100),[-8,8]);axis('xy'); colormap('jet'); colorbar;
title('n=300')
xlabel('y-Axis in mm')
ylabel('x-Axis in mm')
figure();
imagesc(Hz(:,:,300),[-8,8]);axis('xy'); colormap('jet'); colorbar;
title('n=300')
xlabel('y-Axis in mm')
ylabel('x-Axis in mm')
figure();
imagesc(Hz(:,:,500),[-8,8]);axis('xy'); colormap('jet'); colorbar;
title('n=500')
xlabel('y-Axis in mm')
ylabel('x-Axis in mm')
toc