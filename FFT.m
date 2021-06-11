plot_ez=ez;
Fs_x = 1/dx;       % pixels per centimeter
Fs_y = 1/dx;

Dx = 1/Fs_x;     % centimeters per pixel
Dy = 1/Fs_y;

[ M, N ] = size(plot_ez);      % pixels
xx = Dx*(0:N-1)';               % centimeters
yy = Dy*(0:M-1)';
  
dFx = Fs_x/N;              % cycles per centimeter
dFy = Fs_y/M;
  
Fx = (-Fs_x/2:dFx:Fs_x/2-dFx)';     % cycles per centimeter
Fy = (-Fs_y/2:dFy:Fs_y/2-dFy)';
  
ezt=fft2(abs(plot_ez));
[XX,ZZ]=meshgrid(Fx,Fy);
fft_ezt=fftshift(abs(ezt));
s2=log(1+fft_ezt);
mesh(XX,ZZ,s2)
xlabel('Fx coordinate'); ylabel('Fy coordinate');