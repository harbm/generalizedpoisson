%test poisson3D
clear;
opts.method='\';
opts.tol=1E-6;

%% set up mesh
domain=[1 12 1 12 1 12];
n=20*[1 1 1];
h=[(domain(2)-domain(1))/(n(1)-1) (domain(4)-domain(3))/(n(2)-1) (domain(6)-domain(5))/(n(3)-1)];
mesh={(domain(1):h(1):domain(2)),(domain(3):h(2):domain(4)),(domain(5):h(3):domain(6))};
[mesh{1},mesh{2},mesh{3}]=ndgrid(mesh{1},mesh{2},mesh{3});

%% trial function
Bx=3*pi/(domain(2)-domain(1));By=3*pi/(domain(4)-domain(3));Bz=3*pi/(domain(6)-domain(5));
Vref=cos(Bx*(mesh{1}-domain(1))).*cos(By*(mesh{2}-domain(3))).*cos(Bz*(mesh{3}-domain(5)));

%% const epsilon
epsln=mesh{1}*0+1;
rho=-(Bx^2+By^2+Bz^2)*cos(Bx*(mesh{1}-domain(1))).*cos(By*(mesh{2}-domain(3))).*cos(Bz*(mesh{3}-domain(5)));

%% non-const epsilon
epsln=mesh{1}.*mesh{2}.*mesh{3};
rho=-Bx^2.*mesh{1}.*mesh{2}.*mesh{3}.*cos(Bx.*(-domain(1)+mesh{1})).*cos(By.*(-domain(3)+mesh{2})).*cos(Bz.*(-domain(5)+mesh{3})) - ...
    By^2.*mesh{1}.*mesh{2}.*mesh{3}.*cos(Bx.*(-domain(1)+mesh{1})).*cos(By.*(-domain(3) + mesh{2})).*cos(Bz.*(-domain(5) + mesh{3})) - ...
    Bz^2.*mesh{1}.*mesh{2}.*mesh{3}.*cos(Bx.*(-domain(1)+mesh{1})).*cos(By.*(-domain(3) + mesh{2})).*cos(Bz.*(-domain(5) + mesh{3})) - ...
    Bx.*mesh{2}.*mesh{3}.*cos(By.*(-domain(3)+mesh{2})).*cos(Bz.*(-domain(5)+mesh{3})).*sin(Bx.*(-domain(1) + mesh{1})) - ...
    By.*mesh{1}.*mesh{3}.*cos(Bx.*(-domain(1)+mesh{1})).*cos(Bz.*(-domain(5)+mesh{3})).*sin(By.*(-domain(3) + mesh{2})) - ...
    Bz.*mesh{1}.*mesh{2}.*cos(Bx.*(-domain(1)+mesh{1})).*cos(By.*(-domain(3)+mesh{2})).*sin(Bz.*(-domain(5) + mesh{3}));

%% impose BC on rhs
rho(1,:,:)=Vref(1,:,:);
rho(:,1,:)=Vref(:,1,:);
rho(:,:,1)=Vref(:,:,1);
rho(end,:,:)=Vref(end,:,:);
rho(:,end,:)=Vref(:,end,:);
rho(:,:,end)=Vref(:,:,end);

%% solve
V=poissonNH3D(rho,epsln,h,opts);
V([1 end],[1 end],:)=Vref([1 end],[1 end],:);
V([1 end],:,[1 end])=Vref([1 end],:,[1 end]);
V(:,[1 end],[1 end])=Vref(:,[1 end],[1 end]);
err=norm(V(:)-Vref(:),inf)
surf(log10(abs(V(:,:,2)-Vref(:,:,2))))