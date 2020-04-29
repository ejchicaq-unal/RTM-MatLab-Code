%% Forward modeling high order SG FDM
%
% By Emiro Chica
% 20.02.2014
%

%%Copyright 2020 Emiro Chica Quiñones - echica@terralica.com
%
%Licensed under the Apache License, Version 2.0 (the "License");
%you may not use this file except in compliance with the License.
%You may obtain a copy of the License at
%
%  http://www.apache.org/licenses/LICENSE-2.0
%
%Unless required by applicable law or agreed to in writing, software
%distributed under the License is distributed on an "AS IS" BASIS,
%WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%See the License for the specific language governing permissions and
%limitations under the License.
%

%% Load velocity model

load 'velocityModel'

[nz,nx] = size(velocityModel);

%% Defines the grid size

dx = 10; % grid size in surface direction
dz = 10; % grid size in depth direction

x = (1:nx)*dx; % offset vector
z = (1:nz)*dz; % depth vectors

%% Depth of source point in meters

prof=10;

izs=prof/dz;

%% Create shot gathers
% Use the velocity model to simulate a seismic survey.  The wave equations
% is solved using finite differences for a defined initial wavefield.

% Time sampling interval in miliseconds

dt=0.001;

nt=21; %

t  = (0:nt-1).*dt;

M=3; %(2M)th order SDF

N=10; % Transition area width

% Add region around velocity model for applying absorbing
% boundary conditions (2*M nodes wide)

velocityModel_2=zeros(nz+2*M+N,nx);

for x1=1:nx
    for z1=1:nz+2*M+N
        velocityModel_2(z1,x1)=0; % velocityModel_2(z1,x1)=300;
    end
end

for x1=1:nx
    for z1=1:nz
        velocityModel_2(z1+2*M+N,x1)=velocityModel(z1,x1);
    end
end


V = [repmat(velocityModel_2(:,1),1,2*M+N) velocityModel_2 repmat(velocityModel_2(:,end),1,2*M+N)];
V(end+1:end+2*M+N,:) = repmat(V(end,:),2*M+N,1);

%% Definition of the source waveleth

%% uncomment this function to use a Ricker waveleth
% Ricker waveleth
f  = 60;
t0=1/f;
T = dt*(nt-1);
tt = 0:dt:T;
tau = tt-t0;
rw = (1-tau.*tau*f^2*pi^2).*exp(-tau.^2*pi^2*f^2);

%% uncomment this function to use a sinusoidal pulse as waveleth
% sinusoidal pulse
%f  = 20;    % frequency
%t0=1/f;     % period
%T = dt*(nt-1);
%rw=zeros(nt);
%for i=0:t0/dt
%    rw(i+1)=sin(2*pi*f*i*dt);
%end

%  Algorithm to solve the equations to find the coefficients of the time-space
%  domain for SDF based on the paper:
%  Prestack RTM with a time-space domain adaptive high-order staggered-grid
%  Finite Difference Method
%  Yan et all, Exploration Geophysics, 2013, 44, 77-86
%  http://dx.doi.org/10.1071/EG12047

vmin = min(velocityModel(:));
r=vmin*dt/dx;

% Plane wave propagation direction angle is posible to use teta=pi/8 too
teta=pi/16;

K=zeros(M,M);

for i=1:M
    for j=1:M
        K(j,i)=(2*i-1)^(2*j-2);
    end
end

xx=K(2,:);

bb=zeros(1,M);
bn=zeros(1,M);
beta=zeros(1,M);
c=zeros(1,M);
d=zeros(1,M);
a=zeros(1,M);

for n=1:M
    beta(n)=((-1)^(n-1))/(factorial(2*n-1));
    c(n)=((cos(teta))^(2*n-1))*beta(n);
    d(n)=(sin(teta)^(2*n-1))*beta(n);
end

bn(1)=1; % Ecuacion 7a

for n=2:M
    b1=0;
    for j=1:n
        b1=beta(j)*beta(n+1-j)+b1;
    end
    b2=0;
    for j=2:n-1
        b2=(bb(j)*bb(n+1-j))*(c(j)*c(n+1-j)+d(j)*d(n+1-j))+b2;
    end
    bn(n)=(b1*r^(2*n-2)-b2)/(2*(c(1)*c(n)+d(1)*d(n))); % Ecuacion 7b
end

% The matrix 'bb' is an almost singular matrix, therefore it is not possible to
% solve equation (6) by inverting the matrix.
% 'bb' is a Vandermonde matriz

for i=1:M
    a(i)=2*i-1;
end

bb=bn;

n=length(bb);

% Vandermonde matrix solution based on the paper:
% Solution of Vandermonde Systems of Equations
% Bjorck et al
% Mathematics of Computation (1970)
% Volume 24, number 112

% Appendix
for k=1:n-1
    for i=n:-1:k+1
        bb(i)=bb(i)-xx(k)*bb(i-1);
    end
end
for k=n-1:-1:1
    for i=k+1:n
        bb(i)=bb(i)/(xx(i)-xx(i-k));
    end
    for i=k:n-1
        bb(i)=bb(i)-bb(i+1);
    end
end

am=bb./a;

%% Wavefield modeling

data = zeros(size(nt,nx));
figure(gcf)

shot_1 = 339;    % Number of the first record - this can be changed
shot_2 = 339;    % Number of the last record - this can be changed
is = 5;          % Source interval
I =  zeros(nz,nx);


for ixs =  shot_1:is:shot_2 % Records loop

    % Wavefield modeling in the direction of time (Forward Modeling)

    texto = ['Record No.: ',num2str(ixs)];
    disp(' ')
    disp('===============================================================')
    disp(texto)

    disp('Calculating wavefield modeling in the direction of time (Forward Modeling)')

    tic

    [data, FwdSnapshot] = FwdModelingFunction(V,rw,nz,dz,nx,dx,nt,dt,ixs,izs,M,am,N);

    FwdSnapshot = FwdSnapshot(2*M+N+1:end-(2*M+N),2*M+N+1:end-(2*M+N),:);

    % save(['/out/Fwd/FwdSnapshot',num2str(ixs-(2*M+N)),'.mat'],'FwdSnapshot','-v7.3');
    % save(['/out/Fwd/FwdSnapshot',num2str(ixs),'.mat'],'FwdSnapshot','-v7.3');

    data = data(2*M+N+1:end-(2*M+N),:);

    %save(['/out/shot/shot',num2str(ixs-(2*M+N)),'.mat'],'data')
    %save(['/out/shot/shot',num2str(ixs),'.mat'],'data','-v7.3')

    toc

    % Modelamiento del campo de ondas en la direccion contraria del tiempo
    % (Backward % Modeling)

    disp('Calculating modeling in the opposite direction of time (Backward Modeling')

    tic

    [~, BckwdSnapshot] = BackModelingFunction(V,data,nz,nx,dx,nt,dt,M,am,N);

    BckwdSnapshot = BckwdSnapshot(2*M+N+1:end-(2*M+N),2*M+N+1:end-(2*M+N),:);

    %save(['/out/Bkwd/BckwdSnapshot',num2str(ixs),'.mat'],'BckwdSnapshot','-v7.3');
    %save(['/out/Bkwd/BckwdSnapshot',num2str(ixs),'.mat'],'BckwdSnapshot','-v7.3');

    toc

%% Image condition with decomposition of fields in the z direction
    Imz = zeros(nz,nx);
    Imx = zeros(nz,nx);
    SFTz = zeros(nz,1); %Fourier transform of the source
    RFTz = zeros(nz,1); %Fourier transform of receivers
    SFTx = zeros(nx,1); %Fourier transform of the source
    RFTx = zeros(nx,1); %Fourier transform of receivers

    disp('Calculating image condition in z direction')

    tic

    for xx = 1:nx
        Iz = 0;
        for tt = 1:nt

            for zz = 1:nz
                SFTz(zz) = FwdSnapshot(zz,xx,tt);
                RFTz(zz) = BckwdSnapshot(zz,xx,tt);
            end

            Sz = fft(SFTz);
            Rz = fft(RFTz);

            for zz = 1:round((nz+1)/2)

                Sz(zz) = 0;
                Rz(zz) = 0;
            end

            sz = ifft(Sz);
            rz= ifft(Rz);

            Iz = sz.*rz+Iz;

        end

        for zz = 1:nz

        Imz(zz,xx) = 2*real(Iz(zz));

        end

    %save(['/out/ImaCond/Imz',num2str(ixs),'.mat'],'Imz');
    %save(['/out/ImaCond/Imz',num2str(ixs),'.mat'],'Imz');
    %save(['/out/ImaCond/Imz',num2str(ixs),'.mat'],'Imz');

    end

    toc

%% Image condition with decomposition of fields in the x direction

    disp('Calculating image condition in z direction')

    tic

    for zz = 1:nz
        Ix = 0;
        for tt = 1:nt

            for xx = 1:nx
                SFTx(xx) = FwdSnapshot(zz,xx,tt);
                RFTx(xx) = BckwdSnapshot(zz,xx,tt);
            end

            Sx = fft(SFTx);
            Rx = fft(RFTx);

            for xx = 1:round((nx+1)/2)

                Sx(xx) = 0;
                Rx(xx) = 0;
            end

            sx = ifft(Sx);
            rx = ifft(Rx);

            Ix = sx.*rx+Ix;

        end

        for xx = 1:nx

        Imx(zz,xx) = 2*real(Ix(xx));

        end
    %save(['/out/ImaCond/Imx',num2str(ixs),'.mat'],'Imx');
    %save(['/out/ImaCond/Imx',num2str(ixs),'.mat'],'Imx');
    %save(['/out/ImaCond/Imx',num2str(ixs),'.mat'],'Imx');

    end

    toc

    disp('===============================================================')

    I = Imx+Imz+I;

    subplot(2,1,1)
    imagesc(x,z,Imx(:,:)+Imz(:,:))
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title(['Current record migrated ',num2str(ixs)]);
    colormap(gray)
    caxis ([-5 5])

    subplot(2,1,2)
    imagesc(x,z,I(:,:))
    xlabel('Distance (m)'); ylabel('Depth (m)');
    title(['Accumulated migrated image until registration ',num2str(ixs)]);
    colormap(gray)
    caxis ([-350 350])

    drawnow;

end

%save(('/out/ImaCond/I.mat'),'I');
%save(('/out/ImaCond/I.mat'),'I');
%save(('/out/ImaCond/I.mat'),'I');
