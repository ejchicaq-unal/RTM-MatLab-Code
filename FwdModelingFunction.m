function [data, FwdSnapshot] = FwdModelingFunction(V,rw,nz,~,nx,dx,nt,dt,ixs,izs,M,am,N,surec)

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


%% Definition of imput variables
% rw(nt)            Source wavelet
% v(nz,nx)          velocity model
% nx                number of horizontal samples
% nz                number of depth samples
% nt                numer of time samples
% dx                horizontal distance per sample
% dz                depth distance per sample
% dt                time difference per sample
% ixs               x shot position
% izs               z shot position
% M                 order of space derivative
% am                SDF coefficients
% N                 grids number in transition area

% Initialization of data storage arrays
[nzz,nxx] = size(V);
% this matrix stores the values of the wave field in surface
data = zeros(nxx,nt);
% this matrix stores the values of the wave field
fdm  = zeros(nzz,nxx,3);
% this matrix stores the values of the wave field y zone II
fdmII  = zeros(nzz,nxx,3);

%% Wave modeling in time direction

% Inclusion of the source function 'rw (t)' in the wavefield
fdm(izs+2*M+N,ixs+2*M+N,2) = rw(1);

for g = 1:nx
    % The array 'data' stores the values to form the surface records
    data(g+2*M+N,1) = fdm(2*M+N+surec(g),2*M+N+g,2);
end

FwdSnapshot = zeros(nzz,nxx,nt);   % 'FwdSnapshot' matrix stores the values of
                                   % the wave field para for each time increment

%% ========== matrix indices ==========

% Zone I & II for the calculation of Finite Difference with higher order
% wave equation 2M

iz = 2*M+1:nzz - (2*M);
ix = 2*M+1:nxx - (2*M);

%% Finite Difference Wave Field Calculation with Higher Order Wave Equation 2M
 % Taken from: Yan, Liu y Zhang. "Prestack reverse-time migration with a time space
 % domain adaptive high-order staggered-grid finite difference method"
 % Exploration Geophysics, 2013, 44.
 % DOI: https://doi.org/10.1071/EG12047

 % coefficient of the wave equation
load('velocityModel')
aa = (V*dt/dx).^2;
rr=V*dt/dx;

for it = 2:nt % % time iteration cycle

 % Finite differences in zones I and II
 % Wave field calculation by Finite Difference
 % with second order wave equation in one direction
 % (Clayton and Engquist, 1977)
 % Taken from: Liu, Ding y Sen,
 % "Comparisions between the hybrid ABC and the PML method"
 % SEG Annual Meeting, 2011 | https://doi.org/10.1190/1.3627807
 % Equation (4), pag 78 Emiro Chica RTM thesis

    fdm_a = zeros(nzz,nxx);
    fdm_b = zeros(nzz,nxx);
    fdm_c = zeros(nzz,nxx);
    p1 = zeros(nzz,nxx);
    p2 = zeros(nzz,nxx);

    fdm_a(iz,ix)=2*fdm(iz,ix,2)-fdm(iz,ix,1);

    for i=1:M % indice en la direccion 'x'
        for j=1:M % indice en la direccion 'z'
            p1(iz,ix)=(fdm(iz+i+j-1,ix,2)-fdm(iz+i-j,ix,2))-(fdm(iz-i+j,ix,2)-fdm(iz-i-j+1,ix,2));
            fdm_b(iz,ix)=fdm_b(iz,ix)+am(i)*am(j)*p1(iz,ix);
        end
    end

    for i=1:M % indice en la direccion 'x'
        for j=1:M % indice en la direccion 'z'
            p2(iz,ix)=(fdm(iz,ix+i+j-1,2)-fdm(iz,ix+i-j,2))-(fdm(iz,ix-i+j,2)-fdm(iz,ix-i-j+1,2));
            fdm_c(iz,ix)=fdm_c(iz,ix)+am(i)*am(j)*p2(iz,ix);
        end
    end

    fdm(iz,ix,3) = fdm_a(iz,ix)+aa(iz,ix).*fdm_b(iz,ix)+aa(iz,ix).*fdm_c(iz,ix);

    % Finite differences in zones II and III
    % Wave field calculation by Finite Difference
    % with second order wave equation in one direction
    % (Clayton and Engquist, 1977)
    % Taken from: Liu, Ding y Sen,
    % "Comparisions between the hybrid ABC and the PML method"
    % SEG Annual Meeting, 2011 | https://doi.org/10.1190/1.3627807
    % Equation (2) and (4), pag 78 Emiro Chica RTM thesis

    fdmII(2*M+N+1:2*M+N+nz,2*M+N+1:2*M+N+nx,3) = fdm(2*M+N+1:2*M+N+nz,2*M+N+1:2*M+N+nx,3);

    % Zone II

    % left strip

    for i = 2*M+N:-1:2*M+2 % x

        for j = 2*M+N+1:2*M+N+nz  % z

        fdmII(j,i,3) = 1/(1 + rr(j,i))*(rr(j,i)*(fdmII(j,i+1,3) - (fdmII(j,i+1,1) - fdmII(j,i,1))) - ...
        ((fdmII(j,i+1,1)+fdmII(j,i+1,3)-2*fdmII(j,i+1,2))+(fdmII(j,i,1)-2*fdmII(j,i,2))) + ...
        (0.5*rr(j,i)^2)*((fdmII(j+1,i+1,2)+fdmII(j-1,i+1,2)-2*fdmII(j,i+1,2)) + ...
        (fdmII(j+1,i,2)+fdmII(j-1,i,2)-2*fdmII(j,i,2))));

        end

    end

    % right strip
    for i = 2*M+N+nx+1:2*M+N+nx+N-1 % x

        for j = 2*M+N+1:2*M+N+nz % z

        fdmII(j,i,3) = 1/(1 + rr(j,i))*(rr(j,i)*(fdmII(j,i-1,3) - (fdmII(j,i-1,1) - fdmII(j,i,1))) - ...
        ((fdmII(j,i-1,1)+fdmII(j,i-1,3)-2*fdmII(j,i-1,2))+(fdmII(j,i,1)-2*fdmII(j,i,2))) + ...
        (0.5*rr(j,i)^2)*((fdmII(j+1,i-1,2)+fdmII(j-1,i-1,2)-2*fdmII(j,i-1,2)) + ...
        (fdmII(j+1,i,2)+fdmII(j-1,i,2)-2*fdmII(j,i,2))));

        end

    end

    % top strip

    for j = 2*M+N:-1:2*M+2 % z |for j = nzz - (2*M+N-1):nzz-(2*M+1) % z | for j = nzz - (2*M+N-1):nzz-(2*M+1)

         for i = 2*M+1:2*M+N+nx+N %2*M+N+1-f:2*M+N+nx+f % x | for i = 2*M+2:nxx - (2*M+1) % x 1 for i = 2*M+N+1:nxx - (2*M+N)

        fdmII(j,i,3) = 1/(1 + rr(j,i))*(rr(j,i)*(fdmII(j+1,i,3) - (fdmII(j+1,i,1) - fdmII(j,i,1))) - ...
        ((fdmII(j+1,i,1)+fdmII(j+1,i,3)-2*fdmII(j+1,i,2))+(fdmII(j,i,1)-2*fdmII(j,i,2))) + ...
        (0.5*rr(j,i)^2)*((fdmII(j+1,i-1,2)+fdmII(j+1,i+1,2)-2*fdmII(j+1,i,2)) + ...
        (fdmII(j,i-1,2)+fdmII(j,i+1,2)-2*fdmII(j,i,2))));

         end

    end

    % bottom strip

    for j = 2*M+N+nz+1:2*M+N+nz+N-1 % z |for j = nzz - (2*M+N-1):nzz-(2*M+1) % z | for j = nzz - (2*M+N-1):nzz-(2*M+1)

         for i = 2*M+1:2*M+N+nx+N %2*M+N+1-f:2*M+N+nx+f % x | for i = 2*M+2:nxx - (2*M+1) % x 1 for i = 2*M+N+1:nxx - (2*M+N)

        fdmII(j,i,3) = 1/(1 + rr(j,i))*(rr(j,i)*(fdmII(j-1,i,3) - (fdmII(j-1,i,1) - fdmII(j,i,1))) - ...
        ((fdmII(j-1,i,1)+fdmII(j-1,i,3)-2*fdmII(j-1,i,2))+(fdmII(j,i,1)-2*fdmII(j,i,2))) + ...
        (0.5*rr(j,i)^2)*((fdmII(j-1,i+1,2)+fdmII(j-1,i-1,2)-2*fdmII(j-1,i,2)) + ...
        (fdmII(j,i+1,2)+fdmII(j,i-1,2)-2*fdmII(j,i,2))));

         end

    end

    % Lower left corner

    for i = 2*M+N:-1:2*M+2 % x

        for j = nzz-(2*M+N-1):(nzz -(2*M+1)); % z

        fdmII(j,i,3) = 1/(1 + rr(j,i))*(rr(j,i)*((fdmII(j-1,i+1,3) + fdmII(j-1,i-1,2)) - fdmII(j,i,2)) + ...
        (fdmII(j-1,i+1,2)+fdmII(j,i,2)-fdmII(j-1,i+1,3)));


        end
    end

        % Upper left corner

    for i = 2*M+N:-1:2*M+2 % x

        for j = 2*M+N:-1:2*M+2; % z

        fdmII(j,i,3) = 1/(1 + rr(j,i))*(rr(j,i)*((fdmII(j+1,i+1,3) + fdmII(j+1,i-1,2)) - fdmII(j,i,2)) + ...
        (fdmII(j+1,i+1,2)+fdmII(j,i,2)-fdmII(j+1,i+1,3)));
        end
    end

    % Lower right corner

     for i = 2*M+N+nx+1:2*M+N+nx+N-1 % x

        for j = nzz-(2*M+N-1):(nzz -(2*M+1)); % z

    fdmII(j,i,3) = 1/(1 + rr(j,i))*(rr(j,i)*((fdmII(j-1,i-1,3) + fdmII(j-1,i+1,2)) - fdmII(j,i,2)) + ...
        (fdmII(j-1,i-1,2)+fdmII(j,i,2)-fdmII(j-1,i-1,3)));

        end
     end

     % Rigth upper corner

     for i = nxx - (2*M+N-1):nxx - (2*M+1) % x

        for j = 2*M+N:-1:2*M+2; % z

    fdmII(j,i,3) = 1/(1 + rr(j,i))*(rr(j,i)*((fdmII(j+1,i-1,3) + fdmII(j+1,i+1,2)) - fdmII(j,i,2)) + ...
        (fdmII(j+1,i-1,2)+fdmII(j,i,2)-fdmII(j+1,i-1,3)));

        end
     end

    % Zone III

    % Left contour

    i = 2*M+1;% x

        for j = 2*M+1:2*M+N+nz+N; % z

        fdmII(j,i,3) = 1/(1 + rr(j,i))*(rr(j,i)*(fdmII(j,i+1,3) + fdmII(j,i+1,2) - fdmII(j,i,2)) + ...
        (fdmII(j,i+1,2)+fdmII(j,i,2)-fdmII(j,i+1,3)));

        end

   % Right contour

    i = 2*M+N+nx+N; % x

        for j = 2*M+1:2*M+N+nz+N; % z

        fdmII(j,i,3) = 1/(1 + rr(j,i))*(rr(j,i)*(fdmII(j,i-1,3) + fdmII(j,i-1,2) - fdmII(j,i,2)) + ...
        (fdmII(j,i-1,2)+fdmII(j,i,2)-fdmII(j,i-1,3)));

        end

   % Bottom contour

    j = 2*M+N+nz+N; % z

        for i = 2*M+2:2*M+N+nx+N-1 % x

     fdmII(j,i,3) = 1/(1 + rr(j,i))*(rr(j,i)*(fdmII(j-1,i,3) + fdmII(j-1,i,2) - fdmII(j,i,2)) + ...
        (fdmII(j-1,i,2)+fdmII(j,i,2)-fdmII(j-1,i,3)));

        end

    % Top contour

    j = 2*M+1; % z

        for i = 2*M+2:2*M+N+nx+N-1 % x

     fdmII(j,i,3) = 1/(1 + rr(j,i))*(rr(j,i)*(fdmII(j+1,i,3) + fdmII(j+1,i,2) - fdmII(j,i,2)) + ...
        (fdmII(j+1,i,2)+fdmII(j,i,2)-fdmII(j+1,i,3)));

        end

    w = zeros(N,1); % w(i) son pesos que varian de 0 a 1

    w(N) = 0;

    for i=N-1:-1:1
        w(i) = w(i+1)+1/(N-1);
    end

    % Left strip

    for i=N:-1:1
        for j=2*M+1:2*M+N+nz+N
            fdm(j,2*M+i,3)=(1-w(i))*fdm(j,2*M+i,3)+w(i)*fdmII(j,2*M+i,3);
        end
    end

    % Right strip

    for i=N:-1:1
        for j=2*M+1:2*M+N+nz+N
            fdm(j,2*M+N+nx+i,3)=(1-w(i))*fdm(j,2*M+N+nx+i,3)+w(i)*fdmII(j,2*M+N+nx+i,3);
        end
    end

    % Bottom strip

    for i=N:-1:1
        for j=2*M+N+1:2*M+N+nx
            fdm(2*M+N+nz+i,j,3)=(1-w(i))*fdm(2*M+N+nz+i,j,3)+w(i)*fdmII(2*M+N+nz+i,j,3);
        end
    end

     % Top strip

    for i=N:-1:1
        for j=2*M+N+1:2*M+N+nx
            fdm(2*M+i,j,3)=(1-w(i))*fdm(2*M+i,j,3)+w(i)*fdmII(2*M+i,j,3);
        end
    end

    % update fdm (wave field) for the next iteration in time

    fdm(:,:,1) = fdm(:,:,2);
    fdm(:,:,2) = fdm(:,:,3);
    fdmII(:,:,1) = fdmII(:,:,2);
    fdmII(:,:,2) = fdmII(:,:,3);
    update the source function 'rw (it)'
    fdm(izs+2*M+N,ixs+2*M+N,2) = fdm(izs+2*M+N,ixs+2*M+N,2)+rw(it);
    fdm(izs+2*M+N,ixs+2*M+N+1,2) = fdm(izs+2*M+N,ixs+2*M+N+1,2)+rw(it);

    % update data
for g = 1:nx
    data(g+2*M+N,it) = fdm(2*M+N+surec(g),2*M+N+g,2);
end
    FwdSnapshot(:,:,it) = fdm(:,:,2);

end % ends time loop
