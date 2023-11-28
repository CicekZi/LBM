
clear all
clc

%% inputting from user
Uinitial = 0.1;	    % Uinitial in lattice Units
Vinitial = 0.0 ;	% Vinitial in lattice Units
Rhoinitial = 1.0;	% Rhoinitial in lattice Units


Tau = 0.65 ;		    % Tau - Relaxation time
R = 20;
nx = 600;
ny = 150;
Nu_physical= 1e-3;
channel_height=0.01;
Nu = (Tau-0.5)/3 ;
Re=Uinitial*2*ny/Nu;
Re_cylinder=Uinitial*2*R/Nu;
FD=R*Rhoinitial*Uinitial^2/105.6430/Re_cylinder;
t_lattice=channel_height^2/ny^2/3*(Tau-1/2)/Nu_physical;
dt=t_lattice;
dx=channel_height/ny;

%% initializing


%% lattice D2Q9
ex = [ 1  0 -1  0  1 -1 -1  1  0 ] ;
ey = [ 0  1  0 -1  1  1 -1 -1  0 ] ;
w  = [ 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36 4/9 ] ;
cs = 1/sqrt(3) ;	% Cs sound speed
matrix_a_b=[7 4 8;3 9 1;6 2 5]; % matrix used in force calculation in convergance for finding a =0,1,2,...9
%% 
Gx =  0;         % Body force in lattice units in x direction gravity-
Gy =  0;            % Body force in lattice units in y direction

% dx = 2e-4 ;         % dx=channel height/number oflattices if it was a real situation
%dt is calculated dt = 1e-5 ;         % dt

Utol = 1e-3 ;		% Toelarnce of U
Vtol = 1e-3 ;		% Toelarnce of V
Uresidual = 100 ;
Vresidual = 100 ;

%% define the parameters to be used
f = zeros (nx,ny,9);
ftemp = zeros (nx,ny,9);
feq = zeros (nx,ny,10);
fstar = zeros (nx,ny,9);
U = zeros (nx,ny);
V = zeros (nx,ny);
UW = zeros (nx,ny);
VW = zeros (nx,ny);
Ubf= zeros (nx,ny);
Vbf= zeros (nx,ny);
Rho = zeros (nx,ny);
delta= zeros (nx,ny);


%% Initialization
U(:,:) = 0 ;
V(:,:) = 0 ;
Rho(:,:) = Rhoinitial ;

Uold = U ;
Vold = V ;

ttt=200000000;

%% Cylinder nodes

isfluid = ones(nx,ny); % all nodes are part of fluid unless we note
isfluid2 = ones(nx,ny); % all nodes are part of fluid unless we note


% isfluid(:,1)  = 0; % south wall/ down wall
% isfluid(:,ny) = 0;% north wall / upper wall

ii=floor(nx/5+1);
jj=floor(105);

ii2=floor(nx/5 +1 );
jj2=floor(35);

for i=1:nx
    for j=1:ny
        if sqrt(  (i-ii)^2  +  (j-jj)^2  ) < R
            isfluid(i,j)=0;
            U(i,j)= 0.0;
            V(i,j)= 0.0;
        end
    end
end


for k=1:nx
    for m=1:ny
        if sqrt(  (k-ii2)^2  +  (m-jj2)^2  ) < R
            isfluid2(k,m)=0;
            U(k,m)= 0.0;
            V(k,m)= 0.0;
        end
    end
end

%% velocity_calibration

U(1,:)=Uinitial;
V(1,:)=0;
% V(nx,:)=0;
rho0_out=1;
rho0_in=1;
% U(:,1)= Uinitial;
% V(:,1)= 0;
% U(:,ny)= Uinitial;
% V(:,ny)= 0;

for i=1:nx
  for j=1:ny
    if isfluid(i,j)==0
      U(i,j)= 0.0;
      V(i,j)= 0.0;
    end
  end
end

for i=1:nx
  for j=1:ny
    if isfluid2(i,j)==0
      U(i,j)= 0.0;
      V(i,j)= 0.0;
    end
  end
end


%% feq calculation
for j=1:ny
  for i=1:nx
    for a=1:9
      feq(i,j,a)=Rho(i,j) * w(a) * ...
      (1 + 3 * ( ex(a)*U(i,j) + ey(a) * V(i,j) ) + ...
      9.0/2.0 * ( ex(a)*U(i,j) + ey(a) * V(i,j) )^2 -...
      3.0/2.0 *(U(i,j)^2+V(i,j)^2));
    end
  end
end


%% initializing f

f=feq;
loop_counter = 1;

while Uresidual > Utol || Vresidual > Vtol

    g = w .* 1 .* (Gx .* ex + Gy .* ey ) ./ cs^2;

   


  %% collision

for j=1:ny
  for i=1:nx
    if isfluid(i,j) == 1
      for a=1:9
        temp = f(i,j,a) - (f(i,j,a) - feq(i,j,a))/Tau + g(a) ;
        f(i,j,a) = temp;
      end
    end
  end
end

  %% C_D calculation

% Calculating fluid force on the cylinder; momentum exchange approach


Force_x=0;
Force_y=0;
Force_x_second=0;
Force_y_second=0;

for j=1:ny
  for i=1:nx
    if isfluid(i,j)==0
      ff_x=0;
      ff_y=0;
      ff_x_second=0;
      ff_y_second=0;
      for a=-1:1
        for b=-1:1
          if isfluid(i+a,j+b)==1 % for b and a =+1 -1 0 different cases a+2 and b+3 return same results the desired ex and ey
            ff_x=-ex(a+2)*( 2*f(i+a,j+b,matrix_a_b(b+2,a+2)) ) + ff_x; %%*dx/dt first method suggested by J. Gotz et al, Large scale simulation of fluid structure using Lattice Boltzmann Methods and the 'physics engine'
            ff_y=-ey(b+3)*( 2*f(i+a,j+b,matrix_a_b(b+2,a+2)) ) + ff_y; %%*dx/dt same as above
            ff_x_second=-ex(a+2)*(f(i,j,matrix_a_b(-b+2,-a+2))   + f(i+a,j+b,matrix_a_b(b+2,a+2)) ) + ff_x_second;%*dx/dt
            ff_y_second=-ey(b+3)*(f(i,j,matrix_a_b(-b+2,-a+2))   + f(i+a,j+b,matrix_a_b(b+2,a+2)) ) + ff_y_second;%*dx/dt
          end
        end
      end

      Force_x=Force_x+ ff_x;
      Force_y=Force_y+ ff_y;
      Force_x_second=Force_x_second+ ff_x_second;
      Force_y_second=Force_y_second+ ff_y_second;
    end
  end
end

C_d_lattice_second= abs(Force_x_second)/(Uinitial)^2/(R);%based on second momentum method
% C_d_keeper(loop_counter)=C_d_lattice_second; % when converging, average
% cd can be calculated from values obtained by cd_keeper

fprintf(' %i Force_x= %f, cd second lattice=%f, dx/dt=%f, Uresidual = %f\n',loop_counter,  Force_x, C_d_lattice_second,dx/dt, Uresidual)

%% streaming
for i=1:nx

  in=i-1;
  ip=i+1;

  for j=1:ny
    if j>1
      jn=j-1;
    else
      jn=ny;
    end
    if j<ny
      jp=j+1;
    else
      jp=1;
    end

    % ftemp id distribution at time = t+1
    % f is distribution at current time
    if jn~=0
      ftemp(i , jn  , 4) = f(i,j,4);
    end
    ftemp(i , j   , 9) = f(i,j,9);
    if jp~=ny+1
      ftemp(i , jp  , 2) = f(i,j,2);
    end
    if in~=0
      ftemp(in, j   , 3) = f(i,j,3);
      if jp~=ny+1
        ftemp(in, jp  , 6) = f(i,j,6);
      end
      if jn~=0
        ftemp(in, jn  , 7) = f(i,j,7);
      end
    end
    if ip~=nx+1
      if jp~=ny+1
        ftemp(ip, jp  , 5) = f(i,j,5);
      end
      if jn~=0
        ftemp(ip, jn  , 8) = f(i,j,8);
      end
      ftemp(ip, j   , 1) = f(i,j,1);
    end
  end
end

%% velocity_calibration

U(1,:)=Uinitial;
V(1,:)=0;
% V(nx,:)=0;
rho0_out=1;
rho0_in=1;
% U(:,1)= Uinitial;
% V(:,1)= 0;
% U(:,ny)= Uinitial;
% V(:,ny)= 0;

for i=1:nx
  for j=1:ny
    if isfluid(i,j)==0
      U(i,j)= 0.0;
      V(i,j)= 0.0;
    end
  end
end

for i=1:nx
  for j=1:ny
    if isfluid2(i,j)==0
      U(i,j)= 0.0;
      V(i,j)= 0.0;
    end
  end
end


%% Velocity and Pressure boundary, He & Zou

for j=1:ny

  rho0_in= (ftemp(1,j,9)+ftemp(1,j,2)+ftemp(1,j,4)+2.*(ftemp(1,j,3)+ftemp(1,j,6)+ftemp(1,j,7)))/(1-U(1,j));
  ru=rho0_in*U(1,j);
  rv=rho0_in*V(1,j);
  ftemp(1,j,1)= ftemp(1,j,3)+(2/3)*ru;
  ftemp(1,j,5)= ftemp(1,j,7)+  1/6*ru + 1/2*rv + 1/2*(ftemp(1,j,4)-ftemp(1,j,2));
  ftemp(1,j,8)= ftemp(1,j,6)+  1/6*ru - 1/2*rv + 1/2*(ftemp(1,j,2)-ftemp(1,j,4));



  
  U(nx,j)= -1+(ftemp(nx,j,9)+ftemp(nx,j,2)+ftemp(nx,j,4)+2.*(ftemp(nx,j,1)+ftemp(nx,j,5)+ftemp(nx,j,8)))/rho0_out;
  ru=rho0_out*U(nx,j);
  V(nx,j)=0;
  rv=0;
  ftemp(nx,j,3)=ftemp(nx,j,1) - 2/3*ru;
  ftemp(nx,j,6)=ftemp(nx,j,8) + 1/2*(ftemp(nx,j,4)-ftemp(nx,j,2)) - 1/6*ru;
  ftemp(nx,j,7)=ftemp(nx,j,5) + 1/2*(ftemp(nx,j,2)-ftemp(nx,j,4)) - 1/6*ru;
  

end

  f = ftemp;
%% bounce-back
  for j=1:ny
  for i=1:nx
    if isfluid(i,j) == 0  %%on down and upper boundaries

      temp = f(i,j,1);
      f(i,j,1) = f(i,j,3);
      f(i,j,3) = temp;

      temp = f(i,j,2);
      f(i,j,2) = f(i,j,4);
      f(i,j,4) = temp;

      temp = f(i,j,5);
      f(i,j,5) = f(i,j,7);
      f(i,j,7) = temp;

      temp = f(i,j,6);
      f(i,j,6) = f(i,j,8);
      f(i,j,8) = temp;
    end
  end
  end

for j=1:ny
  for i=1:nx
    if isfluid2(i,j) == 0  %%on down and upper boundaries

      temp = f(i,j,1);
      f(i,j,1) = f(i,j,3);
      f(i,j,3) = temp;

      temp = f(i,j,2);
      f(i,j,2) = f(i,j,4);
      f(i,j,4) = temp;

      temp = f(i,j,5);
      f(i,j,5) = f(i,j,7);
      f(i,j,7) = temp;

      temp = f(i,j,6);
      f(i,j,6) = f(i,j,8);
      f(i,j,8) = temp;
    end
  end
end

%% velocity_calibration

U(1,:)=Uinitial;
V(1,:)=0;
% V(nx,:)=0;
rho0_out=1;
rho0_in=1;
% U(:,1)= Uinitial;
% V(:,1)= 0;
% U(:,ny)= Uinitial;
% V(:,ny)= 0;

for i=1:nx
  for j=1:ny
    if isfluid(i,j)==0
      U(i,j)= 0.0;
      V(i,j)= 0.0;
    end
  end
end

for i=1:nx
  for j=1:ny
    if isfluid2(i,j)==0
      U(i,j)= 0.0;
      V(i,j)= 0.0;
    end
  end
end

%% calculating macroscopic velocities


Rhotemp = 0 ;
MomtempX = 0 ;
MomtempY = 0 ;
for I = 1 : nx
  for J = 1 : ny
    if isfluid(I,J) == 1
      for a = 1 : 9
        Rhotemp = Rhotemp + f(I,J,a) ;
        MomtempX = MomtempX + f(I,J,a)*ex(a) ;
        MomtempY = MomtempY + f(I,J,a)*ey(a) ;
      end
      Rho(I,J) = Rhotemp ;
      U(I,J) = MomtempX /Rho(I,J) ;
      V(I,J) = MomtempY /Rho(I,J) ;
      Rhotemp = 0 ;
      MomtempX = 0 ;
      MomtempY = 0 ;
    end
  end
end


%% velocity_calibration

U(1,:)=Uinitial;
V(1,:)=0;
% V(nx,:)=0;
rho0_out=1;
rho0_in=1;
% U(:,1)= Uinitial;
% V(:,1)= 0;
% U(:,ny)= Uinitial;
% V(:,ny)= 0;

for i=1:nx
  for j=1:ny
    if isfluid(i,j)==0
      U(i,j)= 0.0;
      V(i,j)= 0.0;
    end
  end
end

for i=1:nx
  for j=1:ny
    if isfluid2(i,j)==0
      U(i,j)= 0.0;
      V(i,j)= 0.0;
    end
  end
end


%% feq calculation
for j=1:ny
  for i=1:nx
    for a=1:9
      feq(i,j,a)=Rho(i,j) * w(a) * ...
      (1 + 3 * ( ex(a)*U(i,j) + ey(a) * V(i,j) ) + ...
      9.0/2.0 * ( ex(a)*U(i,j) + ey(a) * V(i,j) )^2 -...
      3.0/2.0 *(U(i,j)^2+V(i,j)^2));
    end
  end
end


%%  checking the convergance

Uresidual = 0 ;
Vresidual = 0 ;

for I = 1 : 2*ii
  for J = 1 : ny

    if U(I,J)> 0
      Uresidual = Uresidual + abs(U(I,J)-Uold(I,J))/ U(I,J) ;
    end

    if V(I,J) > 0
      Vresidual = Vresidual + abs(V(I,J)-Vold(I,J))/ V(I,J) ;
    end

  end
end

Uold = U ;
Vold = V ;


vorticity = ((circshift(V, [0, -1]) - circshift(V, [0, 1])) ./ (2.*dx)) - ((circshift(U, [-1, 0]) - circshift(U, [1, 0]))./ (2.*dx));

Vorticity= reshape(vorticity, [nx, ny]);



% for p=1:nx
%     for n=1:ny
% 
%         if Vorticity(nx,ny) < mean(Vorticity)
% 
%             Vorticity(p,n) = nan;
% 
%         end
% 
%     end
% 
% end



%% loop counter
loop_counter = loop_counter +1;
%% animate

if mod(loop_counter,1)==0

  u = (sqrt(U.^2+V.^2));
        quiver(Vorticity',Vorticity', 'Color', 'r', 'LineWidth',2)      

  % imagesc(u');
  axis equal off;
  drawnow
  fprintf(' Iteration = %i, Uresidual = %f, Vresidual = %f \n', loop_counter, Uresidual, Vresidual)
end
end

