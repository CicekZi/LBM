%% Cylinder nodes

isfluid = ones(nx,ny); % all nodes are part of fluid unless we note
isfluid2 = ones(nx,ny); % all nodes are part of fluid unless we note


% isfluid(:,1)  = 0; % south wall/ down wall
% isfluid(:,ny) = 0;% north wall / upper wall

ii=floor(nx/5+1);
jj=floor(175);

ii2=floor(nx/5 +1 );
jj2=floor(75);

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