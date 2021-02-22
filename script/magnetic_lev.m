

m=0.02;
km=2.058E-4;
Km=km/m;
R=0.92;
L=10^-5;
g=9.81;

hmax=0.07;
hmin=0.03;
hstart=0.03;

imax=3;
imin=0.1;
istart=0;

hdot_start=0;

X1=[0.03;0.05;0.08]; x2=0; X3=[0.1;1.5;3];
B=[0;0;1/L];

Astar=zeros(3,3,9);
dstar=zeros(3,9);
k=1;

for i=1:3
    x1=X1(i);
 
    for j=1:3
        x3=X3(j);
        
        xstar=[x1;x2;x3];
        f=[x2;g-Km*(x3/x1)^2;-R/L*x3];

        astar=2*Km*x3^2/x1^3;
        bstar=-2*Km*x3/x1^2;

        Astar(:,:,k)=[0 1 0;astar 0 bstar; 0 0 -R/L];
        dstar(:,k)=f-Astar(:,:,k)*xstar;
        k=k+1;
    end
end
%% Riccati Equation
%finding the gains Ki
B=[B;0];
Q=[eye(3), zeros(3,1); 0 0 0 10^4];
K=zeros(9,4);

for i=1:9
A=[Astar(:,:,i), zeros(3,1);1 0 0 0];
[X,~,G]=care(A,B,Q);
K(i,:)=-G;
end

astar=Astar(2,1,1);
bstar=Astar(2,3,1);
