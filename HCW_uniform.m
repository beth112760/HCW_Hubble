format long
clear all; close all; clc;

epsilon = 10;
N = 45;% Number of nodes
M = 1; % Number of elements
n = 1.1313658 * 10^(-3);

% Initial Parameteres
Neq = 12;
mu    = 398600.4;
a     = 6778.14+515;%;    % radius of the orbit
n_var = sqrt(mu/(a^3));
tf = pi/n_var;
% Initial conditions for x_bar
% x1 x2 y1 y2 z1 z2
x0 = [50,0.025,50,0.025,50,0.025];
%Final conditions for x_bar
% x1 x2 y1 y2 z1 z2
xf = [0,0,0,0,0,0];

BCtype = 'fixed';%'fixed','free','P0-Pf','P0-Vf','V0-Pf'
BC = [x0,xf]';

initial_guess = [ones(1*N*M,1); ones(1*N*M,1); ones(1*N*M,1); ones(1*N*M,1);ones(1*N*M,1); ones(1*N*M,1); ones(1*N*M,1); ones(1*N*M,1);ones(1*N*M,1); ones(1*N*M,1); ones(1*N*M,1); ones(1*N*M,1)]; 

%% Approximation matrix
D = zeros(N,N,M);
phi = zeros(N,N,M);
phid = zeros(N,N,M);

for k = 1 : M
    t0e            = (k-1)*tf/M; % Initial time of each segment (element)
    tfe            = k*tf/M;     % Final time of each segment (element)
    te(:,:,k)      = linspace(t0e,tfe,N); 
    for i = 1 : N
        phi(i,:,k)        = rbf0(epsilon,te(:,i,k),te(:,:,k));
        phid(i,:,k)       = rbf1(epsilon,te(:,i,k),te(:,:,k));
    end
        D(:,:,k) = phid(:,:,k)/(phi(:,:,k));
end


 
%% CRBF Solution using fsolve
ip=0;
Jlocal_pattern = ones(Neq*N);
J_pattern = kron(eye(M),Jlocal_pattern);

J_local = zeros(Neq*N,Neq*N,M);
J    = zeros(Neq*N*M, Neq*N*M); % Local matrices are to be concatenated sequentially

ns = 0;
Elb   = zeros(1,M);
Erb   = Elb;
SLb   = zeros(Neq,1,M);

SRb   = SLb;
R = zeros(N,M);
sizeJ_local = size(J_local);
LJ_local    = sizeJ_local(1);

f = @(X) HCWNAE_local_unconstrained(X,BC,n_var,D,N,Neq,ns,M,ip,J,J_local,Elb,Erb,SLb,SRb,R,LJ_local,BCtype);
options = optimoptions(@fsolve,'Display','iter','Algorithm','levenberg-marquardt','SpecifyObjectiveGradient',false,'JacobPattern',J_pattern,'StepTolerance',1e-20,'FunctionTolerance',1e-20,'UseParallel',true,'FiniteDifferenceType','central');

[xx,fval1,exitflag1,output1] = fsolve(f,initial_guess,options);

x = reshape(xx,[Neq*N,M]);
X1 = x(1:N,:);
X2 = x(N+1:2*N,:);
Y1 = x(2*N+1:3*N,:);
Y2 = x(3*N+1:4*N,:);
Z1 = x(4*N+1:5*N,:);
Z2 = x(5*N+1:6*N,:);
L1 = x(6*N+1:7*N,:);
L2 = x(7*N+1:8*N,:);
L3 = x(8*N+1:9*N,:);
L4 = x(9*N+1:10*N,:);
L5 = x(10*N+1:11*N,:);
L6 = x(11*N+1:12*N,:);

%% Validation and ERROR
T = te(:);

syms C_1 C_2 C_3 C_4 C_5 C_6

eqn1 =  x0(4)*(2/n_var - (2*cos(n_var*tf))/n_var) - x0(1)*(3*cos(n_var*tf)-4) + (x0(2)*sin(n_var*tf))/n_var + (4*C_2)/(n_var^2) + (16*C_3)/(n_var^3) - (4*C_2*cos(n_var*tf))/(n_var^2) - (16*C_3*cos(n_var*tf))/(n_var^3) + (13*C_1*sin(n_var*tf))/(2*n_var^3) - (11*C_4*sin(n_var*tf))/(n_var^2) - (3*C_3*tf^2)/n_var - (4*C_1*tf)/(n_var^2) + (6*C_4*tf)/n_var - (5*C_1*tf*cos(n_var*tf))/(2*n_var^2) + (5*C_4*tf*cos(n_var*tf))/n_var - (5*C_2*tf*sin(n_var*tf))/(2*n_var) - (5*C_3*tf*sin(n_var*tf))/(n_var^2) == 0;
eqn2 =  x0(2)*cos(n_var*tf) + 2*x0(4)*sin(n_var*tf) + 3*n_var*x0(1)*sin(n_var*tf) - (1/(2*n_var^2))*(8*C_1 - 12*C_4*n_var - 8*C_1*cos(n_var*tf) - 22*C_3*sin(n_var*tf) + 12*C_3*n_var*tf + 12*C_4*n_var*cos(n_var*tf) - 3*C_2*n_var*sin(n_var*tf) + 10*C_3*n_var*tf*cos(n_var*tf) - 5*C_1*n_var*tf*sin(n_var*tf) + 5*C_2*n_var^2*tf*cos(n_var*tf) + 10*C_4*n_var^2*tf*sin(n_var*tf)) == 0;
eqn3 =  x0(3) - x0(2)*(2/n_var - (2*cos(n_var*tf))/n_var) + x0(1)*(6*sin(n_var*tf) - 6*n_var*tf) - x0(4)*(3*tf - (4*sin(n_var*tf))/n_var) + (28*C_4)/(n_var^2) - (16*C_1)/(n_var^3) + (3*C_3*tf^3)/2 - (9*C_4*tf^2)/2 + (16*C_1*cos(n_var*tf))/(n_var^3) - (28*C_4*cos(n_var*tf))/(n_var^2) + (11*C_2*sin(n_var*tf))/(n_var^2) + (38*C_3*sin(n_var*tf))/(n_var^3) + (3*C_1*tf^2)/n_var - (6*C_2*tf)/n_var - (28*C_3*tf)/(n_var^2) - (5*C_2*tf*cos(n_var*tf))/n_var - (10*C_3*tf*cos(n_var*tf))/(n_var^2) + (5*C_1*tf*sin(n_var*tf))/(n_var^2) - (10*C_4*tf*sin(n_var*tf))/n_var ==0;
eqn4 =  x0(4)*(4*cos(n_var*tf)-3) - x0(1)*(6*n_var - 6*n_var*cos(n_var*tf)) - 2*x0(2)*sin(n_var*tf) + (9*C_3*tf^2)/2 - (6*C_2)/n_var - (28*C_3)/(n_var^2) - 9*C_4*tf + (6*C_2*cos(n_var*tf))/n_var + (28*C_3*cos(n_var*tf))/(n_var^2) - (11*C_1*sin(n_var*tf))/(n_var^2) + (18*C_4*sin(n_var*tf))/n_var - 10*C_4*tf*cos(n_var*tf) + 5*C_2*tf*sin(n_var*tf) + (6*C_1*tf)/n_var + (5*C_1*tf*cos(n_var*tf))/n_var + (10*C_3*tf*sin(n_var*tf))/n_var ==0;
eqn5 =  x0(5)*cos(n_var*tf) + (x0(6)*sin(n_var*tf))/n_var - (C_5*(tf*cos(n_var*tf)/2 - sin(n_var*tf)/(2*n_var)))/(n_var^2) - (C_6*tf*sin(n_var*tf))/(2*n_var) ==0;
eqn6 =  x0(6)*cos(n_var*tf) - n_var*x0(5)*sin(n_var*tf) + (C_5*tf*sin(n_var*tf))/(2*n_var) - C_6*((tf*cos(n_var*tf))/2 + (sin(n_var*tf))/(2*n_var)) ==0;

eqns = [eqn1 eqn2 eqn3 eqn4 eqn5 eqn6];
S = solve(eqns);

Svert1 = vertcat(S.C_1);
Svert2 = vertcat(S.C_2);
Svert3 = vertcat(S.C_3);
Svert4 = vertcat(S.C_4);
Svert5 = vertcat(S.C_5);
Svert6 = vertcat(S.C_6);

C1 = sym2poly(Svert1);
C2 = sym2poly(Svert2);
C3 = sym2poly(Svert3);
C4 = sym2poly(Svert4);
C5 = sym2poly(Svert5);
C6 = sym2poly(Svert6);

C_1 = C1;
C_2 = C2;
C_3 = C3;
C_4 = C4;
C_5 = C5;
C_6 = C6;

for i = 1:length(T)
x1(i) =  x0(4)*(2/n_var - (2*cos(n_var*te(i)))/n_var) - x0(1)*(3*cos(n_var*te(i))-4) + (x0(2)*sin(n_var*te(i)))/n_var + (4*C_2)/(n_var^2) + (16*C_3)/(n_var^3) - (4*C_2*cos(n_var*te(i)))/(n_var^2) - (16*C_3*cos(n_var*te(i)))/(n_var^3) + (13*C_1*sin(n_var*te(i)))/(2*n_var^3) - (11*C_4*sin(n_var*te(i)))/(n_var^2) - (3*C_3*te(i)^2)/n_var - (4*C_1*te(i))/(n_var^2) + (6*C_4*te(i))/n_var - (5*C_1*te(i)*cos(n_var*te(i)))/(2*n_var^2) + (5*C_4*te(i)*cos(n_var*te(i)))/n_var - (5*C_2*te(i)*sin(n_var*te(i)))/(2*n_var) - (5*C_3*te(i)*sin(n_var*te(i)))/(n_var^2);
x2(i) =  x0(2)*cos(n_var*te(i)) + 2*x0(4)*sin(n_var*te(i)) + 3*n_var*x0(1)*sin(n_var*te(i)) - (1/(2*n_var^2))*(8*C_1 - 12*C_4*n_var - 8*C_1*cos(n_var*te(i)) - 22*C_3*sin(n_var*te(i)) + 12*C_3*n_var*te(i) + 12*C_4*n_var*cos(n_var*te(i)) - 3*C_2*n_var*sin(n_var*te(i)) + 10*C_3*n_var*te(i)*cos(n_var*te(i)) - 5*C_1*n_var*te(i)*sin(n_var*te(i)) + 5*C_2*n_var^2*te(i)*cos(n_var*te(i)) + 10*C_4*n_var^2*te(i)*sin(n_var*te(i)));
x3(i) =  x0(3) - x0(2)*(2/n_var - (2*cos(n_var*te(i)))/n_var) + x0(1)*(6*sin(n_var*te(i)) - 6*n_var*te(i)) - x0(4)*(3*te(i) - (4*sin(n_var*te(i)))/n_var) + (28*C_4)/(n_var^2) - (16*C_1)/(n_var^3) + (3*C_3*te(i)^3)/2 - (9*C_4*te(i)^2)/2 + (16*C_1*cos(n_var*te(i)))/(n_var^3) - (28*C_4*cos(n_var*te(i)))/(n_var^2) + (11*C_2*sin(n_var*te(i)))/(n_var^2) + (38*C_3*sin(n_var*te(i)))/(n_var^3) + (3*C_1*te(i)^2)/n_var - (6*C_2*te(i))/n_var - (28*C_3*te(i))/(n_var^2) - (5*C_2*te(i)*cos(n_var*te(i)))/n_var - (10*C_3*te(i)*cos(n_var*te(i)))/(n_var^2) + (5*C_1*te(i)*sin(n_var*te(i)))/(n_var^2) - (10*C_4*te(i)*sin(n_var*te(i)))/n_var;
x4(i) =  x0(4)*(4*cos(n_var*te(i))-3) - x0(1)*(6*n_var - 6*n_var*cos(n_var*te(i))) - 2*x0(2)*sin(n_var*te(i)) + (9*C_3*te(i)^2)/2 - (6*C_2)/n_var - (28*C_3)/(n_var^2) - 9*C_4*te(i) + (6*C_2*cos(n_var*te(i)))/n_var + (28*C_3*cos(n_var*te(i)))/(n_var^2) - (11*C_1*sin(n_var*te(i)))/(n_var^2) + (18*C_4*sin(n_var*te(i)))/n_var - 10*C_4*te(i)*cos(n_var*te(i)) + 5*C_2*te(i)*sin(n_var*te(i)) + (6*C_1*te(i))/n_var + (5*C_1*te(i)*cos(n_var*te(i)))/n_var + (10*C_3*te(i)*sin(n_var*te(i)))/n_var;
x5(i) =  x0(5)*cos(n_var*te(i)) + (x0(6)*sin(n_var*te(i)))/n_var - (C_5*(te(i)*cos(n_var*te(i))/2 - sin(n_var*te(i))/(2*n_var)))/(n_var^2) - (C_6*te(i)*sin(n_var*te(i)))/(2*n_var);
x6(i) =  x0(6)*cos(n_var*te(i)) - n_var*x0(5)*sin(n_var*te(i)) + (C_5*te(i)*sin(n_var*te(i)))/(2*n_var) - C_6*((te(i)*cos(n_var*te(i)))/2 + (sin(n_var*te(i)))/(2*n_var));
end

xx1 = vertcat(x1);
xx2 = vertcat(x2);
yy1 = vertcat(x3)';
yy2 = vertcat(x4)';
zz1 = vertcat(x5)';
zz2 = vertcat(x6)';

load("uncon_HCW_results_GPOPS_Hubble_v2.mat")
timevec = [solution.phase(1).time];
te = timevec';
x_1 = X1(:);
x_2 = X2(:);
y_1 = Y1(:);
y_2 = Y2(:);
z_1 = Z1(:);
z_2 = Z2(:);


x1vecc = [solution.phase(1).state(:,[1])];
 y1vecc= [solution.phase(1).state(:,[2])];
z1vecc = [solution.phase(1).state(:,[3])];
 x2vecc= [solution.phase(1).state(:,[4])];
 y2vecc = [solution.phase(1).state(:,[5])];
 z2vecc= [solution.phase(1).state(:,[6])];

x1Vec = interp1(timevec,x1vecc,T,'linear');
x2Vec = interp1(timevec,x2vecc,T,'linear');
y1Vec = interp1(timevec,y1vecc,T,'linear');
y2Vec = interp1(timevec,y2vecc,T,'linear');
z1Vec = interp1(timevec,z1vecc,T,'linear');
z2Vec = interp1(timevec,z2vecc,T,'linear');



x1vec = x1Vec(:);
x2vec = x2Vec(:);
y1vec = y1Vec(:);
y2vec = y2Vec(:);
z1vec = z1Vec(:);
z2vec = z2Vec(:);

abserror_x1 = (x_1 - x1vec);
abserror_x2 = (x_2 - x2vec);
abserror_y1 = (y_1 - y1vec);
abserror_y2 = (y_2 - y2vec);
abserror_z1 = (z_1 - z1vec);
abserror_z2 = (z_2 - z2vec);


figure(1)
plot(T,x_1,'*')
hold on
plot(T,y_1,'*')
plot(T,z_1,'*')
plot(T,x1vec,'o-')
plot(T,y1vec,'o-')
plot(T,z1vec,'o-')
plot(T,xx1,'LineWidth',1.25)
plot(T,yy1,'LineWidth',1.25)
plot(T,zz1,'LineWidth',1.25)
xlabel('Time')
ylabel('Position')
title('Unconstrained HCW Equations')
legend('X Position (CRBF)','Y Position (CRBF)','Z Position (CRBF)','X Position (GPOPS-II)','Y Position (GPOPS-II)','Z Position (GPOPS-II)','X Position (Closed-Form)','Y Position (Closed-Form)','Z Position (Closed-Form)')

figure(2)
plot(T,x_2,'*')
hold on
plot(T,y_2,'*')
plot(T,z_2,'*')
plot(T,x2vec,'o-')
plot(T,y2vec,'o-')
plot(T,z2vec,'o-')
plot(T,xx2,'LineWidth',1.25)
plot(T,yy2,'LineWidth',1.25)
plot(T,zz2,'LineWidth',1.25)
xlabel('Time')
ylabel('Velocity')
title('Unconstrained HCW Equations')
legend('X Velocity (CRBF)','Y Velocity (CRBF)','Z Velocity (CRBF)','X Velocity (GPOPS-II)','Y Velocity (GPOPS-II)','Z Velocity (GPOPS-II)','X Position (Closed-Form)','Y Position (Closed-Form)','Z Position (Closed-Form)')


figure(3)
plot3(x_1,y_1,z_1,'*')
hold on

plot3(x1vec,y1vec,z1vec,'o-')

plot3(xx1,yy1,zz1,'LineWidth',1.25)
xlabel('X Position')
ylabel('Y Position')
zlabel('Z Position')

plot3(x_1(1,1),y_1(1,1),z_1(1,1),'.r','MarkerSize',30)
plot3(x_1(end,1),y_1(end,1),z_1(end,1),'.g','MarkerSize',30)
legend('Optimal Trajectory (CRBF)','Optimal Trajectory (GPOPS-II)','Optimal Trajectory (Closed-Form)','Start','End')

figure;
plot(T(1:end-2),abs(abserror_x1(1:end-2)))
hold on
plot(T(1:end-1),abs(abserror_y1(1:end-1)))
plot(T(1:end-1),abs(abserror_z1(1:end-1)))
xlabel('Time')
ylabel('Absolute Error')
title('Position Absolute Error between CRBF-Local Collocation and GPOPS-II')
legend('X','Y','Z')


figure;
plot(T(1:end-1),abs(abserror_x2(1:end-1)))
hold on
plot(T(1:end-1),abs(abserror_y2(1:end-1)))
plot(T(1:end-1),abs(abserror_z2(1:end-1)))
xlabel('Time')
ylabel('Absolute Error')
title('Velocity Absolute Error between CRBF-Local Collocation and GPOPS-II')
legend('X','Y','Z')

abserror_x1 = (x_1 - xx1');
abserror_x2 = (x_2 - xx2');
abserror_y1 = (y_1 - yy1);
abserror_y2 = (y_2 - yy2);
abserror_z1 = (z_1 - zz1);
abserror_z2 = (z_2 - zz2);

figure;
plot(T(1:end-2),abs(abserror_x1(1:end-2)))
hold on
plot(T(1:end-1),abs(abserror_y1(1:end-1)))
plot(T(1:end-1),abs(abserror_z1(1:end-1)))
xlabel('Time')
ylabel('Absolute Error')
title('Position Absolute Error between CRBF-Local Collocation and Closed-Form Solution')
legend('X','Y','Z')

 
figure;
plot(T(1:end-1),abs(abserror_x2(1:end-1)))
hold on
plot(T(1:end-1),abs(abserror_y2(1:end-1)))
plot(T(1:end-1),abs(abserror_z2(1:end-1)))
xlabel('Time')
ylabel('Absolute Error')
title('Velocity Absolute Error between CRBF-Local Collocation and Closed-Form Solution')
legend('X','Y','Z')

abserror_x1 = (x1vec - xx1');
abserror_x2 = (x2vec - xx2');
abserror_y1 = (y1vec - yy1);
abserror_y2 = (y2vec - yy2);
abserror_z1 = (z1vec - zz1);
abserror_z2 = (z2vec - zz2);

figure;
plot(T(1:end-2),abs(abserror_x1(1:end-2)))
hold on
plot(T(1:end-1),abs(abserror_y1(1:end-1)))
plot(T(1:end-1),abs(abserror_z1(1:end-1)))
xlabel('Time')
ylabel('Absolute Error')
title('Position Absolute Error between GPOPS-II and Closed-Form Solution')
legend('X','Y','Z')


figure;
plot(T(1:end-1),abs(abserror_x2(1:end-1)))
hold on
plot(T(1:end-1),abs(abserror_y2(1:end-1)))
plot(T(1:end-1),abs(abserror_z2(1:end-1)))
xlabel('Time')
ylabel('Absolute Error')
title('Velocity Absolute Error between GPOPS-II and Closed-Form Solution')
legend('X','Y','Z')