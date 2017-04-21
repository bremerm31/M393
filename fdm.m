%This is a reproduction of 
%    Numerical simulations of the energy-supercritical nonlinear
%    schroedinger Equation
%To begin we consider the radial domain from (0,Rmax)
Rmax = 200.;
%we discretize the domain with N points
N = 1000;

h = Rmax/(N-1);

R = linspace(0,Rmax,N);
U = 1i*ones(N+4,2); %add a little extra for BCs
K = 1i*ones(N,4);

%we now will write the ODE as dUdt=F(U)
dt = 0.001;
tend = 1.;
t = 0.;
it = 1;

while ( tend - t > 100*eps(tend) )
  
  K(:,1) = assembleRHS(U(:,it),h,R);
  K(:,2) = assembleRHS(U(:,it)+dt/2*K(:,1),h,R);
  K(:,3) = assembleRHS(U(:,it)+dt/2*K(:,2),h,R);
  K(:,4) = assembleRHS(U(:,it)+dt*K(:,2),h,R);
  
  it = 3 - it;
  U(:,it) = U(:,3-it) - 1i*dt/6*sum(K,2);
  
  t = t + dt;
end