function [] = laminar_flow()

% Constants:
visc = 1e-3; %viscosity, Pa*s
rho = 1e3; %density, kg/m3
V_up = 5e-5; %velocity, m/s
B = 1e-3; %distance, m
dp_dx = -1; %Pa/m



%Grid Parameters:
N = 1000; %Number of grid points
dy = B/(N+1); %grid spacing, m
G = -dy^2/visc*dp_dx;


%Define matrix

A = spalloc(N,N,3*N-2);
b = zeros(N,1);

j = 1;
A(1,1) = 2;
A(1,2) = -1;
b(1) = G;

for j=2:N-1
    A(j,j-1)=-1;
    A(j,j)=2;
    A(j,j+1)=-1;
    b(j)=G;
end

j=N;
A(N,N-1)=-1;
A(N,N)=2;
b(N)=G+V_up;

v=A\b;

y = linspace(dy,B-dy,N)

v_anal = V_up*y/B + 1/2/visc*dp_dx*(y.^2-B*y);

plot(y,v_anal,y,v,'.')
xlabel('y [m]')
ylabel('v [m/s]')
legend('v_a_n_a_l','v_f_d')
end
