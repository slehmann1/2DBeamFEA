%%% Geometric Characteristics %%%
r = 11.25/1000; % Beam Radius (m)
A = pi * (r)^2;
I = pi/4*(r)^4;
E = 200 * 10^9; % Pa

% Nodal coordinates: xi yi pi
nodes = [0 0 0;
    1 0 0;
    2 0 0;
    0.5 sin(pi/4) 0;
    1.5 sin(pi/4) 0];

n = size(nodes,1); % number of nodes

% Elements: nodeid, nodeid
elem = [1 2;
    2 3;
    1 4;
    4 2;
    4 5;
    2 5;
    5 3];

m = size(elem,1); % number of elements

% Stiffness matrix
K = zeros(3*n,3*n);

% Define stiffness matrix
for i = 1:m
    a = elem(i,1);
    xa = nodes(a,1);
    ya = nodes(a,2);
    
    b = elem(i,2);
    xb = nodes(b,1);
    yb = nodes(b,2);
    
    L = sqrt((yb-ya)^2+(xb-xa)^2);
    
    Kh = E*A/L;
    Kv = 12*E*I/L^3;
    Kv_bar = 6*E*I/L^2;
    Km=4*E*I/L;
    
    theta = atan2(yb-ya,xb-xa);
    
    % Define rotation matrix
    R = [cos(theta) sin(theta) 0 0 0 0;
        -sin(theta) cos(theta), 0 0 0 0;
        0, 0, 1 0 0 0;
        0 0 0 cos(theta) sin(theta) 0;
        0 0 0 -sin(theta) cos(theta) 0
        0 0 0 0 0 1];
    
    % Define unrotated elemental stiffness matrix
    k_q_el = [Kh 0 0 -Kh 0 0;
        0 Kv Kv_bar 0 -Kv Kv_bar;
        0 Kv_bar Km 0 -Kv_bar Km/2;
        -Kh 0 0 Kh 0 0;
        0 -Kv -Kv_bar 0 Kv -Kv_bar;
        0 -Kv_bar Km/2 0 -Kv_bar Km];
    
    % Define rotated elemental stiffness matrix
    k_el = R * k_q_el * R';
    
    % Update overall stiffness matrix
    a = (a-1)*3+1;
    b = (b-1)*3+1;
    K(a:a+2, a:a+2) = K(a:a+2, a:a+2) + k_el(1:3, 1:3);
    K(a:a+2, b:b+2) = K(a:a+2, b:b+2) + k_el(1:3, 4:6);
    K(b:b+2, a:a+2) = K(b:b+2, a:a+2) + k_el(4:6, 1:3);
    K(b:b+2, b:b+2) = K(b:b+2, b:b+2) + k_el(4:6, 4:6);
end

% Boundary conditions
fixed_dofs = [1 2 3 4 5 6 7 8 9]; % dofs constrained to be 0
free_dofs = setdiff([1:3*n],fixed_dofs);

P_app = -600000;
f = [0 0 0 0 0 0 0 0 0 0 P_app 0 0 P_app 0];

% Solution
u = zeros(3*n,1);
u(free_dofs) = f(free_dofs)/K(free_dofs,free_dofs);

disp('Stiffnes Matrix: ');
disp(K)

% Print nodal results
for i =1:n
    % undeformed configuration
    x = nodes(i,1);
    y = nodes(i,2);
    p=nodes(i,3);
    
    % deformed configuration
    xd = x + u(3*i-2);
    yd = y + u(3*i-1);
    pd = p + u(3*i);
    disp(['Undeformed: ',num2str(x),' ', num2str(y), ' ', num2str(p) ,' Deformed: ', num2str(xd),' ', num2str(yd), ' ', num2str(pd)])
end

% Plot results
figure(1)
hold on

for q =1:m
    i = elem(q,1);
    j = elem(q,2);
    
    % undeformed configuration
    xi = nodes(i,1);
    yi = nodes(i,2);
    xj = nodes(j,1);
    yj = nodes(j,2);
    
    % deformed configuration
    xid = xi+u(3*i-2)*10;
    yid = yi+u(3*i-1)*10;
    xjd = xj+u(3*j-2)*10;
    yjd = yj+u(3*j-1)*10;
    
    figure(1)
    % undeformed; deformed;
    plot([xi,xj],[yi,yj],'k-o',[xid,xjd],[yid,yjd],'b-o')
    xlim([0,2])
    ylim([0,1])
end

title('Structural displacements','fontsize',15)
legend({'Undeformed cofiguration','Deformed cofiguration (10 X Scale)'},'Location','southeast')