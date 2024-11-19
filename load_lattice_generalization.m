clc;
clear;
 
% First of all we have to initialize a system, in our case
% we have to create a matrix of masses, which are connected by a spring

% lattice constants:
M=1;
k=100;
% distance and size:
lattice_width=5;
lattice_height=5;
l0=1;
xdis=1;
ydis=1;

%initial Coordinates
x_matrix=zeros(lattice_height,lattice_width);
y_matrix=zeros(lattice_height,lattice_width);

%initial velocities
vx_matrix=zeros(lattice_height,lattice_width);
vy_matrix=zeros(lattice_height,lattice_width);

%masses on the lattice
mass_matrix=ones(lattice_height,lattice_width);

%precalculation
t(1)=0;
N=10000;
dt=0.008;

for i=1:lattice_height
    for j=1:lattice_width
        x_matrix(i,j)=(j-1)*xdis;
        y_matrix(i,j)=(i-1)*ydis;
    end
end

%with this plot we will see what is the condition for lattice
% before perturbation

subplot(2,2,1)
plot(x_matrix,y_matrix,'o','MarkerFaceColor','black');
axis([ -1 lattice_width*xdis -1 lattice_height*ydis])
title('Initial Lattice')
%___________________________________________________

%perturbation
%___________________________________________________
vx_matrix(:,1)=vx_matrix(:,1)+0.15;
vy_matrix(1,2)=vy_matrix(1,2);
%_____________________
subplot(2,2,2)
plot(x_matrix,y_matrix,'o','MarkerFaceColor','black');
axis([ -1 lattice_width*xdis -1 lattice_height*ydis])
title('Perturbed Lattice')
%___________________________________________________

%Now we will implement RK4 algorithm manually.
%we  have to break down interactions of row loads and column loads
%thats why we will have to write 2 nested for loops
for i=1:N
%-------------
for a = 1:lattice_height
    for b = 1:lattice_width-1
        % Horizontal RK4 calculations

        % Initial orientation
        orient = [x_matrix(a, b+1) - x_matrix(a, b), y_matrix(a, b+1) - y_matrix(a, b)];
        norm_orient = norm(orient);

        % Calculate intermediate accelerations k1, k2, k3, and k4 for both points
        % Step k1
        ac1_k1 = (1 - l0 / norm_orient) * orient * k / mass_matrix(a, b);
        ac2_k1 = -(1 - l0 / norm_orient) * orient * k / mass_matrix(a, b+1);

        % Half step for k2
        orient_k2 = orient + 0.5 * dt * [ac1_k1(1), ac1_k1(2)];
        norm_orient_k2 = norm(orient_k2);
        ac1_k2 = (1 - l0 / norm_orient_k2) * orient_k2 * k / mass_matrix(a, b);
        ac2_k2 = -(1 - l0 / norm_orient_k2) * orient_k2 * k / mass_matrix(a, b+1);

        % Half step for k3
        orient_k3 = orient + 0.5 * dt * [ac1_k2(1), ac1_k2(2)];
        norm_orient_k3 = norm(orient_k3);
        ac1_k3 = (1 - l0 / norm_orient_k3) * orient_k3 * k / mass_matrix(a, b);
        ac2_k3 = -(1 - l0 / norm_orient_k3) * orient_k3 * k / mass_matrix(a, b+1);

        % Full step for k4
        orient_k4 = orient + dt * [ac1_k3(1), ac1_k3(2)];
        norm_orient_k4 = norm(orient_k4);
        ac1_k4 = (1 - l0 / norm_orient_k4) * orient_k4 * k / mass_matrix(a, b);
        ac2_k4 = -(1 - l0 / norm_orient_k4) * orient_k4 * k / mass_matrix(a, b+1);

        % Update velocity matrices using RK4 weighted sum
        vx_matrix(a, b) = vx_matrix(a, b) + dt * (ac1_k1(1) + 2*ac1_k2(1) + 2*ac1_k3(1) + ac1_k4(1)) / 6;
        vy_matrix(a, b) = vy_matrix(a, b) + dt * (ac1_k1(2) + 2*ac1_k2(2) + 2*ac1_k3(2) + ac1_k4(2)) / 6;

        vx_matrix(a, b+1) = vx_matrix(a, b+1) + dt * (ac2_k1(1) + 2*ac2_k2(1) + 2*ac2_k3(1) + ac2_k4(1)) / 6;
        vy_matrix(a, b+1) = vy_matrix(a, b+1) + dt * (ac2_k1(2) + 2*ac2_k2(2) + 2*ac2_k3(2) + ac2_k4(2)) / 6;

        % Update position matrices
        x_matrix(a, b) = x_matrix(a, b) + vx_matrix(a, b) * dt;
        y_matrix(a, b) = y_matrix(a, b) + vy_matrix(a, b) * dt;

        x_matrix(a, b+1) = x_matrix(a, b+1) + vx_matrix(a, b+1) * dt;
        y_matrix(a, b+1) = y_matrix(a, b+1) + vy_matrix(a, b+1) * dt;
    end
end

for c = 1:lattice_width
    for d = 1:lattice_height-1
        % Vertical RK4 calculations

        % Initial orientation
        orient = [x_matrix(d+1, c) - x_matrix(d, c), y_matrix(d+1, c) - y_matrix(d, c)];
        norm_orient = norm(orient);

        % Step k1
        ac1_k1 = (1 - l0 / norm_orient) * orient * k / mass_matrix(d, c);
        ac2_k1 = -(1 - l0 / norm_orient) * orient * k / mass_matrix(d+1, c);

        % Half step for k2
        orient_k2 = orient + 0.5 * dt * [ac1_k1(1), ac1_k1(2)];
        norm_orient_k2 = norm(orient_k2);
        ac1_k2 = (1 - l0 / norm_orient_k2) * orient_k2 * k / mass_matrix(d, c);
        ac2_k2 = -(1 - l0 / norm_orient_k2) * orient_k2 * k / mass_matrix(d+1, c);

        % Half step for k3
        orient_k3 = orient + 0.5 * dt * [ac1_k2(1), ac1_k2(2)];
        norm_orient_k3 = norm(orient_k3);
        ac1_k3 = (1 - l0 / norm_orient_k3) * orient_k3 * k / mass_matrix(d, c);
        ac2_k3 = -(1 - l0 / norm_orient_k3) * orient_k3 * k / mass_matrix(d+1, c);

        % Full step for k4
        orient_k4 = orient + dt * [ac1_k3(1), ac1_k3(2)];
        norm_orient_k4 = norm(orient_k4);
        ac1_k4 = (1 - l0 / norm_orient_k4) * orient_k4 * k / mass_matrix(d, c);
        ac2_k4 = -(1 - l0 / norm_orient_k4) * orient_k4 * k / mass_matrix(d+1, c);

        % Update velocity matrices using RK4 weighted sum
        vx_matrix(d, c) = vx_matrix(d, c) + dt * (ac1_k1(1) + 2*ac1_k2(1) + 2*ac1_k3(1) + ac1_k4(1)) / 6;
        vy_matrix(d, c) = vy_matrix(d, c) + dt * (ac1_k1(2) + 2*ac1_k2(2) + 2*ac1_k3(2) + ac1_k4(2)) / 6;

        vx_matrix(d+1, c) = vx_matrix(d+1, c) + dt * (ac2_k1(1) + 2*ac2_k2(1) + 2*ac2_k3(1) + ac2_k4(1)) / 6;
        vy_matrix(d+1, c) = vy_matrix(d+1, c) + dt * (ac2_k1(2) + 2*ac2_k2(2) + 2*ac2_k3(2) + ac2_k4(2)) / 6;

        % Update position matrices
        x_matrix(d, c) = x_matrix(d, c) + vx_matrix(d, c) * dt;
        y_matrix(d, c) = y_matrix(d, c) + vy_matrix(d, c) * dt;

        x_matrix(d+1, c) = x_matrix(d+1, c) + vx_matrix(d+1, c) * dt;
        y_matrix(d+1, c) = y_matrix(d+1, c) + vy_matrix(d+1, c) * dt;
    end
end

%finally to draw the lattice in motion live 
    subplot(2,2,3)
    plot(x_matrix,y_matrix,'o','MarkerFaceColor','black');
    axis([ -1 lattice_width*xdis -1 lattice_height*ydis])
    title('Lattice in motion')
    drawnow limitrate
    
end  

%___________________________________________________
