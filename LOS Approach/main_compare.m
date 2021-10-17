c = 10; %sliding mode coefficient
M = 0.25; % gain

%Defining initial values of all parameters
%Constant velocities
V_a = 100; 
V_d = 200;
V_m = 200;
V = [V_a; V_d; V_m];

%Initial position coordinates
x_a = 0;
y_a = 0;
x_d = 0;
y_d = 0;
x_m = 5000;
y_m = 0;

%Initial heading angles
g_a = 45*pi/180;
g_d = 0;
g_m = 150*pi/180;

p = [x_a, y_a ; x_d, y_d ; x_m, y_m]; %2D position array
g = [g_a; g_d; g_m]; %Heading angle array

%LOS angles
l_am = atan((y_m-y_a)/(x_m-x_a));
l_dm = atan((y_m-y_d)/(x_m-x_d));
l = [l_am; l_dm]; %lambda array

%LOS distances
r_am = ((x_a-x_m)^2 + (y_a-y_m)^2)^0.5;
r_dm = ((x_d-x_m)^2 + (y_d-y_m)^2)^0.5;
r = [r_am; r_dm]; %Los distance array

% Rate of change of LOS distance
r_dot_am = V_m*cos(g_m-l_am) - V_a*cos(g_a-l_am);
r_dot_dm = V_m*cos(g_m - l_dm) - V_d*cos(g_d - l_dm);
r_dot = [r_dot_am; r_dot_dm];

% Rate of change of LOS angle
l_dot_am =( V_m*sin(g_m - l_am) - V_a*sin(g_a - l_am) )/r_am;
l_dot_dm =( V_m*sin(g_m - l_dm) - V_d*sin(g_d - l_dm) )/r_dm;
l_dot = [l_dot_am ; l_dot_dm];

%Defining the angle between two LOS
phi = l_am - l_dm;
phi_dot = l_dot_am - l_dot_dm;
x2 = phi_dot;
x1 = phi;

%Defining the sliding surface
S = x2 + c*x1;

%Lateral Accelerations
a_a = -3*r_dot_am*l_dot_am; %Aircraft accelrarion is PN stratergy


% begin propagating
t_f = 15;
h = 0.01;

xm_PN(1) = x_m;
ym_PN(1) = y_m;


n = floor(t_f/h);
time_PN(1) = 0;
time_PN(2) = h;


for ii = 2:n        % Run Integrator, step-size = h_RK4, final time = t_f
    [p, g, l, l_dot, r, r_dot, aa_new, ad_new, am_new, S, phi, phi_dot] = RK4_PN_target(p, g, l, l_dot, r, r_dot, a_a, h, c, M, S, V);
    
    if abs(ad_new) >= 20*9.81
        ad_new = 20*9.81 * sign(ad_new);
    end
    
    aaa_PN(ii) = aa_new;
    amm_PN(ii) = am_new;
    add_PN(ii) = ad_new;
    
    Sarray_PN(ii) = S;
    
    phi_array_PN(ii) = phi;
    phi_dot_array_PN(ii) = phi_dot;
    
    xm_PN(ii) = p(3,1);
    ym_PN(ii) = p(3,2);
    xa_PN(ii) = p(1,1);
    ya_PN(ii) = p(1,2);
    xd_PN(ii) = p(2,1);
    yd_PN(ii) = p(2,2);
    
    
    if r(2)< 1
        break;
    end
    
    time_PN(ii+1) = time_PN(ii) + h;
    
end



%%%%%%%%%%%%%%%%%%%%%%


%Resetting all values for 0 acceleartion target simulation

%Defining initial values of all parameters
%Constant velocities
V_a = 100; 
V_d = 200;
V_m = 200;
V = [V_a; V_d; V_m];

%Initial position coordinates
x_a = 0;
y_a = 0;
x_d = 0;
y_d = 0;
x_m = 5000;
y_m = 0;

%Initial heading angles
g_a = 45*pi/180;
g_d = 0;
g_m = 150*pi/180;

p = [x_a, y_a ; x_d, y_d ; x_m, y_m]; %2D position array
g = [g_a; g_d; g_m]; %Heading angle array

%LOS angles
l_am = atan((y_m-y_a)/(x_m-x_a));
l_dm = atan((y_m-y_d)/(x_m-x_d));
l = [l_am; l_dm]; %lambda array

%LOS distances
r_am = ((x_a-x_m)^2 + (y_a-y_m)^2)^0.5;
r_dm = ((x_d-x_m)^2 + (y_d-y_m)^2)^0.5;
r = [r_am; r_dm]; %Los distance array

% Rate of change of LOS distance
r_dot_am = V_m*cos(g_m-l_am) - V_a*cos(g_a-l_am);
r_dot_dm = V_m*cos(g_m - l_dm) - V_d*cos(g_d - l_dm);
r_dot = [r_dot_am; r_dot_dm];

% Rate of change of LOS angle
l_dot_am =( V_m*sin(g_m - l_am) - V_a*sin(g_a - l_am) )/r_am;
l_dot_dm =( V_m*sin(g_m - l_dm) - V_d*sin(g_d - l_dm) )/r_dm;
l_dot = [l_dot_am ; l_dot_dm];

%Defining the angle between two LOS
phi = l_am - l_dm;
phi_dot = l_dot_am - l_dot_dm;
x2 = phi_dot;
x1 = phi;

%Defining the sliding surface
S = x2 + c*x1;

%Lateral Accelerations
a_a = 0; %Aircraft accelrarion considered to be 0


% begin propagating
t_f = 15;
h = 0.01;

xm(1) = x_m;
ym(1) = y_m;


n = floor(t_f/h);
time(1) = 0;
time(2) = h;


for ii = 2:n        % Run Integrator, step-size = h_RK4, final time = t_f
    [p, g, l, l_dot, r, r_dot, aa_new, ad_new, am_new, S, phi, phi_dot] = RK4_new(p, g, l, l_dot, r, r_dot, a_a, h, c, M, S, V);
    
    if abs(ad_new) >= 20*9.81
        ad_new = 20*9.81 * sign(ad_new);
    end
    
    aaa(ii) = aa_new;
    amm(ii) = am_new;
    add(ii) = ad_new;
    
    Sarray(ii) = S;
    
    phi_array(ii) = phi;
    phi_dot_array(ii) = phi_dot;
    
    xm(ii) = p(3,1);
    ym(ii) = p(3,2);
    xa(ii) = p(1,1);
    ya(ii) = p(1,2);
    xd(ii) = p(2,1);
    yd(ii) = p(2,2);
    
    
    if r(2)< 1
        break;
    end
    
    time(ii+1) = time(ii) + h;
    
end



figure();
plot(xm, ym, 'b');
hold on
plot(xa, ya, 'g');
plot(xd, yd, 'r');
plot(xm_PN, ym_PN, 'Color','b', 'LineStyle','--');
plot(xa_PN, ya_PN, 'Color','g','LineStyle', '--');
plot(xd_PN, yd_PN, 'Color','r','LineStyle', '--');
hold off
xlabel('x');
ylabel('y');


figure();
plot(time,Sarray);
hold on
plot(time_PN,Sarray_PN, 'LineStyle','--');
hold off
xlabel('t');
ylabel('S');

figure();
plot(time,phi_dot_array);
hold on
plot(time_PN,phi_dot_array_PN, 'LineStyle','--');
hold off
xlabel('t');
ylabel('phi_dot');

figure();
plot(phi_array, phi_dot_array);
hold on
plot(phi_array_PN, phi_dot_array_PN, 'LineStyle','--');
hold off
xlabel('phi');
ylabel('phi_dot');

figure();
plot(time,amm);
hold on
plot(time_PN,amm_PN,'LineStyle', '--');
hold off
xlabel('t');
ylabel('am');

figure();
plot(time,aaa);
hold on
plot(time_PN,aaa_PN , 'LineStyle','--');
hold off
xlabel('t');
ylabel('aa');

figure();
plot(time,add);
hold on
plot(time_PN,add_PN,'LineStyle', '--');
hold off
xlabel('t');
ylabel('ad');
