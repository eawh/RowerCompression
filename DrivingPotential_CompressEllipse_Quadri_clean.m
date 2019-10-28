%% Quadriflagellates 
% The ELLIPTICAL method, where
% both the force direction and the displacement rate are used to determine
% where the switch points occur and the potentials are fitted
% you find all that you need in

%% Changing the directory containing the clicked data and loading the data
cd /scratch/Rotors/QuadriTracking
force = load('quadri_gallop_3000fps_leftclclk_int.clclk_force','-mat');
% force = load('quadri_trot_3000fps_leftclclk_int.clclk_force','-mat');
% force = load('quadri_trot_3000fps_rightclclk_int.clclk_force','-mat');

% INFORMATION on what the object loaded actually includes
%%% force.Stroke.mean_cilium_length
%%%   force.Force.displ.cod_cy centre of drag x
%%% force.Force.displ.cod_cy centre of drag y
%%% force.Force.displ.commoncyls_F_y_pN   force y
%%% force.Force.displ.commoncyls_F_x_pN   force x
%%%  force.Force.displ.commoncyls_F_para_pN  force para to the cell, to the
%%%  surface of the cell
%%%  force.Force.displ.commoncyls_F_perp_pN  force perp to the cell

%% STEP 1: defining variables for the position and force data from the object

pix2mic = force.px2mum; % the  pixel to micron ratio from the video. 

% The centre of drag for x (cod_x) and y (cod_y) measured in pixels
cod_x = arrayfun(@(x) x.cod_cx, force.Force.displ);
cod_y = arrayfun(@(x) x.cod_cy, force.Force.displ);

% Converting the centre of drag position to pixels and redefining the centre to being (0,0)
cod_x = (cod_x-mean(cod_x))*pix2mic;
cod_y = (cod_y - mean(cod_y))*pix2mic;

% The force along x and y measured from the clicking. The average for three neighbouring measurements is used for each point
Fx = smooth(arrayfun(@(x) x.commoncyls_F_x_pN, force.Force.displ),3);
Fy = smooth(arrayfun(@(x) x.commoncyls_F_y_pN, force.Force.displ),3);

% The time for each frame clicked, not always clicked at equal time points apart
timestamp = arrayfun(@(x) x.timestamp_s, force.Force.displ);


%% STEP 2: Fitting an ellipse to the data
cd /scratch % changing to the directory with the fit_ellipse function. This function was downloaded off the matlab function exchange
traj_ellipse = fit_ellipse(cod_x, cod_y); 

%% STEP 2b: to look at the fitted ellipse x and y values are found describing the fit
ellipsephi = linspace(0,2*pi);

% Calculating the untilted trajectory
xellipse_w = traj_ellipse.a*cos(ellipsephi);
yellipse_w = traj_ellipse.b*sin(ellipsephi);

% Then rotating the ellipse 
%rotmat = rotx(traj_ellipse.phi);
%rotmat = rotmat(2:end,2:end);
xellipse = xellipse_w*cos(-traj_ellipse.phi) - yellipse_w*sin(-traj_ellipse.phi); 
yellipse = xellipse_w*sin(-traj_ellipse.phi) + yellipse_w*cos(-traj_ellipse.phi); 

 % Plotting the linked points, the points with the force and the fitted ellipse.
 figure(6)
   plot(smooth(cod_x,1), smooth(cod_y,1),'k')
hold on
  scatter(cod_x, cod_y,[],1:length(cod_y),'filled')
 quiver(cod_x,cod_y, Fx', Fy')
plot(xellipse,yellipse ,'k--') 
hold off

axis equal
 
%% STEP 3: Having checked the elliptical fit works, rotate the points

% define the major axis angle of the major axis
% When the y direction is longer, that rotates the ellipse by 90 degrees
major_axis = pi/2*(traj_ellipse.b>traj_ellipse.a)-traj_ellipse.phi;

% Then rotate all the x and y coordinates so the major axis aligns with the
% x axis. Only interesting in the resulting y values for the sign of the
% force
relative_maj = -cod_x*sin(major_axis) + cos(major_axis)*cod_y; % this the displacement above or below the major axis, which groups the points into the two different strokes

% the position is the distance parallel to the major axis
dist_maj = cod_x*cos(major_axis) + sin(major_axis)*cod_y;

%% STEP 3a: Plotting to distinguish between the power and recovery sections

% Plotting the trajectory including the newly ascribed sign through colour
figure(7)
   scatter(cod_x(relative_maj<0), cod_y(relative_maj<0),'b','filled')
hold on
   scatter(cod_x(relative_maj>0), cod_y(relative_maj>0),'r','filled')
   set(gca,'ColorOrderIndex',1)
 quiver(cod_x(relative_maj<0),cod_y(relative_maj<0), Fx(relative_maj<0)', Fy(relative_maj<0)')
 quiver(cod_x(relative_maj>0),cod_y(relative_maj>0), Fx(relative_maj>0)', Fy(relative_maj>0)')
plot(xellipse,yellipse ,'k--') 
hold off

axis equal
box on
xlabel('xpos')
ylabel('ypos')

% Now plotting the distance and force for the rower potential

figure(8)
scatter(dist_maj(relative_maj<0), hypot(Fx(relative_maj<0), Fy(relative_maj<0)),'r','filled')
hold on
scatter(dist_maj(relative_maj>0), -hypot(Fx(relative_maj>0), Fy(relative_maj>0)),'b','filled')
hold off

ylabel('Force pN')
xlabel('$d_{maj}\,\mu$m','interpreter','latex','Fontsize',16)

%% STEP 4: Then fitting the profiles
% then separating into sides
% The choice of which is postive or negative depends on whether the
% trajectory is going cw or ccw

% Defining the force as positive when moving right, which depends on the direction of the orbit, i.e. clockwise (cw) or counterclockwise(caw)

% % This is for the trajectory labelled left which goes cw
% dmaj_neg = dist_maj(relative_maj>0);
% dmaj_pos = dist_maj(relative_maj<0);
% F_neg = -hypot(Fx(relative_maj>0), Fy(relative_maj>0));
% F_pos = hypot(Fx(relative_maj<0), Fy(relative_maj<0));

% This is for the trajectory labelled right which goes ccw
dmaj_neg = dist_maj(relative_maj<0);
dmaj_pos = dist_maj(relative_maj>0);
F_neg = -hypot(Fx(relative_maj<0), Fy(relative_maj<0));
F_pos = hypot(Fx(relative_maj>0), Fy(relative_maj>0));


% For the distance on the positive side, I need to switch the directions.
% Since relative to the trap the distance should be decreasing over the
% simulation
xplot = linspace(0, mean([range(dmaj_pos),range(dmaj_neg)]));

dmaj_pos = max(xplot) -(dmaj_pos+max(xplot)/2);

% fitting the polynomial to the distance and forces
fit_maj_pos = polyfit(dmaj_pos, F_pos',3);
fit_maj_neg = polyfit(dmaj_neg+max(xplot)/2 , F_neg',3);

%% STEP 4a: Looking at the fits

figure(9)
plot(xplot, polyval(fit_maj_neg,xplot),'LineWidth',1.5)
hold on
plot(xplot, polyval(fit_maj_pos,xplot),'LineWidth',1.5)
set(gca,'ColorOrderIndex',1)
scatter(dmaj_pos, F_pos,'r','filled')
scatter(dmaj_neg+max(xplot)/2, F_neg,'b','filled')
hold off

ylabel('Force pN')
xlabel('position microns')

% Another figure than colour codes them but plots them in order of
% appearance. To help me know if I need to reflect the fit in the negative
% case
redind = find(relative_maj<0);
blueind = find(relative_maj>0);

figure(10)
plot(redind, hypot(Fx(relative_maj<0), Fy(relative_maj<0)),'x')
hold on
plot(blueind, hypot(Fx(relative_maj>0), Fy(relative_maj>0)),'x')
hold off

xlabel('index')
ylabel('force')



% Plotting the trajectory including the newly ascribed sign through colour
% This is for the right flagella where I swapped the colour because the
% direction of the oscillation switched
figure(11)
   scatter(cod_x(relative_maj>0), cod_y(relative_maj>0),'b','filled')
hold on
   scatter(cod_x(relative_maj<0), cod_y(relative_maj<0),'r','filled')
   set(gca,'ColorOrderIndex',1)
 quiver(cod_x(relative_maj>0),cod_y(relative_maj>0), Fx(relative_maj>0)', Fy(relative_maj>0)')
 quiver(cod_x(relative_maj<0),cod_y(relative_maj<0), Fx(relative_maj<0)', Fy(relative_maj<0)')
plot(xellipse,yellipse ,'k--') 
hold off

axis equal
box on
xlabel('xpos')
ylabel('ypos')

% 
% 
% set(gcf,'PaperPositionMode','auto')
% set(gcf,'InvertHardCopy','off')
% set(gcf,'Color','none')
% 
% 
