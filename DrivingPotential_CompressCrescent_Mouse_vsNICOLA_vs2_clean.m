%% Mouse
% The CRESCENT method is used for the mouse data and airways data, where
% both the force direction and the displacement rate are used to determine
% where the switch points occur and the potentials are fitted
% Version 2 includes the time stamp when fitting. Also the mouse values
% should use parallel and perpendicular coordinates for not the NICOLA data

% This particular code should replicate the results for the video from Pelliciotta2020 in PNAS for the mouse brain cell. It can be applied to other crescent shaped trajectories but force and position should be changed to coordinates relative to the cell wall, see the commented section before STEP 1a. At various stages this script requires human decision making and is not automated, which will require tailoring for another trajectory. Key points are: 
 - data input for the fitted trajectory, the fitting protocol used to find the trajectory can only take complete number of cycles in the input. Appears at the beginning of STEP 2
 - whether it is appropriate to take the average frequency of the two fits. Appears at the end of STEP 2
 - The correct indices are chosen for the endpoints. Appears in STEP 3
 - Worth double checking at STEP 8 that the distance is always decreasing towards the trap.

% This defines a variable that is later defined in a self-dependent if statement. This prevents errors from an undefined variable.
fitrange = 1;

% Changing to the directory to load the structure created by clicking the video
cd /scratch/Rotors/MouseCilium

% %you find all that you need in
force = load('brain_cilium_v2.clclk_force','-mat');

%% Some INFORMATION about the loaded structure
%%% force.Stroke.mean_cilium_length
%%%   force.Force.displ.cod_cy centre of drag x
%%% force.Force.displ.cod_cy centre of drag y
%%% force.Force.displ.commoncyls_F_y_pN   force y
%%% force.Force.displ.commoncyls_F_x_pN   force x
%%%  force.Force.displ.commoncyls_F_para_pN  force para to the cell, to the
%%%  surface of the cell
%%%  force.Force.displ.commoncyls_F_perp_pN  force perp to the cell

%% STEP 1: Defining appropriate variables for the position, force and time
% don't use the video value, instead this is what Nicola says is correct for the pixel to micron value
pix2mic =  0.13*40/90;
% If structure value is correct use; pix2mic = force.px2mum;

% This is the option to use when using NICOLA data. The cell wall couldn't be found properly for this cilium so use the x and y data instead. The video was oriented so that x is the direction of the power stroke
% The position of the centre of drag in x (cod_x) and in y (cod_y) measured in pixels
cod_x = arrayfun(@(x) x.cod_cx, force.Force.displ);
cod_y = arrayfun(@(x) x.cod_cy, force.Force.displ);

% % the force from the cilium in x and y
% Fx = smooth(arrayfun(@(x) x.commoncyls_F_x_pN, force.Force.displ),3);
% Fy = smooth(arrayfun(@(x) x.commoncyls_F_y_pN, force.Force.displ),3);

% This is the option to use when NOT using NICOLA data
cod_x = arrayfun(@(x) x.cod_cpara, force.Force.displ);
cod_y = arrayfun(@(x) x.cod_cperp, force.Force.displ);

% Fx = smooth(arrayfun(@(x) x.commoncyls_F_para_pN, force.Force.displ),3);
% Fy = smooth(arrayfun(@(x) x.commoncyls_F_perp_pN, force.Force.displ),3);


% the centre of drag position is translated to the origin and converted too microns.
cod_x = (cod_x-mean(cod_x))*pix2mic;
cod_y = (cod_y - mean(cod_y))*pix2mic;

% the time stamp for each clicked frame, this is translated to measure from time zero
timestamp = arrayfun(@(x) x.timestamp_s, force.Force.displ);
timestamp = timestamp-timestamp(1);

%% STEP 1a: The fitting function used assumes there is a complete number of cycles. The code below plots the all points with the ones included in the fit range a different colour
% It is used to find the appropriate interval of points to fit the crescent.

% fitrange is the vector of the point indices used in the fit. The if statement prevents an error occurring if the fitrange has already been defined in another fit but there are more points than this case can handle
if(length(fitrange)>length(cod_x))
    fitrange = 1:length(cod_x);
end

% plotting the ellipse of the points
figure(1)
subplot(2,1,1)
plot(cod_x(fitrange), cod_y(fitrange))
hold on
plot(cod_x([fitrange(end) setdiff(1:length(cod_x),fitrange)]),...
    cod_y([fitrange(end) setdiff(1:length(cod_x),fitrange)] ))
hold off

% plotting as a sequence of points, where it is easier to determine the index of each point
subplot(2,1,2)
plot(cod_x)
hold on
plot(cod_y)
scatter(fitrange,cod_x(fitrange))
scatter(fitrange, cod_y(fitrange))
hold off


%% STEP 2: Fitting second order fourier series to the x and y coordinates
% OPTION 1: This is in the simple case that period is similar for both. Otherwise use
% option 2. Option 2 can't be compressed but can be useful when assessing the fitrange


% Defining the fit range as determined by eye in the previous section.
% % For the Nicola data
fitrange = 1:67; 
% An approximation for the period of the stroke. In this case the indices cover two cycles so hard the fitrange is used to define the period
Tapprox = timestamp(fitrange(32))-timestamp(fitrange(1));

% the zero of the time stamp is defined as the first point used in the fit
timestamp = timestamp- timestamp(fitrange(1));

% fitting a second order Fourier series to the x and y data
xfit = fit(timestamp(fitrange)', cod_x(fitrange)', 'fourier2',...
    'Lower',[-Inf -Inf -Inf -Inf -Inf 1.2*pi/Tapprox],'Upper', [Inf Inf Inf Inf Inf 3*pi/Tapprox]);
yfit = fit(timestamp(fitrange)', cod_y(fitrange)', 'fourier2',...
    'Lower',[-Inf -Inf -Inf -Inf -Inf 1.2*pi/Tapprox],'Upper', [Inf Inf Inf Inf Inf 3*pi/Tapprox]);

% tcycles = linspace(0,4*pi/xfit.w,400);
% plotting the fit with the original data to check the fit isn't doing anything crazy.
figure(1)
% % OPTION 1:
subplot(1,2,1)
plot(xfit, timestamp, cod_x)

xlabel('index')
ylabel('x pos')

subplot(1,2,2)
plot(yfit, timestamp, cod_y)

xlabel('index')
ylabel('y pos')
% OPTION 2. Not valid for the compression process and may no longer align with the variable names. Kept for completeness.
% subplot(1,2,1)
% plot(1:length(cod_x), cod_x)
% hold on
% plot(tcycles, xmdl(parameters,tcycles))
% hold off
% 
% xlabel('index')
% ylabel('x pos')
% 
% subplot(1,2,2)
% plot(1:length(cod_y), cod_y)
% hold on
% plot(tcycles, ymdl(parameters,tcycles))
% hold off
% 
% xlabel('index')
% ylabel('y pos')

%% STEP 2a: The final model needs to have the same period.

% define the angular frequency as the average of the two fits
wmean = (xfit.w+yfit.w)/2;

phaseshift = timestamp(1); % in the event that the original timepoint wasn't 0. Kept in case a phase shift needs to be introduced
% Creating models with the average period instead. A new name is used to preserve the confidence intervals and other information that is lost when overwriting a fit object.
xfitmean = xfit;
xfitmean.w = wmean;
yfitmean = yfit;
yfitmean.w = wmean;
tcycles = linspace(0,(4*pi)/wmean,400)+phaseshift;
xtraj = xfitmean(tcycles);
ytraj = yfitmean(tcycles);

% % OPTIONS 2: already have the same period from fitting
% xtraj = xmdl(parameters,tcycles);
% ytraj = ymdl(parameters,tcycles);

% Plotting the new fit with the average cycle with the original data as a sanity check.
figure(2)
plot(xtraj, ytraj,'k','LineWidth',1.5)
hold on
scatter(cod_x, cod_y,25,'b','filled')
plot(cod_x,cod_y,'r')
hold off

%% STEP 3: Defining the end points of the trajectory

% Define the endpoints using the minimum speed. midrd_w is the value and loc_rd the indices of the endpoints 
[minrd_w, loc_rd] = findpeaks(-hypot(xtraj(1:end-1)-xtraj(2:end), ytraj(1:end-1)-ytraj(2:end)));

% Quick plot to check where the endpoints have been measured
%figure(3)
% This uses the distance from the origin
%scatter(xtraj,ytraj,[], hypot(xtraj,ytraj),'filled')
%hold on
%scatter(xtraj(loc_rd), ytraj(loc_rd),'kx')
%hold off



% This section selects the indices for the inside and outside sections of trajectory. Originally I attempted to create if statements to automatically select the indices, but it was faster to select them by hand. This section will require editing if applied to new trajectories
if (length(loc_rd) ==4) 
    inner_ind = linspace(tcycles(loc_rd(1)), tcycles(loc_rd(2)),150);
    outer_ind = linspace(tcycles(loc_rd(2)), tcycles(loc_rd(3)),150);
elseif (length(loc_rd) >=6)
%       THIS IS FOR THE ATNO CASE
    [maxval, maxind] = max(abs(diff(xtraj(loc_rd(1:4)))));

% 
    inner_ind = linspace(tcycles(loc_rd(maxind)), ...
        tcycles(loc_rd(maxind+2)) ,150);
    outer_ind = linspace(tcycles(loc_rd(maxind+2)),...
        tcycles(loc_rd(maxind+4)),150); 
%     %     
else
    fprintf('ERROR: too many minimums return to by hand')
end

% Actually considering the trajectory in two parts, inside and outside and
% before or after twist
cd /scratch
% the coordinates for both cases
xinner = xfit.a0 + xfit.a1*cos(inner_ind*wmean) + xfit.a2*cos(2*wmean*inner_ind) + ...
    xfit.b1*sin(inner_ind*wmean) + xfit.b2*sin(2*inner_ind*wmean);  % the x values for the inside of the trajectory
xouter = xfit.a0 + xfit.a1*cos(outer_ind*wmean) + xfit.a2*cos(2*wmean*outer_ind) + ...
    xfit.b1*sin(outer_ind*wmean) + xfit.b2*sin(2*outer_ind*wmean);	% the x values for the outside

% the same but for the y
yinner =  yfit.a0 + yfit.a1*cos(inner_ind*wmean) + yfit.a2*cos(2*wmean*inner_ind) + ...
    yfit.b1*sin(inner_ind*wmean) + yfit.b2*sin(2*inner_ind*wmean);
youter = yfit.a0 + yfit.a1*cos(outer_ind*wmean) + yfit.a2*cos(2*wmean*outer_ind) + ...
    yfit.b1*sin(outer_ind*wmean) + yfit.b2*sin(2*outer_ind*wmean);

% another plot to check what's going on
figure(4)
scatter(xtraj(1:end-1),ytraj(1:end-1),[], ...
    hypot(xtraj(1:end-1)-xtraj(2:end), ytraj(1:end-1)-ytraj(2:end)),'filled')
hold on
plot(xinner,yinner,'b','Linewidth',2.5)
plot(xouter,youter,'r','Linewidth',0.5)
plot(0.5*(xinner+fliplr(xouter)),0.5*(yinner+fliplr(youter)),'Color',[0.9100 0.4100 0.1700],'Linewidth',1.5)
scatter(xtraj(loc_rd), ytraj(loc_rd),'kx')

hold off

box on
axis equal

xlabel('x pos')
ylabel('y pos')
%% STEP 4: Calculating lines that are perpendicular to the central curve at each point
% This requires a gradient and a point for each line.

% finding the midpoint between the inside and outside trajectories for  x (centralx) and y (centraly) 
centralx = 0.5*(xinner+fliplr(xouter));
centraly = 0.5*(yinner+fliplr(youter));

% Defining the gradient between each point, to then define the gradient of the perpendicular line
mcurve = (centraly(1:end-1)-centraly(2:end))./(centralx(1:end-1)-centralx(2:end));
mperp = -1./mcurve; % perpendicular gradient

% The points near the curve can be problematic, to avoid this, they have been set equal to the the gradient slightly further in
mperp(1:4) = mperp(5);
mperp(end-7:end) = mperp(end-8);
% Doing the same with the m curvature for later when classifying the sign
% of the points
mcurve(1:4) = mcurve(5);
mcurve(end-7:end) = mcurve(end-8);

% Defining the midpoint of each sequential pair of points which had defined the gradient
cperpx = 0.5*(centralx(1:end-1)+centralx(2:end));
cperpy = 0.5*(centraly(1:end-1)+centraly(2:end));

% The midpoint and gradient then calculate the y intercept for each perpendicular ray, i.e. y = perp * x + cperp
cperp = cperpy - mperp.*cperpx;

% A quick plot to look at the trajectory, central axis and the perpendicular rays from the central line.
figure(5)
plot(xtraj,ytraj,'r')
axis equal

for crl = 1:length(mperp)
    refline(mperp(crl), cperp(crl))
    axis([-5 5 -3 4.2])
%     axis([-2.7 2.5 -1.7 2.2])
    drawnow
%      pause
end
hold on
plot(xtraj,ytraj,'r')
plot(0.5*(xinner+fliplr(xouter)),0.5*(yinner+fliplr(youter)),'Color',[0.9100 0.4100 0.1700],'Linewidth',1.5)
hold off
% axis([-2.7 2.5 -1.7 2.2])
axis([-5 5 -3 4.2])

%% STEP 6: Using the perpendicular rays to project the data onto the central axis. Now projecting onto the curve using these lines

% Defining the length of the central line: centralL
centralL = cumsum(hypot(centralx(1:end-1)-centralx(2:end), centraly(1:end-1) - centraly(2:end)));
centralL = [0 centralL];

% preallocating the memory for: distance along the central line (dcentral), the distance from the data point to the central line to use a check (dmincheck), and the index of the point considered (dind)
dcentral = zeros(length(cod_x),2);
dmincheck = zeros(size(cod_x));
dind = zeros(size(cod_x));

% The projection for each point
for ccp = 1:length(cod_x)
    % STEP 1: Calculate the distance of the point [cod_x(i), cod_y(i)] from each line
    % This is just the formula from the wikipedia page
    dlines = abs(mperp*cod_x(ccp) - cod_y(ccp) + cperp)./sqrt(mperp.^2+1);

    % as a safety check, any points more than half the centralL from the
    % test points have an artifically high distance to exclude them
    distCL = hypot(cperpx-cod_x(ccp), cperpy-cod_y(ccp));
    dlines(distCL>centralL(end)/2.5) = 20;
    
    % STEP 2: find the minimum distance and save it as a check and use the
    % index
    [dmincheck(ccp), minind] = min(dlines);
    
    % STEP 3: The projected point is the central x and central y points
    % associated with the index of the closest line
    dcentral(ccp, :) = [centralx(minind), centraly(minind)];
    dind(ccp) = minind;

end


%% Plotting the projected points from a given point ontop of the trajctory

figure(6)
scatter(xtraj(1:end-1),ytraj(1:end-1),25, ...
    hypot(xtraj(1:end-1)-xtraj(2:end), ytraj(1:end-1)-ytraj(2:end)),'filled')
hold on
plot(centralx, centraly, 'Color',[0.9100 0.4100 0.1700],'Linewidth',1.5)
scatter(cod_x, cod_y,[],'b','filled')
line([cod_x', dcentral(:,1)]', [cod_y', dcentral(:,2)]','Color','k')
% scatter(dcentral(:,1), dcentral(:,2),[],'rx')
hold off

box on
axis equal

xlabel('x pos')
ylabel('y pos') 

%% STEP 6a: Now that the points are mapped to the central line, the distance from the trap
% is evaluated as the distance along the arclength

centralL = cumsum(hypot(centralx(1:end-1)-centralx(2:end), centraly(1:end-1) - centraly(2:end)));
centralL = [0 centralL];

% Then the index of of each point is known from the projection section
% This is used to measure the distance along the central axis for each data point centralP
centralP = centralL(dind);

%% STEP 7: Categorising the forces as positive or negative depending ig the centralP
% is increasing or decreasing
% The sign of the force is defined using the change in position of the centre of drag and the orientation of the force vector

% The first bit defines Fsign which is the sign of the force defined using the change in the position of the centre of drag to set the sign
if (centralx(1)>0)
    Fsign = smooth(-sign(diff(centralP,1)),6);
else
    Fsign = smooth(sign(diff(centralP,1)),6);
end

% make the vector the same length by doubling the end value. Can't use circular wrap because a complete set of cycles haven't always been tracked.
Fsign = [Fsign; Fsign(end)];

% setting the zero values to the direction of the previous value. This
% should only be a problem if I use an even smoothing value
Fsign(Fsign==0) = Fsign(find(Fsign==0)+1);

% Then just taking the sign, removing the effect of smoothing.
Fsign = sign(Fsign);

% VERSION 1:
% % Then note the cases where the sign is different if the angle of the force
% % is below 45 degrees and the sign is opposite to whats above
% FsignAngle = (atan(Fy./Fx)<1).*sign(Fx);
% Ferr = find(Fsign.*FsignAngle==-1);

% VERSION 2: There is an error if the direction of the force vector
% measured relative to the gradient of the curve indicates another sign
% than what was set by the change in the position
Frelangle = atan2(Fy, Fx) - atan(mcurve(dind))';	% The angle of the force measured relative to the angle of the central curve at the projection point

% the sign/direction of the force assigned using the force vector. Angles in the first and fourth quadrants give +1, angles in the 2nd and 3rd give -1
FsignAngle = (abs((angle(exp(1i*Frelangle))))<pi/2) -(abs((angle(exp(1i*Frelangle))))>pi/2);

% In the event the gradient of the central curve is close to vertical, the F sign values is unreliable and instead is set equal to the Fsineangle
Fsign((abs(atan(mcurve(dind)))>pi/2.5) & (abs(atan(Fy./Fx))<pi/2.5)') =...
    FsignAngle((abs(atan(mcurve(dind)))>pi/2.5) & (abs(atan(Fy./Fx))<pi/2.5)');

% In all other cases, if Fsign and Fsignangle indicate different directions, the point is marked as a potential error and can be excluded from the analysis later
Ferr = find(Fsign.*FsignAngle == -1); % indices of the points that are inconsistent about direction using the two methods.

%% Another quick sanity check plot. Then plot the assigned force value over the trajectory
% The colour will distinguish between + and - cases.

figure(7)
scatter(xtraj(1:end-1),ytraj(1:end-1),25, ...
    hypot(xtraj(1:end-1)-xtraj(2:end), ytraj(1:end-1)-ytraj(2:end)),'filled')
hold on
plot(centralx, centraly, 'Color',[0.9100 0.4100 0.1700],'Linewidth',1.5)
scatter(cod_x(Fsign<0), cod_y(Fsign<0),[],'b','filled')
scatter(cod_x(Fsign>0), cod_y(Fsign>0),[],'r','filled')
set(gca,'ColorOrderIndex',1)
quiver(cod_x(Fsign<0), cod_y(Fsign<0), Fx(Fsign<0)', Fy(Fsign<0)')
quiver(cod_x(Fsign>0), cod_y(Fsign>0), Fx(Fsign>0)', Fy(Fsign>0)')
% line([cod_x', dcentral(:,1)]', [cod_y', dcentral(:,2)]','Color','k')
% scatter(dcentral(:,1), dcentral(:,2),[],'rx')
plot(cod_x(Ferr), cod_y(Ferr),'*w')

plot(cod_x(Ferr), cod_y(Ferr),'xk')
hold off

box on
axis equal

xlabel('x pos')
ylabel('y pos') 

%% STEP 8: Now the part where I plot the force against the central distance and fit the force polynomial
% remember I will need to reverse the direction of the positive force
% results since the relative distance is always decreasing
% % % % % % THIS CODE HERE IS WHAT WAS ADDED FOR THE VERSION 2
% The amplitude is the average range of the two data ranges for positive and negative forces. This prevents problems from extrapolation, which can occur when the two strokes have different domains or can't be tracked well in one phase.
dom_pos = range(centralP(setdiff(find(Fsign>0),Ferr)));
dom_neg = range( centralP(setdiff(find(Fsign<0),Ferr)));
amp = 0.5*(dom_pos + dom_neg); % the amplitude of the potential, i.e. the distance eventually travelled by the rower within the traps.

% It is also useful to define the minimum and maximum value for both domains to be used later.
mm_pos = minmax(centralP(setdiff(find(Fsign>0),Ferr)));
mm_neg = minmax(centralP(setdiff(find(Fsign<0),Ferr)));

% Calculating the distance along the central axis but this time the trap is placed at zero and the distance should always be decreasing towards the trap. Consistent with the usual implementation of the traps. The points are centred with respect to the amplitude, i.e. if the amp is larger than the range there will be equal amounts of extrapolation at both sides of the trap, or if the amp is smaller an equal distance will overshoot both sides of the profile.
% If the first value of cetnralx is positive, then measuring from the wrong
% end and need to flip everything. In this case the error points have been excluded
xvalneg =  centralP(setdiff(find(Fsign<0),Ferr))-mm_neg(1)+0.5*(amp-dom_neg);
xvalpos =mm_pos(2) - 0.5*(dom_pos-amp) - centralP(setdiff(find(Fsign>0),Ferr));

% the same but the error points are included
xvalposWerr =  mm_pos(2) - 0.5*(dom_pos-amp) -centralP(Fsign>0);
xvalnegWerr = centralP(Fsign<0)-mm_neg(1)+0.5*(amp-dom_neg);

% Just the error points results
xvalposERR = mm_pos(2) - 0.5*(dom_pos-amp) -centralP(Ferr);
xvalnegERR = centralP(Ferr)-mm_neg(1)+0.5*(amp-dom_neg);


% Fitting polynomials to force, the error points have not been included.
% At this stage not including the error points in the fit
fit_pos = polyfit(xvalpos, ...
    hypot(Fx(setdiff(find(Fsign>0),Ferr)), Fy(setdiff(find(Fsign>0),Ferr)))',3);
fit_neg = polyfit(xvalneg, ...
    -hypot(Fx(setdiff(find(Fsign<0),Ferr)), Fy(setdiff(find(Fsign<0),Ferr)))',3);

% A plot of the fit and the points. The colour indicates positive or negative. A black and white star is added to each point that has been excluded because of error. The error points have been plotted in both the positive and negative sides, as it is unclear which is the correct sign for the force
figure(8)
scatter(xvalnegWerr,-hypot(Fx(Fsign<0), Fy(Fsign<0)),[],'b','filled')
hold on
scatter(xvalposWerr, hypot(Fx(Fsign>0), Fy(Fsign>0)),[],'r','filled')
plot(centralL, polyval(fit_neg, centralL))
plot(centralL, polyval(fit_pos, centralL))
plot(xvalnegERR, -hypot(Fx(Ferr), Fy(Ferr)),'*w')
plot(xvalnegERR, -hypot(Fx(Ferr), Fy(Ferr)),'xk')
plot(xvalposERR, hypot(Fx(Ferr), Fy(Ferr)),'*w')
plot(xvalposERR, hypot(Fx(Ferr), Fy(Ferr)),'xk')
plot(xvalposERR, hypot(Fx(Ferr), Fy(Ferr)),'*w')
plot(xvalposERR, hypot(Fx(Ferr), Fy(Ferr)),'xk')
plot(xvalnegERR, -hypot(Fx(Ferr), Fy(Ferr)),'*w')
plot(xvalnegERR, -hypot(Fx(Ferr), Fy(Ferr)),'xk')
hold off
xlabel('d central')
ylabel('force pN')
 % % % % % % END OF THE NEW CODE


%% STEP 9: Creating a summary panel of the different plots that can be saved easily

figure(20)
subplot(2,2,1)
plot(cod_x, cod_y,'k')
hold on
scatter(cod_x,cod_y,[], 1:length(cod_x),'filled')
quiver(cod_x, cod_y, Fx', Fy')
hold off

axis equal
box on
xlabel('x pos')
ylabel('y pos')

subplot(2,2,2)

scatter(xtraj(1:end-1),ytraj(1:end-1),25, ...
    hypot(xtraj(1:end-1)-xtraj(2:end), ytraj(1:end-1)-ytraj(2:end)),'filled')
hold on
plot(centralx, centraly, 'Color',[0.9100 0.4100 0.1700],'Linewidth',1.5)
scatter(cod_x, cod_y,[],'b','filled')
line([cod_x', dcentral(:,1)]', [cod_y', dcentral(:,2)]','Color','k')
% scatter(dcentral(:,1), dcentral(:,2),[],'rx')
hold off

box on
axis equal

xlabel('x pos')
ylabel('y pos')

subplot(2,2,3)
scatter(xtraj(1:end-1),ytraj(1:end-1),25, ...
    hypot(xtraj(1:end-1)-xtraj(2:end), ytraj(1:end-1)-ytraj(2:end)),'filled')
hold on
plot(centralx, centraly, 'Color',[0.9100 0.4100 0.1700],'Linewidth',1.5)
scatter(cod_x(Fsign<0), cod_y(Fsign<0),[],'b','filled')
scatter(cod_x(Fsign>0), cod_y(Fsign>0),[],'r','filled')
set(gca,'ColorOrderIndex',1)
quiver(cod_x(Fsign<0), cod_y(Fsign<0), Fx(Fsign<0)', Fy(Fsign<0)')
quiver(cod_x(Fsign>0), cod_y(Fsign>0), Fx(Fsign>0)', Fy(Fsign>0)')
% line([cod_x', dcentral(:,1)]', [cod_y', dcentral(:,2)]','Color','k')
% scatter(dcentral(:,1), dcentral(:,2),[],'rx')
plot(cod_x(Ferr), cod_y(Ferr),'*w')

plot(cod_x(Ferr), cod_y(Ferr),'xk')
hold off

box on
axis equal

xlabel('x pos')
ylabel('y pos') 

subplot(2,2,4)
scatter(xvalnegWerr,-hypot(Fx(Fsign<0), Fy(Fsign<0)),[],'b','filled')
hold on
scatter(xvalposWerr, hypot(Fx(Fsign>0), Fy(Fsign>0)),[],'r','filled')
plot(linspace(0,amp), polyval(fit_neg, linspace(0,amp)))
plot(linspace(0,amp), polyval(fit_pos, linspace(0,amp)))
plot(xvalnegERR, -hypot(Fx(Ferr), Fy(Ferr)),'*w')
plot(xvalnegERR, -hypot(Fx(Ferr), Fy(Ferr)),'xk')
plot(xvalposERR, hypot(Fx(Ferr), Fy(Ferr)),'*w')
plot(xvalposERR, hypot(Fx(Ferr), Fy(Ferr)),'xk')
plot(xvalposERR, hypot(Fx(Ferr), Fy(Ferr)),'*w')
plot(xvalposERR, hypot(Fx(Ferr), Fy(Ferr)),'xk')
plot(xvalnegERR, -hypot(Fx(Ferr), Fy(Ferr)),'*w')
plot(xvalnegERR, -hypot(Fx(Ferr), Fy(Ferr)),'xk')
hold off
xlabel('d central')
ylabel('force pN')
box on
xlabel('d central')
ylabel('force pN')

% Calculating the area under each curve, which can be used to assess how wide each trajectory is compared with the amplitude.

trajArea = 0.5*diff(xinner)*abs(yinner(1:end-1)+yinner(2:end))' - ...
    0.5*diff(xouter)*abs(youter(1:end-1)+youter(2:end))'

