% MATLAB code for EAS 4200, Homework #1, Problem #2, Fall 2024
% Jeffrey L. Kauffman  <JLKauffman@ucf.edu>  31 AUG 2024

% Do not try to run this file immediately -- there are many partial
% equations that you need to complete.  Do read the comments carefully,
% and follow along with the homework assignment.

% Given Values and Constants 
m = 1800;
S = 20;
cbar = 1.65; % keep going with numbers given in problem statement
CLalpha = 2*pi;
CLmax_pos = 1.88;%
CLmax_neg = -1.63;%
nmax_pos = 5.5;%
nmax_neg = -4;%

VC = 125;%
VD = 190;%
% should also define density (sea level) and gravity -- will need these
rho = 1.227;%
g = 9.81;%

% Calculate the critical flight speed values
VA = sqrt(2*nmax_pos*m*g/rho/S/CLmax_pos);
VG = sqrt(2*nmax_neg*m*g/rho/S/CLmax_neg);% same equation as VA, but with the max negative n, CLmax values
VS_pos = sqrt(2*m*g/rho/S/CLmax_pos);% same equation as VA, but with n = 1, positive CLmax
VS_neg = sqrt(-1*2*m*g/rho/S/CLmax_neg);% I'll let you figure out this one

% Part 2.2
% Create an array of velocity values to use in plotting
Vpos = linspace(0,VA,50);       % 50 points linearly spaced from 0 to VA
npos = ((Vpos.^2)/(m*g))*((rho*S*CLmax_pos)/2);% insert equation for the positive load factor n from Problem 1
       % Tip: it should have a V^2 term -- type it in Matlab as Vpos.^2
       % Reason: Vpos is a vector, so Vpos^2 will create an error -- 
       % a period before an operator performs "element-wise" operations
       % e.g., [1 2 3 4].^2 = [1 4 9 16], but [1 2 3 4]^2 = error

%Now create the points that will actually make up the line...
% Start by creating an array of velocity values to use in plotting
%Vpos = linspace(0,VA,50); % 50 points linearly spaced from 0 to VA
%npos = % insert equation for the positive load factor n from Problem 1
        % Tip: it should have a V^2 term -- type it in Matlab as Vpos.^2
        % Reason: Vpos is a vector, so Vpos^2 will create an error --
        % a period before an operator performs "element-wise" operations
        % e.g., [1 2 3 4].^2 = [1 4 9 16], but [1 2 3 4]^2 = error

Vneg = linspace(0,VG,50);
nneg = ((Vneg.^2)/(m*g))*((rho*S*CLmax_neg)/2);% insert equation; again use the .^ for Vneg.^2

% Part 2.3
% Create a new figure and turn on hold -- this way successive plot calls
% will just add the curve to the graph rather than clearing out the
% graph and starting over
figure;  hold on;


% Actually plot the OA stall line
plot(Vpos,npos,'b--')  % 'b--' means blue dashed line
% and now the 0G stall line
plot(Vneg,nneg,'b--')
% Add the AD line
plot([VA VD],[nmax_pos nmax_pos],'b--')  % syntax: plot([x1 x2],[y1 y2])
% and the GF line
plot([VG VC],[nmax_neg nmax_neg], 'b--')% fill in the rest
% Plot the right side of the envelope
plot([VC VD],[nmax_neg 0], 'b--')% fill in the rest
% and the left side (easiest to do this in two calls)
plot([VS_pos VS_pos],[1 0], 'b')% fill in the rest
plot([VS_neg VS_neg],[-1 0], 'b')% fill in the rest
plot([VD VD],[0 nmax_pos], 'b--')

plot([VC VC],[nmax_neg 0], 'b')

plot([0 200],[0 0], 'k')

% Add some labels to the plot
xlabel('Velocity V (m/s)');  ylabel('Load Factor n');


% Part 2.4
mug = (2*m)/(rho*S*cbar*(CLalpha));% equation for airplane mass ratio
Kg = (0.88*mug)/(5.3+mug);% equation for gust alleviation factor

% Gust of 7.62 m/s (25 ft/s) at VD
Ugust = 7.62;
% Create an array of velocities from 0 to VD
VDgust = linspace(0,VD,50);
% Calculate an array of values along the positive load factor line
nDgust_pos = 1 + rho*VDgust*Kg*Ugust*S*CLalpha/2/m/g;
% Calculate an array of values along the negative load factor line
nDgust_neg = 1 - rho*VDgust*Kg*Ugust*S*CLalpha/2/m/g;% insert equation

% Same thing, but for the stronger gust of 15.24 m/s (50 ft/s) at VC
Ugust = 15.24;%
% Create an array of velocities from 0 to VC
VCgust = linspace(0,VC,50);%
% Calculate an array of values along the positive load factor line
nCgust_pos = 1 + rho*VCgust*Kg*Ugust*S*CLalpha/2/m/g;%
% Calculate an array of values along the negative load factor line
nCgust_neg = 1 - rho*VCgust*Kg*Ugust*S*CLalpha/2/m/g;%


% Part 2.5
% Now plot the gust load factors
plot(VDgust,nDgust_pos,'r-.')  % 'r-.' is for red dash-dot lines
plot(VDgust,nDgust_neg,'r-.')
% and more for VC curves...
plot(VCgust,nCgust_pos,'r-.')%
plot(VCgust,nCgust_neg,'r-.')%
% Connect the points where +V_D and -V_D lines reach V_D... it might be
% easier to calculate the maximum values of the +V_C and +V_D lines (and
% similar for -V_C and -V_D
maxnDgust_pos = max(nDgust_pos);
maxnDgust_neg = min(nDgust_neg);
maxnCgust_pos = max(nCgust_pos);
maxnCgust_neg = min(nCgust_neg);
% similar for _neg and then for nCgust
plot([VD VD], [maxnDgust_pos maxnDgust_neg],'r-.')
plot([VC VD], [maxnCgust_pos maxnDgust_pos],'r-.')
% Now connect those last points, from where the +V_C line reaches V_C to
% where the +V_D line reaches V_D
plot([VC VD], [maxnCgust_neg maxnDgust_neg],'r-.')

% General cleanup; the envelope should fit nicely in 0<v<200 and -5<n<7:
axis([0 200 -5 7])
