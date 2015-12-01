% plotQuenched.m

% Created by Krista Longnecker, klongnecker@whoi.edu, 9/19/2013
% Revised by Jamie Collins, james.r.collins@aya.yale.edu, on 10/31/2013 to
% use a user-defined function for parameter optimization, rather than
% MATLAB's interactive fitting tool
% Revised by J.R.C., 12/1/2015, to make suitable for upload to GitHub

% Purpose: Obtains best-fit parameters for a cubic function
% quenchcurve_cubicfitfunc.m, for conversion of H-number and cpm to dpm,
% based on liquid scintillation counter efficiency. The function
% coefficients are written to an output file Quench_curve_fit_params.txt,
% which are then reimported by a companion script, riBPdata.m, for use in
% rate calculations. plotQuenched.m allows the user to control the manner
% in which dpm values are calculated, rather than rely on the LSC's own
% internal calculation routine.

% For details on the 3H-leucine microcentrifige method, see:
%
%   Kirchman, D., E. Knees, and R. Hodson. 1985. Leucine Incorporation and
%   Its Potential as a Measure of Protein-Synthesis by Bacteria in Natural
%   Aquatic Systems. Appl. Environ. Microbiol. 49: 599-607.
%
%   and
%
%   Kirchman, D. 2001. Measuring bacterial biomass production and growth
%   rates from leucine incorporation in natural aquatic environments.
%   In Methods in Microbiology. J. H. Paul, editor. Academic Press. 227-
%   237.

% Dependencies/assumptions:
%
%  1. The m-file quenchcurve_cubicfitfunc.m, the function for which
%  parameters are to be obtained
%
%  2. nlleasqr.m, for computing non-linear least-squares regression, and
%  the sub-dependency dfdp.m, which computes numerical partial derviatives
%
%  3. A file containing quench standard data from your liquid scintillation
%  counter (example data are contained in a file
%  "Quench_standards_20130924.xlsx," which is included in the repo where
%  this script is maintained.
%
%  4. User specification of the activity of the quench standard series
%  used to generate the quench standard data, specified below as sDPM

%% User specify file paths, other required inputs

% location of quench standard data
[Quench_num Quench_txt Quench_raw ] =xlsread('Quench_standards_20130924.xlsx');

sDPM = 258000; % DPM of the quenched standard series; user should update appropriately

%% Meat and potatoes

Hnumber = Quench_num(:,8); % extract H-numbers
CPM = num(:,9); % extract CPMs

for a = 1
    %calculate the efficiency for the standard curve:
    ce = CPM/sDPM;
    plot(Hnumber,ce,'o','markeredgecolor','r')
    hold on
end
clear a j

% pass things off to cubic function of form p(1)*x^3 + p(2)*x^2 + p(3)*x + p(4)

% some initial guesses for nlleasqr

pin = [1e-9 -1e-6 -2e-2 1];

% now, use nlleasqr to obtain best-fit values of our parameters

[f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2]=nlleasqr(Hnumber,ce,pin,'quenchcurve_cubicfitfunc');
yf = quenchcurve_cubicfitfunc(Hnumber,p);
plot(Hnumber,ce,'o',Hnumber,yf);
xlabel('H number');
ylabel('CPM/sDPM');

% return values of our parameters
disp([char(10) 'Values of best-fit parameters for a cubic function of form' char(10)...
    'y = p(1)*x^3 + p(2)*x^2 + p(3)*x + p(4)' char(10) char(10)...
    '   p(1)         p(2)         p(3)            p(4)' char(10)])
format SHORTG;
disp(p');

%write the values to a file

dlmwrite('Quench_curve_fit_params.txt',p);
