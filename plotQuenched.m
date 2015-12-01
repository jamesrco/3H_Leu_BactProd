%plotQuenched_forJamie
%KL 9/19/2013
%revised JC 10/31/2013:
%  presumes you have the least-squares fitting routine nlleasqr in an accessible place
%  for MATLAB, and a function to be fit (I am using a cubic function)

[num txt raw ] =xlsread('Quench_standards_20130924.xlsx');
Hnumber = num(:,8);
CPM = num(:,9);
sDPM = 258000; %DPM of the quenched standard series, ref date 6/7/2011

for a = 1
    %calculate the efficiency for the standard curve:
    ce = CPM/sDPM;
    plot(Hnumber,ce,'o','markeredgecolor','r')
    hold on
end
clear a j

%%Jamie: at this point, I cheated and used the interactive fitting tool in
%%MATLAB. In the figure that comes up, Tools/Basic Fitting. I decided that
%%a cubic fit was best, and the tool gave me the formula and coefficients.
%%That is what gets used in the riBPdata m-file. 

%%Below is the type of output you get from the Basic Fitting:
% 
% y = p1*x^3 + p2*x^2 +
%       p3*x + p4 
% 
% Coefficients:
%   p1 = 8.1248e-009
%   p2 = -1.5438e-006
%   p3 = -0.0021008
%   p4 = 0.59562

% good advice, but going to pass things off to my own cubic function of
% form p(1)*x^3 + p(2)*x^2 + p(3)*x + p(4)

% some initial guesses

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

dlmwrite('/Users/jrcollins/Dropbox/KN207-3 NA VICE/Data/BP workups/Quench_curve_fit_params.txt',p);
