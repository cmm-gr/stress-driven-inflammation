%
%**** This script plots Figure 5 (Example 3(c)) in:
%
%	  Latorre M, Spronck B, Humphrey JD (2020) Complementary Roles of
%	  Mechanotransduction and Inflammation in Vascular Homeostasis.
%	  Proceedings of the Royal Society A
%
%**** marcos.latorre@yale.edu

clearvars
close all
%
kinf   = 2/7;						% rate parameter for inflammatory response (days^-1)
Kinf   = 2.5;						% gain parameter for adaptive inflammatory response (-)
svmax  = 170;						% stress threshold for maladaptative inflammatory response (kPa)
mudotp = 0.102;						% rate of change of inflammatory cell infiltration + (days^-2)
mudotm = 8e-04;						% rate of change of inflammatory cell infiltration - (days^-2)
%
parInf = [kinf,Kinf,svmax,mudotp,mudotm];
%
days = 224;							% total simulation time (days)
%
SPs = [7 7];						% periods for pressure elevation (days)
KPs = [0 0];						% gains for hypertension-induced increase in active tone (-)
PRs = [1 0];						% flags for passive properties recovery (1 = yes | 0 = no)
%
lines = {'-','-.'};					% line styles
%
fign = 5;							% figure number
%
for plotcase = 1:length(SPs)		% case [1,2] -> PR = [1,0]
	%
	SP = SPs(plotcase);
	KP = KPs(plotcase);
	PR = PRs(plotcase);
	%
	line = char(lines(plotcase-floor((plotcase-1)/length(lines))*length(lines)));
	%
	BiThinStressInflam('DTA_pas.mat',parInf,days,SP,KP,PR,fign,line) % compute/plot case (passive data)
	%
end
%
subplot(341)
hl = legend('Adaptive','Maladaptive');
set(hl,'Location','NorthEast','Box','Off')
%