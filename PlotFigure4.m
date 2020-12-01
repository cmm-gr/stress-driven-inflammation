%
%**** This script plots Figure 4 (Example 3(b)) in:
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
svmax  = 150.25;					% stress threshold for maladaptative inflammatory response (kPa)
mudotp = 0.102;						% rate of change of inflammatory cell infiltration + (days^-2)
mudotm = 0;							% rate of change of inflammatory cell infiltration - (days^-2)
%
parInf = [kinf,Kinf,svmax,mudotp,mudotm];
%
days = 28;							% total simulation time (days)
%
SPs = [7  7   7  7];				% periods for pressure elevation (days)
KPs = [0 1.5 3.5 4];				% gains for hypertension-induced increase in active tone (-)
PRs = [1  1   1  1];				% flags for passive properties recovery (1 = yes | 0 = no)
%
lines = {'-','--','-.',':'};		% line styles
%
fign = 4;							% figure number
%
for plotcase = 1:length(SPs)		% case [1,2,3,4] -> KP = [0,1.5,3.5,4]
	%
	SP = SPs(plotcase);
	KP = KPs(plotcase);
	PR = PRs(plotcase);
	%
	line = char(lines(plotcase-floor((plotcase-1)/length(lines))*length(lines)));
	%
	BiThinStressInflam('DTA_pas_act.mat',parInf,days,SP,KP,PR,fign,line) % compute/plot case (passive + active data)
	%
end
%
subplot(341)
hl = legend(['K_P = ',num2str(KPs(1))],['K_P = ',num2str(KPs(2))],['K_P = ',num2str(KPs(3))],['K_P = ',num2str(KPs(4))]);
set(hl,'Location','East','Box','Off')
%