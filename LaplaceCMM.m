function Res = LaplaceCMM(ltFZ,par,mR,q,kmco,kmc,lTauInv,act,qact,c1m,c2m,c1c,c2c,Gm,Gc)
%
%**** it computes circumferential (med/adv) stretches and axial force at
%	  current G&R time from equilibrium and compatibility equations
%
%**** marcos.latorre@yale.edu

%
%** UNKNOWNS
%
ltM = ltFZ(1);								% medial circum stretch
ltA = ltFZ(2);								% advent circum stretch
fz  = ltFZ(3);								% vessel axial force
%
%** RETRIEVE KNOWN VALUES
%
P      = par(1);							% inner pressure
Eps    = par(2);							% relative change in flow rate
lz     = par(3);							% axial stretch
sgmIo  = par(4);							% tr(sgm) at o
JM     = par(5);							% medial Jacobian
JA     = par(6);							% advent Jacobian
rio    = par(7);							% innner radius at o
hMo    = par(8);							% medial thickness at o
hAo    = par(9);							% advent thickness at o
ce     = par(10);							% elastin shear modulus
Ge     = par(11:13);						% deposition stretch of elastin
rhoReM = par(14);							% ref mass density of medial elastin
rhoReA = par(15);							% ref mass density of advent elastin
rhoRmM = par(16);							% ref mass density of medial smc
alp    = par(17);							% diagonal collagen orientation
Tmax   = par(18);							% baseline active tone
lM     = par(19);							% stretch for maximum contractile strength
l0     = par(20);							% stretch for minimum contractile strength
CB     = par(21);							% baseline parameter for flow-induced contraction
CS     = par(22);							% gain parameter for flow-induced contraction
kact   = par(23);							% rate parameter for active remodeling
rho    = par(24);							% artery mass density
t      = par(25:26);						% hereditary integral limits (passive)
ta	   = par(27:28);						% hereditary integral limits (active)
Ds     = par(29);							% time step increment
%
%** KINEMATICS / GEOMETRY
%
lTauInv(1,:) = [1/ltM	1/ltA	1/lz];		% update current stretches
%
ltTauInvM = lTauInv(:,1);					% medial circum stretches (history)
ltTauInvA = lTauInv(:,2);					% advent circum stretches (history)
lzTauInv  = lTauInv(:,3);					% axial stretches (history)
%
lrM = JM/ltM/lz;							% medial radial stretch
lrA = JA/ltA/lz;							% advent radial stretch
%
FM = [lrM,ltM,lz];							% medial deformation grad
FA = [lrA,ltA,lz];							% advent deformation grad
%
hM = lrM*hMo;								% medial thickness
hA = lrA*hAo;								% advent thickness
ri = (ltM*(2*rio+hMo)-hM)/2;				% inner radius computed from ltM
%
%** ELASTIN
%
FeM = FM.*Ge;								% medial elastin deformation grad
FeA = FA.*Ge;								% advent elastin deformation grad
SeMod = ce*[1,1,1];							% second P-K stresses at constituent level
sgmWeM = rhoReM/(JM*rho)*FeM.*SeMod.*FeM;	% associated Cauchy stress [r,t,z] in media
sgmWeA = rhoReA/(JA*rho)*FeA.*SeMod.*FeA;	% associated Cauchy stress [r,t,z] in advent
pM = sgmWeM(1);								% medial Lagrange multiplier
pA = sgmWeA(1);								% advent Lagrange multiplier
%
%** SMC and COLLAGEN
%
%* circum, diag, and axial stretches (including G's)
%
lftM = ltM*ltTauInvM.*Gm;
lfdM = sqrt((ltM*ltTauInvM*sin(alp)).^2+(lz*lzTauInv*cos(alp)).^2).*Gc;
lftA = ltA*ltTauInvA.*Gc;
lfdA = sqrt((ltA*ltTauInvA*sin(alp)).^2+(lz*lzTauInv*cos(alp)).^2).*Gc;
lfz  = lz*lzTauInv.*Gc;
%
lftM(lftM<1) = 1; lfdM(lfdM<1) = 1; lftA(lftA<1) = 1; lfdA(lfdA<1) = 1; lfz(lfz<1) = 1;
%
lfM = [ltM*ltTauInvM lz*lzTauInv ltM*ltTauInvM lz*lzTauInv];
lfA = [ltA*ltTauInvA lz*lzTauInv ltA*ltTauInvA lz*lzTauInv];
%
%* auxiliary second-P-K-like stresses
%
SfModM(:,1) = c1m.*(lftM.^2-1).*exp(c2m.*(lftM.^2-1).^2).*Gm.^2;
SfModM(:,2) = c1c.*(lfz.^2 -1).*exp(c2c.*(lfz.^2 -1).^2).*Gc.^2;
SfModM(:,3) = c1c.*(lfdM.^2-1).*exp(c2c.*(lfdM.^2-1).^2).*Gc.^2*sin(alp)^2;
SfModM(:,4) = c1c.*(lfdM.^2-1).*exp(c2c.*(lfdM.^2-1).^2).*Gc.^2*cos(alp)^2;
SfModA(:,1) = c1c.*(lftA.^2-1).*exp(c2c.*(lftA.^2-1).^2).*Gc.^2;
SfModA(:,2) = c1c.*(lfz.^2 -1).*exp(c2c.*(lfz.^2 -1).^2).*Gc.^2;
SfModA(:,3) = c1c.*(lfdA.^2-1).*exp(c2c.*(lfdA.^2-1).^2).*Gc.^2*sin(alp)^2;
SfModA(:,4) = c1c.*(lfdA.^2-1).*exp(c2c.*(lfdA.^2-1).^2).*Gc.^2*cos(alp)^2;
%
%* Cauchy stresses at constituent level
%
sgmWfInt = [ 1/(JM*rho)*lfM.*SfModM.*lfM	1/(JA*rho)*lfA.*SfModA.*lfA ];
%
%* intramural stress stimulus, rate parameters, and survival functions
%
Dsgm       = ( P*ri/(hM+hA) + fz/(pi*(hM+hA)*(2*ri+hM+hA)) )/sgmIo - 1;
breaks     = size(q,1);
kmc(1,:)   = kmco*(1+Dsgm^2);
q(2:end,:) = q(2:end,:).*( ones(breaks-1,1) * exp(-Ds*mean(kmc)) );
%
%* Cauchy stresses at mixture level
%
sgmWf  = simpsons(mR.*q.*sgmWfInt,t(1),t(end),[]);	% all of them
sgmWfM = sgmWf(1:4);								% media
sgmWfA = sgmWf(5:8);								% adventitia
%
%* active contribution to circumferential stress
%
aact = simpsons(kact*act.*qact,ta(1),ta(end),[]);	% referential active length
lact = ri/aact;										% active stretch
flact = 1-((lM-lact)/(lM-l0))^2;					% parabolic term (Rachev's model)
if flact < 0, flact = 0; end
%
C = CB - CS*(Eps/(ri/rio)^3-1);						% ratio vaso-constrictors/dilators
%
sgmAmM = Tmax * (rhoRmM/JM)/rho * lact * (1-exp(-C^2)) * flact; % active stress
%
%** MIXTURE
%
sgmtM = sgmWeM(2) + sgmWfM(1) + sgmWfM(3) - pM + sgmAmM;	% medial circum stress
sgmtA = sgmWeA(2) + sgmWfA(1) + sgmWfA(3) - pA;				% advent circum stress
%
sgmzM = sgmWeM(3) + sgmWfM(2) + sgmWfM(4) - pM;				% medial axial stress
sgmzA = sgmWeA(3) + sgmWfA(2) + sgmWfA(4) - pA;				% advent axial stress
%
%** EQUATIONS TO SOLVE IN RESIDUAL FORM
%
Res = [  sgmtM*hM+sgmtA*hA - P*ri								% circum mech equilibrium
		 
		 sgmzM*pi*hM*(2*ri+hM)+sgmzA*pi*hA*(2*ri+2*hM+hA) - fz	%  axial mech equilibrium
		
		(2*ri+2*hM+hA)/(2*rio+2*hMo+hAo) - ltA ];				% geometrical compatibility
end