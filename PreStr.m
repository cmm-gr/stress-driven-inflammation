function [Pfo,so] = PreStr(PAR,parT0,parA,riio,emco)
%
%**** it computes mechanobiologically equilibrated stresses,
%	  pressure, and axial force at original homeostic state (o)
%
%**** marcos.latorre@yale.edu

%
%** PAR
%
c   = PAR(1);						% c elastin
Get = PAR(2);						% circumferential deposition stretch elastin
Gez = PAR(3);						% axial deposition stretch elastin
Bt  = PAR(4);						% fraction of circumferential collagen within the adventitia
Bz  = PAR(5);						% fraction of axial collagen within the adventitia
alp = PAR(6);						% orientation of diagonal collagen wrt axial direction
%
betaM = [Bz 1-Bz];					% medial betas [bzM 2*bdM]
betaA = [Bt Bz 1-Bt-Bz];			% adventitial betas [btA bzA 2*bdA]
%
%** parT0
%
c1m = parT0(1);						% c1t muscle
c2m = parT0(2);						% c2t muscle
c1c = parT0(3);						% c1t collagen
c2c = parT0(4);						% c2t collagen
Gm  = parT0(5);						% circumferential deposition stretch (combined medial collagen and smc)
Gc  = parT0(6);						% deposition stretch (collagen)
%
%** parA
%
Tmax = parA(1);						% baseline active tone
lM   = parA(2);						% stretch for maximum contractile strength
l0   = parA(3);						% stretch for minimum contractile strength
CB   = parA(4);						% baseline parameter for flow-induced contraction
%
%** emco and riio
%
phiMo = [emco(1:2) emco(3)*betaM];	% local mass fractions of medial [e mt cz 2*cd]
phiAo = [emco(4)   emco(5)*betaA];	% local mass fractions of adventitial [e ct cz 2*cd]
%
rio  = riio(1);						% inner radius
rMAo = riio(2);						% M-A radius
roo  = riio(3);						% outer radius
%
lzoo = 1;							% axial stretch from o
%
hMo = rMAo-rio;						% medial thickness
hAo = roo-rMAo;						% adventitial thickness
%
SMo = pi/lzoo*(rMAo^2-rio^2);		% medial cross-sectional area
SAo = pi/lzoo*(roo^2-rMAo^2);		% adventitial cross-sectional area
%
Ge = [1/Get/Gez 1/Get/Gez Get Get Gez Gez];	% [GerM GerA GetM GetA GezM GezA]
%
stMo = phiMo(1)*c*(Ge(3)^2-Ge(1)^2) + ...	% circ. stress media
	   phiMo(2)*c1m*(Gm^2-1)*exp(c2m*(Gm^2-1)^2)*Gm^2 + ...
	   phiMo(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*sin(alp)^2 + ...
	   phiMo(2)*Tmax*(1-exp(-CB^2))*1*(1-((lM-1)/(lM-l0))^2);
%
stAo = phiAo(1)*c*(Ge(4)^2-Ge(2)^2) + ...	% circ. stress adventitia
	   phiAo(2)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2 + ...
	   phiAo(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*sin(alp)^2;
%
szMo = phiMo(1)*c*(Ge(5)^2-Ge(1)^2) + ...	% axial stress media
	   phiMo(3)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2 + ...
	   phiMo(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*cos(alp)^2;
%
szAo = phiAo(1)*c*(Ge(6)^2-Ge(2)^2) + ...	% axial stress adventitia
	   phiAo(3)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2 + ...
	   phiAo(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*cos(alp)^2;
%
sto = (stMo*hMo+stAo*hAo)/(hMo+hAo);		% mean circ. stress
szo = (szMo*SMo+szAo*SAo)/(SMo+SAo);		% mean axial stress
%
Pfo(1) = (stMo*hMo + stAo*hAo)/rio;			% Pressure at o
Pfo(2) = szMo*SMo + szAo*SAo;				% axial force at o
%
so = [stMo stAo szMo szAo sto szo sto+szo];	% stresses at o
%
end