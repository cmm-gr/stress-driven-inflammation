%
%**** This function computes G&R for a bilayered thin-walled cylindrical
%	  artery with stress-dependent inflammatory effects
%
%**** marcos.latorre@yale.edu

function BiThinStressInflam(filename,parInf,days,SP,KP,PR,fign,line)
%% ***** INITIALIZATION AT ORIGINAL HOMEOSTATIC STATE (o) (s = 0)
%
%** read/load input variables
%
load(filename)						% load parameter values needed by CMM for G&R
%
kinf   = parInf(1);					% rate parameter for inflammatory response (days^-1)
Kinf   = parInf(2);					% gain parameter for adaptive inflammatory response (-)
svmax  = parInf(3);					% stress threshold for maladaptative inflammatory response (kPa)
mudotp = parInf(4);					% rate of change of inflammatory cell infiltration + (days^-2)
mudotm = parInf(5);					% rate of change of inflammatory cell infiltration - (days^-2)
%
%** in-vivo bilayered wall geometry (mm)
%
rio  = riio(1);						% inner radius
rMAo = riio(2);						% media-adventitia (M-A) interface radius
roo  = riio(3);						% outer radius
%
hMo = rMAo-rio;						% medial thickness
hAo = roo-rMAo;						% adventitial thickness
%
%** PAR (inflammation-independent parameters)
%
ce  = PAR(1);						% elastin modulus (kPa)
Get = PAR(2);						% circumferential deposition stretch elastin (-)
Gez = PAR(3);						% axial deposition stretch elastin (-)
Bt  = PAR(4);						% fraction of circumferential collagen (for adventitia) (-)
Bz  = PAR(5);						% fraction of axial collagen (for both media and adventitia) (-)
alp = PAR(6);						% orientation of diagonal collagen wrt axial direction (rad)
%
Ge = [1/Get/Gez Get Gez];			% elastin deposition stretches [r t z] in media and adventitia
%
betaM = [Bz 1-Bz];					% medial betas [BzM 2*BdM]
betaA = [Bt Bz 1-Bt-Bz];			% adventitial betas [BtA BzA 2*BdA]
%
phioM = [emco(1:2) emco(3)*betaM];	% local mass fractions of medial [e mt cz 2*cd]
phioA = [emco(4)   emco(5)*betaA];	% local mass fractions of adventitial [e ct cz 2*cd]
%
%** parT0 (inflammation-dependent parameters at day 0)
%
c1m0 = parT0(1);					% c1 muscle (kPa)
c2m0 = parT0(2);					% c2 muscle (-)
c1c0 = parT0(3);					% c1 collagen (kPa)
c2c0 = parT0(4);					% c2 collagen (-)
Gm0  = parT0(5);					% circumferential deposition stretch (combined medial collagen and smc) (-)
Gc0  = parT0(6);					% deposition stretch (collagen) (-)
%
%** parT4 (inflammation-dependent parameters at day 28 = week 4)
%
c1m4 = parT4(1);					% c1 muscle (kPa)
c2m4 = parT4(2);					% c2 muscle (-)
c1c4 = parT4(3);					% c1 collagen (kPa)
c2c4 = parT4(4);					% c2 collagen (-)
Gm4  = parT4(5);					% circumferential deposition stretch (combined medial collagen and smc) (-)
Gc4  = parT4(6);					% deposition stretch (collagen) (-)
%
%** parA (active stress parameters)
%
Tmaxo = parA(1);					% baseline active tone
lM    = parA(2);					% stretch for maximum contractile strength
l0    = parA(3);					% stretch for minimum contractile strength
CB    = parA(4);					% baseline parameter for flow-induced contraction
CS    = parA(5);					% gain parameter for flow-induced contraction
%
%** G&R parameters
%
Tqoc = 7;							% collagen half-life (hypertensive) (days)
Tqom = (1/etaq)*Tqoc;				% smooth muscle half-life (days)
kmco = [1/Tqom,1/Tqoc,1/Tqoc,1/Tqoc,...	% ko circ m M, ax c M, diag c M, diag c M
		1/Tqoc,1/Tqoc,1/Tqoc,1/Tqoc];	% ko circ c A, ax c A, diag c A, diag c A
%
Kim = etaUm*KicM;					% smooth muscle gains [Kim Ksm Kfm] (-)
Ksm = etaUm*KscM;
Kfm = etaUm*KfcM;
%
KicA = etaUc*KicM;					% adventitial collagen gains [KicA KscA KfcA] (-)
KscA = etaUc*KscM;
KfcA = etaUc*KfcM;
%
Tact = Tqom;						% characteristic time for active remodeling (days)
kact = 1/Tact;						% associated rate parameter
%
Tinf = 1/kinf;						% characteristic time for inflammatory cell clearance (days)
%
%** parameters needed for temporal integrations
%
Ds = min(Tqom,Tqoc)/20;				% G&R time step
%
numHl = 20;							% number of half-lives (Tqo) to consider for heredity integrals
%
domain = numHl*max(Tqom,Tqoc);		% total integration domain
t = (0:Ds:domain)';					% integration times for past histories
breaks = length(t);					% breaks within integration domain
%
domact = numHl*Tact;				% total integration domain for active length
tact = (0:Ds:domact)';				% integration times for past history
breact = length(tact);				% breaks within integration domain
%
dominf = numHl*Tinf;				% total integration domain for inflammation
tinf = (0:Ds:dominf)';				% integration times for past history
breinf = length(tinf);				% breaks within integration domain
%
q = exp(-t*kmco);					% survival (decay) functions
SimpsErr = simpsons(q,t(1),t(end),[])./(1./kmco) - ones(size(kmco)); %* -> 0 if properly integrated
%
qact = exp(-tact*kact);				% active remodeling (decay) function
SimpsAct = simpsons(qact,tact(1),tact(end),[])/(1/kact) - 1; %* -> 0 if properly integrated
%
qinf = exp(-tinf*kinf);				% inflammatory cell clearance (decay) function
SimpsInf = simpsons(qinf,tinf(1),tinf(end),[])/(1/kinf) - 1; %* -> 0 if properly integrated
%
%** initialize past history arrays (with original homeostatic values)
%
rho = 1050;							% arterial mass density (kg/m^3)
%
rhooM = phioM*rho;					% medial [e mt cz 2*cd] mass densities
rhooA = phioA*rho;					% advent [e ct cz 2*cd] mass densities
%
rhoReM = rhooM(1);					% referential medial elastin density (constant)
rhoReA = rhooA(1);					% referential advent elastin density (constant)
%
rhoRf = ones(breaks,1)*[rhooM([2,3,4,4]),...	% circ m M, ax c M, diag c M, diag c M
						rhooA([2,3,4,4])];		% circ c A, ax c A, diag c A, diag c A
%
kmc = [kmco; kmco];					% rate parameters
%
mN  = (ones(breaks,1)*kmco).*rhoRf;	% nominal mass density production rates
Ups = ones(size(rhoRf));			% stimulus functions
mR  = mN.*Ups;						% referential mass density production rates
%
ltTauInvM = ones(breaks,1);			% medial circum stretches
ltTauInvA = ones(breaks,1);			% advent circum stretches
lzTauInv  = ones(breaks,1);			% axial stretches
%
lTauInv = [ltTauInvM ltTauInvA lzTauInv];	% stretches
%
c1m = c1m0*ones(breaks,1);			% material parameters
c2m = c2m0*ones(breaks,1);
c1c = c1c0*ones(breaks,1);
c2c = c2c0*ones(breaks,1);
Gm  =  Gm0*ones(breaks,1);
Gc  =  Gc0*ones(breaks,1);
%
act = ones(breact,1)*rio;			% inner radius for active length integral
%
muinf = zeros(breinf,1);			% inflammatory cell infiltration rate (normalized)
%
%** loading
%
[Pfo,so] = PreStr(PAR,parT0,parA,riio,emco);	% P, f, stresses at o
%
Po    = Pfo(1);						% inner pressure at o
fo    = Pfo(2);						% axial force at o
sgmIo = so(7);						% tr(sgm) at o
%
fctr  = 1.3635;						% relative increase in pressure at day 28 (h)
fctrf = 1.1500;						% relative increase in pressure at day 224
%
Epso = 1;							% relative change in flow rate at o
%
%** initialize solution variables/arrays
%
s = 0:Ds:days;						% G&R times
%
steps = length(s);					% total G&R steps
%
P    = Po*ones(steps,1);			% pressure
ri   = rio*ones(steps,1);			% inner radius
hM   = hMo*ones(steps,1);			% medial thickness
hA   = hAo*ones(steps,1);			% adventitial thickness
Upsi = ones(steps,8);				% stimulus functions
Eps  = Epso*ones(steps,1);			% relative cardiac output Q/Qo
lz   = ones(steps,1);				% axial stretch relative to state o
f    = fo*ones(steps,1);			% axial force
rhoR = ones(steps,1)*[rhooM([2,3,4,4]),...	% circ m M, ax c M, diag c M, diag c M
					  rhooA([2,3,4,4])];	% circ c A, ax c A, diag c A, diag c A
Dinf = zeros(steps,1);				% inflammatory cell density relative to maximum
st = so(5)*ones(steps,1);			% circum. stress
sz = so(6)*ones(steps,1);			%  axial  stress
sv = so(7)/3*ones(steps,1);			%  volum. stress
%
ltM = 1; ltA = 1;					% circum stretches at o
JM  = 1; JA  = 1;					% Jacobians at o
%
%% ***** COMPUTE GROWTH & REMODELING FOR s > 0
%
options = optimoptions(@fsolve,'Display','none','FunctionTolerance',1e-9,'StepTolerance',1e-9,'OptimalityTolerance',1e-9);
%
j = 0;			% maximum number of local iterations performed (upon convergence)
%
for k = 1:steps
	%
	%** profile for pressure elevation -> maintenance -> recovery -> plateau 
	%
	if s(k) <= SP			% increase P over SP days (SP input)
		P(k) = (1 + (fctr-1)*1/2*(1-cos(pi*s(k)/SP)))*Po;
	elseif s(k) <= 28		% maintain P until day 28
		P(k) = fctr*Po;
	elseif s(k) <= 224		% decrease P until day 224
		P(k) = (fctr-(fctr-fctrf)*sin(pi/2*(s(k)-28)/(224-28)))*Po;
	else					% maintain P afterwards (?)
		P(k) = fctrf*Po;
	end
	%
	%** hypertension-induced increase in active tone (if any)
	%
	Tmax = (1 + KP*(P(k)/Po-1))*Tmaxo;		% gain KP (input)
	%
	%** shift past history arrays so that current state is located first
	%
	qaux = q;
	qaux(2:end,:) = qaux(1:end-1,:);
	kmc(2,:) = kmc(1,:);
	%
	c1m(2:end,:) = c1m(1:end-1,:);	c2m(2:end,:) = c2m(1:end-1,:);
	c1c(2:end,:) = c1c(1:end-1,:);	c2c(2:end,:) = c2c(1:end-1,:);
	Gm(2:end,:)  = Gm(1:end-1,:);	Gc(2:end,:)  = Gc(1:end-1,:);
	%
	rhoRf(2:end,:)   = rhoRf(1:end-1,:);
	mN(2:end,:)      = mN(1:end-1,:);
	mR(2:end,:)      = mR(1:end-1,:);
	lTauInv(2:end,:) = lTauInv(1:end-1,:);
	Ups(2:end,:)     = Ups(1:end-1,:);
	%
	act(2:end) = act(1:end-1);
	%
	muinf(2:end) = muinf(1:end-1);
	%
	mufi1 = muinf(1);						% updated iteratively
	%
	rhoRfiter = zeros(size(rhoRf(1,:)));	% set iterative densities to zero
	%
	i = 0;									% local iteration counter
	%
	if k > 1, f(k) = f(k-1); end			% axial force guess
	%
	%** perform local iterations until desired tolerance is attained
	%
	while norm(rhoRf(1,:)-rhoRfiter)/norm(rhoRf(1,:)) > 1e-9 && i < 100
		%
		%* compute circumferential stretches and axial force
		%
		parL = [P(k),Eps(k),lz(k),sgmIo,JM,JA,rio,hMo,hAo,ce,Ge,rhoReM,rhoReA,rhoRf(1,1),...
				alp,Tmax,lM,l0,CB,CS,kact,rho,t(1),t(end),tact(1),tact(end),Ds];
		%
		[ltf,fval,eF] = fsolve(@(ltf) LaplaceCMM(ltf,parL,mR,qaux,kmco,kmc,lTauInv,act,qact,c1m,c2m,c1c,c2c,Gm,Gc),[ltM ltA f(k)],options);
		%
		ltM  = ltf(1);									% medial circum stretch
		ltA  = ltf(2);									% advent circum stretch
		f(k) = ltf(3);									% vessel axial force
		%
		%* additional variables from converged solution
		%
		lrM = JM/ltM/lz(k);								% medial radial strech
		lrA = JA/ltA/lz(k);								% advent radial strech
		%
		hM(k) = lrM*hMo;								% medial thickness
		hA(k) = lrA*hAo;								% advent thickness
		ri(k) = (ltM*(2*rio+hMo)-hM(k))/2;				% inner radius
		%
		lTauInv(1,:) = [1/ltM	1/ltA	1/lz(k)];		% update current stretches
		%
		%* shear and wall stress (stimuli)
		%
		Dtau = Eps(k)*rio^3/ri(k)^3-1;					% relative change in shear stress
		%
		st(k) = P(k)*ri(k)/(hM(k)+hA(k));						% circum. stress
		sz(k) = f(k)/(pi*(hM(k)+hA(k))*(2*ri(k)+hM(k)+hA(k)));	%  axial  stress
		sv(k) = (st(k)+sz(k))/3;								%  volum. stress
		%
		Dsv = sv(k)/(sgmIo/3)-1;						% relative change in intramural stress
		%
		%* inflammatory cell density (stimulus)
		%
		if ~sum(sv(1:k)>=svmax)
			muinf(1) = kinf*Kinf*Dsv;					% adaptive inflammation (reversible phase)
			if muinf(1) < 0, muinf(1) = 0; end
		elseif sv(k) >= svmax
			muinf(1) = mufi1 + mudotp*Ds;				% maladaptive inflammation (chronic phase)
			if muinf(1) > kinf, muinf(1) = kinf; end
		elseif sv(k) <= svmax && sum(sv(1:k)>=svmax)
			muinf(1) = mufi1 - mudotm*Ds;				% maladaptive inflammation (recovery phase)
			if muinf(1) < 0, muinf(1) = 0; end
		end
		%
		Dinf(k) = simpsons(muinf.*qinf,tinf(1),tinf(end),[]); % inflammatory cell density
		%
		%* update inflammation-dependent parameters (PR input)
		%
		if PR == 1										% passive properties evolve with inflammation
			DinfPas = Dinf(k);
		else											% passive properties remain maladaptive
			DinfPas = max(Dinf(1:k));
		end
		%
		c1m(1) = c1m0 + (DinfPas^fctInf)*(c1m4-c1m0);
		c2m(1) = c2m0 + (DinfPas^fctInf)*(c2m4-c2m0);
		c1c(1) = c1c0 + (DinfPas^fctInf)*(c1c4-c1c0);
		c2c(1) = c2c0 + (DinfPas^fctInf)*(c2c4-c2c0);
		Gm(1)  =  Gm0 + (DinfPas^fctInf)*( Gm4-Gm0 );
		Gc(1)  =  Gc0 + (DinfPas^fctInf)*( Gc4-Gc0 );
		%
		%* rate parameters, survival functions, and stimulus functions
		%
		kmc(1,:) = kmco*(1+Dsv^2);
		%
		q(2:end,:) = qaux(2:end,:).*( ones(breaks-1,1) * exp(-Ds*mean(kmc)) );
		%
		Ups(1,:) = [ 1 +  Kim*Dsv -  Ksm*Dtau +  Kfm*Dinf(k) ...
					(1 + KicM*Dsv - KscM*Dtau + KfcM*Dinf(k))*ones(1,3) ...
					(1 + KicA*Dsv - KscA*Dtau + KfcA*Dinf(k))*ones(1,4)];    %* a row array
		%
		%* update mass production and integrate mass densities
		%
		rhoRfiter = rhoRf(1,:);							% densities from previous time step (i=0) or iteration (i>0)
		mN(1,:)   = kmc(1,:).*rhoRfiter;				% update nominal rates
		mR(1,:)   = mN(1,:).*Ups(1,:);					% update mass production
		%
		rhoRf(1,:) = simpsons(mR.*q,t(1),t(end),[]);	% update referential mass densities
		%
		JM = ( rhoReM + sum(rhoRf(1,1:3)) ) / rho;		% update medial [e mt cz 2*cd] volume ratios
		JA = ( rhoReA + sum(rhoRf(1,5:7)) ) / rho;		% update advent [e ct cz 2*cd] volume ratios
		%
		i = i+1;										% update iteration counter
		%
	end
	%
	j = max(i,j);
	%
	%** other converged values (for plots)
	%
	act(1)    = ri(k);
	rhoR(k,:) = rhoRf(1,:);
	Upsi(k,:) = Ups(1,:);
	%
end
%
%% ***** PLOT FIGURE
%
scrsz = get(0,'ScreenSize');
figure(fign)
set(gcf,'position',[0.05*scrsz(3) 0.05*scrsz(4) 0.8*scrsz(3) 0.8*scrsz(4)])
%
%** inner pressure
%
subplot(341)
hold on
plot(s,P/P(1),'k','lines',line,'linew',1)
ylabel('Systolic pressure $P/P_o$ (-)','interpreter','latex')
set(gca,'xlim',[0 days],'ylim',[0.95 1.4],'XTick',[0 days/4 days/2 3*days/4 days])
set(gca,'fontsize',12)
%
%** stimulus function for medial smooth muscle
%
subplot(342)
hold on
plot(s,Upsi(:,1),'k','lines',line,'linew',1)
set(gca,'xlim',[0 days],'ylim',[0 3],'XTick',[0 days/4 days/2 3*days/4 days])
ylabel('Med. SM stimulus $\Upsilon^m_M$ (-)','interpreter','latex')
set(gca,'fontsize',12)
%
%** referential mass density of medial smc relative to homeostatic
%
subplot(343)
hold on
plot(s,rhoR(:,1)/rhoR(1,1),'k','lines',line,'linew',1)
set(gca,'xlim',[0 days],'ylim',[0 6],'XTick',[0 days/4 days/2 3*days/4 days])
ylabel('SM density $\rho^m_{MR}/\rho^m_{Mo}$ (-)','interpreter','latex')
set(gca,'fontsize',12)
%
%** medial wall thickness
%
subplot(344)
hold on
plot(s,hM/hMo,'k','lines',line,'linew',1)
set(gca,'xlim',[0 days],'ylim',[0 6],'XTick',[0 days/4 days/2 3*days/4 days])
ylabel('Med. thickness $h_M/h_{Mo}$ (-)','interpreter','latex')
set(gca,'fontsize',12)
%
%** luminal radius
%
subplot(345)
hold on
plot(s,ri/rio,'k','lines',line,'linew',1)
set(gca,'xlim',[0 days],'ylim',[0.9 1.05],'XTick',[0 days/4 days/2 3*days/4 days])
ylabel('Inner radius $a/a_o$ (-)','interpreter','latex')
set(gca,'fontsize',12)
%
%** stimulus function Upsilon for adventitial collagen
%
subplot(346)
hold on
plot(s,Upsi(:,5),'k','lines',line,'linew',1)
set(gca,'xlim',[0 days],'ylim',[0 3],'XTick',[0 days/4 days/2 3*days/4 days])
ylabel('Adv. coll. stimulus $\Upsilon^c_A$ (-)','interpreter','latex')
set(gca,'fontsize',12)
%
%** referential mass density of adventitial collagen relative to homeostatic
%
subplot(347)
hold on
plot(s,rhoR(:,5)/rhoR(1,5),'k','lines',line,'linew',1)
set(gca,'xlim',[0 days],'ylim',[0 6],'XTick',[0 days/4 days/2 3*days/4 days])
ylabel('Adv. coll. density $\rho^c_{AR}/\rho^c_{Ao}$ (-)','interpreter','latex')
set(gca,'fontsize',12)
%
%** adventitial wall thickness
%
subplot(348)
hold on
plot(s,hA/hAo,'k','lines',line,'linew',1)
set(gca,'xlim',[0 days],'ylim',[0 6],'XTick',[0 days/4 days/2 3*days/4 days])
ylabel('Adv. thickness $h_A/h_{Ao}$ (-)','interpreter','latex')
set(gca,'fontsize',12)
%
%** circumferential stress
%
subplot(349)
hold on
plot(s,P.*ri./(hM+hA),'k','lines',line,'linew',1)
set(gca,'xlim',[0 days],'ylim',[0 300],'XTick',[0 days/4 days/2 3*days/4 days])
xlabel('Time (days)','interpreter','latex')
ylabel('Circ. stress $\sigma_{\theta\theta}$ (kPa)','interpreter','latex')
set(gca,'fontsize',12)
%
%** axial stress
%
subplot(3,4,10)
hold on
plot(s,f./(pi*(hM+hA).*(2*ri+hM+hA)),'k','lines',line,'linew',1)
set(gca,'xlim',[0 days],'ylim',[0 300],'XTick',[0 days/4 days/2 3*days/4 days])
xlabel('Time (days)','interpreter','latex')
ylabel('Axial stress $\sigma_{zz}$ (kPa)','interpreter','latex')
set(gca,'fontsize',12)
%
%** volumetric stress and threshold
%
subplot(3,4,11)
hold on
if isempty(get(gca,'children'))
	plot(s,svmax*ones(size(s)),'--k','linew',0.5)
	text(14,svmax+15,strcat('$\sigma^*_v=',int2str(svmax),'$ kPa'),'interpreter','latex','fontsize',12)
end
plot(s,1/3*(P.*ri./(hM+hA)+f./(pi*(hM+hA).*(2*ri+hM+hA))),'k','lines',line,'linew',1)
set(gca,'xlim',[0 days],'ylim',[0 210],'XTick',[0 days/4 days/2 3*days/4 days],'YTick',[0 70 140 210])
xlabel('Time (days)','interpreter','latex')
ylabel('Vol. stress $\sigma_{v}$ (kPa)','interpreter','latex')
set(gca,'fontsize',12)
%
%** inflammatory cell density
%
subplot(3,4,12)
hold on
plot(s,Dinf,'k','lines',line,'linew',1)
set(gca,'xlim',[0 days],'ylim',[-0.05 1.05],'XTick',[0 days/4 days/2 3*days/4 days],'YTick',[0 0.3 0.7 1.0])
xlabel('Time (days)','interpreter','latex')
ylabel('Inflammation $\Delta\varrho_{\varphi}$ (-)','interpreter','latex')
set(gca,'fontsize',12)
%
drawnow
%
end