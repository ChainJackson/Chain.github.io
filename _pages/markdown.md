---
permalink: /markdown/
title: "Codes"
author_profile: true
redirect_from: 
  - /md/
  - /markdown.html
---
<!-- TOC -->
- [MATLAB](#matlab)
  - [1.Symplectic Analysis](#1symplectic-analysis)
  - [2.Zeros \& Chebfun](#2zeros--chebfun)
<!-- TOC -->

# MATLAB

## 1.Symplectic Analysis
```matlab
%%%============================Copyright============================%%%
%%% Version Oct. 17, 2024
%%%
%%% Lizichen Chen <lzcchen@zju.edu.cn>
%%% Department of Engineering Mechanics, Zhejiang University
%%%
%%%===========================Description===========================%%%
%%% Symplectic contact analysis
%%%
%%%=================================================================%%%

clear; clc; format long;

nu= 0.25;
E = 1;
beta = 0.1;
a = -0.5;
b = 0.5;
l = 10;
h = 10;
d = 0.02; %maximum indentation depth

xi = l*cosh(beta*l)/sinh(beta*l)-1/beta;

syms x z mu real

vec = [1 1 exp(beta*x) exp(beta*x)];
Mt = diag(vec);

% bar = waitbar(0,'Start...');

%% eigenvectors and eigen-solutions for zero eigenvalue

% Saint-Venant solutions only require this part

phi0s_0 = [1 0 0 0].';
phi0a_0 = [0 1 0 0].';
phi0s_1 = vpa([0 -nu*x E 0].');
phi0a_1 = [-x 0 0 0].';
phi0a_2 = vpa([0 (nu*x.^2)/2 -E*x 0].');
phi0_3 = vpa([(2*(1+nu)/E*(-E*l*exp(-beta*x)/(beta^2*sinh(beta*l))+(E*x.^2)/(2*beta)-(xi*E/beta+E/beta^2)*x)-nu*x^3/6+(xi*nu*x^2)/2) 0 0 (E*l*exp(-beta*x)/(beta*sinh(beta*l))+E*x/beta-(xi*E/beta+E/beta^2))].');

f0s_0 = phi0s_0;
f0a_0 = phi0a_0;
f0s_1 = phi0s_1 + z*phi0s_0;
f0a_1 = phi0a_1 + z*phi0a_0;
f0a_2 = phi0a_2 + z*phi0a_1 + (z^2)*phi0a_0/2;
f0_3 = vpa(phi0_3 + z*phi0a_2 + (z^2)*phi0a_1/2 + (z^3)*phi0a_0/factorial(3) + xi*(z*phi0s_1 + (z^2)*phi0s_0/2));

%% eigenvectors and eigen-solutions for general eigenvalues

% the corresponding eigenvalues are derived through chebfun[1]
mu_num1 = [0.216770703063582 + 0.107139620004182i
0.378795349308418 + 0.134884827323222i
0.538497796737552 + 0.152222782289058i
0.697264749719046 + 0.164961556015475i
0.855540207813688 + 0.175063152082181i
1.01351974051428 + 0.183444051676859i
1.17130518761091 + 0.190610170120139i
1.32895553408938 + 0.196871496153662i
1.48650759360047 + 0.202432161302456i
1.64398566053789 + 0.207433908702813i
1.80140647972292 + 0.211979256700566i
1.95878200270734 + 0.216144794762847i
2.11612100859733 + 0.219989271520295i
2.27343010345277 + 0.223558750075335i
2.43071436122255 + 0.226890022828818i
2.58797774865504 + 0.230012947581372i
2.74522341511208 + 0.232952090178354i
2.90245389516800 + 0.235727907232762i
3.05967125332764 + 0.238357615465374i
3.21687718938876 + 0.240855842412500i
3.37407311646433 + 0.243235121396738i
3.53126021964467 + 0.245506273491593i
3.68843950071283 + 0.247678706117715i
3.84561181265540 + 0.249760649212011i
4.00277788660121 + 0.251759344015567i
4.15993835306842 + 0.253681195455886i
4.31709375888278 + 0.255531896240977i
4.47424458076788 + 0.257316528745475i
4.63139123635102 + 0.259039649296082i
4.78853409314339 + 0.260705358385119i
4.94567347591917 + 0.262317359541845i
5.10280967281849 + 0.263879008992454i
5.25994294042609 + 0.265393357786577i
5.41707350802170 + 0.266863187721901i
5.57420158115649 + 0.268291042131689i
5.73132734467771 + 0.269679252392493i
5.88845096529876 + 0.271029960846979i
6.04557259379307 + 0.272345140708481i
6.20269236687458 + 0.273626613412135i
6.35981040881639 + 0.274876063795931i
6.51692683284908 + 0.276095053429523i
6.67404174237337 + 0.277285032355608i
6.83115523201516 + 0.278447349465525i
6.98826738854668 + 0.279583261695419i
7.14537829169310 + 0.280693942200320i
7.30248801484104 + 0.281780487639536i
7.45959662566274 + 0.282843924686832i
7.61670418666728 + 0.283885215862355i
7.77381075568887 + 0.284905264769371i
7.93091638632029 + 0.285904920807295i
8.08802112829864 + 0.286884983422585i
8.24512502784955 + 0.287846205950889i
8.40222812799476 + 0.288789299096715i
8.55933046882782 + 0.289714934090879i
8.71643208776158 + 0.290623745560942i
8.87353301975074 + 0.291516334145325i
9.03063329749260 + 0.292393268878157i
9.18773295160813 + 0.293255089368516i
9.34483201080592 + 0.294102307795008i
9.50193050203066 + 0.294935410734140i
9.65902845059789 + 0.295754860838849i
9.81612588031655 + 0.296561098381653i
9.97322281360052 + 0.297354542675364i
10.1303192715703 + 0.298135593382812i
10.2874152741460 + 0.298904631725828i
10.4445108401322 + 0.299662021602654i
10.6016059872956 + 0.300408110621975i
10.7587007324367 + 0.301143231060922i
10.9157950914548 + 0.301867700753675i
11.0728890794087 + 0.302581823916593i
11.2299827105717 + 0.303285891915252i
11.3870759984827 + 0.303980183978222i
11.5441689559935 + 0.304664967861964i
11.7012615953120 + 0.305340500470808i
11.8583539280428 + 0.306007028435598i
12.0154459652238 + 0.306664788654295i
12.1725377173614 + 0.307314008797460i
12.3296291944618 + 0.307954907781353i
12.4867204060611 + 0.308587696211090i
12.6438113612526 + 0.309212576796101i
12.8009020687127 + 0.309829744739967i
12.9579925367247 + 0.310439388106467i
13.1150827732008 + 0.311041688163590i
13.2721727857030 + 0.311636819707074i
13.4292625814627 + 0.312224951364918i
13.5863521673981 + 0.312806245884183i
13.7434415501319 + 0.313380860401332i
13.9005307360067 + 0.313948946697183i
14.0576197310994 + 0.314510651437577i
14.2147085412360 + 0.315066116400656i
14.3717971720037 + 0.315615478691656i
14.5288856287636 + 0.316158870946031i
14.6859739166620 + 0.316696421521670i
14.8430620406411 + 0.317228254680875i
15.0001500054491 + 0.317754490762767i
15.1572378156500 + 0.318275246346727i
15.3143254756321 + 0.318790634407391i
15.4714129896168 + 0.319300764461766i
15.6285003616662 + 0.319805742708898i
15.7855875956911 + 0.320305672162549i
15.9426746954575 + 0.320800652777339i
16.0997616645938 + 0.321290781568682i
16.2568485065966 + 0.321776152726898i
16.4139352248371 + 0.322256857725865i
16.5710218225665 + 0.322732985426466i
16.7281083029215 + 0.323204622175180i
16.8851946689291 + 0.323671851898058i
17.0422809235116 + 0.324134756190349i
17.1993670694912 + 0.324593414402023i
17.3564531095938 + 0.325047903719392i
17.5135390464540 + 0.325498299243074i
17.6706248826179 + 0.325944674062460i
17.8277106205477 + 0.326387099326902i
17.9847962626245 + 0.326825644313766i
18.1418818111521 + 0.327260376493525i
18.2989672683599 + 0.327691361592067i
18.4560526364059 + 0.328118663650328i
18.6131379173798 + 0.328542345081382i
18.7702231133055 + 0.328962466725174i
18.9273082261435 + 0.329379087900919i];

mu_num = [-conj(mu_num1); mu_num1];

eta = [(-beta/2 + sqrt(beta^2-4*mu.^2-4*mu*beta*sqrt(nu))/2) (-beta/2 - sqrt(beta^2-4*mu.^2-4*mu*beta*sqrt(nu))/2) (-beta/2 + sqrt(beta^2-4*mu.^2+4*mu*beta*sqrt(nu))/2) (-beta/2 - sqrt(beta^2-4*mu.^2+4*mu*beta*sqrt(nu))/2)].';

CharaM = [exp(eta(2)*l) exp(eta(3)*l) exp(eta(4)*l); exp(-eta(2)*l) exp(-eta(3)*l) exp(-eta(4)*l); exp(eta(2)*l)/(eta(2)+beta) exp(eta(3)*l)/(eta(3)+beta) exp(eta(4)*l)/(eta(4)+beta)];
D1 = E*mu;
D1vecctor = -D1*[exp(eta(1)*l) exp(-eta(1)*l) exp(eta(1)*l)/(eta(1)+beta)].';
Dtemp = inv(CharaM)*D1vecctor;
Dvec = [D1;Dtemp];

Avec = vpa((mu*nu-(eta+beta).^2/mu)./(E*mu*(eta+beta)).*Dvec);
Bvec = vpa(-(mu-nu*(eta+beta).^2/mu)./(E*eta.*(eta+beta)).*Dvec);
Cvec = vpa(-(eta+beta).*Dvec/mu);

phii = vpa(simplify([sum(Avec.*exp(eta*x)) sum(Bvec.*exp(eta*x)) sum(Cvec.*exp(eta*x)) sum(Dvec.*exp(eta*x))].'));

% f_trivali = exp(mu*z).*phii;

f_trival = [];

for k_ = 1:length(mu_num)
    if k_ > length(mu_num)/2
        f_trival = [f_trival vpa(real(Mt*subs(exp(-mu*h).*((real(subs(exp(mu*z),mu,mu_num(k_)))+imag(subs(exp(mu*z),mu,mu_num(k_)))*1i).*(real(simplify(subs(phii,mu,mu_num(k_)),'Criterion','preferReal','Steps',100))+imag(simplify(subs(phii,mu,mu_num(k_)),'Criterion','preferReal','Steps',100))*1i)),mu,mu_num(k_)))) vpa(imag(Mt*subs(exp(-mu*h).*((real(subs(exp(mu*z),mu,mu_num(k_)))+imag(subs(exp(mu*z),mu,mu_num(k_)))*1i).*(real(simplify(subs(phii,mu,mu_num(k_)),'Criterion','preferReal','Steps',100))+imag(simplify(subs(phii,mu,mu_num(k_)),'Criterion','preferReal','Steps',100))*1i)),mu,mu_num(k_))))];
    else
        f_trival = [f_trival vpa(real(Mt*subs((real(subs(exp(mu*z),mu,mu_num(k_)))+imag(subs(exp(mu*z),mu,mu_num(k_)))*1i).*(real(simplify(subs(phii,mu,mu_num(k_)),'Criterion','preferReal','Steps',100))+imag(simplify(subs(phii,mu,mu_num(k_)),'Criterion','preferReal','Steps',100))*1i),mu,mu_num(k_)))) vpa(imag(Mt*subs((real(subs(exp(mu*z),mu,mu_num(k_)))+imag(subs(exp(mu*z),mu,mu_num(k_)))*1i).*(real(simplify(subs(phii,mu,mu_num(k_)),'Criterion','preferReal','Steps',100))+imag(simplify(subs(phii,mu,mu_num(k_)),'Criterion','preferReal','Steps',100))*1i),mu,mu_num(k_))))];
    end
end

f_special0 = vpa([Mt*f0s_0 Mt*f0a_0 Mt*f0s_1 Mt*f0a_1 Mt*f0a_2 Mt*f0_3]);

f = [f_special0 f_trival];

%% derivation of coefficients via calculus of variations

%parallel computing
pa=parpool(20);

Em = [];
Fm = [];

num_len = length(f(1,:));

parfor i_ = 1:num_len
    for j_ = 1:num_len
        temp1 = matlabFunction(transpose(f(3:4,i_))*f(1:2,j_),'vars',{x,z});
        temp2 = matlabFunction(f(3,i_)*f(1,j_),'vars',{x,z});
        temp3 = matlabFunction(f(2,i_)*f(4,j_),'vars',{x,z});
        temp4 = matlabFunction(f(1,i_)*f(3,j_),'vars',{x,z});
        Em(i_,j_) = -integral(@(x) temp1(x,h),-l,l,'ArrayValued',true,'RelTol',0,'AbsTol',1e-20) + integral(@(x) temp2(x,0),a,b,'ArrayValued',true,'RelTol',0,'AbsTol',1e-20) - integral(@(x) temp3(x,0),-l,l,'ArrayValued',true,'RelTol',0,'AbsTol',1e-20) - integral(@(x) temp4(x,0),-l,a,'ArrayValued',true,'RelTol',0,'AbsTol',1e-20) - integral(@(x) temp4(x,0),b,l,'ArrayValued',true,'RelTol',0,'AbsTol',1e-20);
    end
    temp5 = matlabFunction(f(3,i_),'vars',{x,z});
    Fm(i_) = integral(@(x) d*temp5(x,0),a,b,'ArrayValued',true,'RelTol',0,'AbsTol',1e-20);
    % str=['Waiting...',num2str(100*i_/num_len),'%'];
    % waitbar(i_/num_len,bar,str)
end
delete(pa);

Mm = inv(Em)*(Fm.');
ansvec = vpa(f*Mm)

%% plot results on the surface
w = ansvec(1)
ww = matlabFunction(w);
fplot(@(x) ww(x,0),[-l,l])
set(gca,'YDir','reverse');

figure;
u = ansvec(2);
uu = matlabFunction(u);
fplot(@(x) uu(x,0),[-l,l])
set(gca,'YDir','reverse');

figure;
sigma = ansvec(3);
ss = matlabFunction(sigma);
fplot(@(x) ss(x,0),[-l,l])
set(gca,'YDir','reverse');

figure;
tau = ansvec(4);
tt = matlabFunction(tau);
fplot(@(x) tt(x,0),[-l,l])
set(gca,'YDir','reverse');

%% Reference
% [1] Chebfun Team (2024). Chebfun - current version (https://github.com/chebfun/chebfun), GitHub. Retrieved October 16, 2024.
```

## 2.Zeros & Chebfun
```matlab
clear; clc; format long;

nu= 0.25;
E = 1;
beta = 0.1;
l = 10;
h = 10;

syms mu

eta1 = -beta/2 + sqrt(beta^2-4*mu.^2-4*mu*beta*sqrt(nu))/2;
eta2 = -beta/2 - sqrt(beta^2-4*mu.^2-4*mu*beta*sqrt(nu))/2;
eta3 = -beta/2 + sqrt(beta^2-4*mu.^2+4*mu*beta*sqrt(nu))/2;
eta4 = -beta/2 - sqrt(beta^2-4*mu.^2+4*mu*beta*sqrt(nu))/2;

equa =  (eta1-eta2).*(eta3-eta4)+(eta1-eta4).*(eta2-eta3).*cosh((eta1-eta2-eta3+eta4)*l)+(eta1-eta3).*(eta4-eta2).*cosh((eta1-eta2+eta3-eta4)*l);
funct = matlabFunction(equa,'vars',mu);

Zerostemp = chebfun(funct,[0.1 5],'splitting','on')
mu_zerotemp = roots(Zerostemp,'complex','norecursion');
Verifications1 = double(vpa(subs(equa,mu,mu_zerotemp)));
index = find(abs(real(Verifications1))<=1e-6 & abs(imag(Verifications1))<=1e-6);
mu_zeros1 = mu_zerotemp(index);
AN1 = Verifications1(index);
index1 = find(real(mu_zeros1)>=0 & imag(mu_zeros1)>=0);
mu_zeros1_fin = mu_zeros1(index1);
VerificationsX = double(abs(vpa(subs(equa,mu,mu_zeros1_fin))));
```
