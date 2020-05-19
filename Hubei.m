%%% Solution of ModelAntibody
%par0 = [0.8,0.8];
%lower = [0,0];
%upper = [1,1];
%p_estimate = fmincon(@(par)ODE_fit(exp_C, par),par0,[],[],[],[],lower,upper);
par0 = zeros(1,11);
lower = [1/30,0,0,0.1,0,0,1/10,0,1/10,0,1/30];
upper = [1/10,1,1,1,1,1,1,1,1,1/5,1];
options = optimset('MaxFunEvals',100000,'MaxIter',50000, "TolFun",1e-30);
p_estimate = fmincon(@(par)ODE_fit(exp_C,d_Hospitalised,d_infection,c_Quarantine, par),par0,[],[],[],[],lower,...
    upper,[],options);
%p_estimate = fminsearch(@(par)ODE_fit(exp_C,exp_h, par),par0);
diff = ODE_fit(exp_C,d_Hospitalised,d_infection,c_Quarantine, p_estimate);
diff2 = ODE_fit(exp_C,d_Hospitalised,d_infection,c_Quarantine, par0);
S_0 = 75560000; %3.34e+6;
L_0 = 10; % Assumed
I_0 = 0;
A_0 = 2; %Assumed double of I_0
Sq_0 = 0;
T_0 = 1;
Ii_0 = 0;
H_0 = 1;
Hm_0 = 1000;
C_0 = 1;
y0 = [S_0 L_0 I_0 A_0 Sq_0 T_0 Ii_0 H_0 Hm_0 C_0];
[t, y] = ode45(@(t,y)DeqnModel(t,y,p_estimate),[1:1:139],y0);
Sq = y(:,5).';
Sqcum = double(139);
for i=1:139
    if i == 1
        Sqcum(i) = Sq(i);
    else
        Sqcum(i)=Sq(i)+Sqcum(i-1);
    end
end
tiledlayout(3,3)
nexttile
plot(y(:,10))
hold on
plot(exp_C)
title("cum confirmed")
nexttile
plot(y(:,8)+y(:,9))
hold on
plot(d_Hospitalised)
title("d_H")
nexttile
plot(y(:,3))
hold on 
plot(d_infection)
title("I")
nexttile
plot(Sqcum)
hold on
plot(c_Quarantine)
title("Sq cum")
nexttile
plot(Sq)
title("Sq")
nexttile
plot(y(:,8))
hold on
plot(y(:,9))
title("H and Hm")
nexttile
plot(y(:,1))
title("S")
% == ODE ==
function ydot = DeqnModel(t,y, par)

S = y(1);
L = y(2);
I = y(3);
A = y(4);
Sq = y(5);
T = y(6);
Ii = y(7);
H = y(8);
Hm = y(9);
C = y(10);

ydot = zeros(10,1);

r = 1/14;%fixed
gamma = 0.154;%fixed
zeta = 1/14;%fixed
betas = 5e-9;%fixed
betaa = 5e-9*0.5; %fixed
betai = 5e-9*0.1;%fixed
lambda = 1/4; %fixed
eta = par(1); %fixed

theta = 0.817;
rhoh = par(2);
if t <= 55
    rhom = 0;%fixed
    alpham = 0;%fixed
    alphai = par(4);
    alphah = par(11);%fixed
    alphaq = par(5);
    eta2 =0;%fixed
    tau = par(6);
    q = 0;%fixed
elseif (t <=68) 
    rhom = 0;
    alphai = 0.1;%fixed
    alphah = par(11);%fixed
    alphaq = par(8);
    eta2 = 0;
    tau = par(9);
    alpham = 0;
    q = par(3);%fixed
elseif (t<=101)
    rhom = par(7);
    alphai = 0;%fixed
    alphah = 1/5;%fixed
    alphaq = par(8);
    tau = par(9);
    alpham = 1/5;
    eta2 = 1/10;
    q = par(3);%fixed
else
    rhom=0;
    alpham = 0;
    eta2 = 1;
    alphai = 0;%fixed
    alphah = 1/5;%fixed
    alphaq = par(8);
    tau = par(9);
    q = par(3);%fixed
end

ydot(1) = -(betas*I+betaa*A+betai*Ii)*S-q*S+r*Sq;
ydot(2) = (betas*I+betaa*A+betai*Ii)*S - lambda*L;
ydot(3) = theta*lambda*L - tau*I;
ydot(4) = (1-theta)*lambda*L - gamma*A;
ydot(5) = q*S+alphaq*T-(r+tau)*Sq;
ydot(6) = tau*Sq+tau*I-(alphaq+alphai+alphah+alpham)*T;
ydot(7) = alphai*T-(zeta+rhoh+rhom)*Ii;
ydot(8) = alphah*T +rhoh*Ii-eta*H;
ydot(9) = alpham*T +rhom*Ii-eta2*Hm;
ydot(10)= (alphai+alpham+alphah)*T;

end
% == Error Function ==
function err = ODE_fit(exp_C,d_Hospitalised,d_infection,c_Quarantine,par)
S_0 = 75560000; %3.34e+6;
L_0 = 10; % Assumed
I_0 = 0;
A_0 = 2; %Assumed double of I_0
Sq_0 = 0;
T_0 = 1;
Ii_0 = 0;
H_0 = 1;
Hm_0 = 1000;
C_0 = 1;
y0 = [S_0 L_0 I_0 A_0 Sq_0 T_0 Ii_0 H_0 Hm_0 C_0];
[t, y] = ode45(@(t,y)DeqnModel(t,y,par),[1:1:139],y0);
C = y(:,10).';
H = y(:,8).';
Hm = y(:,9).';
d_H = zeros(1,139);
d_H(1:68) = H(1:68);
d_H(69:end)=H(69:end)+Hm(69:end);
I = y(:,3).';
Sq = y(:,5).';
Sqcum = double(139);
for i=1:139
    if i == 1
        Sqcum(i) = Sq(i);
    else
        Sqcum(i)=Sq(i)+Sqcum(i-1);
    end
end
err = sum(((C-exp_C)/max(exp_C)).^2)...
    +sum(((d_H(56:139)-d_Hospitalised(56:139))/max(d_Hospitalised(56:139))).^2)...
    +0*sum(((I-d_infection)/(max(d_infection)).^2))...
    +sum(((Sqcum(56:139)-c_Quarantine(56:139))/max(c_Quarantine(56:139))).^2);
end
