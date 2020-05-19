par0 = zeros(1,9);
lower = [0,0,0,0,1/10,1/10,0,0,1/10];
upper = [1,1/2,1,1/5,1/2.6,1,1,1,1];
options = optimset('MaxFunEvals',100000,'MaxIter',50000, "TolFun",1e-10);
p_estimate = fmincon(@(par)ODE_fit(exp_C,d_hospitalised,active_infection,par),par0,[],[],[],[],lower,upper,[],options);
diff = ODE_fit(exp_C,d_hospitalised,active_infection, p_estimate);
diff2 = ODE_fit(exp_C,d_hospitalised,active_infection, par0);
S_0 = 113460000; %3.34e+6;
L_0 = 10; % Assumed
I_0 = 1;
A_0 = 2; %Assumed double of I_0
Sq_0 = 0;
T_0 = 1;
H_0 = 1;
Hm_0 = 0;
Ii_0 = 0;
C_0 = 1;
y0 = [S_0 L_0 I_0 A_0 Sq_0 T_0 Ii_0 H_0 Hm_0 C_0];
[t, y] = ode45(@(t,y)DeqnModel(t,y,p_estimate),[1:1:90],y0);
C = y(:,10).';
H = y(:,8).';
tiledlayout(2,2)
nexttile
plot(y(:,2))
hold on;
plot(y(:,4),'--og')
hold on;
plot(y(:,5))
nexttile
plot(C)
hold on
plot(exp_C)
title('cum confirmed')
nexttile
plot(H);
hold on
plot(d_hospitalised)
title('H')
nexttile
plot(y(:,3))
hold on 
plot(active_infection)
title("I")
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
lambda = 1/5.2; %fixed
gamma=0.154; %fixed
theta = 0.5; %fixed
betas = 1e-9;%fixed
betaa = 1e-9*0.5; %fixed
betai = 0;
alpham=0;
rhom=0;
eta=par(4);
alphai =0;
zeta = 0;
rhoh =0;


if t <= 5
    q = par(1);
    tau = par(2);
    alphaq = par(3); %fixed
    alphah = par(5);
else 
    alphaq = par(6);
    q = par(7);
    alphah = par(8);
    tau = par(9);
end

ydot(1) = -(betas*I+betaa*A+betai*Ii)*S-q*S+r*Sq;
ydot(2) = (betas*I+betaa*A+betai*Ii)*S - lambda*L;
ydot(3) = theta*lambda*L - tau*I;
ydot(4) = (1-theta)*lambda*L - gamma*A;
ydot(5) = q*S+alphaq*T-(r+tau)*Sq;
ydot(6) = tau*Sq+tau*I-(alphaq+alphai+alphah+alpham)*T;
ydot(7) = alphai*T-(zeta+rhoh+rhom)*Ii;
ydot(8) = alphah*T +rhoh*Ii-eta*H;
ydot(9) = alpham*T +rhom*Ii-eta*Hm;
ydot(10)= (alphai+alpham+alphah)*T;

end
function err = ODE_fit(exp_C,d_hospitalised,active_infection,par)
S_0 = 113460000; %3.34e+6;
L_0 = 10; % Assumed
I_0 = 1;
A_0 = 2; %Assumed double of I_0
Sq_0 = 0;
T_0 = 1;
H_0 = 1;
Hm_0 = 0;
Ii_0 = 0;
C_0 = 1;
y0 = [S_0 L_0 I_0 A_0 Sq_0 T_0 Ii_0 H_0 Hm_0 C_0];
[t, y] = ode45(@(t,y)DeqnModel(t,y,par),[1:1:90],y0);
C = y(:,10).';
H = y(35:90,8).';
err = sum(((C-exp_C)/max(exp_C)).^2)+0*sum(((H-d_hospitalised(35:90))/max(d_hospitalised(35:90))).^2)+...
0*sum(((y(:,3).'-active_infection)/max(active_infection)).^2);
end

