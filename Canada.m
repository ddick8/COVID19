%%% Solution of ModelAntibody
par0 = zeros(1,10);
lower = [0,0,0,0,0,0,0,0,1/10,0];
upper = [1,1,1,1,1,1,1,1,1,1,1,1,1];
options = optimset('MaxFunEvals',100000,'MaxIter',50000, "TolFun",1e-30);
p_estimate = fmincon(@(par)ODE_fit(exp_C,daily_testing,d_H, par),par0,[],[],[],[],lower,...
upper,[],options);
%p_estimate = fminsearch(@(par)ODE_fit(exp_C,exp_h, par),par0);
diff = ODE_fit(exp_C,daily_testing,d_H, p_estimate);
diff2 = ODE_fit(exp_C,daily_testing,d_H, par0);

S_0 = 37671158; %3.34e+6;
L_0 = 10; % Assumed
I_0 = 1;
A_0 = 10; %Assumed double of I_0
Sq_0 = 0;
T_0 = 10;
Ii_0 = 0;
H_0 = 0;
Hm_0 = 0;
C_0 = 1;
y0 = [S_0 L_0 I_0 A_0 Sq_0 T_0 Ii_0 H_0 Hm_0 C_0];
[t, y] = ode45(@(t,y)DeqnModel(t,y,p_estimate),[1:1:82],y0);
T = y(:,6).';
tiledlayout(2,3)
nexttile
plot(y(:,2))
hold on;
plot(y(:,3),'--og')
hold on;
plot(y(:,4),'--og')
hold on;
plot(y(:,5),'+')
hold on;
plot(y(:,7))
nexttile
plot(y(:,10));
hold on;
plot(exp_C);
title('cum confirmed')
nexttile
plot(T);
hold on;
plot(daily_testing)
title('Testing')
nexttile
plot(y(:,1))
title('S')
nexttile
plot(y(:,8))
title('d_H')
hold on
plot(d_H);
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
betas = 1e-9;%fixed
betaa = 1e-9*0.5; %fixed
betai = 1e-9*0.1;%fixed
lambda = 1/5.2; %fixed
rhom = 0;%fixed
alpham = 0;%fixed
theta = 0.5;

eta = par(1); %fixed
rhoh = par(2);
alphah = par(3);%fixed
alphai = par(5);
if t <= 51
    tau = par(4);
    q = par(6);
    alphaq = 0;
else
    q = par(8);
    tau = par(9);
    alphaq = par(10);
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

% == Error Function ==
function err = ODE_fit(exp_C,daily_testing,d_H, par)
S_0 = 37671158; %3.34e+6;
L_0 = 10; % Assumed
I_0 = 1;
A_0 = 10; %Assumed double of I_0
Sq_0 = 0;
T_0 = 10;
Ii_0 = 0;
H_0 = 0;
Hm_0 = 0;
C_0 = 1;
y0 = [S_0 L_0 I_0 A_0 Sq_0 T_0 Ii_0 H_0 Hm_0 C_0];
[t, y] = ode45(@(t,y)DeqnModel(t,y,par),[1:1:82],y0);
C = y(:,10).';
T = y(:,6).';
err = sum(((C-exp_C)/max(exp_C)).^2)+sum(((T(46:82)-daily_testing(46:82))/max(daily_testing(46:82))).^2)+...
sum(((y(45:82,8).'-d_H(45:82))/max(d_H(45:82))).^2);
%err = sum(((H_t.'-exp_h)/max(exp_h)).^2);
end
