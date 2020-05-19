%%% Solution of ModelAntibody
par0 = zeros(1,12);
lower = [1/30,0,1/10,0,0,0,0,0,0,0,1/5,0];
upper = [1/5,1,1,1/3,1/5,1,1,1,1/5,1e-4,1,1e-5];
options = optimset('MaxFunEvals',100000,'MaxIter',50000, "TolFun",1e-30);
p_estimate = fmincon(@(par)ODE_fit(exp_C,c_testing, par),par0,[],[],[],[],lower,...
    upper,[],options);
%p_estimate = fminsearch(@(par)ODE_fit(exp_C,exp_h, par),par0);
diff = ODE_fit(exp_C,c_testing, p_estimate);
diff2 = ODE_fit(exp_C,c_testing, par0);

S_0 = 60480209; %3.34e+6;
L_0 = 100; % Assumed
I_0 = 100;
A_0 = 50; %Assumed double of I_0
Sq_0 = 0;
T_0 = 4324;
Ii_0 = 100;
H_0 = 10;
Hm_0 = 0;
C_0 = 132;
y0 = [S_0 L_0 I_0 A_0 Sq_0 T_0 Ii_0 H_0 Hm_0 C_0];
[t, y] = ode45(@(t,y)DeqnModel(t,y,p_estimate),[1:1:78],y0);
T = y(:,6).';
T_cum = double(78);
for i=1:78
    if i == 1
        T_cum(i) = T(i);
    else
        T_cum(i)=T(i)+T_cum(i-1);
    end
end
tiledlayout(2,2)
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
plot(T_cum);
hold on;
plot(c_testing)
title('Testing')
nexttile
plot(y(:,1))
title('S')
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
betas = 1e-8;%fixed
betaa = 1e-8*0.5; %fixed
betai = 1e-8*0.1;%fixed
lambda = 1/4.8; %fixed
rhom = 0;%fixed
alpham = 0;%fixed
theta = 0.5;
eta = par(1); %fixed
rhoh = par(2);
if t <= 38
    alphai = par(5);
    alphaq = 0*par(9);
    alphah = par(4);%fixed
    q = par(12);
    tau = par(3);
else
    alphai = par(6);
    alphaq = par(7);
    alphah = par(8);%fixed
    q = par(10);
    tau = par(11);
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
function err = ODE_fit(exp_C,c_testing, par)
S_0 = 60480209; %3.34e+6;
L_0 = 100; % Assumed
I_0 = 100;
A_0 = 50; %Assumed double of I_0
Sq_0 = 0;
T_0 = 4324;
Ii_0 = 100;
H_0 = 10;
Hm_0 = 0;
C_0 = 132;
y0 = [S_0 L_0 I_0 A_0 Sq_0 T_0 Ii_0 H_0 Hm_0 C_0];
[t, y] = ode45(@(t,y)DeqnModel(t,y,par),[1:1:78],y0);
C = y(:,10).';
T = y(:,6).';
T_cum = double(78);
for i=1:78
    if i == 1
        T_cum(i) = T(i);
    else
        T_cum(i)=T(i)+T_cum(i-1);
    end
end
err = sum(((C-exp_C)/max(exp_C)).^2)+sum(((T_cum(26:78)-c_testing(26:78))/max(c_testing(26:78))).^2);
%err = sum(((H_t.'-exp_h)/max(exp_h)).^2);
end
