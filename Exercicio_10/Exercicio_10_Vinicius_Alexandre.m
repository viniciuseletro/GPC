% Dados iniciais
clear all;
close all;
clc;

% Dados iniciais
ttotal = 1800;
Ts = 60;
N = 3;
Nu = 3;
niter = round(ttotal/Ts);

% Dados de controle
ref = 5;
lambda = 1;
Ql = lambda*eye(Nu);
W = ref*ones(N,1);
MU0 = [1;zeros(Nu-1,1)];
Mt      = tril(ones(Nu));
M       = inv(Mt);


% Caracteristicas do tanque
hant = 1;   % Altura inicial
qi = 0.1392; % Entrada do sistema (vazao de entrada)
A = 10;     % Area do tanque
a = 0.01*pi;    % Area do tubo de saida
k1 = a*sqrt(2*9.8);     % Constante

% Parametros da funcao
prmt_tanque.A = A;
prmt_tanque.k1 = k1;
prmt_tanque.h0 = hant;

% Altura inicial
h0 = hant;
imp(1) = hant;
impd(1) = hant;

% Calcula G
delta = 0.01;
Gnu = calcula_G_EPSAC(N,Nu,delta,qi,Ts,prmt_tanque);

% Condicoes iniciais
yt = 1;
xt = 0;
u = qi;
pert = zeros(1,niter);
pert(round(niter/2):end) = 0.1;
t = 0;

for i=1:niter
    
    % Saida y(t)
    yt(i+1) = yt(i) + (Ts/A)*(u(i)+pert(i)) - (Ts*k1/A)*sqrt(yt(i));
    
    % x(t)
    xt(i+1) = yt(i) + (Ts/A)*u(i) - (Ts*k1/A)*sqrt(yt(i));
    
    % n(t)
    nt(i) = yt(i+1) - xt(i+1);
    
    % Ubase
    ub = u(i)*ones(1,N);
   
    for j=1:N
        
        % Calculo de Ybase
        if j ==1
            yb(j+1) = yt(i+1) + (Ts/A)*(ub(j)) - (Ts*k1/A)*sqrt(yt(i+1)) + nt(i);
        else
            yb(j+1) = yb(j) + (Ts/A)*(ub(j)) - (Ts*k1/A)*sqrt(yb(j)) + nt(i);
        end
        
    end
    
    yb = yb(2:end);
    
    % Calculo de G dinamico
    Gnu = calcula_G_EPSAC(N,Nu,delta,u(i),Ts,prmt_tanque);
    
    % Calculo dos ganhos
    K1      = pinv(Gnu'*Gnu + M'*Ql*M)*Gnu';
    K1      = K1(1,:);
    K2      = pinv(Gnu'*Gnu + M'*Ql*M)*M'*Ql;
    K2      = K2(1,:);
    
    % Vetor de U(t-1)
    U0 = [u(i);zeros(Nu-1,1)];
    
    % Calculo de U_otimo
    uopt = K1*(W-yb') - K2*(M*ub(1:Nu)'-U0);
    
    % Calculo de U
    u(i+1) = uopt(1) + ub(1);
    
    % Atualiza tempo
    t(i+1) = t(i) + 1;
    
end

% Graficos
subplot(2,1,1)
plot(t,yt,'LineWidth',2)
title('Altura do fluido no tanque','FontSize',14)
xlabel('Amostras')
ylabel('h (m)')
grid on

subplot(2,1,2)
plot(t,u,'r','LineWidth',2)
title('Sinal de controle')
xlabel('Amostras')
grid on

