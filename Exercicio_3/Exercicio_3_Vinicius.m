%   UNIVERSIDADE FEDERAL DO CEARÁ 
%   PROGRAMA DE PÓS-GRADUAÇÃO EM ENGENHARIA ELÉTRICA
%   CONTROLE ADAPTATIVO PREDITIVO
%   Vinícius Alexandre de Mesquita
%   Exercício 3

clc, clear;
close all;
Ts = 1;        % Tempo de amostragem
N = 10;        % Horizonte de predição
Nu = 10;       % Fator
d = 0;        % Atraso
A = [1 -0.8];     
B = [zeros(1,d) 0.4 0.6];
Ce = [1 0];
lambda = 0.01;     % Lambda (Ponderação)
W=ones(N+d,1);

[G,F,E,At,Bt,D] = calcula_G_Vinicius(A,B,Ce,d,N);

Q=lambda*eye(N);

dt = 1;
for i=1:4
    tsim = 21;
    H = [1 0];
    tempo = 0:dt:tsim;
    x=[0;0];
    ref=ones(1,length(tempo));
    du=zeros(1,length(tempo));
    U=zeros(1,length(tempo));

    ng = length(G);
    M=tril(ones(ng));

    Hfc = 2*(G'*G + Q);

    Ymax = [1.1 3 3 1.2];
    Ymin = [0 0 0 0];
    Umax = [5 1.1 5 1.2];
    Umin = [-5 -0.2 -5 -0.4];
    dUmax = [5 5 0.3 0.4];
    dUmin = [-5 -5 -0.3 -0.4];
    [y,yi] = filter(B,A,0);
    U0 = zeros(ng,1);
    umax = Umax(i)*ones(ng,1);
    umin = Umin(i)*ones(ng,1);
    dumax = dUmax(i)*ones(ng,1);
    dumin = dUmin(i)*ones(ng,1);
    ymax = Ymax(i)*ones(ng,1);
    ymin = Ymin(i)*ones(ng,1);
    I = eye(ng);
    x = [0;0];
    e=0;

    for k=1:tsim
        % Cálculo do erro
        e = y(k) - H*x;

        % Resposta livre
        f = F(1,:)*x + E*e;

        % Cálculo de delta U
    %     dU = kw*ref(k) - kf*x - ke*e;

        bfc = 2*G'*(f-W);
        U0(1) = U(k);
        B1 = umax - M*U0;
        B2 = M*U0 - umin;
        B5 = ymax - f;
        B6 = f - ymin;

        Ar = [M;-M;G;-G;I;-I];
        Br = [B1;B2;B5;B6;dumax;-dumin];

        dU = quadprog(Hfc,bfc,Ar,Br);
        dU = dU(1,1);
        du(k) = dU;
        % Cálculo de U
        U(k+1) = dU + U(k);

        % Cálculo dos estados
        x = Bt*dU + D*e + At*x;

        % Cálculo da saída
        [y(k+1),yi] = filter(B,A,U(k+1),yi);

    end

    figure(i)
    subplot(3,1,1)
    stairs(tempo(1:end-1),y(1:end-1),'Linewidth',2)
    title (['Gráfico y(t)      \lambda = ' num2str(lambda) '     y_{max} = '...
        num2str(Ymax(i)) '  y_{min} = ' num2str(Ymin(i))])

    subplot(3,1,2)
    stairs(tempo(1:end-1),U(2:end),'Linewidth',2)
    title(['Gráfico u(t)  ' '  u_{max} = ' num2str(Umax(i)) ...
        '  u_{min} = ' num2str(Umin(i))])

    subplot(3,1,3)
    stairs(tempo(1:end-1),du(2:end),'Linewidth',2)
    title(['Gráfico du(t)    \Delta u_{max} = ' num2str(dUmax(i)) ...
        '  \Delta u_{min} = ' num2str(dUmin(i))])
end