%   UNIVERSIDADE FEDERAL DO CEARÁ 
%   PROGRAMA DE PÓS-GRADUAÇÃO EM ENGENHARIA ELÉTRICA
%   CONTROLE ADAPTATIVO PREDITIVO
%   Vinícius Alexandre de Mesquita
%   Reprodução exercício - Livro

clear,clc;
close all;
Ts=1;        % Tempo de amostragem
N=3;        % Horizonte de predição
Nu=3;       % Fator
d=0;        % Atraso
A=[1 -0.8];     
B=[zeros(1,d) 0.4 0.6];
C=1;
lambda=0.8;     % Lambda (Ponderação)
lbd0=0;            % Índice do vetor lambda que será zerado
[Hcl1,Ur1,dU1]=GPCgen(B,A,d,Ts,N,Nu,lambda,lbd0);

%% Exercício
[yy,xx]=step(Hcl1,12);
duu=step(dU1,12);
plot(xx,yy,'b','LineWidth',2);
hold on
plot(xx,duu,'r','LineWidth',2);
title({'Exemplo Livro'},'FontSize',16);
xlabel('Discrete time');
legend({'y(t)','\Delta u(t)'},'Location','best','FontSize',14)

%% Outras abordagens

lambda2=0.1;
N2=3;
Nu2=3;
lbd02=0;
[Hcl2,Ur2,dU2]=GPCgen(B,A,d,Ts,N2,Nu2,lambda2,lbd02);

lambda3=100;
N3=10;
Nu3=10;
lbd03=0;
[Hcl3,Ur3,dU3]=GPCgen(B,A,d,Ts,N3,Nu3,lambda3,lbd03);

lambda4=2;
N4=3;
Nu4=3;
lbd04=1;
[Hcl4,Ur4,dU4]=GPCgen(B,A,d,Ts,N4,Nu4,lambda4,lbd04);

figure
step(Hcl2,Hcl3,Hcl4)
title({'Outras abordagens'},'FontSize',14);
xlabel('Discrete time');
legend({...
    ['\lambda = ' num2str(lambda2) ' N=' num2str(N2) ' Nu=' num2str(Nu2)]...
    ,['\lambda = ' num2str(lambda3) ' N=' num2str(N3) ' Nu=' num2str(Nu3)]...
    ,['\lambda = ' num2str(lambda4) ' N=' num2str(N4) ' Nu=' num2str(Nu4) ' \lambda_{1,1} = 0']}...
    ,'Location','southeast','FontSize',14)
