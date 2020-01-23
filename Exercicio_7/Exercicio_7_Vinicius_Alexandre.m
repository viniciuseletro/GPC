%   UNIVERSIDADE FEDERAL DO CEARÁ
%   PROGRAMA DE PÓS-GRADUAÇÃO EM ENGENHARIA ELÉTRICA
%   CONTROLE ADAPTATIVO PREDITIVO
%   Vinícius Alexandre de Mesquita

close all;
clear,clc;

% Valores iniciais
Ts      = 1;
N       = 5;
Nu     = 5;

% Matriz do sistema
A                = [1 -0.8];
B                = [0.4 0.6];
Processo   = filt([0 B],A,Ts);

% Matriz do polinômio C
Ce         = [1 0];

% Ponderação Lambda
lambda      = 0.8;
Qlambda   = lambda*eye(Nu);

% Ponderação Delta 
delta      = 1;
Qdelta   = delta*eye(N);

% Cálculo do polinômio Ã e B~
Ac         = conv([1 -1],A);
Btil        = conv([1 -1],B);

nac       = length(Ac)-1;
nbtil      = length(Btil);

Atil        = [Ac(2:end) zeros(1,nbtil-nac)];

natil      = length(Atil);

C           = [Ce(2:end) zeros(1,natil-length(Ce)+1)];

% Cálculo da matriz Abarra
AA        = [-Atil' [eye(natil-1);zeros(1,natil-1)]];

% Cálculo da matriz Bbarra
BB         = Btil';

% Cálculo da matriz Dbarra
DD        = C' - Atil'; 

% Cálculo da matriz H
H           = [1 zeros(1,nbtil-1)];

% Pré-alocação das matrizes G, F e E
G           = zeros(N,N);
F            = zeros(N,natil);
E            = zeros(N,1);

% Matriz S armazena a soma dos termos h
S            = 0;

for i=1:N
    
    if i==1
        G(1,1)    = H*BB;
    else
        G(i,:)      = [H*(AA^(i-1))*BB G(i-1,1:end-1)];
    end
    F(i,:)          = H*(AA^i);

    E(i,:)          = H*(AA^(i-1))*DD;  
    
    S(i+1)       = S(i)+G(i,1);              
    
end

S          = S(2:end);

% Matriz Gnu
Gnu                           = G(:,1:Nu);
Gnu(Nu:end,end)     = S(1:N-Nu+1);

% Matriz M
Mt      = tril(ones(Nu));
M       = inv(Mt);

% Cálculo dos ganhos
K1      = (inv(Gnu'*Qdelta*Gnu+M'*Qlambda*M))*Gnu'*Qdelta;
K2      = (inv(Gnu'*Qdelta*Gnu+M'*Qlambda*M))*M'*Qlambda;

% Primeira linha dos ganhos
k1      = K1(1,:);
k2      = K2(1,:);

% Cálculo de KF
KF      = k1*F;

% Cálculo de KE
KE      = k1*E;

% Cálculo de KW
KW    = sum(k1);



% Comparacao com anterior

[ny,nent] = size(Processo);

[AAant,BBant,DDant,Gant,Fant,Eant,Hant] = CalculaG(N,Processo,Ce);

[nnn,nnnn]=size(BBant);
nest = nnn/nnnn;

Gnuant = Gant(:,1:nent*Nu);

% Multiplica a entrada de referencias
W = repmat(eye(nent),N,1);

%Tamanho de G
ng = length(Gant);

% Restricoes

Ymax = 10;
Ymin = 0;
Umax = 0.3;
Umin = -0.2;
dUmax = 10;
dUmin = -10;

%Matriz A de restriï¿½ï¿½es
Ar=[M;-M;eye(ng);-eye(ng);Gant;-Gant];

% Matriz Hessiana H
Hfc = (Gnuant'*Qdelta*Gnuant+Qlambda);

% Parametros s-function

Parametros.nin          = length(W);
Parametros.nent        = nent;
Parametros.Nu           = Nu;
Parametros.N             = N;
Parametros.G             = Gant;
Parametros.M            = M;
Parametros.Hfc          = Hfc;
Parametros.Ar            = Ar;
Parametros.Umax      = Umax;
Parametros.Umin       = Umin;
Parametros.dUmax    = dUmax;
Parametros.dUmin     = dUmin;
Parametros.ymax       = Ymax;
Parametros.ymin        = Ymin;
Parametros.Qdelta     = Qdelta;
Parametros.Ts            = Ts;


% Atualiza simulink
sim('Exercicio_7_Vinicius_Comparacao',30)

[yer,u] = simul.signals.values;
tempo = simul.time;
y1 = yer(:,1);
y2 = yer(:,2);
u1 = u(:,1);
u2 = u(:,2);

% Geracao de graficos
figure

subplot(2,1,1)
stairs(tempo,y1,'-','LineWidth',2)
hold on
stairs(tempo,y2,'m--o','LineWidth',1.0)
grid on;
ylabel('Outputs and References')
legend('y_{Btil}','y_{quadprog}',...
    'Location','southeast')
hold off
title ({[' \lambda = ' num2str(lambda)],...
    ['Grafico y(t)']})


subplot(2,1,2)
stairs(tempo,u1,'-','LineWidth',2)
hold on
stairs(tempo,u2,'m--o','LineWidth',1.0)
grid on;
ylabel('Inputs')
legend('u_{Btil}','u_{quadprog}')
title(['Grafico u(t) ' ...
    '    u_{max} = ' num2str(Umax) ...
    '    u_{min} = ' num2str(Umin)]);



