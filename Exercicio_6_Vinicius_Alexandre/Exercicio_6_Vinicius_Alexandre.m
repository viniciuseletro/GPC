%   UNIVERSIDADE FEDERAL DO CEARÁ
%   PROGRAMA DE PÓS-GRADUAÇÃO EM ENGENHARIA ELÉTRICA
%   CONTROLE ADAPTATIVO PREDITIVO
%   Vinícius Alexandre de Mesquita


clear,close all, clc;
Ts = 0.03;
N = 5;
Nu = 5;
lambda = 1;     % Lambda (Ponderação)
Kp = 1.0;                % Ganho na planta
% dd = 15;              % Atraso de transporte na planta

Ct11 = tf(1,[0.7 1]);
Ct12 = tf(5,[0.3 1]);
Ct21 = tf(1,[0.5 1]);
Ct22 = tf(2,[0.4 1]);

Ct = [Ct11,Ct12;Ct21,Ct22];
[ny,nent] = size(Ct);
Cd = c2d(Ct,Ts);
Cd.TimeUnit = 'minutes';
Cc=[1 0 0 0];
[AA,BB,DD,G,F,E,H] = CalculaG(N,Cd,Cc);

[nnn,nnnn]=size(BB);
nest = nnn/nnnn;

% Cálculo de M
Minv = eye(N*nent);
for j=nent+1:nent*N
    Minv(j,j-nent) = -1;
end

M = inv(Minv);
%%%%%%%%%%%%%%%%%%%%
Qdelta = eye(N*nent);

Gnu = G(:,1:nent*Nu);
Qlambda = lambda*eye(Nu*nent);
Qlambda(1,1) = 0;
Qlambda(2,2) = 0;

%Multiplica a entrada de referências
W = repmat(eye(nent),N,1);

%Tamanho de G
ng = length(G);

% Restrições
Ymax = [ 1,1,0.61; 1,1,0.305;];
Ymin = [ 0 0 0; 0 0 0];
Umax = [ 0.3 1 1; 0.2 1 1];
Umin = [ -.1 -1 -1; -.1 -1 -1];
dUmax = [ 1 0.3 1;1 0.1 1];
dUmin = [ -1 -.05 -1;-1 -.05 -1];

%Matriz A de restrições
Ar=[M;-M;eye(ng);-eye(ng);G;-G];

% Matriz Hessiana H
Hfc = (Gnu'*Qdelta*Gnu+Qlambda);

for ii=1:3
    % Parâmetros da S-function
    Parametros.nin          = length(W);
    Parametros.nent        = nent;
    Parametros.Nu           = Nu;
    Parametros.N             = N;
    Parametros.G             = G;
    Parametros.M            = M;
    Parametros.Hfc          = Hfc;
    Parametros.Ar            = Ar;
    Parametros.Umax      = Umax(:,ii);
    Parametros.Umin       = Umin(:,ii);
    Parametros.dUmax    = dUmax(:,ii);
    Parametros.dUmin    = dUmin(:,ii);
    Parametros.ymax      = Ymax(:,ii);
    Parametros.ymin       = Ymin(:,ii);
    Parametros.Qdelta    = Qdelta;
    Parametros.Ts            = Ts;
    
    % Simulação no simulink e geração de gráficos
    
    sim('Exercicio_6_Vinicius_2',4)
    
    [yer,u,du] = simul.signals.values;
    tempo = simul.time;
    y1 = yer(:,1);
    y2 = yer(:,2);
    ref1 = yer(:,3);
    ref2 = yer(:,4);
    u1 = u(:,1);
    u2 = u(:,2);
    du1 = du(:,1);
    du2 = du(:,2);
    
    figure(ii)
    
    subplot(3,1,1)
    plot(tempo,y1,'LineWidth',1.5)
    hold on
    plot(tempo,y2,'LineWidth',1.5)
    plot(tempo,ref1,'k--','LineWidth',0.1)
    plot(tempo,ref2,'m--','LineWidth',0.1)
    ylabel('Outputs and References')
    legend('y_1','y_2','Referência 1','Referência 2',...
        'Location','southeast')
    hold off
    title ({[' \lambda = ' num2str(lambda) ...
        '    Q_{\lambda 1,1} = ' num2str(Qlambda(1,1))...
        '    Q_{\lambda 2,2} = ' num2str(Qlambda(2,2))],...
        [    'Gráfico y(t) ' ...
        '    y_1_{max} = ' num2str(Ymax(1,ii)) ...
        '  y_2_{max} = ' num2str(Ymax(2,ii)) ...
        '    y_1_{min} = ' num2str(Ymin(1,ii)) ...
        '  y_2_{min} = ' num2str(Ymin(2,ii))]})
    
    
    subplot(3,1,2)
    plot(tempo,u1,'LineWidth',1.5)
    hold on
    plot(tempo,u2,'LineWidth',1.5)
    ylabel('Inputs')
    legend('u_1','u_2')
    title(['Gráfico u(t) ' ...
        '    u_1_{max} = ' num2str(Umax(1,ii)) ...
        '  u_2_{max} = ' num2str(Umax(2,ii)) ...
        '    u_1_{min} = ' num2str(Umin(1,ii)) ...
        '  u_2_{min} = ' num2str(Umin(2,ii))])
    
    subplot(3,1,3)
    plot(tempo,du1,'LineWidth',1.5)
    hold on
    plot(tempo,du2,'LineWidth',1.5)
    ylabel('dU')
    legend('\Delta u_1','\Delta u_2')
    title(['Gráfico \Delta u(t) ' ...
        '    \Delta u_1_{max} = ' num2str(dUmax(1,ii)) ...
        '  \Delta u_2_{max} = ' num2str(dUmax(2,ii)) ...
        '    \Delta u_1_{min} = ' num2str(dUmin(1,ii)) ...
        '  \Delta u_2_{min} = ' num2str(dUmin(2,ii))])

end





