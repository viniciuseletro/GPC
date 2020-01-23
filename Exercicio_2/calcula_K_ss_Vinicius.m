function [kw,kf,ke,G,F,E,At,Bt,D,H] = calcula_K_ss_Vinicius(A,B,Ce,d,N,Nu,lambda,lbd0)

    W=ones(N+d,1);
    Ac = conv(A,[1 -1]);
    Atio = Ac(2:end);
    C = [Ce(2:end) zeros(1,length(Atio)-length(Ce)+1)];

    na = length(A);
    nb = length(B);
    nc = length(C);

    At = [-Atio' [eye(na-1);zeros(1,na-1)]];
    Bt = B';
    C = [C zeros(1,(na-nc))];
    D = C'-[Atio';zeros(nc-na,1)];
    H = [1 zeros(1,na-1)];

    G = zeros(N,N);
    F = zeros(N,na);
    E = zeros(N,1);

    for i=1:N
        if i==1
            G(1,1) = H*Bt;
        else
            G(i,:) = [H*(At^(i-1))*Bt G(i-1,1:end-1)];
        end
        F(i,:) = H*(At^i);
        E(i,:) = H*(At^(i-1))*D;
    end

    Gnu=G(:,1:Nu);
    Q=lambda*eye(Nu);
    if lbd0>0
        for i=1:lbd0
            Q(i,i)=0;
        end
    end
    K=(Gnu'*Gnu+Q)\Gnu';
    k1=K(1,:);
    
    %Ganho KW
    kw=k1*W;
    
    %Ganho KF
    [mf,nf]=size(F);
    for i=1:nf
        kf(i)=k1*F(:,i);
    end
    
    %Ganho KE
    ke=(k1*E)';

end

