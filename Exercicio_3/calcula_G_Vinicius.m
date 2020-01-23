function [G,F,E,At,Bt,D,H] = calcula_G_Vinicius(A,B,Ce,d,N)

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

end

