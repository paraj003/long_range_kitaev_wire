function avg_magt = avg_G_per_turbo(m_even, m_odd,narr,x,w)
    %calculates avg_mag at different times denoted by x(narr) and returns an array avg_magt 
    %using vectorized formula (See
    %Simplifying_Joe's_Equations.lyx for details)
    %cell array : 1-> even , 2-> odd
    m{1}=m_even;
    m{2}=m_odd;
    etaarr=[1,-1];
    mag=zeros(2,length(x)); %array to store for different bc, for different times
    avg_mag=zeros(length(narr));
    %I will avoid doing both bc because i know that it is the same value.
    %and just take twice the value
    for bc=1:1
        eta=etaarr(bc);
        [Udag, D] = eig(m{bc}); U = ctranspose(Udag);
        fprintf("Diagonalized eta=%d\n",eta)
        N = length(D) / 2;
        %define different correlation matrices (i have flipped signs of
        %eta from Joe's notes. Seems to make the answer better agree.
        AiAj=1/8*(zeros(N)+diag(-ones(N-1,1),1)+diag(ones(N-1,1),-1)+diag(eta,N-1)+diag(-eta,-N+1));
        AiAjd=1/8*(zeros(N)+diag(2*ones(N,1))+diag(-ones(N-1,1),1)+diag(-ones(N-1,1),-1)+diag(-eta,N-1)+diag(-eta,-N+1));
        AidAj=1/8*(zeros(N)+diag(2*ones(N,1))+diag(ones(N-1,1),1)+diag(ones(N-1,1),-1)+diag(eta,N-1)+diag(eta,-N+1));
        AidAjd=1/8*(zeros(N)+diag(ones(N-1,1),1)+diag(-ones(N-1,1),-1)+diag(-eta,N-1)+diag(eta,-N+1));
        %time dependence
        odd=[1:2:2*N];
        even=[2:2:2*N];
        for n=1:length(x)
            t=x(n);
            if(rem(t,100)==0)
                fprintf("t=%d\n",t)
            end
            W= Udag * diag(exp(-2i*t*diag(D))) * U;
            %define some matrices for faster multiplication/less number of
            %operations
            WOminus=W(odd,odd)-W(even,odd);
            WoTplus=transpose(W(odd,odd)+W(even,odd));
            WeTplus=transpose(W(odd,even)+W(even,even));
            Weminus=W(odd,even)-W(even,even);
            WOminusAiAj=WOminus*AiAj;
            WOminusAiAjd=WOminus*AiAjd;
            WeminusAidAj=Weminus*AidAj;
            WeminusAidAjd=Weminus*AidAjd;
            %Calculate the first off diagonal of the matrix multiplication
            %and sum it. Also calculate the (N,1) matrix 
            A1=sum(sum(WOminusAiAj(1:N-1,:).*transpose(WoTplus(:,2:N)),2));
            A1b=eta*(WOminusAiAj(N,:)*WoTplus(:,1));
            B1=sum(sum(WOminusAiAjd(1:N-1,:).*transpose(WeTplus(:,2:N)),2));
            B1b=eta*(WOminusAiAjd(N,:)*WeTplus(:,1));
            C1=sum(sum(WeminusAidAj(1:N-1,:).*transpose(WoTplus(:,2:N)),2));
            C1b=eta*(WeminusAidAj(N,:)*WoTplus(:,1));
            D1=sum(sum(WeminusAidAjd(1:N-1,:).*transpose(WeTplus(:,2:N)),2));
            D1b=eta*(WeminusAidAjd(N,:)*WeTplus(:,1));

%             A1=sum(diag(WOminus*AiAj*WoTplus,1));
%             A1b=eta*diag(WOminus*AiAj*WoTplus,-N+1);
%             B1=sum(diag(WOminus*AiAjd*WeTplus,1));
%             B1b=eta*diag(WOminus*AiAjd*WeTplus,-N+1);
%             C1=sum(diag(Weminus*AidAj*WoTplus,1));
%             C1b=eta*diag(Weminus*AidAj*WoTplus,-N+1);
%             D1=sum(diag(Weminus*AidAjd*WeTplus,1));
%             D1b=eta*diag(Weminus*AidAjd*WeTplus,-N+1);

%             A1=sum(diag((W(odd,odd)-W(even,odd))*AiAj*transpose(W(odd,odd)+W(even,odd)),1));
%             A1b=eta*diag((W(odd,odd)-W(even,odd))*AiAj*transpose(W(odd,odd)+W(even,odd)),-N+1);
%             B1=sum(diag((W(odd,odd)-W(even,odd))*AiAjd*transpose(W(odd,even)+W(even,even)),1));
%             B1b=eta*diag((W(odd,odd)-W(even,odd))*AiAjd*transpose(W(odd,even)+W(even,even)),-N+1);
%             C1=sum(diag((W(odd,even)-W(even,even))*AidAj*transpose(W(odd,odd)+W(even,odd)),1));
%             C1b=eta*diag((W(odd,even)-W(even,even))*AidAj*transpose(W(odd,odd)+W(even,odd)),-N+1);
%             D1=sum(diag((W(odd,even)-W(even,even))*AidAjd*transpose(W(odd,even)+W(even,even)),1));
%             D1b=eta*diag((W(odd,even)-W(even,even))*AidAjd*transpose(W(odd,even)+W(even,even)),-N+1);
            mag(bc,n)=-1/N*(A1+A1b+B1+B1b+C1+C1b+D1+D1b);%-1/N*(sum(diag(A1,1))+eta*A1(N,1)+sum(diag(B1,1))+eta*B1(N,1)+sum(diag(C1,1))+eta*C1(N,1)+sum(diag(D1,1))+eta*D1(N,1));

            %fprintf("mag(t)=%d\n",mag(bc,n))
        end
    end
    mag(2,:)=mag(1,:); % use this line to avoid doing both odd and even parts.
    %figure()
    %plot(x, squeeze(mag(bc,:)))

    for n=1:length(narr)
        avg_magt(n)=sum(w(1:narr(n)).*squeeze(mag(1,1:narr(n)))+w(1:narr(n)).*squeeze(mag(2,1:narr(n))))/(x(narr(n))-x(1));
    end
end
