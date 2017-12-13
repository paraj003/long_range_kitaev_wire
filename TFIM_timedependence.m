%Generates data and saves to file.
clear

Delta=0.001; %define a time-step 
Tf=10^3; %Final time
T0=0;%Initial time
[x, w] = GetTrapRule((Tf-T0)/Delta+1, T0, Tf);
Ns = [100];% array of system sizes
Bs =2;% 1.8:0.05:2.2;
J = 1;
narr=floor(logspace(1,log10((Tf-T0)/Delta+1),25));%generate log-space array to get times.
%plot_avg(Ns, Bs, J, x, w);
ms=calculate_avgt(Ns, Bs, J,narr, x, w); 





function [x, w] = GetTrapRule(N, t0, tf)
    % Computes the points and weights for the Trapazoidal quadrature
    % rule from t0 to tf with N function evaluations.
    % Should converge >~ N^-2.
    N = N - 1; delta = (tf - t0) / N;
    x = zeros(1, N+1); w = zeros(1, N+1);
    x(1) = t0; x(end) = tf; 
    w(1) = delta / 2.0; w(end) = delta / 2.0;
    for n=1:N-1
        x(n+1) = t0 + n * delta;
        w(n+1) = delta;
    end
end

function i = Quad(func, x, w)
    % Approximates the integral of func using the quadrature scheme given
    % by the points x and the weights w.
    % x and w must be the same length.
    i = 0;
    for n=1:length(x)
        fprintf("Evaluating t=%d\n",x(n))
        i = i + w(n) * func(x(n));
    end
end

function m = create_m(N, B, J, eta)
    % create the matrix m where H_eta = A^dag m A (+ const) for
    % nearest neighbor interaction.
    % eta = 0: open boundary conditions.
    % eta = 1: periodic boundary conditions.
    % eta = -1: antiperiodic boundary conditions.
    % see appendix C.1.
    
    m = zeros(2*N, 2*N); a = J / 2;
    for i=1:N
        m(2*i-1, 2*i-1) = B / 2.0;
        m(2*i, 2*i) = -B / 2.0;
        j = i + 1;
        if j <= N
            [m(2*i-1, 2*j-1), m(2*j-1, 2*i-1)] = deal(-a);
            [m(2*i-1, 2*j), m(2*j, 2*i-1)] = deal(-a);
            [m(2*i, 2*j-1), m(2*j-1, 2*i)] = deal(a);
            [m(2*i, 2*j), m(2*j, 2*i)] = deal(a);
        end
    end
    i = N; j = 1;
    [m(2*i-1, 2*j-1), m(2*j-1, 2*i-1)] = deal(-eta*a);
    [m(2*i-1, 2*j), m(2*j, 2*i-1)] = deal(-eta*a);
    [m(2*i, 2*j-1), m(2*j-1, 2*i)] = deal(eta*a);
    [m(2*i, 2*j), m(2*j, 2*i)] = deal(eta*a);
end

function m = update_m(m, newB)
    % Rather than create a new m everytime, just update what has changed.
    % in this case only the B field.
    for k=1:2:length(m)
        m(k, k) = newB / 2.0;
        m(k+1, k+1) = -newB / 2.0;
    end
end

function mag = G_per_noboundary(W, eta)
    % helper function for the function G_per below. 
    % Computes the nonboundary portion of G_per(t) for just one 
    % sector (even or odd) depending on the W and eta provided.
    % W is for A(t) = W(t)A
    
    N = length(W) / 2;

    A = @(i,j,k) (W(2*i-1,2*j-1)-W(2*i,2*j-1))*(W(2*i+1,2*k-1)+W(2*i+2,2*k-1));
    B = @(i,j,k) (W(2*i-1,2*j-1)-W(2*i,2*j-1))*(W(2*i+1,2*k)+W(2*i+2,2*k));
    C = @(i,j,k) (W(2*i-1,2*j)-W(2*i,2*j))*(W(2*i+1,2*k-1)+W(2*i+2,2*k-1));
    D = @(i,j,k) (W(2*i-1,2*j)-W(2*i,2*j))*(W(2*i+1,2*k)+W(2*i+2,2*k));
    
    mag = 0;
    for i=1:N-1
        mag = mag + eta*A(i,N,1) - eta*A(i,1,N);
        mag = mag + eta*B(i,1,N) + eta*B(i,N,1);
        mag = mag - eta*C(i,1,N) - eta*C(i,N,1);
        mag = mag - eta*D(i,N,1) + eta*D(i,1,N);
        for j=1:N-1
            mag = mag + A(i,j+1,j) - A(i,j,j+1) ...
                      - B(i,j+1,j) - B(i,j,j+1) + 2*B(i,j,j) ...
                      + C(i,j,j+1) + C(i,j+1,j) + 2*C(i,j,j) ...
                      + D(i,j,j+1) - D(i,j+1,j);
        end
        %For j=N
        mag = mag + 2*(B(i,N,N) + C(i,N,N));
    end
    mag = -mag / (8*N);
end

function mag = G_per_justboundary(W, eta)
    % helper function for the function G_per below. 
    % Computes the boundary portion of G_per(t) for just one 
    % sector (even or odd) depending on the W and eta provided.
    % W is for A(t) = W(t)A
    
    N = length(W) / 2;
    A = @(j,k) (W(2*N-1,2*j-1)-W(2*N,2*j-1))*(W(1,2*k-1)+W(2,2*k-1));
    B = @(j,k) (W(2*N-1,2*j-1)-W(2*N,2*j-1))*(W(1,2*k)+W(2,2*k));
    C = @(j,k) (W(2*N-1,2*j)-W(2*N,2*j))*(W(1,2*k-1)+W(2,2*k-1));
    D = @(j,k) (W(2*N-1,2*j)-W(2*N,2*j))*(W(1,2*k)+W(2,2*k));
    
    mag = 0;
    mag = mag + eta*A(N,1) - eta*A(1,N);
    mag = mag + eta*B(1,N) + eta*B(N,1);
    mag = mag - eta*C(1,N) - eta*C(N,1);
    mag = mag - eta*D(N,1) + eta*D(1,N);
    for j=1:N-1
        mag = mag + A(j+1,j) - A(j,j+1) ...
                  - B(j+1,j) - B(j,j+1) + 2*B(j,j) ...
                  + C(j,j+1) + C(j+1,j) + 2*C(j,j) ...
                  + D(j,j+1) - D(j+1,j);
    end
    %For j=N
    mag = mag + 2*(B(N,N) + C(N,N));

    mag = -mag / (8*N);
end

function mag = G_per(W_even, W_odd)
    % Computes G_per(t) from equation 19.
    % W_even and W_odd are the time evolution A(t) = W(t)A for
    % the eta = 1 and -1 sectors respectively.
    mag = G_per_noboundary(W_even, 1) + G_per_noboundary(W_odd, -1) + ...
          G_per_justboundary(W_even, 1) - G_per_justboundary(W_odd, -1);
end

function avg_mag = avg_G_per(m_even, m_odd, x, w)
    % time average <G_per(t)>_t. Equation 19.
    % Approximates the average using the quadrature
    % scheme given by the points x and the weights w, where x is time.
    
    % Diagonlize only once, then use it in the function below.
    
    [Udag_even, D_even] = eig(m_even); U_even = ctranspose(Udag_even);
    [Udag_odd, D_odd] = eig(m_odd); U_odd = ctranspose(Udag_odd);
    N = length(D_even) / 2;
    fprintf("Diagonalized");

    function r = f(t)
        % Could replace all this with just
        %   W_even = expm(-2i*t*m_even)
        %   W_odd  = expm(-2i*t*m_odd)
        % but I think this is faster because I'd imagine that expm 
        % diagonalizes the matrix in order to exponentiate it. This way
        
        % we only diagonalize it once and just exponentiate the diagonal
        % matrix D for M = U^dag D U.
        
        E_even = zeros(2*N, 2*N); E_odd = zeros(2*N, 2*N);
        for j=1:2*N
            E_even(j, j) = exp(-2i*t*D_even(j, j));
            E_odd(j, j) = exp(-2i*t*D_odd(j, j));
        end
        W_even = Udag_even * E_even * U_even;
        W_odd = Udag_odd * E_odd * U_odd;
        
        r = G_per(W_even, W_odd);
    end

    % time average.
    avg_mag = Quad(@f, x, w) / (x(end) - x(1));
end
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


function plot_avg(Ns, Bs, J, x, w)
    % plot <G_per>_t = the time average of G_per(t) vs B. See equation 19. 
    % Take the time average with the quadrature rule defined by the points 
    % x and weights w. 
    % Ns and Bs are lists of values of N and B to find <G_per>_t for.
    
    ms = zeros(length(Ns), length(Bs));
    for j=1:length(Ns)
        m_even = create_m(Ns(j), 0, J, 1);
        m_odd = create_m(Ns(j), 0, J, -1);
        fprintf("N = %d\n", Ns(j));
        for i=1:length(Bs)
            fprintf("Bs = %d\n", Bs(i));
            % Update m with new B value. Don't recreate m. Waste of time.
            m_even = update_m(m_even, Bs(i));
            m_odd = update_m(m_odd, Bs(i));
            ms(j, i) = avg_G_per(m_even, m_odd, x, w);
        end
    end
    
    %Make the plot
    figure();
    plot(Bs, ms, '.-');
    set(0, 'defaulttextinterpreter', 'latex');
    legend("N = " + string(Ns));
    xlabel('$$B$$', 'FontSize', 14);
    ylabel('$$\langle G_{per}(t)\rangle_t$$', 'FontSize', 14);
    title(sprintf( ...
      "$J$ = %g; %d-point trapazoidal rule average for $t$ = [%g, %g]", ...
      J, length(x), x(1), x(end)), 'FontSize', 12);
end

function ms=calculate_avgt(Ns, Bs, J,narr, x, w)
    % plot <G_per>_t = the time average of G_per(t) vs B.
    %Plot avg for different times given by x(narr)
    
    ms = zeros(length(Ns), length(Bs),length(narr));
    for j=1:length(Ns)
        m_even = create_m(Ns(j), 0, J, 1);
        m_odd = create_m(Ns(j), 0, J, -1);
        fprintf("N = %d\n", Ns(j));
        for i=1:length(Bs)
            fprintf("Bs = %d\n", Bs(i));
            % Update m with new B value. Don't recreate m. Waste of time.
            m_even = update_m(m_even, Bs(i));
            m_odd = update_m(m_odd, Bs(i));
            ms(j, i,:) = avg_G_per_turbo(m_even, m_odd,narr, x, w);
        end
    mdata=ms(j,:,:);
    N=Ns(j);
    filename=sprintf('data/Avg_mag_N_%G_t_%G:%G:%G_Bs%G:%G.mat',Ns(j),x(1),x(2)-x(1),x(end),Bs(1),Bs(end))
    save(filename,'x','w','N','Bs','J','narr','mdata')
    end
end
    
