function Si = efast(Y,OMi,MI,time_points,var)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Input:
%                   Y:  synthetic data from main function
%                   OMi:    Frequency vectors
%                   MI: Maximum Fourier coefficient
%                   time_points:    time instnace of interest for SA
%                   var:    Dimension of state variables of systems
%
%       Output:
%                   Si:     first-order sensitivity index
%                   Sti:    Overall sensitivity index
%                   rangeSi:    Variance of first-order sensitivity index
%                   rabgeSti:   Variance of overall sensitivity index
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[a b c k NR]=size(Y);   
% a: number of samples in each search curve here is 65
% b: number of time points that the SA is required
% here the time points are 50,1000 and 5000, so here b
% is 3.
% c: the dimension of states in dynamics, here the
% c is 4, because four states include in the
% system.
% k: the number of parameters of interest, included
% the dummy parameter.
% NR: the numer of search curve in SA, here the NR
% is 5.

for u=var
    for t=1:length(time_points)%1:length(Y(1,:,1,1,1))  %loop through time, here it is 2
        [t u] %t=[1000];
        ttest_sti=zeros(k,k);
        ttest_si=ttest_sti;
        for i=1:k %loop through parameters
            % Initialize AV,AVi,AVci to zero.
            AV = 0;
            AVi = 0;
            AVci = 0;
            for L=1:NR%length(Y(1,t,u,i,:))
                Y(:,t,u,i,L) = (Y(:,t,u,i,L)-mean(Y(:,t,u,i,L)))';
                % Fourier coeff. at [1:OMi/2].
                N=length(Y(:,t,u,i,L));
                NQ = (N-1)/2;
                N0 = NQ+1;
                COMPL = 0;
                Y_VECP = Y(N0+(1:NQ),t,u,i,L)+Y(N0-(1:NQ),t,u,i,L);
                Y_VECM = Y(N0+(1:NQ),t,u,i,L)-Y(N0-(1:NQ),t,u,i,L);
                for j=1:OMi/2
                    ANGLE = j*2*(1:NQ)*pi/N;
                    C_VEC = cos(ANGLE);
                    S_VEC = sin(ANGLE);
                    AC(j) = (Y(N0,t,u,i,L)+Y_VECP'*C_VEC')/N;
                    BC(j) = Y_VECM'*S_VEC'/N;
                    COMPL = COMPL+AC(j)^2+BC(j)^2;
                end
                % Computation of V_{(ci)}.
                Vci(L) = 2*COMPL;
                %AVci = AVci+Vci;
                % Fourier coeff. at [P*OMi, for P=1:MI].
                COMPL = 0;
                Y_VECP = Y(N0+(1:NQ),t,u,i,L)+Y(N0-(1:NQ),t,u,i,L);
                Y_VECM = Y(N0+(1:NQ),t,u,i,L)-Y(N0-(1:NQ),t,u,i,L);
                for j=OMi:OMi:OMi*MI
                    ANGLE = j*2*(1:NQ)*pi/N;
                    C_VEC = cos(ANGLE');
                    S_VEC = sin(ANGLE');
                    AC(j) = (Y(N0,t,u,i,L)+Y_VECP'*C_VEC)/N;
                    BC(j) = Y_VECM'*S_VEC/N;
                    COMPL = COMPL+AC(j)^2+BC(j)^2;
                end
                % Computation of V_i.
                Vi(L) = 2*COMPL;
                % Computation of the total variance
                % in the time domain.
                V(L) = Y(:,t,u,i,L)'*Y(:,t,u,i,L)/N;
            end %L
            % Computation of sensitivity indexes.
            Si(i,t,u) = mean(Vi)/mean(V);
        end %i
    end %t
end