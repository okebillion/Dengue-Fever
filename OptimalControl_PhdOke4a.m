function [T, X, U, Lambda, flag] = OptimalControl_PhdOke4a(control_switch)
    % This function computes the optimal control and the corresponding solution using forward-backward sweep
    % Always run compute the output of this file by running
    % OptimallControlMainFile.
    test = -1;
    delta = 0.001; %set tolerance
    N = 1000; %number of subdivisions
    Time = 100; %timespan
    T =linspace(0,Time,N+1);
    k = Time/N; %step
    %%

     % Assigning bounds to each of the control variables
     u1max = 1;  u2max = 1;  u3max = 1;
    
    % Pre-allocating memory space for control variables
    U = zeros(3, N + 1);
    
    % Pre-allocating memory space for states variables
    X  = zeros(8, N+1);
    
    % Setting the initial coditions for the state variable
    P1=218541216;  
    X(:, 1) = P1*[0.605389443, 0.000147656,	0.125507567, 0.059062385,...
                  0.199335548, 0.002657807, 0.000516796, 0.007382798]; 
    
    % Pre-allocating memory space for adjoint variables
    Lambda = zeros(8, N+1);
    
    Pi = 13920; gamma_1= 0.04715; gamma_2= 0.04715;  psi_1  = 0.1429; 
    psi_2  = 0.1429; 
    varepsilon = 0.02726; alpha = 0.8973; d_1 =   0.002843; 
    d_2= 0.003310; d_12= 0.004131; tau_1= 0.2954; tau_2= 0.0096;
    d_21= 0.004131; sigma  = 0.4685;  beta_1  = 4.6935;  beta_2  = 4.1196; 
    epsilon_1=  0.04750; epsilon_2=  0.3715; mu= 0.0015; rho=0.2;  psi_3  = 0.0429;
   % beta_1  = 10.6935;  beta_2  = 10.1196; % this is just for visualization
    
   ko1 = gamma_1 + mu; k22 = tau_1 + d_12 + psi_1 + mu; k3 = gamma_2 + mu;
   k4 = tau_2 + d_21 + psi_2 + mu ; %k5 = psi_3 + delta_21 + delta_12 + mu;
    
    R_10 = beta_1*gamma_1*(rho*sigma + alpha + mu)*(tau_1.*epsilon_1 +...
         (psi_3  + d_12 + mu))/(ko1*k22*(psi_3 + d_12 + mu)*(alpha + mu + rho));
    R_20 = beta_2.*gamma_2.*(rho.*sigma + alpha + mu).*(tau_2.*epsilon_2 +...
         (psi_3  + d_21 + mu))./(k3.*k4.*(psi_3 + d_21 + mu).*(alpha + mu + rho));

       
    %Define your Weight
   
     %w1=10; w2=20; w3=20; w4=40; kappa1=6000; kappa2=6000; kappa3=6000;% kappa4 = 5000; %kappa5 = 6000; % kappa are Cost function
     w1=10; w2=20; w3=20; w4=40; kappa1=1200; kappa2=9000; kappa3=12000;% kappa4 = 5000; %kappa5 = 6000; % kappa are Cost function
    function dy = dxdt(t, x, T, U)
        u = interp1(T, U', t)'; 
        u1 = u(1); u2 = u(2); u3 = u(3);% u4 = u(4); u5 = u(5);
        S = x(1); V = x(2); E_1 = x(3); I_1 = x(4); E_2 = x(5); I_2 = x(6); H = x(7); R = x(8);
        Total = S + V + E_1 + I_1 + E_2 + I_2 + H + R;
        dy = [Pi + varepsilon*R - ((beta_1*(I_1 + epsilon_1*H))./Total + (beta_2*(I_2 + epsilon_2*H)./Total))*S*(1-u2) + alpha*V - (u1 + mu)*S;
              u1*S - ((beta_1*(I_1 + epsilon_1*H))./Total + (beta_2*(I_2 + epsilon_2*H)./Total))*(V*sigma)*(1-u2) - (alpha + mu)*V;
              (S + sigma*V)*(1-u2)*(beta_1*(I_1 + epsilon_1*H)./Total) - (gamma_1 + mu)*E_1;
              gamma_1*E_1 - (tau_1 + d_1 + psi_1 + mu)*I_1;
              (S + sigma*V)*(beta_2*(I_2 + epsilon_2*H)./Total)*(1-u2) - (gamma_2 + mu)*E_2;
              gamma_2*E_2 - (tau_2 + d_2 + psi_2 + mu)*I_2;
              tau_1*I_1 + tau_2*I_2 - (u3 + mu + d_12 + d_21)*H;
              psi_1*I_1 + psi_2*I_2 + u3*H - (varepsilon + mu)*R];
    end

    function dlambda = dldt(t, lambda, T, U, X)
        x = interp1(T, X', t)';
        S = x(1); V = x(2); E_1 = x(3); I_1 = x(4); E_2 = x(5); I_2 = x(6); H = x(7); R = x(8);
        Total = S + V + E_1 + I_1 + E_2 + I_2 + H + R;
        u = interp1(T, U', t)'; 
        u1 = u(1); u2 = u(2); u3 = u(3); % u4 = u(4); u5 = u(5);
        lambda1 = lambda(1); lambda2 = lambda(2); lambda3 = lambda(3); lambda4 = lambda(4); 
        lambda5 = lambda(5); lambda6 = lambda(6); lambda7 = lambda(7); lambda8 = lambda(8);
        dlambda = [-(w1 + lambda1*(-(-beta_1*(H*epsilon_1 + I_1)./Total^2 -beta_2*(H*epsilon_2 + I_2)./Total^2)*S*(1-u2) ... 
                   - (beta_1*(H*epsilon_1+I_1)./Total + beta_2*(H*epsilon_2 + I_2)./Total)*(1-u2)  - u1 - mu) ... 
                   + lambda2*(u1-(1-u2)*sigma*(-beta_1*(H*epsilon_1 + I_1)./Total^2 - beta_2*(H*epsilon_2 + I_2)./Total^2)*V)...
                   + lambda3*((1-u2)*beta_1*(H*epsilon_1+I_1)./Total - (V*sigma + S)*(1-u2)*beta_1*(H*epsilon_1 + I_1)./Total^2)...
                   + lambda5*((1-u2)*beta_2*(H*epsilon_2 + I_2)./Total - (V*sigma + S)*(1-u2)*beta_2*(H*epsilon_2 + I_2)./Total^2)); 

                   -(lambda1*(-(-beta_1*(H*epsilon_1 + I_1)./Total^2 - beta_2*(H*epsilon_2+I_2)./Total^2)*S*(1-u2) + alpha)...
                   + lambda2*(-sigma*(1-u2)*(-beta_1*(H*epsilon_1+I_1)./Total^2 - beta_2*(H*epsilon_2+I_2)./Total^2)*V ...
                   - sigma*(1-u2)*(beta_1*(H*epsilon_1+I_1)./Total + beta_2*(H*epsilon_2+I_2)./Total) - alpha - mu) ...
                   + lambda3*(sigma*(1-u2)*beta_1*(H*epsilon_1+I_1)./Total - (V*sigma +S)*(1-u2)*beta_1*(H*epsilon_1+I_1)./Total^2) ...
                   + lambda5*((1-u2)*sigma*beta_2*(H*epsilon_2 + I_2)./Total - (V*sigma + S)*(1-u2)*beta_2*(H*epsilon_2+I_2)./Total^2));

                   -(-lambda1*(-beta_1*(H*epsilon_1+I_1)./Total^2 - beta_2*(H*epsilon_2+I_2)./Total^2)*S*(1-u2)...
                   - lambda2*sigma*(1-u2)*(-beta_1*(H*epsilon_1+I_1)./Total^2 - beta_2*(H*epsilon_2+I_2)./Total^2)*V ...
                   + lambda3*(-(V*sigma + S)*(1-u2)*beta_1*(H*epsilon_1+I_1)./Total^2 - gamma_1- mu) + lambda4*gamma_1...
                   - lambda5*(V*sigma + S)*(1-u2)*beta_2*(H*epsilon_2 + I_2)./Total^2);

                   -(w2 - lambda1*(beta_1./Total - beta_1*(H*epsilon_1 + I_1)./Total^2 - beta_2*(H*epsilon_2 + I_2)./Total^2)*S*(1-u2) ...
                   - lambda2*sigma*(1-u2)*(beta_1./Total - beta_1*(H*epsilon_1+I_1)./Total^2 -beta_2*(H*epsilon_2 + I_2)./Total^2)*V ...
                   + lambda3*((V*sigma + S)*(1-u2)*beta_1./Total - (V*sigma + S)*(1-u2)*beta_1*(H*epsilon_1 + I_1)./Total^2) ...
                   + lambda4*(- tau_1 - d_1 - psi_1- mu) - lambda5*(V*sigma + S)*(1-u2)*beta_2*(H*epsilon_2 + I_2)./Total^2 ...
                   + lambda7*tau_1 + lambda8*psi_1);

                   -(-lambda1*(-beta_1*(H*epsilon_1+I_1)./Total^2 - beta_2*(H*epsilon_2+I_2)./Total^2)*S*(1-u2) ...
                   - lambda2*sigma*(1-u2)*(-beta_1*(H*epsilon_1+I_1)./Total^2 - beta_2*(H*epsilon_2+I_2)./Total^2)*V ...
                   - lambda3*(V*sigma + S)*(1-u2)*beta_1*(H*epsilon_1+I_1)./Total^2 ...
                   + lambda5*(-(V*sigma + S)*(1-u2)*beta_2*(H*epsilon_2+I_2)./Total^2 - gamma_2-mu) + lambda6*gamma_2);

                   -(w3-lambda1*(-beta_1*(H*epsilon_1+I_1)./Total^2 + beta_2./Total - beta_2*(H*epsilon_2+I_2)./Total^2)*S*(1-u2) ...
                   - lambda2*sigma*(1-u2)*(-beta_1*(H*epsilon_1+I_1)./Total^2 + beta_2./Total-beta_2*(H*epsilon_2+I_2)./Total^2)*V ...
                   -lambda3*(V*sigma + S)*(1-u2)*beta_1*(H*epsilon_1+I_1)./Total^2 + lambda5*((V*sigma + S)*(1-u2)*beta_2./Total ...
                   - (V*sigma + S)*(1-u2)*beta_2*(H*epsilon_2 + I_2)./Total^2) + lambda6*(-tau_2 - d_2 - psi_2 - mu) + lambda7*tau_2 + lambda8*psi_2);

                   -(w4-lambda1*(beta_1*epsilon_1./Total - beta_1*(H*epsilon_1+I_1)./Total^2 + beta_2*epsilon_2./Total ...
                   - beta_2*(H*epsilon_2+I_2)./Total^2)*S*(1-u2)-lambda2*sigma*(1-u2)*(beta_1*epsilon_1./Total - beta_1*(H*epsilon_1+I_1)./Total^2 ...
                   + beta_2*epsilon_2./Total - beta_2*(H*epsilon_2 + I_2)./Total^2)*V ... 
                   + lambda3*((V*sigma + S)*(1-u2)*beta_1*epsilon_1./Total-(V*sigma + S)*(1-u2)*beta_1*(H*epsilon_1+I_1)./Total^2) ...
                   + lambda5*((V*sigma + S)*(1-u2)*beta_2*epsilon_2./Total - (V*sigma + S)*(1-u2)*beta_2*(H*epsilon_2+I_2)./Total^2)...
                   + lambda7*(- u3 - mu -d_12 - d_21) + lambda8*u3);

                   -(lambda1*(varepsilon -(-beta_1*(H*epsilon_1+I_1)./Total^2 - beta_2*(H*epsilon_2+I_2)./Total^2)*S*(1-u2)) ...
                   - lambda2*sigma*(1-u2)*(-beta_1*(H*epsilon_1+I_1)./Total^2-beta_2*(H*epsilon_2+I_2)./Total^2)*V ...
                   -lambda3*(V*sigma + S)*(1-u2)*beta_1*(H*epsilon_1+I_1)./Total^2 -lambda5*(V*sigma + S)*(1-u2)*beta_2*(H*epsilon_2+I_2)./Total^2 ...
                   + lambda8*(-varepsilon -mu))];
    
    end

    iter = 0;
    while(test<0 && iter < 20) % while the tolerance is reached, repeat
        % This step stores current values  as previous ones
        % Storing the preceeding values of the control variables
        Old_U = U;

        % Storing the preceeding values of the state variables
        Old_X = X;
    
        % Storing the preceeding values of the adjoint variables
        Old_Lambda = Lambda;

        % forward sweep with RK order 4  to evaluate the state variables
        for i = 1:N
            X(:, i+1) = rk4(@(t, x)dxdt(t, x, T, U), (i-1)*k , X(:, i), k);
        end
    
        % Backward sweep with RK order 4
        for j = N:-1:1
            Lambda(:, j) =  rk4(@(t, y)dldt(t, y, T, U, X), j*k , Lambda(:, j + 1), -k);
        end

       SS = X(1,:); VV= X(2,:); II_1 = X(4,:); II_2 = X(6,:); HH = X(7,:);
       EE_1 = X(3,:); EE_2 = X(5,:); RR = X(8, :);
       L1= VV.*(Lambda(2,:)- Lambda(3,:))*sigma - SS.*(Lambda(3,:) - Lambda(1,:));
       L2= VV.*(Lambda(2,:) - Lambda(5,:))*sigma - SS.*(Lambda(5,:) - Lambda(1,:));
       Zeta= (HH.*epsilon_1 + II_1)*beta_1.*L1 + (HH.*epsilon_2 + II_2)*beta_2.*L2;
       NN = SS + VV + EE_1 + II_1 + II_2 + EE_2 + HH + RR;

        U(1, :) = max(0, min(SS.*(Lambda(1,:) - Lambda(2,:))./kappa1, u1max));
        U(2, :) = max(0, min(Zeta./(kappa2*NN), u2max));
        U(3, :) = max(0, min(HH.*(Lambda(7,:) - Lambda(8,:))./kappa3, u3max));

        U = 0.5*(control_switch.*(U + Old_U));
        tempX      = delta*sum(abs(X), 2) - sum(abs(Old_X - X), 2);
        tempLambda = delta*sum(abs(Lambda), 2) - sum(abs(Old_Lambda - Lambda), 2);
        tempU      = delta*sum(abs(U), 2) - sum(abs(Old_U - U), 2);
        test       = min([min(tempX),  min(tempLambda), min(tempU)]);
        iter = iter + 1;
    end 
    flag =  iter < 20;
end

function ynp1 = rk4(dydt, tn, yn, dt)
    k1 = dydt(tn, yn);
    k2 = dydt(tn + 0.5*dt, yn + 0.5 * k1 * dt);
    k3 = dydt(tn + 0.5*dt, yn + 0.5 * k2 * dt);
    k4 = dydt(tn + dt, yn + k3 * dt);
    ynp1 = yn + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end


