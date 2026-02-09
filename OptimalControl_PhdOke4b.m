[T1, X1, U1, Lambda1] = OptimalControl_PhdOke4a([0;0;0]);
[T2, X2, U2, Lambda2] = OptimalControl_PhdOke4a([1;1;1]);


% Label = ["S","V","E_1","I_1","E_2","I_2","H","R"];
% for i = 1:8
%    subplot(2,2,i)
%    plot(T1, X1(i,:),"r", T2, X2(i,:),"b",'LineWidth', 2);
%   % ylabel(Label(i))
%    title(['Vaccination Only: ', title_suffix]);
%     xlabel('Time')
%     ylabel(Label(i))
%     legend('No Control', 'Control'); grid on;
% end
beta_1  = 0.3935;  beta_2  = 0.1196; epsilon_1=  0.04750; epsilon_2=  0.3715;

% figure;
% plot(T1, (beta_1* (X1(4,:) + epsilon_1* X1(7,:))) .* X1(1,:) ./ ...
%     (X1(1,:) + X1(2,:) + X1(3,:) + X1(4,:) + X1(7,:) + X1(8,:)), 'r', ...
%      T1, (beta_2* (X1(6,:) + epsilon_2 * X1(7,:))) .* X1(1,:) ./ ...
%     (X1(1,:) + X1(2,:) + X1(5,:) + X1(6,:) + X1(7,:) + X1(8,:)),'k', 'LineWidth', 2);
% xlabel('Time');
% ylabel('Incedence');
% grid on;
% 
% figure;
% plot(T2, (beta_1* (X2(4,:) + epsilon_1* X2(7,:))) .* X2(1,:) ./ ...
%     (X2(1,:) + X2(2,:) + X2(3,:) + X2(4,:) + X2(7,:) + X2(8,:)), '--r', ...
%      T2, (beta_2* (X2(6,:) + epsilon_2 * X2(7,:))) .* X2(1,:) ./ ...
%     (X2(1,:) + X2(2,:) + X2(5,:) + X2(6,:) + X2(7,:) + X2(8,:)),'--k', 'LineWidth', 2);
% xlabel('Time');
% ylabel('Incedence');
% grid on;
% 
% 
% 
Label = ["S","V","E_1","I_1","E_2","I_2","H","R"]; 

% First figure: Plot from S to I_1
figure;
for i = 1:4
    subplot(2,2,i)
    plot(T1, X1(i,:), "r", T2, X2(i,:), "b", 'LineWidth', 2);
    %title('All Control: ');
    xlabel('Time');
    ylabel(Label(i));
    legend('No Control', 'Control');
    grid on;
end

% Second figure: Plot from E_2 to R
figure;
for i = 5:8
    subplot(2,2,i-4) % Adjust index for subplot
    plot(T1, X1(i,:), "r", T2, X2(i,:), "b", 'LineWidth', 2);
    %title('All Control:')
    xlabel('Time');
    ylabel(Label(i));
    legend('No Control', 'Control');
    grid on;
end

% % Plot of prevalence
%  figure
%  plot(T1, (X1(4,:) + X1(6,:) + X1(7,:)),"-k", T2, (X2(4,:) + X2(6,:) + X2(7,:)), '-r',LineWidth=4)
%  ylabel('Disease Prevalence')
%  xlabel('Time (Days)')
%  legend('$\bf{u_3 =  u_2 = u_3 = 0}$','interpreter','latex')
% 
% figure
% % for single control
%  plot(T2, U2(1,:),"-b" ,LineWidth=4)
%  ylabel('Control Profile')
%  xlabel('Time (Days)')
%  legend('$u_2$','interpreter','latex')
% % 
% % % 
%  figure
% %for double contol
% plot(T2, U2(2,:),"-b" ,T2, U2(3,:), "-r", LineWidth=4)
%  ylabel('Control Profile')
%  xlabel('Time (Days)')
%  legend('$u_3$','interpreter','latex')

%Strategy 3 (Triple  Control)
%  figure
% plot(T2, U2(1,:),"-b" ,T2, U2(2,:), "-r", T2, U2(3,:),"-k" ,LineWidth=4)
%  ylabel('Control Profile')
%  xlabel('Time (Days)')
%  legend('$u_3$','interpreter','latex')


%w1=10; w2=10; w3=10; w4=10; kappa1=6000; kappa2=6000; kappa3=6000;% kappa4 = 5000; %kappa5 = 6000; % kappa are Cost function
w1=20; w2=50; w3=50; w4=20; kappa1=12000; kappa2=90000; kappa3=10000;% kappa4 = 5000; %kappa5 = 6000; % kappa are Cost function
A= sum(X1(4,:) + X1(6,:) + X1(7,:))
B= sum(X2(4,:) + X2(6,:) + X2(7,:))
TIA= A-B
EI= (TIA/A)*100
TC= 1/2*sum(kappa1*U2(1,:)) +  1/2*sum(kappa2*U2(2,:))  +  1/2*sum(kappa3*U2(3,:)) % +  1/2*sum(kappa4*U2(4,:)) % + 1/2*sum(kappa5*U2(5,:))
ICER = TC/TIA