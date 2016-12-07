function CME

clear; clc;

k = [40, 40/6.67, 2, 1, 0.4, 0];
P_0 = [1; 0; 0];
n_states = length(P_0);
state_names = {'*','A*','B*'};
r_states = [-k(6), 0, k(5)];

Q = [ -(k(1) + k(6)), k(2), k(5)
    k(1), -(k(2) + k(3)), k(4)
    k(6), k(3), -(k(4) + k(5))];

dQd1 = [ -k(1), k(2), 0
    k(1), -k(2), 0
    0, 0, 0];
dQd2 = [ 0, 0, 0
    0, -k(3), k(4)
    0, k(3), -k(4)];
dQd3 = [ -k(6), 0, k(5)
    0, 0, 0
    k(6), 0, -k(5)];

% Solve for time scales
[~, vals] = eigs(Q);
vals = diag(vals);
vals(abs(vals) < 1e-8) = [];
taus = -1 ./ vals;
tau = max(taus)

% Solve time-dependent problem
t_points = 100;
t_vec = linspace(0,6,t_points);
c_t = zeros(1,t_points);
P_t = zeros(n_states, t_points);
r_t = zeros(1,t_points);
var_r_t = zeros(1,t_points);
NSC = zeros(3,t_points);
for i = 1:length(t_vec)
    
    P_t(:,i) = expm(Q * t_vec(i)) * P_0;
    c_t(i) = exp(-t_vec(i)/tau);
    r_t(i) = r_states * P_t(:,i);
    var_r_t(i) = r_states .^2 * P_t(:,i) - r_t(i)^2;
    
    dP_dtheta = deriv_mat_exp(dQd1, t_vec(i));
    NSC(1,i) = r_states * dP_dtheta / r_t(i); 
    dP_dtheta = deriv_mat_exp(dQd2, t_vec(i));
    NSC(2,i) = r_states * dP_dtheta / r_t(i);
    dP_dtheta = deriv_mat_exp(dQd3, t_vec(i));
    NSC(3,i) = r_states * dP_dtheta / r_t(i) + 1;
    
end

    function dP_dtheta = deriv_mat_exp(dQ, t)       % Computes the derivative of exp(Q)
        
        [alph_vec, Z] = ode15s(@mat_exp_diff, [0, 1], zeros(1,n_states) );
        dP_dtheta = Z(end,:)';
        
        function dzdt = mat_exp_diff(alpha,z)
            dzdt = expm( alpha * Q * t) * (dQ * t) * expm((1-alpha) * Q * t) * P_0;
        end
    end

%% Plot species populations
set(0,'defaultlinelinewidth',1.5)
set(0,'defaultaxeslinewidth',1.5)
figure
hold on
for state = 1:n_states
    plot(t_vec, P_t(state,:));
end
hold off
box('on')
ax = gca;
ax.FontSize = 18;
xlabel('Time (s)','FontSize',18)
ylabel('State prob.','FontSize',18)
legend(state_names);
legend('boxoff');

figure
plot(t_vec, c_t);
% errorbar(t_vec, r_t, sqrt(var_r_t))
box('on')
ax = gca;
ax.FontSize = 18;
xlabel('Time (s)','FontSize',18)
ylabel('Correlation','FontSize',18)

figure
plot(t_vec, r_t);
% errorbar(t_vec, r_t, sqrt(var_r_t))
box('on')
ax = gca;
ax.FontSize = 18;
xlabel('Time (s)','FontSize',18)
ylabel('Rate (1/s)','FontSize',18)

figure
hold on
for state = 1:n_states
    plot(t_vec, NSC(state,:));
end
hold off
box('on')
ax = gca;
ax.FontSize = 18;
xlabel('Time (s)','FontSize',18)
ylabel('NSC','FontSize',18)
legend({'1','2','3'});
legend('boxoff');
disp(NSC(:,end))

end