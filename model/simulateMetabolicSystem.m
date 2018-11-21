function [T, F] = simulateMetabolicSystem
% This script defines a kinetic model for the simple network used in the
% CompSysBio course and allows the network to be simulated using standard numerical solver routines in Matlab.
% The model does not contain any regulatory interactions.

% HdJ 16/3/17

    clear all;
    close all;
    format long;

    % Define model
    model = @metabolicModel;
    
    % Define parameters
    K1 = 50;
    K2 = 50;
    k1 = 10;
    k2 = 5;
    k3 = 10;
    k4 = 5;
    k5 = 5;
    beta = 1; 
    
    p = [K1, K2, k1, k2, k3, k4, k5, beta];

    
    % Define time interval
    tspan = [0:0.1:6];

    
    % Define initial conditions
    s1_0 = 100;
    s2_0 = 80;
    x1_0 = 0;
    x2_0 = 0;
    m_0 = 0;
    b_0 = 0.1;
    
    x0 = [s1_0, s2_0, x1_0, x2_0, m_0, b_0];
 
    % Simulation options. 
    %odeOptions = odeset('AbsTol',absMatrix,'RelTol',.001);
    odeOptions = odeset('RelTol',.01);

    % Run simulations
    [T, X] = ode15s(@(t, x) model(t, x, p), tspan, x0, odeOptions);

    
    % Plot results
    
    figure; 
    subplot(2,2,1); hold on;
    title('Biomass');
    plot(T,X(:,6),'-b');
    subplot(2,2,2); hold on;
    title('Log biomass');
    plot(T,X(:,6),'-b');
    set(gca,'YScale','log');
    subplot(2,2,3); hold on;
    title('Substrates');
    plot(T,X(:,1),'-b'); 
    plot(T,X(:,2),'-r');
    legend('s_1','s_2');
    subplot(2,2,4); hold on;
    title('Metabolites');
    plot(T,X(:,3),'-b');
    plot(T,X(:,4),'-r');
    plot(T,X(:,5),'-g');
    legend('x_1','x_2','m');


    
    

end