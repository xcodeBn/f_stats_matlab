
clear;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Field Dimensions - Maximum x and y dimensions of the network area (in meters)
xm = 100; % x dimension
ym = 100; % y dimension

% Coordinates of the Sink - Position of the base station
sink.x = 0.5 * xm;
sink.y = 0.5 * ym;

% Number of Nodes in the field
n = 100;

% Optimal Election Probability of a node to become cluster head
p = 0.1;

% Energy Model (all values in Joules)
% Initial Energy of nodes
Eo = 0.5;
% Energy consumption parameters for various operations
ETX = 50 * 0.000000001; % Energy consumed during transmission
ERX = 50 * 0.000000001; % Energy consumed during reception
Efs = 10 * 0.000000000001; % Energy for free space model
Emp = 0.0013 * 0.000000000001; % Energy for multi-path fading model
EDA = 5 * 0.000000001; % Energy consumed during data aggregation

% Values for Heterogeneity 
% Percentage of nodes that are advanced
m = 0.1;
% Parameter alpha for energy model
a = 1;

% Maximum number of rounds for the simulation
rmax = 200*5;

%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%

% Computation of do - Calculating the optimal distance for energy calculations
do = sqrt(Efs / Emp);

% Create figure windows
figure(1); % Network topology
figure(2); % Energy statistics
figure(3); % Dead nodes statistics
figure(4); % Throughput statistics
figure(5); % Network lifetime metrics

% Creation of the random Sensor Network
figure(1);
title('Wireless Sensor Network Model');
xlabel('X-coordinate (m)');
ylabel('Y-coordinate (m)');
grid on;

% Arrays to store node energy levels for each round
NODE_ENERGY = zeros(n, rmax+1);

% Arrays to store statistics
DEAD = zeros(1, rmax+1);
DEAD_N = zeros(1, rmax+1);
DEAD_A = zeros(1, rmax+1);
CLUSTERHS = zeros(1, rmax+1);
PACKETS_TO_BS = zeros(1, rmax+1);
PACKETS_TO_CH = zeros(1, rmax+1);
AVG_ENERGY = zeros(1, rmax+1);
ENERGY_VARIANCE = zeros(1, rmax+1);
TOTAL_ENERGY = zeros(1, rmax+1);

% Initialize nodes
for i = 1:1:n
    S(i).xd = rand(1, 1) * xm;
    XR(i) = S(i).xd;
    S(i).yd = rand(1, 1) * ym;
    YR(i) = S(i).yd;
    S(i).G = 0;
    % Initially there are no cluster heads only nodes
    S(i).type = 'N';
   
    temp_rnd0 = i;
    % Random Election of Normal Nodes
    if (temp_rnd0 >= m * n + 1) 
        S(i).E = Eo;
        S(i).ENERGY = 0;
        plot(S(i).xd, S(i).yd, 'bo', 'MarkerSize', 6);
        hold on;
    end
    % Random Election of Advanced Nodes
    if (temp_rnd0 < m * n + 1)  
        S(i).E = Eo * (1 + a);
        S(i).ENERGY = 1;
        plot(S(i).xd, S(i).yd, 'g+', 'MarkerSize', 8);
        hold on;
    end
    
    % Initialize node energy for round 0
    NODE_ENERGY(i, 1) = S(i).E;
end

% Plot the base station
S(n + 1).xd = sink.x;
S(n + 1).yd = sink.y;
plot(S(n + 1).xd, S(n + 1).yd, 'rx', 'MarkerSize', 10, 'LineWidth', 3);
legend('Normal Nodes', 'Advanced Nodes', 'Base Station');

% Calculate initial total energy
for i=1:n
    TOTAL_ENERGY(1) = TOTAL_ENERGY(1) + S(i).E;
end
AVG_ENERGY(1) = TOTAL_ENERGY(1)/n;

% First Iteration - Initializing the first iteration
% Counter for CHs - Counts the number of cluster heads
countCHs = 0;
% Counter for CHs per round - Counts the number of cluster heads per round
rcountCHs = 0;
% Initializing cluster counter
cluster = 1;
% Initialize counters
countCHs;
rcountCHs = rcountCHs + countCHs;
% Flag for the first dead node
flag_first_dead = 0;
first_dead = 0;

% Simulation starts here
   for r = 0:1:rmax % Loop for rounds
    % Display current round
    fprintf('Round: %d\n', r);
    
    % Operation for epoch - Updating node status for each epoch
    if (mod(r, round(1/p)) == 0)
        for i = 1:1:n
            S(i).G = 0;
            S(i).cl = 0;
        end
    end
   end

    % Reset figure for network visualization
    figure(1);
    clf;
    title(['Wireless Sensor Network - Round: ', num2str(r)]);
    xlabel('X-coordinate (m)');
    ylabel('Y-coordinate (m)');
    grid on;
    hold on;

    % Number of dead nodes - Counts the total number of dead nodes
    dead = 0;
    % Number of dead Advanced Nodes - Counts the number of dead advanced nodes
    dead_a = 0;
    % Number of dead Normal Nodes - Counts the number of dead normal nodes
    dead_n = 0;

    % Counter for bits transmitted to Base Station and to Cluster Heads
    packets_TO_BS = 0;
    packets_TO_CH = 0;

    % Counter for bits transmitted to Base Station and to Cluster Heads per round
    PACKETS_TO_CH(r + 1) = 0;
    PACKETS_TO_BS(r + 1) = 0;

    % Calculate current energy status and plot nodes
    for i = 1:1:n
        % Checking if there is a dead node
        if (S(i).E <= 0)
            plot(S(i).xd, S(i).yd, 'r.', 'MarkerSize', 10);
            dead = dead + 1;
            if (S(i).ENERGY == 1)
                dead_a = dead_a + 1;
            end
            if (S(i).ENERGY == 0)
                dead_n = dead_n + 1;
            end
        else
            % If node is alive
            S(i).type = 'N';
            if (S(i).ENERGY == 0)  
                plot(S(i).xd, S(i).yd, 'bo', 'MarkerSize', 6);
            end
            if (S(i).ENERGY == 1)  
                plot(S(i).xd, S(i).yd, 'g+', 'MarkerSize', 8);
            end
        end
        
        % Store energy level of node for this round
        NODE_ENERGY(i, r+1) = S(i).E;
    end
    
    % Plot the base station
    plot(S(n + 1).xd, S(n + 1).yd, 'rx', 'MarkerSize', 10, 'LineWidth', 3);
    
    % Recording statistics about dead nodes
    STATISTICS(r + 1).DEAD = dead;
    DEAD(r + 1) = dead;
    DEAD_N(r + 1) = dead_n;
    DEAD_A(r + 1) = dead_a;

    % Calculate average energy and variance
    total_energy = 0;
    sum_squared_diff = 0;
    alive_nodes = 0;
    
    for i = 1:n
        if S(i).E > 0
            total_energy = total_energy + S(i).E;
            alive_nodes = alive_nodes + 1;
        end
    end
    
    if alive_nodes > 0
        avg_energy = total_energy / alive_nodes;
        
        for i = 1:n
            if S(i).E > 0
                sum_squared_diff = sum_squared_diff + (S(i).E - avg_energy)^2;
            end
        end
        
        energy_variance = sum_squared_diff / alive_nodes;
    else
        avg_energy = 0;
        energy_variance = 0;
    end
    
    TOTAL_ENERGY(r+1) = total_energy;
    AVG_ENERGY(r+1) = avg_energy;
    ENERGY_VARIANCE(r+1) = energy_variance;

    % When the first node dies
    if (dead > 0 && flag_first_dead == 0)
        first_dead = r;
        flag_first_dead = 1;
    end

    % Cluster formation
    countCHs = 0;
    cluster = 1;
    
    % Reset cluster heads array
    if r > 0
        clear C;
        clear X;
        clear Y;
    end
    
    % Cluster head selection
    for i = 1:1:n
        if (S(i).E > 0)
            temp_rand = rand;     
            if ((S(i).G) <= 0)
                % Election of Cluster Heads
                if (temp_rand <= (p / (1 - p * mod(r, round(1 / p)))))
                    countCHs = countCHs + 1;
                    packets_TO_BS = packets_TO_BS + 1;
                    PACKETS_TO_BS(r + 1) = packets_TO_BS;
                    
                    S(i).type = 'C';
                    S(i).G = round(1 / p) - 1;
                    C(cluster).xd = S(i).xd;
                    C(cluster).yd = S(i).yd;
                    plot(S(i).xd, S(i).yd, 'k*', 'MarkerSize', 10);
                    
                    distance = sqrt((S(i).xd - (S(n + 1).xd))^2 + (S(i).yd - (S(n + 1).yd))^2);
                    C(cluster).distance = distance;
                    C(cluster).id = i;
                    X(cluster) = S(i).xd;
                    Y(cluster) = S(i).yd;
                    cluster = cluster + 1;
                    
                    % Calculation of Energy dissipated
                    if (distance > do)
                        S(i).E = S(i).E - ((ETX + EDA) * (4000) + Emp * 4000 * (distance^4)); 
                    else
                        S(i).E = S(i).E - ((ETX + EDA) * (4000) + Efs * 4000 * (distance^2)); 
                    end
                end     
            end
        end
    end

    % Update legend if there are cluster heads
    if countCHs > 0
        legend('Normal Nodes', 'Advanced Nodes', 'Dead Nodes', 'Base Station', 'Cluster Heads');
    else
        legend('Normal Nodes', 'Advanced Nodes', 'Dead Nodes', 'Base Station');
    end

    STATISTICS(r + 1).CLUSTERHEADS = cluster - 1;
    CLUSTERHS(r + 1) = cluster - 1;

    % Cluster association - connecting nodes to their nearest cluster head
    for i = 1:1:n
        if (S(i).type == 'N' && S(i).E > 0)
            if (cluster - 1 >= 1)
                min_dis = sqrt((S(i).xd - S(n + 1).xd)^2 + (S(i).yd - S(n + 1).yd)^2);
                min_dis_cluster = 0;
                
                for c = 1:1:cluster - 1
                    temp = sqrt((S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2);
                    if (temp < min_dis)
                        min_dis = temp;
                        min_dis_cluster = c;
                    end
                end

                % If a cluster head was found
                if min_dis_cluster > 0
                    % Draw a line from node to its cluster head
                    line([S(i).xd, C(min_dis_cluster).xd], [S(i).yd, C(min_dis_cluster).yd], 'Color', 'c', 'LineStyle', ':');
                    
                    % Energy dissipated by associated Cluster Head
                    if (min_dis > do)
                        S(i).E = S(i).E - (ETX * (4000) + Emp * 4000 * (min_dis^4)); 
                    else
                        S(i).E = S(i).E - (ETX * (4000) + Efs * 4000 * (min_dis^2));

                  end
                    
                    % Energy dissipated at cluster head
                    S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ((ERX + EDA) * 4000); 
                    PACKETS_TO_CH(r + 1) = PACKETS_TO_CH(r + 1) + 1;
                    
                    S(i).min_dis = min_dis;
                    S(i).min_dis_cluster = min_dis_cluster;
                else
                    % If no cluster head found, send directly to BS
                    dist_to_bs = sqrt((S(i).xd - S(n + 1).xd)^2 + (S(i).yd - S(n + 1).yd)^2);
                    line([S(i).xd, S(n + 1).xd], [S(i).yd, S(n + 1).yd], 'Color', 'm', 'LineStyle', '-.');
                    
                    if (dist_to_bs > do)
                        S(i).E = S(i).E - (ETX * (4000) + Emp * 4000 * (dist_to_bs^4)); 
                    else
                        S(i).E = S(i).E - (ETX * (4000) + Efs * 4000 * (dist_to_bs^2)); 
                    end
                    
                    packets_TO_BS = packets_TO_BS + 1;
                    PACKETS_TO_BS(r + 1) = packets_TO_BS;
                
                
                % No cluster heads, send directly to BS
                dist_to_bs = sqrt((S(i).xd - S(n + 1).xd)^2 + (S(i).yd - S(n + 1).yd)^2);
                line([S(i).xd, S(n + 1).xd], [S(i).yd, S(n + 1).yd], 'Color', 'm', 'LineStyle', '-.');
                
                if (dist_to_bs > do)
                    S(i).E = S(i).E - (ETX * (4000) + Emp * 4000 * (dist_to_bs^4)); 
                else
                    S(i).E = S(i).E - (ETX * (4000) + Efs * 4000 * (dist_to_bs^2)); 
                end
                
                packets_TO_BS = packets_TO_BS + 1;
                PACKETS_TO_BS(r + 1) = packets_TO_BS;
            end
        end
    end
    
    % Update round statistics
    countCHs;
    rcountCHs = rcountCHs + countCHs;
    
    % Plot energy statistics after each round
    figure(2);
    subplot(2,2,1);
    plot(0:r, TOTAL_ENERGY(1:r+1), 'b-', 'LineWidth', 2);
    title('Total Remaining Energy');
    xlabel('Round');
    ylabel('Energy (J)');
    grid on;
    
    subplot(2,2,2);
    plot(0:r, AVG_ENERGY(1:r+1), 'g-', 'LineWidth', 2);
    title('Average Energy per Node');
    xlabel('Round');
    ylabel('Energy (J)');
    grid on;
    
    subplot(2,2,3);
    plot(0:r, ENERGY_VARIANCE(1:r+1), 'r-', 'LineWidth', 2);
    title('Energy Variance');
    xlabel('Round');
    ylabel('Variance');
    grid on;
    
    subplot(2,2,4);
    % Plot energy level of 5 random nodes
    node_indices = min(n, 5);
    for i = 1:node_indices
        plot(0:r, NODE_ENERGY(i, 1:r+1), 'LineWidth', 1.5);
        hold on;
    end
    title('Energy Level of Sample Nodes');
    xlabel('Round');
    ylabel('Energy (J)');
    if node_indices == 5
        legend('Node 1', 'Node 2', 'Node 3', 'Node 4', 'Node 5');
    end
    grid on;
    hold off;
    
    % Plot dead nodes statistics
    figure(3);
    subplot(2,1,1);
    plot(0:r, DEAD(1:r+1), 'r-', 'LineWidth', 2);
    hold on;
    plot(0:r, DEAD_N(1:r+1), 'b--', 'LineWidth', 2);
    plot(0:r, DEAD_A(1:r+1), 'g-.', 'LineWidth', 2);
    title('Number of Dead Nodes');
    xlabel('Round');
    ylabel('Number of Nodes');
    legend('Total Dead', 'Normal Dead', 'Advanced Dead');
    grid on;
    hold off;
    
    subplot(2,1,2);
    alive = n - DEAD(1:r+1);
    plot(0:r, alive, 'g-', 'LineWidth', 2);
    title('Number of Alive Nodes');
    xlabel('Round');
    ylabel('Number of Nodes');
    grid on;
    
    % Plot throughput statistics
    figure(4);
    subplot(2,1,1);
    plot(0:r, PACKETS_TO_BS(1:r+1), 'b-', 'LineWidth', 2);
    title('Packets Sent to Base Station');
    xlabel('Round');
    ylabel('Number of Packets');
    grid on;
    
    subplot(2,1,2);
    plot(0:r, PACKETS_TO_CH(1:r+1), 'r-', 'LineWidth', 2);
    hold on;
    plot(0:r, CLUSTERHS(1:r+1), 'g--', 'LineWidth', 2);
    title('Packets Sent to CHs and Number of CHs');
    xlabel('Round');
    ylabel('Count');
    legend('Packets to CHs', 'Number of CHs');
    grid on;
    hold off;
    
    % Check if all nodes are dead
    if dead >= n
        fprintf('All nodes died at round %d\n', r);
        break;
    end
    
    % Pause briefly to visualize
    if mod(r, 5) == 0
        pause(0.01);
    end
end

% Final Network Lifetime Metrics
figure(5);
bar([first_dead, r], [first_dead, r-first_dead], 'stacked');
title('Network Lifetime Metrics');
xlabel('Rounds');
ylabel('Duration');
grid on;
legend('First Node Death', 'Network Operation After First Death');
text(0.5, first_dead/2, ['FND: ', num2str(first_dead)], 'HorizontalAlignment', 'center');
text(1.5, (r-first_dead)/2 + first_dead, ['LND: ', num2str(r)], 'HorizontalAlignment', 'center');

% Print final statistics
fprintf('\n\nFinal Statistics:\n');
fprintf('First Node Died at Round: %d\n', first_dead);
fprintf('Last Node Died at Round: %d\n', r);
fprintf('Network Lifetime: %d rounds\n', r);
fprintf('Total Packets to BS: %d\n', PACKETS_TO_BS(end));
fprintf('Total Packets to CHs: %d\n', sum(PACKETS_TO_CH));
fprintf('Average Number of Cluster Heads per Round: %.2f\n', mean(CLUSTERHS(1:r+1)));

% Energy Consumption Analysis
figure(6);
% Create a 3D visualization of energy consumption
[X, Y] = meshgrid(1:10:n, 0:ceil(r/10):r);
Z = zeros(size(X));

for i = 1:size(X, 1)
    for j = 1:size(X, 2)
        node_idx = X(i, j);
        round_idx = Y(i, j) + 1;
        if node_idx <= n && round_idx <= r+1
            Z(i, j) = NODE_ENERGY(node_idx, round_idx);
        end
    end
end

surf(X, Y, Z);
title('3D Energy Consumption Map');
xlabel('Node ID');
ylabel('Round');
zlabel('Energy Level (J)');
colorbar;

% Save all figures
for i = 1:6
    figure(i);
    saveas(gcf, sprintf('leach_figure_%d.png', i));
end
