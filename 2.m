%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                       %
% This is the LEACH [1] code.                                                           %
% The same code can be used for FAIR if m=1.                                            %
%                                                                                       %
% [1] W.R.Heinzelman, A.P.Chandrakasan and H.Balakrishnan,                              %
% An application-specific protocol architecture for wireless microsensor networks       %
% IEEE Transactions on Wireless Communications, 1(4):660-670,2002                       %
%                                                                                       %
% Original Author: Mohammad Hossein Homaei                                              %
%                                           Homaei@ieee.org                             %
% Google Scholar: https://scholar.google.com/citations?user=8IGmFIoAAAAJ&hl=en&oi=ao    %
%                                                                                       %
% Modified by: xcodebn@github                                                           %
% Date: May 8, 2025                                                                     %
% Description: Extended LEACH Protocol Simulation with detailed result tracking.        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc; % Clear command window for cleaner output
close all; % Close all existing figures

%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Field Dimensions - Maximum x and y dimensions of the network area (in meters)
xm = 100; % x dimension of the field
ym = 100; % y dimension of the field

% Coordinates of the Sink - Position of the base station
sink.x = 0.5 * xm;
sink.y = 0.5 * ym;

% Number of Nodes in the field
n = 100;

% Optimal Election Probability of a node to become cluster head
p = 0.1;

% Energy Model (all values in Joules)
% Initial Energy of nodes
Eo = 0.5; % Initial energy for normal nodes
% Energy consumption parameters for various operations
ETX = 50 * 0.000000001;       % Energy consumed during transmission (J/bit)
ERX = 50 * 0.000000001;       % Energy consumed during reception (J/bit)
Efs = 10 * 0.000000000001;    % Amplifier energy for free space model (J/bit/m^2)
Emp = 0.0013 * 0.000000000001; % Amplifier energy for multi-path fading model (J/bit/m^4)
EDA = 5 * 0.000000001;        % Energy consumed during data aggregation (J/bit)

% Data packet size
k_bits = 4000; % Packet size in bits

% Values for Heterogeneity 
m = 0.1;  % Percentage of nodes that are advanced
a = 1;    % Energy factor for advanced nodes (Eo*(1+a))

% Maximum number of rounds for the simulation
rmax = 200;

%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%

% Computation of do - Threshold distance for energy model switching
do = sqrt(Efs / Emp);

% Creation of the random Sensor Network
figure(1);
hold on; % Hold on to plot all nodes and sink

for i = 1:1:n
    S(i).xd = rand(1, 1) * xm;
    S(i).yd = rand(1, 1) * ym;
    S(i).G = 0; % Eligibility counter for CH selection (set to 0 if not eligible in this round)
    S(i).cl = 0; % Cluster ID to which node belongs (0 if not clustered or CH)
    S(i).type = 'N'; % Node type: 'N' for Normal, 'C' for CH, 'D' for Dead
    S(i).num_times_ch = 0; % Counter for how many times a node has been CH

    % Assign initial energy based on node type (normal or advanced)
    if (i > m * n) % Normal Nodes (corrected condition based on typical interpretation)
        S(i).E = Eo;
        S(i).is_advanced = 0; % Flag for normal node
        plot(S(i).xd, S(i).yd, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
    else % Advanced Nodes
        S(i).E = Eo * (1 + a);
        S(i).is_advanced = 1; % Flag for advanced node
        plot(S(i).xd, S(i).yd, 's', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
    end
end
S(n + 1).xd = sink.x; % Sink x-coordinate
S(n + 1).yd = sink.y; % Sink y-coordinate
plot(S(n + 1).xd, S(n + 1).yd, 'x', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'LineWidth', 2);
title('Initial Network Deployment');
xlabel('X-coordinate (m)');
ylabel('Y-coordinate (m)');
legend('Normal Node', 'Advanced Node', 'Base Station');
hold off;

% Initialize Statistics Arrays
STATISTICS = struct(); % Using a structure for some original stats
DEAD = zeros(1, rmax + 1);
DEAD_A = zeros(1, rmax + 1); % Dead Advanced Nodes
DEAD_N = zeros(1, rmax + 1); % Dead Normal Nodes
CLUSTERHS = zeros(1, rmax + 1); % Number of CHs per round

% Renaming original PACKETS_TO_BS to reflect its actual meaning (cumulative CH formations)
CUMULATIVE_CH_FORMATIONS_STAT = zeros(1, rmax + 1);
cumulative_ch_formations_count = 0; % Scalar counter

% Packets from member nodes to CHs per round
PACKETS_TO_CH_ROUND = zeros(1, rmax + 1);

% Initialize New Statistics Arrays
ALIVE_NODES_PER_ROUND = zeros(1, rmax + 1);
NETWORK_ENERGY_PER_ROUND = zeros(1, rmax + 1);
AVG_ENERGY_ALIVE_NODES_PER_ROUND = zeros(1, rmax + 1);
ENERGY_VARIANCE_PER_ROUND = zeros(1, rmax + 1);

ENERGY_MEMBER_TX_TO_CH_PER_ROUND = zeros(1, rmax + 1);
ENERGY_CH_RX_FROM_MEMBERS_PER_ROUND = zeros(1, rmax + 1);
ENERGY_CH_TX_TO_BS_PER_ROUND = zeros(1, rmax + 1);

AVG_NODES_PER_CLUSTER_PER_ROUND = zeros(1, rmax + 1);
DATA_PACKETS_FROM_CH_TO_BS_ROUND = zeros(1, rmax + 1);

CONTROL_PACKETS_ADV_PER_ROUND = zeros(1, rmax + 1);
CONTROL_PACKETS_JOIN_PER_ROUND = zeros(1, rmax + 1);
CONTROL_PACKETS_SCHEDULE_PER_ROUND = zeros(1, rmax + 1);
TOTAL_CONTROL_PACKETS_PER_ROUND = zeros(1, rmax + 1);

% Lifetime metrics
first_dead_round = -1; % Round of First Node Death (FND)
half_dead_round = -1;  % Round of Half Nodes Death (HND)
last_dead_round = -1;  % Round of Last Node Death (LND)
flag_first_dead = 0;
flag_half_dead = 0;
flag_last_dead = 0;


%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION LOOP %%%%%%%%%%%%%%%%%%%%%%%%
disp('Simulation Started...');
for r = 0:1:rmax
    fprintf('Round: %d\n', r);

    % Reset round-specific energy accumulators
    current_round_member_tx_energy = 0;
    current_round_ch_rx_energy = 0;
    current_round_ch_tx_energy = 0;

    % Operation for epoch: Reset G for nodes that have been CH for 1/p rounds
    if (mod(r, round(1/p)) == 0)
        for i = 1:1:n
            S(i).G = 0; % Node is eligible to become CH again
            S(i).cl = 0; % Reset cluster association
        end
    end

    % Update plot for current network state (optional for performance on many rounds)
    % Can be conditioned: if mod(r, plot_interval) == 0
    figure(1); 
    clf; % Clear current figure
    hold on;
    
    plot(S(n+1).xd,S(n+1).yd,'x','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',2); % Plot BS
    title(['Network State at Round ', num2str(r)]);
    xlabel('X-coordinate (m)');
    ylabel('Y-coordinate (m)');

    dead_this_round = 0;
    dead_a_this_round = 0;
    dead_n_this_round = 0;
    
    % Update node status and plot (Dead or Alive)
    for i = 1:1:n
        if (S(i).E <= 0) % Node is dead
            if ~strcmp(S(i).type, 'D') % If not already marked dead
                 S(i).type = 'D';
                 dead_this_round = dead_this_round + 1;
                 if (S(i).is_advanced == 1)
                     dead_a_this_round = dead_a_this_round + 1;
                 else
                     dead_n_this_round = dead_n_this_round + 1;
                 end
            end
            plot(S(i).xd,S(i).yd,'.','MarkerSize',10,'MarkerFaceColor','red','MarkerEdgeColor','red');
        else % Node is alive
            S(i).type = 'N'; % Assume normal node initially for the round
            if (S(i).is_advanced == 0)  
                plot(S(i).xd,S(i).yd,'o','MarkerFaceColor','b','MarkerEdgeColor','b');
            else  
                plot(S(i).xd,S(i).yd,'s','MarkerFaceColor','g','MarkerEdgeColor','g');
            end
        end
        hold on;
    end
    
    % Update total dead node counts
    if r == 0
        DEAD(r+1) = dead_this_round;
        DEAD_A(r+1) = dead_a_this_round;
        DEAD_N(r+1) = dead_n_this_round;
    else
        DEAD(r+1) = DEAD(r) + dead_this_round;
        DEAD_A(r+1) = DEAD_A(r) + dead_a_this_round;
        DEAD_N(r+1) = DEAD_N(r) + dead_n_this_round;
    end

    % Store dead node statistics (using original STATISTICS struct for compatibility if needed)
    STATISTICS(r+1).DEAD = DEAD(r+1);
    
    % Check for FND, HND, LND
    if (DEAD(r+1) > 0 && flag_first_dead == 0)
        first_dead_round = r;
        flag_first_dead = 1;
    end
    if (DEAD(r+1) >= n/2 && flag_half_dead == 0)
        half_dead_round = r;
        flag_half_dead = 1;
    end
    if (DEAD(r+1) == n && flag_last_dead == 0)
        last_dead_round = r;
        flag_last_dead = 1;
        fprintf('All nodes are dead. Ending simulation prematurely at round %d.\n', r);
        rmax = r; % Adjust rmax to end simulation if all nodes die
        % Break from the main loop if needed or let it finish calculations for this round
    end

    % Cluster Head Selection Phase
    num_clusterheads_this_round = 0;
    CHs_this_round = struct('xd',[],'yd',[],'distance',[],'id',[]); % Temporary storage for CHs

    for i = 1:1:n
        if (S(i).E > 0 && S(i).G <= 0) % If node is alive and eligible
            temp_rand = rand;     
            % Probability threshold for becoming a CH
            threshold = p / (1 - p * mod(r, round(1/p)));
            if (temp_rand <= threshold)
                num_clusterheads_this_round = num_clusterheads_this_round + 1;
                
                S(i).type = 'C';
                S(i).num_times_ch = S(i).num_times_ch + 1; % Increment CH counter for this node
                S(i).G = round(1/p) -1; % Node cannot be CH for the next 1/p rounds
                
                CHs_this_round(num_clusterheads_this_round).xd = S(i).xd;
                CHs_this_round(num_clusterheads_this_round).yd = S(i).yd;
                CHs_this_round(num_clusterheads_this_round).id = i;
                plot(S(i).xd, S(i).yd, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 8); % Plot CH
                
                % Energy consumed by CH to send aggregated data to Base Station
                distance_to_bs = sqrt((S(i).xd - S(n+1).xd)^2 + (S(i).yd - S(n+1).yd)^2);
                CHs_this_round(num_clusterheads_this_round).distance = distance_to_bs;
                
                energy_ch_tx = 0;
                if (distance_to_bs > do)
                    energy_ch_tx = (ETX + EDA) * k_bits + Emp * k_bits * (distance_to_bs^4); 
                else
                    energy_ch_tx = (ETX + EDA) * k_bits + Efs * k_bits * (distance_to_bs^2); 
                end
                S(i).E = S(i).E - energy_ch_tx;
                current_round_ch_tx_energy = current_round_ch_tx_energy + energy_ch_tx;
                
                cumulative_ch_formations_count = cumulative_ch_formations_count + 1; % Global cumulative count
            end     
        end
    end
    CUMULATIVE_CH_FORMATIONS_STAT(r+1) = cumulative_ch_formations_count;
    CLUSTERHS(r+1) = num_clusterheads_this_round;
    STATISTICS(r+1).CLUSTERHEADS = num_clusterheads_this_round;

    % Member Node Association and Data Transmission to CH
    num_packets_to_ch_this_round = 0;
    if num_clusterheads_this_round > 0
        for i = 1:1:n
            if (strcmp(S(i).type,'N') && S(i).E > 0) % If node is Normal (not CH, not Dead) and Alive
                min_dis_to_ch = inf;
                assigned_ch_id = -1;

                for c_idx = 1:1:num_clusterheads_this_round
                    dist_sq = (S(i).xd - CHs_this_round(c_idx).xd)^2 + (S(i).yd - CHs_this_round(c_idx).yd)^2;
                    if dist_sq < min_dis_to_ch
                        min_dis_to_ch = dist_sq;
                        assigned_ch_id = CHs_this_round(c_idx).id;
                    end
                end
                min_dis_to_ch = sqrt(min_dis_to_ch);

                if assigned_ch_id ~= -1 % If a CH was found
                    S(i).cl = assigned_ch_id; % Assign node to cluster
                    
                    % Energy consumed by member node to transmit data to CH
                    energy_member_tx = 0;
                    if (min_dis_to_ch > do)
                        energy_member_tx = ETX * k_bits + Emp * k_bits * (min_dis_to_ch^4);
                    else
                        energy_member_tx = ETX * k_bits + Efs * k_bits * (min_dis_to_ch^2);
                    end
                    S(i).E = S(i).E - energy_member_tx;
                    current_round_member_tx_energy = current_round_member_tx_energy + energy_member_tx;

                    % Energy consumed by CH to receive this data
                    % (EDA is associated with reception per packet in this model)
                    energy_ch_rx = (ERX + EDA) * k_bits; 
                    S(assigned_ch_id).E = S(assigned_ch_id).E - energy_ch_rx;
                    current_round_ch_rx_energy = current_round_ch_rx_energy + energy_ch_rx;
                    
                    num_packets_to_ch_this_round = num_packets_to_ch_this_round + 1;
                else
                    % Node could not associate (e.g., no CHs or too far - not explicitly modeled here but good to note)
                    % In this basic LEACH, if CHs exist, nodes join the closest one.
                end
            end
        end
    end
    PACKETS_TO_CH_ROUND(r+1) = num_packets_to_ch_this_round;

    % Update Stats for this round
    ALIVE_NODES_PER_ROUND(r+1) = n - DEAD(r+1);
    
    current_total_network_energy = 0;
    alive_energies_vector = [];
    for k_node = 1:n
        if S(k_node).E > 0
            current_total_network_energy = current_total_network_energy + S(k_node).E;
            alive_energies_vector = [alive_energies_vector S(k_node).E];
        end
    end
    NETWORK_ENERGY_PER_ROUND(r+1) = current_total_network_energy;

    if ALIVE_NODES_PER_ROUND(r+1) > 0
        AVG_ENERGY_ALIVE_NODES_PER_ROUND(r+1) = current_total_network_energy / ALIVE_NODES_PER_ROUND(r+1);
        if length(alive_energies_vector) > 1
            ENERGY_VARIANCE_PER_ROUND(r+1) = var(alive_energies_vector);
        else
            ENERGY_VARIANCE_PER_ROUND(r+1) = 0; % Variance is 0 if only one or no alive nodes
        end
    else
        AVG_ENERGY_ALIVE_NODES_PER_ROUND(r+1) = 0;
        ENERGY_VARIANCE_PER_ROUND(r+1) = 0;
    end

    ENERGY_MEMBER_TX_TO_CH_PER_ROUND(r+1) = current_round_member_tx_energy;
    ENERGY_CH_RX_FROM_MEMBERS_PER_ROUND(r+1) = current_round_ch_rx_energy;
    ENERGY_CH_TX_TO_BS_PER_ROUND(r+1) = current_round_ch_tx_energy;

    if CLUSTERHS(r+1) > 0
        num_member_nodes = ALIVE_NODES_PER_ROUND(r+1) - CLUSTERHS(r+1);
        AVG_NODES_PER_CLUSTER_PER_ROUND(r+1) = num_member_nodes / CLUSTERHS(r+1);
    else
        AVG_NODES_PER_CLUSTER_PER_ROUND(r+1) = 0;
    end
    
    DATA_PACKETS_FROM_CH_TO_BS_ROUND(r+1) = CLUSTERHS(r+1); % Each CH sends one aggregated packet

    CONTROL_PACKETS_ADV_PER_ROUND(r+1) = CLUSTERHS(r+1); % CHs send advertisements
    % Active member nodes (alive, not CH) send join requests
    CONTROL_PACKETS_JOIN_PER_ROUND(r+1) = max(0, ALIVE_NODES_PER_ROUND(r+1) - CLUSTERHS(r+1)); 
    CONTROL_PACKETS_SCHEDULE_PER_ROUND(r+1) = CLUSTERHS(r+1); % CHs send TDMA schedules
    TOTAL_CONTROL_PACKETS_PER_ROUND(r+1) = CONTROL_PACKETS_ADV_PER_ROUND(r+1) + ...
                                           CONTROL_PACKETS_JOIN_PER_ROUND(r+1) + ...
                                           CONTROL_PACKETS_SCHEDULE_PER_ROUND(r+1);
    
    hold off;
    pause(0.01); % Pause for visualization, adjust as needed
    
    if DEAD(r+1) == n % All nodes are dead
        break; % Exit simulation loop
    end
end % End of simulation rounds loop
disp('Simulation Ended.');

% Adjust rmax if simulation ended early due to all nodes dying
actual_rounds = find(DEAD == n, 1, 'first');
if isempty(actual_rounds)
    actual_rounds = rmax + 1; % All rounds completed, or not all nodes died
else
    % If LND occurred, then that is the number of rounds for which we have full data
    % But plots should go up to 'actual_rounds - 1' for 0-indexed rounds, or 'actual_rounds' for 1-indexed rounds
    % Since 'r' goes from 0 to rmax, arrays are indexed r+1. So up to 'actual_rounds'.
end
plot_rounds = 0:min(r, rmax); % Rounds for plotting based on actual simulation progress


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   PLOTTING STATISTICS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Plotting results...');

% Figure 2: Network Lifetime Metrics
figure(2);
subplot(2,1,1);
plot(plot_rounds, ALIVE_NODES_PER_ROUND(1:length(plot_rounds)), 'b-', 'LineWidth', 1.5);
title('Number of Alive Nodes per Round');
xlabel('Round');
ylabel('Number of Alive Nodes');
grid on;
hold on;
if flag_first_dead, plot(first_dead_round, n - DEAD(first_dead_round+1), 'rx', 'MarkerSize', 8, 'DisplayName', sprintf('FND (Round %d)', first_dead_round)); end
if flag_half_dead, plot(half_dead_round, n - DEAD(half_dead_round+1), 'gx', 'MarkerSize', 8, 'DisplayName', sprintf('HND (Round %d)', half_dead_round)); end
if flag_last_dead, plot(last_dead_round, n - DEAD(last_dead_round+1), 'kx', 'MarkerSize', 8, 'DisplayName', sprintf('LND (Round %d)', last_dead_round)); end
legend('Alive Nodes', 'FND', 'HND', 'LND', 'Location', 'best');
hold off;

subplot(2,1,2);
plot(plot_rounds, DEAD(1:length(plot_rounds)), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Total Dead Nodes');
hold on;
plot(plot_rounds, DEAD_A(1:length(plot_rounds)), 'm--', 'LineWidth', 1.5, 'DisplayName', 'Dead Advanced Nodes');
plot(plot_rounds, DEAD_N(1:length(plot_rounds)), 'c:', 'LineWidth', 1.5, 'DisplayName', 'Dead Normal Nodes');
title('Number of Dead Nodes per Round');
xlabel('Round');
ylabel('Number of Dead Nodes');
legend('Location', 'best');
grid on;
hold off;

% Figure 3: Energy Metrics
figure(3);
subplot(2,2,1);
plot(plot_rounds, NETWORK_ENERGY_PER_ROUND(1:length(plot_rounds)), 'r-', 'LineWidth', 1.5);
title('Total Network Energy Remaining per Round');
xlabel('Round');
ylabel('Total Energy (J)');
grid on;

subplot(2,2,2);
plot(plot_rounds, AVG_ENERGY_ALIVE_NODES_PER_ROUND(1:length(plot_rounds)), 'g-', 'LineWidth', 1.5);
title('Average Energy per Alive Node per Round');
xlabel('Round');
ylabel('Average Energy (J)');
grid on;

subplot(2,2,3);
plot(plot_rounds, ENERGY_VARIANCE_PER_ROUND(1:length(plot_rounds)), 'm-', 'LineWidth', 1.5);
title('Variance of Residual Energy among Alive Nodes');
xlabel('Round');
ylabel('Energy Variance (J^2)');
grid on;

subplot(2,2,4);
plot(plot_rounds, ENERGY_MEMBER_TX_TO_CH_PER_ROUND(1:length(plot_rounds)), 'b-', 'DisplayName', 'Members TX to CH');
hold on;
plot(plot_rounds, ENERGY_CH_RX_FROM_MEMBERS_PER_ROUND(1:length(plot_rounds)), 'r--', 'DisplayName', 'CHs RX from Members');
plot(plot_rounds, ENERGY_CH_TX_TO_BS_PER_ROUND(1:length(plot_rounds)), 'g:', 'DisplayName', 'CHs TX to BS');
title('Energy Consumption Components per Round');
xlabel('Round');
ylabel('Energy (J)');
legend('Location', 'best');
grid on;
hold off;

% Figure 4: Cluster and Packet Metrics
figure(4);
subplot(2,2,1);
plot(plot_rounds, CLUSTERHS(1:length(plot_rounds)), 'b-', 'LineWidth', 1.5);
title('Number of Active Cluster Heads per Round');
xlabel('Round');
ylabel('Number of CHs');
grid on;

subplot(2,2,2);
plot(plot_rounds, AVG_NODES_PER_CLUSTER_PER_ROUND(1:length(plot_rounds)), 'g-', 'LineWidth', 1.5);
title('Average Number of Nodes per Cluster');
xlabel('Round');
ylabel('Avg Nodes/Cluster');
grid on;

subplot(2,2,3);
plot(plot_rounds, PACKETS_TO_CH_ROUND(1:length(plot_rounds)), 'k-', 'DisplayName', 'Data to CHs');
hold on;
plot(plot_rounds, DATA_PACKETS_FROM_CH_TO_BS_ROUND(1:length(plot_rounds)), 'm--', 'DisplayName', 'Data to BS (from CHs)');
title('Data Packets Transmitted per Round');
xlabel('Round');
ylabel('Number of Packets');
legend('Location', 'best');
grid on;
hold off;

subplot(2,2,4);
plot(plot_rounds, TOTAL_CONTROL_PACKETS_PER_ROUND(1:length(plot_rounds)), 'c-', 'LineWidth', 1.5, 'DisplayName', 'Total Control Pkts');
hold on;
plot(plot_rounds, CONTROL_PACKETS_ADV_PER_ROUND(1:length(plot_rounds)), 'r:', 'DisplayName', 'ADV Pkts');
plot(plot_rounds, CONTROL_PACKETS_JOIN_PER_ROUND(1:length(plot_rounds)), 'b-.', 'DisplayName', 'JOIN Pkts');
plot(plot_rounds, CONTROL_PACKETS_SCHEDULE_PER_ROUND(1:length(plot_rounds)), 'g--', 'DisplayName', 'SCHEDULE Pkts');
title('Control Packet Overhead per Round');
xlabel('Round');
ylabel('Number of Control Packets');
legend('Location', 'best');
grid on;
hold off;

% Figure 5: CH Selection Distribution
figure(5);
ch_selection_counts = zeros(1, n);
for i_node = 1:n
    ch_selection_counts(i_node) = S(i_node).num_times_ch;
end
histogram(ch_selection_counts, 'BinMethod', 'integers');
title('Distribution of CH Selections per Node');
xlabel('Number of Times Node became CH');
ylabel('Number of Nodes');
grid on;

% Print Lifetime Metrics to Command Window
fprintf('\n%%%%%%%%%%%%%%%% NETWORK LIFETIME METRICS %%%%%%%%%%%%%%%%\n');
if flag_first_dead
    fprintf('First Node Death (FND) occurred at round: %d\n', first_dead_round);
else
    fprintf('First Node Death (FND) did not occur within %d rounds.\n', rmax);
end

if flag_half_dead
    fprintf('Half Nodes Death (HND) occurred at round: %d\n', half_dead_round);
else
    fprintf('Half Nodes Death (HND) did not occur within %d rounds (Total dead: %d/%d).\n', rmax, DEAD(end), n);
end

if flag_last_dead
    fprintf('Last Node Death (LND) occurred at round: %d\n', last_dead_round);
else
    fprintf('Last Node Death (LND) did not occur within %d rounds (Total dead: %d/%d).\n', rmax, DEAD(end), n);
end
fprintf('Total CH formations during simulation: %d\n', cumulative_ch_formations_count);
fprintf('Total Data Packets Sent to CHs: %d\n', sum(PACKETS_TO_CH_ROUND));
fprintf('Total Data Packets Sent from CHs to BS: %d\n', sum(DATA_PACKETS_FROM_CH_TO_BS_ROUND));
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');

disp('All results plotted and summarized.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   END OF STATISTICS SECTION    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
