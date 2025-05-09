%% LEACH Protocol Simulation in MATLAB
% This script simulates the LEACH (Low-Energy Adaptive Clustering Hierarchy) protocol
% for wireless sensor networks

clear all;
close all;
clc;

%% Network Parameters
% Field dimensions
xm = 100; % Length of the field in meters
ym = 100; % Width of the field in meters

% Base station position (at the center)
sink.x = 50;
sink.y = 50;

% Nodes
n = 100; % Number of nodes
p = 0.1; % Probability of becoming cluster head
rounds = 20; % Number of rounds to simulate

% Energy parameters (in Joules)
Eo = 0.5; % Initial energy of each node
ETX = 50e-9; % Energy for transmitting one bit (J/bit)
ERX = 50e-9; % Energy for receiving one bit (J/bit)
EDA = 5e-9; % Energy for data aggregation (J/bit)
Efs = 10e-12; % Free space model (d^2 loss)
Emp = 0.0013e-12; % Multi-path model (d^4 loss)
d0 = sqrt(Efs/Emp); % Threshold distance

% Packet sizes
packetLength = 4000; % bits
controlPacketLength = 200; % bits

%% Initialize Nodes
for i = 1:n
    % Random node position
    S(i).xd = rand * xm;
    S(i).yd = rand * ym;
    
    % Initial energy
    S(i).E = Eo;
    
    % Node is alive
    S(i).alive = 1;
    
    % Node has not been a cluster head yet
    S(i).clusterHead = 0;
    S(i).type = 'N'; % N: normal node, C: cluster head
    
    % Calculate distance to base station
    S(i).distToBS = sqrt((S(i).xd - sink.x)^2 + (S(i).yd - sink.y)^2);
    
    % Initialize round that node died
    S(i).diedRound = 0;
end

%% Network Statistics Initialization
stats.dead = zeros(1, rounds);
stats.alive = zeros(1, rounds);
stats.totalEnergy = zeros(1, rounds);
stats.packetsToBS = zeros(1, rounds);
stats.clusterHeads = zeros(1, rounds);

%% Start Simulation
figure(1);
for r = 1:rounds
    fprintf('Round: %d\n', r);
    
    % Reset cluster head count and packets sent
    countCHs = 0;
    packetsSent = 0;
    
    % Cluster head selection
    for i = 1:n
        S(i).type = 'N';
        S(i).CH = 0;
    end
    
    % Algorithm for cluster head selection
    for i = 1:n
        if S(i).alive == 1
            temp_rand = rand;
            if temp_rand <= p
                countCHs = countCHs + 1;
                S(i).type = 'C';
                S(i).CH = 1;
                
                % Announce cluster head status
                % Broadcast to all nodes requires energy
                distanceSum = 0;
                for j = 1:n
                    if S(j).alive == 1
                        distanceSum = distanceSum + sqrt((S(i).xd - S(j).xd)^2 + (S(i).yd - S(j).yd)^2);
                    end
                end
                
                % Energy spent by cluster head for broadcasting
                energy = ETX * controlPacketLength + Efs * controlPacketLength * (distanceSum^2);
                S(i).E = S(i).E - energy;
            end
        end
    end
    
    % Assign nodes to clusters
    for i = 1:n
        if S(i).alive == 1 && S(i).CH == 0 % If node is alive and not a cluster head
            min_dis = Inf;
            min_dis_cluster = 0;
            
            % Find closest cluster head
            for j = 1:n
                if S(j).alive == 1 && S(j).CH == 1
                    temp_dis = sqrt((S(i).xd - S(j).xd)^2 + (S(i).yd - S(j).yd)^2);
                    if temp_dis < min_dis
                        min_dis = temp_dis;
                        min_dis_cluster = j;
                    end
                end
            end
            
            % If a cluster head was found
            if min_dis_cluster ~= 0
                S(i).cluster = min_dis_cluster;
                
                % Energy spent by node to transmit to cluster head
                if min_dis > d0
                    S(i).E = S(i).E - (ETX * packetLength + Emp * packetLength * (min_dis^4));
                else
                    S(i).E = S(i).E - (ETX * packetLength + Efs * packetLength * (min_dis^2));
                end
                
                % Energy spent by cluster head to receive
                S(min_dis_cluster).E = S(min_dis_cluster).E - (ERX * packetLength);
            else
                % No cluster head found, node sends directly to base station
                distance = S(i).distToBS;
                
                % Energy spent to transmit to base station
                if distance > d0
                    S(i).E = S(i).E - (ETX * packetLength + Emp * packetLength * (distance^4));
                else
                    S(i).E = S(i).E - (ETX * packetLength + Efs * packetLength * (distance^2));
                end
                
                packetsSent = packetsSent + 1;
            end
        end
    end
    
    % Cluster heads aggregate data and send to base station
    for i = 1:n
        if S(i).alive == 1 && S(i).CH == 1
            % Count members in this cluster
            membersCount = 0;
            for j = 1:n
                if S(j).alive == 1 && S(j).CH == 0 && S(j).cluster == i
                    membersCount = membersCount + 1;
                end
            end
            
            % Aggregate data (energy cost for data aggregation)
            S(i).E = S(i).E - (EDA * packetLength * membersCount);
            
            % Send to base station
            distance = S(i).distToBS;
            if distance > d0
                S(i).E = S(i).E - (ETX * packetLength + Emp * packetLength * (distance^4));
            else
                S(i).E = S(i).E - (ETX * packetLength + Efs * packetLength * (distance^2));
            end
            
            packetsSent = packetsSent + 1;
        end
    end
    
    % Update node status (alive/dead)
    for i = 1:n
        if S(i).E <= 0 && S(i).alive == 1
            S(i).alive = 0;
            S(i).diedRound = r;
        end
    end
    
    % Calculate statistics for this round
    deadCount = 0;
    for i = 1:n
        if S(i).alive == 0
            deadCount = deadCount + 1;
        end
    end
    stats.dead(r) = deadCount;
    stats.alive(r) = n - deadCount;
    stats.clusterHeads(r) = countCHs;
    stats.packetsToBS(r) = packetsSent;
    
    totalEnergy = 0;
    for i = 1:n
        if S(i).alive == 1
            totalEnergy = totalEnergy + S(i).E;
        end
    end
    stats.totalEnergy(r) = totalEnergy;
    
    % Plot network state
    clf;
    hold on;
    
    % Plot base station
    plot(sink.x, sink.y, 'k^', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    
    % Plot nodes
    for i = 1:n
        if S(i).alive == 1
            if S(i).CH == 1
                plot(S(i).xd, S(i).yd, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
            else
                plot(S(i).xd, S(i).yd, 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
            end
        else
            plot(S(i).xd, S(i).yd, 'ko', 'MarkerSize', 6);
        end
    end
    
    % Plot clusters (lines between nodes and their cluster heads)
    for i = 1:n
        if S(i).alive == 1 && S(i).CH == 0 && S(i).cluster > 0
            if S(S(i).cluster).alive == 1
                line([S(i).xd, S(S(i).cluster).xd], [S(i).yd, S(S(i).cluster).yd], 'Color', 'g');
            end
        end
    end
    
    % Plot cluster head to base station connections
    for i = 1:n
        if S(i).alive == 1 && S(i).CH == 1
            line([S(i).xd, sink.x], [S(i).yd, sink.y], 'Color', 'r', 'LineStyle', '--');
        end
    end
    
    title(sprintf('Round %d: %d alive, %d dead, %d CHs', r, stats.alive(r), stats.dead(r), countCHs));
    xlabel('X (m)');
    ylabel('Y (m)');
    axis([0 xm 0 ym]);
    grid on;
    drawnow;
    pause(0.5);
end

%% Plot Results
figure(2);
subplot(2,2,1);
plot(1:rounds, stats.alive, 'b-', 'LineWidth', 2);
title('Number of Alive Nodes');
xlabel('Round');
ylabel('Nodes');
grid on;

subplot(2,2,2);
plot(1:rounds, stats.dead, 'r-', 'LineWidth', 2);
title('Number of Dead Nodes');
xlabel('Round');
ylabel('Nodes');
grid on;

subplot(2,2,3);
plot(1:rounds, stats.clusterHeads, 'g-', 'LineWidth', 2);
title('Number of Cluster Heads');
xlabel('Round');
ylabel('Cluster Heads');
grid on;

subplot(2,2,4);
plot(1:rounds, stats.totalEnergy, 'm-', 'LineWidth', 2);
title('Total Network Energy');
xlabel('Round');
ylabel('Energy (J)');
grid on;

figure(3);
plot(1:rounds, stats.packetsToBS, 'b-', 'LineWidth', 2);
title('Packets Sent to Base Station');
xlabel('Round');
ylabel('Packets');
grid on;
