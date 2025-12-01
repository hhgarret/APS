% clear
% clc
% rng(1,"twister");

% speed of sound (m/s)
c = 343;

% 2D sensor location [x y]
sensor_num = 100;
% -------- Square 
% This is very close to ideal layout
% Adding randomness with 30 cm range to all points.
n_rows = sqrt(sensor_num);
x = repmat((linspace(0,3*(n_rows-1),n_rows))',n_rows,1);
y = repelem((linspace(0,3*(n_rows-1),n_rows))',n_rows);
P_actual = [x y];

% wind velocity (m/s) [w_x w_y] up to 2, no z value
w_actual = [0 -1];

% Calculate time of flight between nodes
t_actual = zeros(sensor_num,sensor_num);
d_actual = zeros(sensor_num,sensor_num);
for i=1:sensor_num
    for j=1:sensor_num
        if j~=i
            Pij = P_actual(i,:)-P_actual(j,:);
            t_actual(i,j) = (norm(Pij)^2)/(norm(Pij)*c+dot(w_actual,Pij));
            d_actual(i,j) = norm(Pij);
        end
    end
end

% Generate Gaussian noise
% -- Note: a vector of n random values from a normal distribution with mean
% mu and standard deviation sd: mu+sd*randn(1,n)

% ------- 5% of each value
noise = 0.05*(t_actual.*randn(sensor_num,sensor_num));


% ------- 5% of mean t_actual value
% noise_std = 0.4 * mean(abs(t_actual),"all"); 
% noise = noise_std * randn(size(t_actual));

% ------- 5% of distance
% noise = 0.01 * log(1+d_actual).*randn(size(t_actual)); 

% ------- Add the noise 
t_measured = t_actual + noise; 
% Making sure t_measured values are positive
for i = 1:sensor_num
    for j=1:sensor_num
        if t_measured(i,j)<0
            t_measured(i,j) = t_actual(i,j) - noise(i,j);
        end
    end
end
t_measured(boolean(eye(length(t_actual)))) = 0;

d_measured = zeros(sensor_num,sensor_num);
for i = 1:sensor_num
    for j=1:sensor_num
        d_measured(i,j) = 2*c/((1/t_measured(i,j))+1/t_measured(j,i));
    end
end

%% writematrix(P_actual,".\APS_coordinates_original.txt")
%% writematrix(t_measured,".\t_measured.txt")


%% Plot Results
p_original = readmatrix(".\APS_coordinates_original.txt");
p_LM = readmatrix(".\APS_coordinates_optimized.txt");
p_MDS = readmatrix(".\APS_MDS_coordinates_calculated.txt");
p_code = readmatrix(".\APS_RotationVerification_optimized_rotated_coordinates.txt");

d_original2 = pdist2(p_original, p_original);
d_code2 = pdist2(p_code, p_code);
d_LM = pdist2(p_LM, p_LM);
d_MDS = pdist2(p_MDS, p_MDS);


figure
hold on
plot(p_original(:,1),p_original(:,2),'bo','MarkerFaceColor','b','MarkerSize',5)
plot(p_MDS(:,1),p_MDS(:,2),'mo','MarkerFaceColor','m','MarkerSize',5)
plot(p_LM(:,1),p_LM(:,2),'ko','MarkerFaceColor','k','MarkerSize',5)
plot(p_code(:,1),p_code(:,2),'ro','MarkerFaceColor','r','MarkerSize',5)
plot(p_code(1,1),p_code(1,2),'go','MarkerFaceColor','r','MarkerSize',10)
plot(p_LM(1,1),p_LM(1,2),'go','MarkerFaceColor','k','MarkerSize',10)
plot(p_MDS(1,1),p_MDS(1,2),'go','MarkerFaceColor','m','MarkerSize',10)
plot(p_original(1,1),p_original(1,2),'go','MarkerFaceColor','b','MarkerSize',10)

grid minor
title('Comparing LM Output and Rotated')
legend('Actual','MDS output','LM output', 'Rotated Output')

% Calculate error
nodes_difference = p_original-p_code;
mean_node_diff = mean(abs(nodes_difference))

%%
sensor_num = 100;
d_code = zeros(sensor_num,sensor_num);
for i=1:sensor_num
    for j=1:sensor_num
        if j~=i
            Pij = p_original(i,:)-p_original(j,:);
            d_actual(i,j) = norm(Pij);
            Pij2 = p_code(i,:)-p_code(j,:);
            d_code(i,j) = norm(Pij2);
        end
    end
end

err = 100*abs(d_actual-d_code)./d_actual;
mean(err,"all","omitnan")

rmse_original_mds = rmse(d_original2, d_MDS, "all")
rmse_original_lm = rmse(d_original2, d_LM, "all")