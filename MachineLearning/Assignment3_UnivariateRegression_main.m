%#!/bin/octave

%First part is from ex1 AndrewNG coursera (in pdf, no here), second part is here
%This script covers the procedure used by AndrewNG in the coursera data but
%applied to a part of Wickham's tidy data paper case study
%Extra procedures like root mean squares and regularization computation will be
%covered in another script "Assignment3_UnivariateRegression_extra.m"

%Used R to export csv with the data
%Used BASH to format the data for OCTAVE

%There is one section for the unscaled data and one for the log scaled data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load, ready up data and variables
clear; close all; clc
fprintf('Start of linear scale computations...\n')
data = load('deviation_tidy_wickham.txt');
X = data(:, 1);
y = data(:, 2);
m = length(y);

%Plot the training set
figure
plot(X, y, 'rx', 'MarkerSize', 10);
xlabel('Observations(n)'); ylabel('Deviation distribution');
title('no scaling - training set')

%Compute the cost function J (quadratic errors in this case)
%assuming thetas of 0 for testing the computation
theta = zeros(2, 1);
X = [ones(m, 1), X(:, 1)];

function J = CostFun1(X, y, m, theta)
  prediction_products = X * theta;
  squared_differences = (prediction_products - y) .^2;
  J = sum(squared_differences) / (2 * m);
endfunction

CostFun1(X, y, m, theta)
fprintf('The latter is the cost of the error function when thetas are zeros\n')

%Plot the error function
theta0_vals = linspace(-10, 10, 100);
theta1_vals = linspace(-1, 4, 100);

J_vals = zeros(length(theta0_vals), length(theta1_vals));

for i = 1 : length(theta0_vals)
  for j = 1 : length(theta1_vals)
    t = [theta0_vals(i); theta1_vals(j)];
    J_vals(i, j) = CostFun1(X, y, m, t);
  end
end

J_vals = J_vals';
figure
surf(theta0_vals, theta1_vals, J_vals)
xlabel('\theta_0'); ylabel('\theta_1');
title('no scaling - 3D cost function')

figure
contour(theta0_vals, theta1_vals, J_vals, logspace(-2, 3, 20))
xlabel('\theta_0'); ylabel('\theta_1');
title('no scaling - contour cost function plot')

%Compute grade and Descent
iterations = 1500;
alpha = 0.00000001;
J_history = zeros(iterations, 1);

for j = 1 : iterations
  delta = (((X * theta - y)' * X) / m)';
  theta = theta - (alpha * delta);

  J_history(j, 1) = CostFun1(X, y, m, theta);
end;

theta'
fprintf('The latter are the parameters for theta0 & theta1');
fprintf(' respectively that minimize the cost function\n')

%Plot the fitted model
old_X = data(:, 1);
figure
plot(old_X, y, 'rx', 'MarkerSize', 10);
xlabel('Observations(n)'); ylabel('Deviation distribution');
title('no scaling - fitted model')
hold on;

plot(X(:, 2), X * theta, '-')
legend('Training data', 'Linear regression')
hold off;

fprintf('End of linear scale computations\n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Start of log scale computations...\n')
%Load, ready up data and variables
clear
data = load('deviation_tidy_wickham.txt');
X = log10(data(:, 1));
y = log10(data(:, 2));
m = length(y);

%Plot the training set
figure
plot(X, y, 'rx', 'MarkerSize', 10);
xlabel('Observations(n)'); ylabel('Deviation distribution');
title('log scale - training set')

%Compute the cost function J (quadratic errors in this case)
%assuming thetas of 0 for testing the computation
theta = zeros(2, 1);
X = [ones(m, 1), X(:, 1)];

function J = CostFun2(X, y, m, theta)
  prediction_products = X * theta;
  squared_differences = (prediction_products - y) .^2;
  J = sum(squared_differences) / (2 * m);
endfunction

CostFun2(X, y, m, theta)
fprintf('The latter is the cost of the error function when thetas are zeros\n')

%Plot the error function
theta0_vals = linspace(-10, 10, 100);
theta1_vals = linspace(-1, 4, 100);

J_vals = zeros(length(theta0_vals), length(theta1_vals));

for i = 1 : length(theta0_vals)
  for j = 1 : length(theta1_vals)
    t = [theta0_vals(i); theta1_vals(j)];
    J_vals(i, j) = CostFun2(X, y, m, t);
  end
end

J_vals = J_vals';
figure
surf(theta0_vals, theta1_vals, J_vals)
xlabel('\theta_0'); ylabel('\theta_1');
title('log scale - 3D cost function')

figure
contour(theta0_vals, theta1_vals, J_vals, logspace(-2, 3, 20))
xlabel('\theta_0'); ylabel('\theta_1');
title('log scale - contour cost function plot')

%Compute grade and Descent
iterations = 1500;
alpha = 0.01;
J_history = zeros(iterations, 1);

for j = 1 : iterations
  delta = (((X * theta - y)' * X) / m)';
  theta = theta - (alpha * delta);

  J_history(j, 1) = CostFun2(X, y, m, theta);
end;

theta'
fprintf('The latter are the parameters for theta0 & theta1');
fprintf(' respectively that minimize the cost function\n')

%Plot the fitted model
old_X = log10(data(:, 1));
figure
plot(old_X, y, 'rx', 'MarkerSize', 10);
xlabel('Observations(n)'); ylabel('Deviation distribution');
title('log scale - fitted model')
hold on;

plot(X(:, 2), X * theta, '-')
legend('Training data', 'Linear regression')
hold off;

fprintf('End of log scale computations\n')
