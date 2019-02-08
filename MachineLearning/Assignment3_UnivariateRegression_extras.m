%#!/bin/octave

%Root mean squares and regularizers for coursera data
clear
regularizer_flag = input('Would you like to use a regularizer for the coursera data? (y/n): ', 's');
if (regularizer_flag == 'n'),
  data = load('AndrewNGcourseraEx1_foodtruckSales.txt');

  m = length(data);
  X = [ones(m, 1), data(:, 1)];
  y = data(:, 2);
  theta = zeros(2, 1);

  %Computing the cost function to plot root mean squares vs theta later
  function J = CostFun3(X, y, m, theta)
    prediction_products = X * theta;
    squared_differences = (prediction_products - y) .^2;
    J = sum(squared_differences) / (2 * m);
  endfunction

  %Applying Grade and Descent, but slightly modified to calculate
  %root mean squares as well
  iterations = 1500;
  alpha = 0.01;
  J_history = zeros(iterations, 1);
  rmse = zeros(iterations, 2);

  for j = 1 : iterations
    FeatByParameMatrix = [theta(1, 1) * X(:, 1), theta(2, 1) * X(:, 2)];
    predicted_values_y = FeatByParameMatrix(:, 1) + FeatByParameMatrix(:, 2);
    residuals = y - predicted_values_y;
    rmse(j, 1) = sqrt(sum(residuals .^2) / m); rmse(j, 2) = j;
    %Computing the root mean squares for each theta update

    delta = (((X * theta - y)' * X) / m)';
    theta = theta - (alpha * delta);

    J_history(j, 1) = CostFun3(X, y, m, theta);

  end;

  %Plot root mean squares error for every cost function value update
  figure
  plot(J_history, rmse(:, 1));
  xlabel('Cost - J(\theta_1 \theta_2)'); ylabel('Root mean square error');
  title('Coursera RMSE changing with each grade & descent Theta update')

  %Plot root mean squares error for every iteration
  figure
  plot(rmse(:, 2), rmse(:, 1));
  ylim([0, 10])
  xlabel('Iterations'); ylabel('Root mean square error');
  title('Coursera RMSE changing with each grade & descent iteration')

  disp('For the coursera data:')
  disp('The parameters for theta found by grade and descent are ')
  theta
  disp('The cost computed by the cost function with these parameters is')
  J_history(iterations, 1)
  disp('The root mean square error for the model using these parameters is ')
  rmse(iterations, 1)

elseif (regularizer_flag == 'y'),
  lambda = input('Enter the regularizer parameter: ');
  if (isnumeric(lambda) == 0),
    disp('Please run the script again with a numeric value for lambda')
    return;
  else

    DATA = load('AndrewNGcourseraEx1_foodtruckSales.txt');
    y = DATA(:, 2);
    m = length(y);
    x_0 = ones(m, 1);
    x_1 = DATA(:, 1);
    X = [x_0 x_1];

    %Modified cost function for regularization
    function J = CostFun4(X, y, m, THETA, theta_1, lambda)
      J = (sum(((X * THETA) - y) .^2) + (lambda * (theta_1 .^2))) / (2 * m);
    endfunction

    theta_0 = 0;
    theta_1= 0;
    THETA= [theta_0;theta_1];

    iterations = 1500;
    alpha = 0.01;
    J_history = zeros(iterations, 1);
    RMSE = zeros(iterations, 2);

    %Modified grade and descent for regularization
    for j = 1 : iterations
      FeatByParameMatrix = [(theta_0 * x_0) (theta_1 * x_1)];
      predicted_values_y = FeatByParameMatrix(:, 1) + FeatByParameMatrix(:, 2);
      residuals = y - predicted_values_y;
      RMSE(j, 1) = sqrt(sum(residuals .^2) / m); RMSE(j, 2) = j;

      theta_0 = theta_0 - (((((X * THETA) - y)' * x_0) * alpha) / m);
      theta_1 = (theta_1 * (1 - (alpha * (lambda / m)))) - (((((X * THETA) - y)' * x_1) * alpha) / m);
      THETA = [theta_0;theta_1];
      J_history(j, 1) = CostFun4(X, y, m, THETA, theta_1, lambda);

    end;

    %Plot the error function
    theta0_vals = linspace(-10, 10, 100);
    theta1_vals = linspace(-1, 4, 100);

    J_vals = zeros(length(theta0_vals), length(theta1_vals));

    for i = 1 : length(theta0_vals)
      for j = 1 : length(theta1_vals)
        t = [theta0_vals(i); theta1_vals(j)];
        J_vals(i, j) = CostFun4(X, y, m, t, j, lambda);
      end
    end

    J_vals = J_vals';
    figure
    surf(theta0_vals, theta1_vals, J_vals)
    xlabel('\theta_0'); ylabel('\theta_1'); zlabel('Cost')
    title('Cousera 3D cost function regularized')

    %Plot root mean squares error for every iteration
    figure
    plot(RMSE(:, 2), RMSE(:, 1));
    ylim([0, 10])
    xlabel('Iterations'); ylabel('Root mean square error');
    title('Coursera RMSE changing with each grade & descent iteration for regularized model')

    %Plot model fit
    figure
    plot(x_1, y, 'rx', 'MarkerSize', 10);
    xlabel('Population of City in 10,000s'); ylabel('Profit in $10,000s');
    title(lambda)
    text (20, 0, "Coursera data");
    hold on;
    plot(X(:, 2), X * THETA, '-')
    legend('Training data', 'Linear regression')
    hold off;


    disp('For the coursera data:')
    disp('The parameters for theta found by grade and descent are ')
    THETA
    disp('The cost computed by the cost function with these parameters is')
    J_history(iterations, 1)
    disp('The root mean square error for the model using these parameters is ')
    RMSE(iterations, 1)

  end;

else
  disp('No valid answer given, exiting script...');
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tested Wickham's non-scaled data and it exceeds float capacity
%Root mean squares and regularizers for Wickham's log-scaled data
clear
regularizer_flag = input('Would you like to use a regularizer for the Wickhams log-scaled data? (y/n): ', 's');
if (regularizer_flag == 'n'),
  data = load('deviation_tidy_wickham.txt');
  data = log10(data);

  m = length(data);
  X = [ones(m, 1), data(:, 1)];
  y = data(:, 2);
  theta = zeros(2, 1);

  %Computing the cost function to plot root mean squares vs theta later
  function J = CostFun7(X, y, m, theta)
    prediction_products = X * theta;
    squared_differences = (prediction_products - y) .^2;
    J = sum(squared_differences) / (2 * m);
  endfunction

  %Applying Grade and Descent, but slightly modified to calculate
  %root mean squares as well
  iterations = 1500;
  alpha = 0.01;
  J_history = zeros(iterations, 1);
  rmse = zeros(iterations, 2);

  for j = 1 : iterations
    FeatByParameMatrix = [theta(1, 1) * X(:, 1), theta(2, 1) * X(:, 2)];
    predicted_values_y = FeatByParameMatrix(:, 1) + FeatByParameMatrix(:, 2);
    residuals = y - predicted_values_y;
    rmse(j, 1) = sqrt(sum(residuals .^2) / m); rmse(j, 2) = j;
    %Computing the root mean squares for each theta update

    delta = (((X * theta - y)' * X) / m)';
    theta = theta - (alpha * delta);

    J_history(j, 1) = CostFun7(X, y, m, theta);

  end;

  %Plot root mean squares error for every cost function value update
  figure
  plot(J_history, rmse(:, 1));
  xlabel('Cost - J(\theta_1 \theta_2)'); ylabel('Root mean square error');
  title('Wickham log-scaled RMSE changing with each grade & descent Theta update')

  %Plot root mean squares error for every iteration
  figure
  plot(rmse(:, 2), rmse(:, 1));
  ylim([0, 10])
  xlabel('Iterations'); ylabel('Root mean square error');
  title('Wickham log-scaled RMSE changing with each grade & descent iteration')

  disp('For the Wickham data:')
  disp('The parameters for theta found by grade and descent are ')
  theta
  disp('The cost computed by the cost function with these parameters is')
  J_history(iterations, 1)
  disp('The root mean square error for the model using these parameters is ')
  rmse(iterations, 1)

elseif (regularizer_flag == 'y'),
  lambda = input('Enter the regularizer parameter: ');
  if (isnumeric(lambda) == 0),
    disp('Please run the script again with a numeric value for lambda')
    return;
  else

    DATA = load('deviation_tidy_wickham.txt');
    DATA = log10(DATA);
    y = DATA(:, 2);
    m = length(y);
    x_0 = ones(m, 1);
    x_1 = DATA(:, 1);
    X = [x_0 x_1];

    %Modified cost function for regularization
    function J = CostFun8(X, y, m, THETA, theta_1, lambda)
      J = (sum(((X * THETA) - y) .^2) + (lambda * (theta_1 .^2))) / (2 * m);
    endfunction

    theta_0 = 0;
    theta_1= 0;
    THETA= [theta_0;theta_1];

    iterations = 1500;
    alpha = 0.01;
    J_history = zeros(iterations, 1);
    RMSE = zeros(iterations, 2);

    %Modified grade and descent for regularization
    for j = 1 : iterations
      FeatByParameMatrix = [(theta_0 * x_0) (theta_1 * x_1)];
      predicted_values_y = FeatByParameMatrix(:, 1) + FeatByParameMatrix(:, 2);
      residuals = y - predicted_values_y;
      RMSE(j, 1) = sqrt(sum(residuals .^2) / m); RMSE(j, 2) = j;

      theta_0 = theta_0 - (((((X * THETA) - y)' * x_0) * alpha) / m);
      theta_1 = (theta_1 * (1 - (alpha * (lambda / m)))) - (((((X * THETA) - y)' * x_1) * alpha) / m);
      THETA = [theta_0;theta_1];
      J_history(j, 1) = CostFun8(X, y, m, THETA, theta_1, lambda);

    end;

    %Plot the error function
    theta0_vals = linspace(-10, 10, 100);
    theta1_vals = linspace(-1, 4, 100);

    J_vals = zeros(length(theta0_vals), length(theta1_vals));

    for i = 1 : length(theta0_vals)
      for j = 1 : length(theta1_vals)
        t = [theta0_vals(i); theta1_vals(j)];
        J_vals(i, j) = CostFun8(X, y, m, t, j, lambda);
      end
    end

    J_vals = J_vals';
    figure
    surf(theta0_vals, theta1_vals, J_vals)
    xlabel('\theta_0'); ylabel('\theta_1'); zlabel('Cost')
    title('Wickham log-scaled 3D cost function regularized')

    %Plot root mean squares error for every iteration
    figure
    plot(RMSE(:, 2), RMSE(:, 1));
    ylim([0, 10])
    xlabel('Iterations'); ylabel('Root mean square error');
    title('Wickham log-scaled RMSE changing with each grade & descent iteration for regularized model')

    %Plot model fit
    figure
    plot(x_1, y, 'rx', 'MarkerSize', 10);
    xlabel('Observations(n)'); ylabel('Deviation distribution');
    title(lambda)
    text (2, -5, "Wickham log-scale");
    hold on;
    plot(X(:, 2), X * THETA, '-')
    legend('Training data', 'Linear regression')
    hold off;

    disp('For the Wickham data:')
    disp('The parameters for theta found by grade and descent are ')
    THETA
    disp('The cost computed by the cost function with these parameters is')
    J_history(iterations, 1)
    disp('The root mean square error for the model using these parameters is ')
    RMSE(iterations, 1)

  end;

else
  disp('No valid answer given, exiting script...');
end;
