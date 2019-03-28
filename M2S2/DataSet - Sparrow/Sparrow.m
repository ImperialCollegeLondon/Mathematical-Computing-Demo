%% Import Data and Changing Variable Names
% Number of Observation = 13
data = readtable('Sparrow.csv');
data.Properties.VariableNames{'x__ays'} = 'days';

%% Design Matrix
X = [ones(13,1) data.days];

%% Least Square Estimator
Y = data.Length;
beta = (X' * X) \ (X' * Y)
Yhat = beta(1) + beta(2) * data.days;

%% Plot
hold on
scatter(data.days,Y)                      % Real Data
plot(data.days, Yhat)                     % Fitted Data
plot(data.days, mean(Y)*ones(13,1))       % Constant Model
xlim([0 18])
hold off

%% Estimate sigma^2
resvec = Y - Yhat;
SSE = norm(resvec)
MSE = SSE / (13-2)

%% Global F Test (Whether Constant Model is Appropriate)
devvec = Y - mean(Y)*ones(13,1);
SST = norm(devvec)
F = (SST - SSE) / ((2-1) * MSE)
pGlobal = 1 - fcdf(F, 1, 11)

%% QQ Plot Normal Assumption
P = X * ((X' * X) \ X');                    
stdres = resvec ./ sqrt((1 - diag(P))*MSE);
stdres = sort(stdres);                    % Sort
hold on
scatter(norminv((1:13)/(13+1)), stdres)
plot(linspace(-1.5,1.5,20), linspace(-1.5,1.5,20))
hold off