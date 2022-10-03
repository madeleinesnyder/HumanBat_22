%% Notes from Willa on ANOVA, GLM, and Logistic, and Linear Regression

%% Multinomial logistic regression for the case where N = {0,1}
[B,dev,stats] = mnrfit(table1,categorical(E_cat));

%% GLM linear regression for the case where N = {0.....n} integeters
b = glmfit(X,y,distr)



