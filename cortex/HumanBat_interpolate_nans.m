function [interp_trajectory] = HumanBat_interpolate_nans(trajectory)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here 

x = trajectory;
nanx = isnan(x);
t    = 1:numel(x);
x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx));

x_movmed_1 = movmedian(x(1,:),20); x_movmed_2 = movmedian(x(2,:),20); x_movmed_3 = movmedian(x(3,:),20);


interp_trajectory  = [x_movmed_1;x_movmed_2;x_movmed_3];
end