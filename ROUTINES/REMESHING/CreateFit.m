function [Fit_Function, gof] = CreateFit(x, y, z)
% CreateFit: Fit function
% INPUTS:
%   x: Vector of x Coordinates of Nodes	
%   y: Vector of y Coordinates of Nodes	
%   z: Vector of z Coordinates of Nodes	 ( preload)
% OUTPUTS:
%      Fit_Function : a fit object representing the fit.
%      gof : structure with goodness-of fit info.

%% Fit
[xData, yData, zData] = prepareSurfaceData( x, y, z );
% Weight e.g. for boundarys could be included.

% Set up fittype and options.
% % 'poly1' - Linear polynomial curve
% % 'poly11' - Linear polynomial surface
% % 'poly2' - Quadratic polynomial curve
% % 'linearinterp' - Piecewise linear interpolation
% % 'cubicinterp' - Piecewise cubic interpolation
% % 'smoothingspline' - Smoothing spline (curve)
% % 'lowess' - Local linear regression (surface)

ft = 'biharmonicinterp';
% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft);%, 'Normalize', 'on' );

% Create Function Handle of Fit
Fit_Function = @(x,y) fitresult(x,y);

% % Plot fit with data.
% figure
% plot( fitresult, [xData, yData], zData );
% % Label axes
% xlabel x
% ylabel y
% zlabel z
% grid on
% h = get(gca,'DataAspectRatio') ;
% set(gca,'DataAspectRatio',[1 1 h(3)])
end


