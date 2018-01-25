function [HereIsYourHandle HereAreYourAxes] = DoMeAFigure(AxisLimits,BackgroundColour)

if nargin < 2
    BackgroundColour = [1 1 1];
end

HereIsYourHandle = figure('Color',BackgroundColour, ...
                          'Renderer', 'OpenGL', ...
                          'Units', 'inches', ...
                          'PaperUnits', 'inches', ...
                          'PaperSize', [10 10], ...
                          'PaperPositionMode', 'manual', ...
                          'PaperPosition', [0 0 10 10], ...
                          'Visible','off');
                      
HereAreYourAxes =  axes('DataAspectRatio', [1,1,1], ...
                        'Position', [0 0 1 1], ...
                        'Color',BackgroundColour,...
                        'Visible','off');
                    
axis(AxisLimits);
axis square image tight
box('off');
hold on % don't reset anything when plotting this next thing...