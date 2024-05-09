
% matlab formatting properties

% consistent colour scheme
colour = [0 0.5 0.95;
          1 0 1;
          0 0.8 0.7;
          1 0.1 0.1;
          0.9 0.9 0.2];

% have horizontal but not vertical gridlines
set(groot, 'DefaultAxesXGrid', 'off');
set(groot, 'DefaultAxesYGrid', 'on');

% big font - only for when have two graphs side by side in overleaf
set(groot,'DefaultAxesFontSize',16)
set(groot, 'DefaultTextInterpreter', 'latex');

% set default line widths
set(groot, 'DefaultAxesLineWidth', 0.5);
set(groot, 'DefaultLineLineWidth', 1);

% if use yyaxis right for a second axis - change second axis colour
%%% set(gca,'YColor',colour(1,:))

% if want semilog axes
%%% set(gca,'xscale','log')

% autosave the picture with high resolution (300dpi)
%%% print('name', '-dpng', '-r300');

