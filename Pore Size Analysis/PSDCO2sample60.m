
R = readmatrix("PSDsample60vCLC.xlsx");

X_SLC60 = R(1:114, 1);
Y_SLC60 = R(1:114, 2);

X_CLC = R(1:104, 6);
Y_CLC = R(1:104, 7);


%% Isotherm Plot

hfig = figure;
plot(X_SLC60,Y_SLC60,'r');
hold on
plot(X_CLC, Y_CLC, 'Color', '#5F9EA0')


% Enhance the plot for publication quality
xlabel('Pore Width (\r{A})', 'Interpreter', 'latex');
ylabel('$dV/dW$ Pore Volume (cm$^3$/g)', 'Interpreter', 'latex');

% Manually set axis limits
xlim([4.6 20]); % Replace with your chosen limits
%set(gca, 'XScale', 'log');

fname = 'PSDCO260micro';
picturewidth = 8; % set this parameter and keep it forever
hw_ratio = 0.8; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',11) % adjust fontsize to your document
set(findall(hfig,'-property','Box'),'Box','on') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
box on;
grid on;
hl = legend ('SLC60', 'CLC','NorthEast');
set(hl, 'Interpreter', 'latex', 'FontSize', 9)
%print(hfig,fname,'-dpdf','-painters','-fillpage')

%% Isotherm Plot

hfig = figure;
plot(X_SLC60,Y_SLC60,'r');
hold on
plot(X_CLC, Y_CLC, 'Color', '#5F9EA0')

cm = prism(6);
cm(3,:) = [];
cm(2,:) = [];
set(gca,'ColorOrder',cm);


% Enhance the plot for publication quality
xlabel('Pore Width (\r{A})', 'Interpreter', 'latex');
ylabel('$dV/dW$ Pore Volume (cm$^3$/g)', 'Interpreter', 'latex');

% Manually set axis limits
xlim([20 240]); % Replace with your chosen limits

fname = 'PSDCO260meso';
picturewidth = 8; % set this parameter and keep it forever
hw_ratio = 0.8; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',11) % adjust fontsize to your document
set(findall(hfig,'-property','Box'),'Box','on') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
box on;
grid on;
hl = legend ('SLC60', 'CLC','NorthEast');
set(hl, 'Interpreter', 'latex', 'FontSize', 9)
%print(hfig,fname,'-dpdf','-painters','-fillpage')

