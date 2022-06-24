% ----------------------------------------------------------------------------
% getFigures: Plot the CA output
% Copyright 2018 A. Gupta and M. Taylor
% Article: A. Gupta, G. Reint, I. Gozen, and M. Taylor, "A cellular automaton
% for modeling of non-trivial biomembrane ruptures"
% bioRxiv 429548; doi: https://doi.org/10.1101/429548
% ----------------------------------------------------------------------------

function F = getFigures(cellStatus, cellLocation, cellKey, automationSize, pinningProb, clusterProb, clusterThresh, percolationThresh, circleMaxRadius, circleCurrentRadius)


cmap = [...
    0.00 0.00 0.00;  %Background   || 1  (Black)
    0.75 0.00 0.00;  %Bilayer      || 2 (Brush Red)
    0.9 0.00 0.00;  %Pinned BL    || 3  (Red)
    0.50 0.00 0.00]; %Fractured BL || 4 (Dark Red)
colormap(cmap)


X = zeros(automationSize, automationSize);
Y = zeros(automationSize, automationSize);
C = zeros(automationSize, automationSize);

for i = 1:automationSize
    for j = 1:automationSize
        X(i, j) = cellLocation(cellKey(i,j), 1);
        Y(i, j) = cellLocation(cellKey(i,j), 2);
        if cellStatus(cellKey(i,j), 1) == 0    % cell is inactive
            C(i,j) = 0.5;
        elseif ((cellStatus(cellKey(i,j), 10) == 1)&&(cellStatus(cellKey(i,j), 4) == 0) ) ||  ((cellStatus(cellKey(i,j), 10) == 1)&&(cellStatus(cellKey(i,j), 5) == 1) )% cell is broken and not a pin or cell is broken and is a cluster pin
            C(i,j) = 3.5;
        elseif (cellStatus(cellKey(i,j), 4) == 1) % cell is pinned
            C(i,j) = 2.5;
        elseif cellStatus(cellKey(i,j), 10) == 1  % cell is broken
            C(i,j) = 3.5;
        else  % cell is active
            C(i,j) = 1.5;
        end
    end
end

surf(X,Y,C, 'EdgeColor', 'flat');
view(2)
axis square
caxis([0 4])
h = colorbar;
set(h, 'XTick', [0.5, 1.5, 2.5, 3.5])
set(h,'XTickLabel', {'Background',  'Bi-Layer', 'Pinning Sites', 'Fractured Area'})

set(gca,'fontsize', 16);
title({['Pin Probability: ',  num2str(pinningProb*100), '%'];...
        ['Cluster Probability: ', num2str(clusterProb*100), '%' ];...
        ['Cluster Threshold: ', num2str(clusterThresh*100), '%'];...
        ['Percolation Threshold: ', num2str(percolationThresh*100), '%']});
set(gca,'xtick',[]);
set(gca,'YTick',[]);
xlabel({['Maxium Radius: ' num2str(circleMaxRadius)];...
        [' Bi-layer Radius: ', ...
        num2str(round(circleCurrentRadius))]});

axis([-40 40 -40 40])
F = getframe(gcf);
