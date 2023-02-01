function [wettedArea, fractureArea] = computeAreas(cellStatus, numCells)

wettedArea = 0;
fractureArea = 0;

for i = 1:numCells
    if cellStatus(i, 1) == 1
        wettedArea = wettedArea + 1;
    end
    
    if (cellStatus(i, 1) == 1) && (cellStatus(i,10) == 1)
        fractureArea = fractureArea + 1;
    end
end