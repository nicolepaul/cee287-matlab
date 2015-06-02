function plotSquare(xvals,heights,options)
    % Same as line plot except everything is square
    % Assumes that xvals and heights are column vectors
    xvals = reshape(vertcat(xvals',xvals'),1,2*length(xvals));
    heights = reshape(vertcat(heights',heights'),1,2*length(heights));
    plot([0 xvals], [heights 0],options)
end