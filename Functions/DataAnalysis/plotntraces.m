function plotntraces(traces)
    figure; hold on;
    xvals=1:size(traces,2);
    for c=1:size(traces,1)
        line(xvals,traces(c,:));
    end
end