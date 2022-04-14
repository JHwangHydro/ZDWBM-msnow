function [nse] = nash_sutcliffe(modeled,observed)
%Nash Sutcliffe Efficiency measure

diff=modeled-observed; %start after warmup period
diffmeanflow=observed-(mean(observed)); 
nse= 1-((sum(diff.^2))/(sum(diffmeanflow.^2)));

end