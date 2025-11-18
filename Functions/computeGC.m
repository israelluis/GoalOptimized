function [gait_cycle,time_series] = computeGC(Misc_time,extra_frame)
% get time as in one frame equals 0.01 s
time        =Misc_time;
time_series =time(1):0.01:time(end);
data_length =length(time_series);

% conversion to gait cycle
initial=0;   
final=100;

data_GC      = data_length-2*extra_frame;
frames_per_GC= final/(data_GC-1);
extra_times  = frames_per_GC*extra_frame;
gait_cycle   = initial-extra_times:frames_per_GC:final+extra_times; % case without extra frame-> gait_cycle =linspace(0,100,length(time_series));
end