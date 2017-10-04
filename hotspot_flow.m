% Sherief Reda (Brown University) and Adel Belouchrani (ENP)
% "Blind Identification of Power Sources in Processors", in IEEE/ACM Design, Automation & Test in Europe, 2017.
% sherief_reda@Brown.edu and adel.belouchrani@enp.edu.dz
clear
errors=[];
num_cores=4;
n_samples=280;
 
% offline-learning
% sampling rate is 1 second per sample
% estimate steady-state matrix R from steady_state data 
% find matrix A : T[k]=AT[k-1]+Bp[k]
A = find_A('DATA_HOTSPOT_4CORES/nat_trace_', n_samples-1, 0);

% find matrix B : T[k]=AT[k-1]+Bp[k]
R = find_R('DATA_HOTSPOT_4CORES/steady_state.txt', n_samples-1, 0);
B=(eye(size(A))-A)*R;

% prep the data for online evaluation
eval_data=load('DATA_HOTSPOT_4CORES/eval_cores');
p_act=eval_data(:, 1:num_cores);
MS=size(eval_data(:, num_cores+1:2*num_cores));
trace=[sum(eval_data(:, 1:num_cores), 2), eval_data(:, num_cores+1:2*num_cores)];
dlmwrite('DATA_HOTSPOT_4CORES/eval_cores_temp', trace, 'delimiter', '\t', 'precision','%.4f');   

% on-line evaluation: use the model to identify p
p=eval_runtime('DATA_HOTSPOT_4CORES/eval_cores_temp', A, B);
p = p';
err=0;
err_perc=0;
maxp=max(p_act(:));
for j=1:size(p, 2)
    err = err+mean(abs(p_act(2:end, j)-p(:, j)));
    err_perc = err_perc+mean(abs(p_act(2:end, j)-p(:, j)))/maxp;
end
num_cores=size(p,2);
err=err/num_cores;
err_perc=err_perc/num_cores;
errors=[errors; err 100*err_perc]

% plot power and thermal traces
T_amb=25;
n_samples = 280;
X1=eval_data(:, 5:8);
t=1:length(eval_data);
figure
subplot(5, 1, 1);
plot(t,X1+T_amb, 'linewidth', 2);
grid on;
ylabel('Temperature (C)');
xlabel('time (s)');
legend('core 1', 'core 2', 'core 3', 'core 4')
xlim([0 1500]);
t=t(1:end-1);
title('(a) thermal simulation results from HotSpot', 'fontsize', 10);
subplot(5, 1, 2);
plot(t, p(:,1), 'b', 'linewidth', 2);
grid on;
ylabel('Power(W)');
xlabel('time (s)');
title('(b) estimated and actual power consumption of core 1', 'fontsize', 10);
hold on;
ylim([0 30]);
xlim([0 1500]);
plot(t, p_act(2:end,1), '--r', 'linewidth', 2);
grid on;
legend('estimated', 'actual')
subplot(5, 1, 3);
plot(t, p(:,2), 'b', 'linewidth', 2);
grid on;
ylabel('Power(W)');
xlabel('time (s)');
title('(c) estimated and actual power consumption of core 2', 'fontsize', 10);
hold on;
xlim([0 1500]);
ylim([0 30]);
plot(t, p_act(2:end,2), '--r', 'linewidth', 2);
grid on;
legend('estimated', 'actual')
subplot(5, 1, 4);
plot(t, p(:,3), 'b', 'linewidth', 2);
grid on;
legend('estimated', 'actual')
ylim([0 30]);
hold on;
plot(t, p_act(2:end,3), '--r', 'linewidth', 2);
grid on;
title('(d) estimated and actual power consumption of core 3', 'fontsize', 10);
ylim([0 30]);
ylabel('Power(W)');
xlabel('time (s)');
xlim([0 1500]);
legend('estimated', 'actual')
subplot(5, 1, 5);
plot(t, p(:,4), 'b', 'linewidth', 2);
grid on;
hold on;
title('(e) estimated and actual power consumption of core 4', 'fontsize', 10);
plot(t, p_act(2:end,4), '--r', 'linewidth', 2);
ylim([0 30]);
ylabel('Power(W)');
xlabel('time (s)');
xlim([0 1500]);
legend('estimated', 'actual')

