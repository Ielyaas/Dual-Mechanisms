test_length = ceil(length(T)/2);
t_test = T(test_length:end);
ca_test = Y(test_length:end,1);
% 
% figure
% plot(t_test,ca_test)

ca_diff_1 = diff(ca_test);
% figure
% plot(t_test(2:end),ca_diff_1)
max_diff = max(ca_diff_1);
diff_tol = 1e-7;

max_idx_init = find(ca_diff_1>(max_diff-diff_tol));
peak_times = t_test(max_idx_init);
period_est = diff(peak_times);

PEAK = findpeaks(c);


INDEX = find(ismember(c,PEAK));