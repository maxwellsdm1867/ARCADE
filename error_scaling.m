function error_vec = error_scaling(in_range,calcuim_trace,A,V_conj,S_r)
error_vec = zeros(1,length(in_range));
for i = 1:length(in_range)
    r = in_range(i);
    V_r = V_conj(1:r,:);%reduce to r-dimensions
    S_r = S(1:r,1:r);%reduce to r-dimensions
    A = S_r*V_r;%project to r-dimensions
    singular_value_rank = 1;
    [B,FitInfo] = lasso(A',calcuim_trace);%lasso
    Bb= B(:,singular_value_rank)';
    sr = length((Bb*A));
    y_hat = (Bb*A)+FitInfo.Intercept(1)*ones(1,sr);
    er1 = (calcuim_trace-y_hat').*(calcuim_trace-y_hat');
    er2 =sum(sqrt( er1))/sr;
    error_vec(i) = er2;
end
end