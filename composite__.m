function [model]= composite__(xxall, gall_t)
%model.alpha = 2;
k_type = ["pow_exp"];
dim_io=size(xxall,2);

options.kernel_type=repmat(k_type,dim_io,1);
options.alpha=repmat(2.0,dim_io,1);
options.nugget_est = true;
options.zero_mean  = false;
options.trend=[ones(size(xxall,1),1), xxall];
model=ppgasp(xxall, gall_t, options);

end