function [pred_model] =predict_intgasp(model_f,model_g,newu,newy)
% This function predicts the output of a linked-emulator based on 
% two Gaussian process models model_f and model_g, for given inputs newu, newy.

% INPUTS:
% - model_f: Contains input/output data, parameters, and linear fits for model f.
% - model_g: Contains input/output data, parameters, and linear fits for model g.
% - newu: New user input values for the internal model.
% - newy: New user input values for the external model.

% OUTPUT:
% - pred_model: Structure with fields 'mean' and 'var' representing the predictive mean 
%   and variance respectively for the new user inputs.

% Set output format to long for high precision.
format long;
% notice that model_g.alpha(1) is a vector of a single value, so only the
% first element of that is used (i.e., model_g.alpha(1)(1))

m=size(model_f.input,1);
M=size(newu,1);

% Dimensionality and model structure extraction
dimu=size(model_f.input,2);
dimz=size(model_f.output,2);
dimy=size(model_g.input,2)-dimz; % after ppgasp, the independent input to g
% and the dependent input from g are
% combined.

param=model_g.range_par;
paramy=param(1:dimy);
paramz=param(dimy+1:end);

% setup design matrix
xdg=model_g.input;              % again, the independent input to g and the
% dependent input from g are combined after
% using ppgasp.

Hxdg=[ones(size(xdg,1),1) xdg]; % here, note that the order of the design
% matrix is different than what it is in
% Ming and Guillas (2020).
%
% here, Hxdg=[1 y(input direct to g) z(i/o from f)]
% or Hxdg=[H(z),w], where H(z)=[1,y] and w[z];

len=dimu+1;
mu_u=zeros(m,dimz);
mu_newu=zeros(M,dimz);

% Calculate the mean trend of PPLE

for dd=1:dimz
    mu_u(:,dd)=model_f.theta_hat(1,dd);
    mu_newu(:,dd)=model_f.theta_hat(1,dd);
    for kk=1:len-1
        
        mu_u(:,dd)    = mu_u(:,dd)    +   model_f.theta_hat(kk+1,dd)*model_f.input(:,kk);
        % unconditional mean of z(:,dd) at u
        % this corresponds to H(x)*bhat in (EQ 4) in Ming and Guillas (2020)
        
        mu_newu(:,dd) = mu_newu(:,dd) +   model_f.theta_hat(kk+1,dd)*newu(:,kk);
        % unconditional mean of z(:,dd) at newu
        % this corresponds to h(x_0)'*bhat in (EQ 4) in Ming and Guillas (2020)
    end
end

% Compute the predictive mean of the PPLE
mu0_newu=zeros(M,dimz);
sig2_newu=zeros(M,dimz);

Ru=kernel_PE(model_f.input,model_f.input,model_f.alpha(1),model_f.range_par);

if model_f.nugget<1e-7
    invRu=inv(Ru+eye(size(Ru))*(1e-6));
else
    invRu=inv(Ru+eye(size(Ru))*(model_f.nugget));
end

for d=1:dimz
    for n=1:size(newu,1)                                     %model_f.param_f(:,d)
        
        r=kernel_PE(model_f.input,newu(n,:),model_g.alpha(1),model_f.range_par);
        mu0_newu(n,d)=mu_newu(n,d) + r'*invRu*(model_f.output(:,d) - mu_u(:,d));
        
        % \mu in Ming and Guillas (2020)
        
        sig2_newu(n,d)=model_f.sigma2_hat(d).^2 * ...
            abs((1 - r'*invRu*r + (mu_newu(n,d) - mu_u(:,d)'*invRu*r)'*...
            inv(mu_u(:,d)'*invRu*mu_u(:,d))*(mu_newu(n,d) - mu_u(:,d)'*invRu*r)));
        
        % \sigma^2 in Ming and Guillas (2020)
    end
end
% \sigma^2 in Ming and Guillas (2020)
% find conditional mean of g given data

Rz=kernel_PE([model_g.input],[model_g.input],model_g.alpha(1),param);

if model_g.nugget<1e-6
    invRz=inv(Rz+eye(size(Rz))*(1e-6));
else
    invRz=inv(Rz+eye(size(Rz))*(model_g.nugget));
end

%;
thetas_g=inv(Hxdg'*invRz*Hxdg)*(Hxdg'*invRz*model_g.output);

% thetas_g corresponds to [theta_hat, beta_hat] in Ming and Guillas (2020)
%
% NOTE: the order is different than what it is
% in Ming and Guillas (2020)
%%%%%%%%%Correction in the dimension%%%%%%%%%%%%
mug=@(x) thetas_g(1,:)'+thetas_g(2:end,:)'*x';

% calculate I'A

A=invRz*(model_g.output - Hxdg*thetas_g);
c=ones(M,m);
xi=ones(M,m);

for ii=1:m
    for j=1:dimy
        c(:,ii)=c(:,ii) .* kernel_PE(newy(:,j),model_g.input(ii,j),model_g.alpha(j),paramy(j));
    end
    
    for jj=1:dimz
        xi(:,ii)=xi(:,ii) .* (1./sqrt(1+2*sig2_newu(:,jj)/paramz(jj)^2)) ...
            .* exp(-(mu0_newu(:,jj) - model_f.output(ii,jj)).^2 ./ (2*sig2_newu(:,jj) + paramz(jj)^2));
    end
end

I=c.*xi;

% EQ 9 in Ming and Guillas (2020)
pred_model.mean=mug([newy, mu0_newu])' + I*A;

% variance calculation
%
% calculate V1
%
% calculate J
cc=ones(m,m,M);
zeta=ones(m,m,M);

for mtest=1:M % element out of a M x 1 vector (a test point at each time)
    for k=1:dimy
        cc1=kernel_PE(newy(mtest,k),model_g.input(:,k),model_g.alpha(k),paramy(k))';
        cc1=repmat(cc1,1,m);
        cc2=cc1';
        cc(:,:,mtest)=cc(:,:,mtest) .* cc1 .* cc2;
    end
end

count=0;
for i=1:m
    for j=1:m
        count=count+1;
        z_sum(count,:)=model_f.output(i,:)+model_f.output(j,:);
        z_subt(count,:)=model_f.output(i,:)-model_f.output(j,:);
    end
end
z_sum=reshape(z_sum,m,m,dimz);
z_subt=reshape(z_subt,m,m,dimz);

for mtest=1:M
    for kk=1:dimz
        zeta(:,:,mtest)=zeta(:,:,mtest) .* (1./sqrt(1+4*sig2_newu(mtest,kk)/paramz(kk).^2)) ...
            .* exp(-(z_sum(:,:,kk)./2 - mu0_newu(mtest,kk)).^2 ./ (paramz(kk).^2./2 + 2*sig2_newu(mtest,kk)) ...
            - (z_subt(:,:,kk).^2 ./ (2*paramz(kk).^2)));
    end
end

J=cc.*zeta;

% rearrange I
I=I';
%One last nested loop for calculating B
for mtest=1:M
    for l=1:dimz
        for j=1:m
            phi=(1/sqrt(1+2*sig2_newu(mtest,l)/paramz(l)^2)) ...
                * exp(-(mu0_newu(mtest,l) - model_f.output(j,l))^2 / (2*sig2_newu(mtest,l) + paramz(l)^2)) ...
                * (2*sig2_newu(mtest,l)*model_f.output(j,l) + paramz(l)^2*mu0_newu(mtest,l)) ...
                / (2*sig2_newu(mtest,l) + paramz(l)^2);
            
            xi_jk=1;
            for k=1:dimz
                if k~=l
                    xi_jk=xi_jk * (1/sqrt(1+2*sig2_newu(mtest,k)/paramz(k)^2)) ...
                        * exp(-(mu0_newu(mtest,k) - model_f.output(j,k))^2 ...
                        / (2*sig2_newu(mtest,k) + paramz(k)^2));
                end
            end
            
            c_k=1;
            for k=1:dimy
                c_k=c_k * kernel_PE(newy(mtest,k),model_g.input(j,k),model_g.alpha(k),paramy(k));
            end
            B(l,j,mtest)=phi.*xi_jk.*c_k;
        end
    end
end

% calculate V1
omega=zeros(dimz,dimz,M);
clear V1;

%%%%%%%%%%%%%%%%VARIANCE CALCULATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_V1 = [];
for i = 1:size(thetas_g(end-dimz+1:end,:),2)
    %for i = 1:2
    temp_V1(i) = A(:,i)'*(J-I*I')*A(:,i) + 2*(thetas_g(end-dimz+1:end,i))'*(B - mu0_newu'*I')*A(:,i);       
end

for i=1:M
    omega(:,:,i)=diag(sig2_newu(i,:));
    
    V1(:,:,i)=diag(A'*(squeeze(J(:,:,i)-I(:,i)*I(:,i)'))*A ...
        + 2*(thetas_g(end-dimz+1:end,:)' * (B(:,:,i) - mu0_newu(i,:)'*I(:,i)') * A));      
end

for k = 1:size(thetas_g(end-dimz+1:end,:),2)
    %Trace needs to be calculated independently for a vector valued function
    trace_V1(k,:) = trace(thetas_g(end-dimz+1:end,k)*thetas_g(end-dimz+1:end,k)'*omega(:,:,M));
end

Q=invRz*Hxdg*inv(Hxdg'*invRz*Hxdg)*Hxdg'*invRz - invRz;

for i=1:M
    hz=[ones(1,1),newy(i,:)]';
    G(:,:,i)=[ones(1,1),newy(i,:),mu0_newu(i,:)]';      % G=[h(z)', \mu']
    % NOTE: G=[\mu', h(z)']' in Ming and Guillas (2020)
end

% calculate C
C=inv(Hxdg'*invRz*Hxdg);

% calculate K
clear K
for i=1:M
    K(:,:,i)=[I(:,i), (I(:,i)*newy(i,:)),B(:,:,i)'];     % K=[I*h(z)', B']
    % NOTE: K=[B', I*h(z)'] in Ming and Guillas (2020)
end

% build matrix P
sz=size(C,1);
zero=zeros(sz-size(omega,2),sz-size(omega,2));
for i=1:M
    P(:,:,i)=blkdiag(omega(:,:,i),zero);
end

% calculate V2
clear V2

sigma_sqred = diag(sigma_sqr(Hxdg, invRz, model_g.output,thetas_g, size(model_g.output,1),size(thetas_g,1)));

for i=1:M    
    V2(:,:,i) = sigma_sqred.*(1 + model_g.nugget + trace(Q*J(:,:,i)) ...
        + G(:,:,i)'*C*G(:,:,i) ...
        + trace(C*P(:,:,i) - 2*C*Hxdg'*invRz*K(:,:,i)));    
end

V1 = V1 + trace_V1;
pred_model.var= V1 + V2;                                  % EQ 10 in Ming and Guillas (2020)
end