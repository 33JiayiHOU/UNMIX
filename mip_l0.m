function [a,est,t] = mip_l0(H,y,Kb,Kc,tau,nbelements,spec_group,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    min      0.5*x'*H*x+f*x or f*x    => cplexmiqp
%    min      norm(C*x-d)^2 or f*x    => cplexlsqmilp
%    st.      Aineq*x <= bineq 
%             Aeq*x    = beq
%             lb <= x <= ub
%  
%   x = cplexlp(f,Aineq,bineq,Aeq,beq,lb,ub);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input:  H: dictionary
%         y: data / hyperspectral image cube [nbands x nsamples]
%         Kb: number of endmembers
%         Kc: number of groups
%         tau: threshold of endmembers
%         nbelements: number of endmembers in each group
%         spec_group: index of group corresponding to each endmember 
%         type: 'group', 'elitist' or 'gs'(global sparsity) or 'mgs'(modified gs)
% Output: a: abundance matrix
%         est: estimates
%         t: elapsed time of each pixel

[N,Q]=size(H);
[l,p]=size(y);
Nbclass=length(nbelements);

if strcmp(type,'group')
    C=[zeros(N,Q),zeros(N,Nbclass),eye(N)];
    CTC=C'*C;

    Aeq=[ones(1,Q),zeros(1,Nbclass+N);H,zeros(N,Nbclass),-eye(N)];
    beq=[1;zeros(N,1)];

    Aineq1=eye(Q);
    Aineq2=zeros(Q,Nbclass);
    for i=1:Q
        Aineq2(i,spec_group(i))=-1;
    end
    Aineq=[Aineq1,Aineq2,zeros(Q,N);-Aineq1,zeros(Q,Nbclass),zeros(Q,N);zeros(1,Q),ones(1,Nbclass),zeros(1,N)];
    bineq=[zeros(2*Q,1);Kc];

    ctype=blanks(Q+Nbclass+N);
    ctype(1:Q)='C';
    ctype(Q+1:Q+Nbclass)='B';
    ctype(Q+Nbclass+1:end)='C';

%     options = cplexoptimset('Display','off','MaxIter',5000);
    options = cplexoptimset('Display','off','MaxIter',1000);
    
    x=zeros(Q+Nbclass+N,p);
    t=zeros(1,p);
    for j=1:p
        tic;
        [x(:,j), resnorm, residual, exitflag, output] = ...
            cplexlsqmilp(C,y(:,j),Aineq,bineq,Aeq,beq,...
            [],[],[],[],[],ctype,[],options);
%         f=-(C'*y(:,j))';
%         [x(:,j),fval,exitflag,output] = cplexmiqp(CTC,f,Aineq,bineq,Aeq,beq,...
%             [],[],[],[],[],ctype,[],options);
        t(j)=toc;
    end
    
    
    
elseif strcmp(type,'elitist')
    C=[zeros(N,2*Q),eye(N)];
    CTC=C'*C;
    
    Aeq=[ones(1,Q),zeros(1,Q+N);H,zeros(N,Q),-eye(N)];
    beq=[1;zeros(N,1)];
    
    Aineq1=eye(Q);
    Aineq2=zeros(Nbclass,Q);
    cum_part=cumsum(nbelements);
    start_ind=1;
    for i=1:Nbclass
        sel=start_ind:cum_part(i);
        Aineq2(i,sel)=1;
        start_ind=cum_part(i)+1;
    end
    Aineq=[Aineq1,-Aineq1,zeros(Q,N);-Aineq1,zeros(Q,Q+N);...
        zeros(Nbclass,Q),Aineq2,zeros(Nbclass,N);zeros(1,Q),ones(1,Q),zeros(1,N)];
    k=1;  % number of the nonzero elements in each group
    bineq=[zeros(2*Q,1);k*ones(Nbclass,1);Kb];
    
    ctype=blanks(Q*2+N);
    ctype(1:Q)='C';
    ctype(Q+1:2*Q)='B';
    ctype(2*Q+1:end)='C';

%     options = cplexoptimset('Display','off');
    options = cplexoptimset('Display','off','MaxIter',1000);
    
    x=zeros(Q*2+N,p);
    t=zeros(1,p);
    for j=1:p
        tic;
        [x(:,j), resnorm, residual, exitflag, output] = ...
            cplexlsqmilp(C,y(:,j),Aineq,bineq,Aeq,beq,...
            [],[],[],[],[],ctype,[],options);
%             f=-(C'*y(:,j))';
%             [x(:,j),fval,exitflag,output] = cplexmiqp(CTC,f,Aineq,bineq,Aeq,beq,...
%                 [],[],[],[],[],ctype,[],options);
        t(j)=toc;
    end
    
    
elseif strcmp(type,'gs')
    C=[zeros(N,2*Q+Nbclass),eye(N)];
    CTC=C'*C;
    
    Aeq=[ones(1,Q),zeros(1,Q+Nbclass+N);H,zeros(N,Q+Nbclass),-eye(N)];
    beq=[1;zeros(N,1)];
    
    Aineq1=eye(Q);
    Aineq2=zeros(Q,Nbclass);
    for i=1:Q
        Aineq2(i,spec_group(i))=-1;
    end
    Aineq=[Aineq1,-Aineq1,zeros(Q,Nbclass+N);...
        -Aineq1,zeros(Q,Q+Nbclass+N);...
        zeros(1,Q),ones(1,Q),zeros(1,Nbclass+N);...
        zeros(Q,Q),Aineq1,Aineq2,zeros(Q,N);...
        zeros(1,2*Q),ones(1,Nbclass),zeros(1,N)];
    bineq=[zeros(2*Q,1);Kb;zeros(Q,1);Kc];
    
    ctype=blanks(2*Q+Nbclass+N);
    ctype(1:Q)='C';
    ctype(Q+1:2*Q+Nbclass)='B';
    ctype(2*Q+Nbclass+1:end)='C';
    
    %     options = cplexoptimset('Display','off');
    options = cplexoptimset('Display','off','MaxIter',1000);
    
    x=zeros(Q*2+Nbclass+N,p);
    t=zeros(1,p);
    for j=1:p
        tic;
        [x(:,j), resnorm, residual, exitflag, output] = ...
            cplexlsqmilp(C,y(:,j),Aineq,bineq,Aeq,beq,...
            [],[],[],[],[],ctype,[],options);
%         f=-(C'*y(:,j))';
%         [x(:,j),fval,exitflag,output] = cplexmiqp(CTC,f,Aineq,bineq,Aeq,beq,...
%             [],[],[],[],[],ctype,[],options);
        t(j)=toc;
    end
  
    
elseif strcmp(type,'mgs')
    C=[zeros(N,2*Q+Nbclass),eye(N)];
    CTC=C'*C;
    
    Aeq=[ones(1,Q),zeros(1,Q+Nbclass+N);H,zeros(N,Q+Nbclass),-eye(N)];
    beq=[1;zeros(N,1)];
    
    Aineq1=eye(Q);
    Aineq2=zeros(Q,Nbclass);
    for i=1:Q
        Aineq2(i,spec_group(i))=-1;
    end
    Aineq=[Aineq1,-Aineq1,zeros(Q,Nbclass+N);...
        -Aineq1,tau*Aineq1,zeros(Q,Nbclass+N);...
        zeros(Q,Q),Aineq1,Aineq2,zeros(Q,N);...
        zeros(1,2*Q),ones(1,Nbclass),zeros(1,N)];
    bineq=[zeros(3*Q,1);Kc];
    
    ctype=blanks(2*Q+Nbclass+N);
    ctype(1:Q)='C';
    ctype(Q+1:2*Q+Nbclass)='B';
    ctype(2*Q+Nbclass+1:end)='C';
    
    %     options = cplexoptimset('Display','off');
    options = cplexoptimset('Display','off','MaxIter',1000);
    
    x=zeros(Q*2+Nbclass+N,p);
    t=zeros(1,p);
    for j=1:p
        tic;
        [x(:,j), resnorm, residual, exitflag, output] = ...
            cplexlsqmilp(C,y(:,j),Aineq,bineq,Aeq,beq,...
            [],[],[],[],[],ctype,[],options);
%         f=-(C'*y(:,j))';
%         [x(:,j),fval,exitflag,output] = cplexmiqp(CTC,f,Aineq,bineq,Aeq,beq,...
%             [],[],[],[],[],ctype,[],options);
        t(j)=toc;
    end
    
    
elseif strcmp(type,'group')==0 && strcmp(type,'elitist')==0 ...
        && strcmp(type,'gs')==0 && strcmp(type,'mgs')==0
    warning('Non-existent type!')
    return
end

a=x(1:Q,:);
% a=a.*(a>1e-3);
est=H*a;

end

