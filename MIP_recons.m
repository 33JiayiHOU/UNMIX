clc
clear all
close all

load('dictionary.mat');

%% Simulation
% randomly select K groups
prompt='Simulate K groups/elitists.\nLower limit: ';
ll=input(prompt);
prompt='Upper limit: ';
ul=input(prompt);
K=(ll:1:ul);

draw_rand_class=cell([1,length(K)]);
for j=1:length(K)
    rng(42);
    larger_class=find(nbelements>1); % group size larger than 1
    rand_class=randperm(length(larger_class));
    ind_rand_class=rand_class(1:K(j));
    draw_rand_class{j}=larger_class(ind_rand_class);
end

%% Choosing type of simulated data
str='0';
while str~='1'&& str~='2'
    prompt='\nType of simulated data:\n 1.Group;\n 2.Elitist.\nYour choice: ';
    str=input(prompt,'s');
    if str~='1'&& str~='2'
%         warning('\nNon-existent choice!')
        disp('Warning: Non-existent choice!')
%     else if str==''
% %         warning('Non-existent choice!')
%             disp('Warning: Non-existent choice!')
%         end
    end
end

draw_rand_spec=cell([1,length(K)]);
sim_a=zeros(Q,length(K));
Sim=zeros(N,length(K));
if str=='1' % Group: a few groups with highly density within each group 
    for j=1:length(K)
        draw_rand_spec_g=[];
        size_group=zeros(1,K(j));
        for i=1:K(j)
            index_spec=find(spec_group==draw_rand_class{j}(i));
            draw_rand_spec_g=[draw_rand_spec_g,index_spec];
            size_group(i)=length(index_spec); 
        end
        draw_rand_spec{j}=draw_rand_spec_g;
        
        a_g_grouped=rand(K(j),1);
        a_g_grouped=a_g_grouped/sum(a_g_grouped);
        % make sure that abundances of each group exceed 0.1
        a_g_grouped=a_g_grouped*(1-K(j)*0.1);
        a_g_grouped=a_g_grouped+0.1;
        a_g=[];
        for i=1:K(j)
            a_g_intragroup=rand(size_group(i),1);
            a_g_intragroup=a_g_grouped(i)*a_g_intragroup/sum(a_g_intragroup);
            a_g=[a_g;a_g_intragroup];
        end
        sim_a_g=zeros(Q,1);
        sim_a_g(draw_rand_spec_g)=a_g;
        sim_a(:,j)=sim_a_g;
        Sim(:,j)=H*sim_a_g;
    end
else if str=='2' % Elitist: at most one spectrum per group
        for j=1:length(K)
            draw_rand_spec_e=[];
            for i=1:K(j)
                index_spec=find(spec_group==draw_rand_class{j}(i));
                rand_ind=randperm(length(index_spec));
                rand_spec=index_spec(rand_ind(1));
                draw_rand_spec_e=[draw_rand_spec_e,rand_spec];
            end
            draw_rand_spec{j}=draw_rand_spec_e;
            
            a_e=rand(K(j),1);
            a_e=a_e/sum(a_e);
            % make sure that abundance of each spectrum exceeds 0.1
            a_e=a_e*(1-K(j)*0.1);
            a_e=a_e+0.1; 
            sim_a_e=zeros(Q,1);
            sim_a_e(draw_rand_spec_e)=a_e;
            sim_a(:,j)=sim_a_e;
            Sim(:,j)=H*sim_a_e;
        end
   end
end

%% Adding noise
prompt='\nAdd white Gaussian noise to signal: snr = ';
snr=input(prompt);
Sim_n=zeros(N,length(K));
for j=1:length(K)
    Sim_n(:,j)=awgn(Sim(:,j),snr,'measured');
end

%% Solving (Cplex for Matlab)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    min      0.5*x'*H*x+f*x or f*x    => cplexmiqp
%    min      norm(C*x-d)^2 or f*x    => cplexlsqmilp
%    st.      Aineq*x <= bineq 
%             Aeq*x    = beq
%             lb <= x <= ub
%  
%   x = cplexlp(f,Aineq,bineq,Aeq,beq,lb,ub);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kb=K;

if str=='1' % Group
    C=[zeros(N,Q),zeros(N,Nbclass),eye(N)];
    CTC=C'*C;

    Aeq=[ones(1,Q),zeros(1,Nbclass+N);H,zeros(N,Nbclass),-eye(N)];
    beq=[1;zeros(N,1)];

    Aineq1=eye(Q);
    Aineq2=zeros(Q,Nbclass);
    for i=1:Q
        Aineq2(i,spec_group(i))=-1;
    end
%     tau=1e-6;
%     Aineq=[Aineq1,Aineq2,zeros(Q,N);-Aineq1,-tau*Aineq2,zeros(Q,N);zeros(1,Q),ones(1,Nbclass),zeros(1,N)];
    Aineq=[Aineq1,Aineq2,zeros(Q,N);-Aineq1,zeros(Q,Nbclass),zeros(Q,N);zeros(1,Q),ones(1,Nbclass),zeros(1,N)];

    ctype=blanks(Q+Nbclass+N);
    ctype(1:Q)='C';
    ctype(Q+1:Q+Nbclass)='B';
    ctype(Q+Nbclass+1:end)='C';

%     options = cplexoptimset ('cplex');
%     options.optimalitytarget=3;
%     options.Display = 'on';
    options = cplexoptimset('Display','on','MaxIter',5000);
%     options = cplexoptimset('Display','on');

    x=zeros(Q+Nbclass+N,length(Kb));
    t=zeros(1,length(Kb));
    for j=1:length(Kb)
        tic;
        bineq=[zeros(2*Q,1);Kb(j)];

        [x(:,j), resnorm, residual, exitflag, output] = ...
            cplexlsqmilp(C,Sim_n(:,j),Aineq,bineq,Aeq,beq,...
            [],[],[],[],[],ctype,[],options);
%         f=-(C'*Sim_n(:,j))';
%         [x(:,j),fval,exitflag,output] = cplexmiqp(CTC,f,Aineq,bineq,Aeq,beq,...
%             [],[],[],[],[],ctype,[],options);
        t(j)=toc;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else if str=='2' % Elitist  
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

        ctype=blanks(Q*2+N);
        ctype(1:Q)='C';
        ctype(Q+1:2*Q)='B';
        ctype(2*Q+1:end)='C';

%         options = cplexoptimset ('cplex');
%         options.optimalitytarget=3;
%         options.Display = 'on';
        options = cplexoptimset('Display','on');

        x=zeros(Q*2+N,length(Kb));
        t=zeros(1,length(Kb));
        for j=1:length(Kb)
            tic;
            bineq=[zeros(2*Q,1);k*ones(Nbclass,1);Kb(j)];

            [x(:,j), resnorm, residual, exitflag, output] = ...
                cplexlsqmilp(C,Sim_n(:,j),Aineq,bineq,Aeq,beq,...
                [],[],[],[],[],ctype,[],options);
%             f=-(C'*Sim_n(:,j))';           
%             [x(:,j),fval,exitflag,output] = cplexmiqp(CTC,f,Aineq,bineq,Aeq,beq,...
%                 [],[],[],[],[],ctype,[],options);
            t(j)=toc;
        end    
   end
end

%% Estimation
ab=x(1:Q,:);
est=H*ab;

for j=1:length(K)
    figure
    plot(Sim(:,j));
    hold on
    plot(Sim_n(:,j),'g');
    hold on
    plot(est(:,j),'r');
    legend('ground truth','measurement','estimate');
    title(['Estimation ( K = ',num2str(K(j)),')']);
end

%% Abundances
ab_th=ab.*(ab>1e-3);

if str=='1'
    cum_part=cumsum(nbelements);
    for j=1:length(K)   
        axis_y=max([sim_a(:,j);ab_th(:,j)])+0.0005;
        Y=[0,axis_y,axis_y,0];
        start_ind=0.5;
        figure 
        for i=1:Nbclass
            if mod(i,2)==1
                X=[start_ind,start_ind,cum_part(i)+0.5,cum_part(i)+0.5];
                fill(X,Y,'g','EdgeColor','none');
                hold on
            end
            start_ind=cum_part(i)+0.5;
        end
        p1=stem(sim_a(:,j),'r');
        hold on
        p2=stem(ab_th(:,j),'x','DisplayName','estimate');
        legend([p1 p2],{'ground truth ','estimate'});
        title(['Abundances-G ( K = ',num2str(K(j)),')']);
    end
else if str=='2'
        for j=1:length(K)    
            figure 
            stem(sim_a(:,j),'r');
            hold on
            stem(ab_th(:,j),'x');
            legend('ground truth ','estimate');
            title(['Abundances-E ( K = ',num2str(K(j)),')']);
        end
    end
end

%% TP_TN_FP_FN - F1_PPV_FDR
est_class_g=cell([1,length(Kb)]);
est_spec_e=cell([1,length(Kb)]);

TP=cell([1,length(Kb)]);
TN=cell([1,length(Kb)]);
FP=cell([1,length(Kb)]);
FN=cell([1,length(Kb)]);
tp=zeros(1,length(Kb));
tn=zeros(1,length(Kb));
fp=zeros(1,length(Kb));
fn=zeros(1,length(Kb));
F1=zeros(1,length(Kb)); 
PPV=zeros(1,length(Kb));
FDR=zeros(1,length(Kb));

if str=='1'
    est_class_g=cell([1,length(Kb)]);
    for j=1:length(Kb)   
        est_class=spec_group'.*(ab_th(:,j)>0);
        est_class=unique(est_class);
        est_class=est_class(2:end);
        est_class_g{j}=est_class;

        TP_g=intersect(draw_rand_class{j},est_class);
        TN_g=setdiff(spec_group,union(draw_rand_class{j},est_class));
        FP_g=setdiff(est_class,TP_g); 
        FN_g=setdiff(draw_rand_class{j},TP_g); 
        TP{j}=TP_g;
        TN{j}=TN_g;
        FP{j}=FP_g;
        FN{j}=FN_g;

        tp(j)=length(TP_g);
        tn(j)=length(TN_g);
        fp(j)=length(FP_g);
        fn(j)=length(FN_g);

        F1(j)=2*tp(j)/(2*tp(j)+fp(j)+fn(j));    
        PPV(j)=tp(j)/(tp(j)+fp(j)); 
        FDR(j)=fp(j)/(tp(j)+fp(j)); 
    end
else if str=='2'
        est_spec_e=cell([1,length(Kb)]);
        for j=1:length(Kb)   
            est_spec=(1:1:Q)'.*(ab_th(:,j)>0);
            est_spec=unique(est_spec);
            est_spec=est_spec(2:end);
            est_spec_e{j}=est_spec;

            TP_e=intersect(draw_rand_spec{j},est_spec);
            TN_e=setdiff((1:1:Q),union(draw_rand_spec{j},est_spec));
            FP_e=setdiff(est_spec,TP_e); 
            FN_e=setdiff(draw_rand_spec{j},TP_e); 
            TP{j}=TP_e;
            TN{j}=TN_e;
            FP{j}=FP_e;
            FN{j}=FN_e;

            tp(j)=length(TP_e);
            tn(j)=length(TN_e);
            fp(j)=length(FP_e);
            fn(j)=length(FN_e);

            F1(j)=2*tp(j)/(2*tp(j)+fp(j)+fn(j));    
            PPV(j)=tp(j)/(tp(j)+fp(j)); 
            FDR(j)=fp(j)/(tp(j)+fp(j)); 
        end
    end
end

%% Save
if str=='1'
    prompt=['MIP_G_snr',num2str(snr),'_K(',num2str(ll),'-',num2str(ul),').mat'];
else if str=='2'
        prompt=['MIP_E_snr',num2str(snr),'_K(',num2str(ll),'-',num2str(ul),').mat'];
    end
end

save(prompt)


