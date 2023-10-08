%_________________________________________________________________________%
%  Whale Optimization Algorithm (WOA) source codes demo 1.0               %
%                                                                         %
%  Developed in MATLAB R2011b(7.13)                                       %
%                                                                         %
%  Author and programmer: Seyedali Mirjalili                              %
%                                                                         %
%         e-Mail: ali.mirjalili@gmail.com                                 %
%                 seyedali.mirjalili@griffithuni.edu.au                   %
%                                                                         %
%       Homepage: http://www.alimirjalili.com                             %
%                                                                         %
%   Main paper: S. Mirjalili, A. Lewis                                    %
%               The Whale Optimization Algorithm,                         %
%               Advances in Engineering Software , in press,              %
%               DOI: http://dx.doi.org/10.1016/j.advengsoft.2016.01.008   %
%                                                                         %
%_________________________________________________________________________%


% The Whale Optimization Algorithm
function [Leader_score,Leader_pos,Convergence_curve]=WOA(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

% initialize position vector and score for the leader
Leader_pos=zeros(1,dim);
Leader_score=inf; %change this to -inf for maximization problems


%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);

Convergence_curve=zeros(1,Max_iter);

t=1;% Loop counter

% Main loop
while t<Max_iter
    for i=1:size(Positions,1)%1——>30
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;%第一个向量
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        % Calculate objective function for each search agent
        fitness=fobj(Positions(i,:));%计算出当前函数值
        % Update the leader
        if fitness<Leader_score % Change this to > for maximization problem
            Leader_score=fitness; % Update alpha
            Leader_pos=Positions(i,:);
        end   
    end
    a=2.*(1-power((t/Max_iter),2))%a=2-t*((2)/Max_iter); % a decreases linearly fron 2 to 0 in Eq. (2.3)
    % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a2=-1+t*((-1)/Max_iter);
    % Update the Position of search agents 
    for i=1:size(Positions,1)
        r1=rand(); % r1 is a random number in [0,1]
        r2=rand(); % r2 is a random number in [0,1]
        A=2*a*r1-a;  % Eq. (2.3) in the paper
        C=2*r2;      % Eq. (2.4) in the paper
        b=1;               %  parameters in Eq. (2.5)
        l=(a2-1)*rand+1;   %  parameters in Eq. (2.5) 
        p = rand();        % p in Eq. (2.6)
        %针对改进2，加入的——start
        f_f=zeros(1,SearchAgents_no);% 申请一个1*SearchAgents_no的元素均为0的矩阵
        % 注意： 打印数组时1.0e+05 *表示矩阵里的数都要乘上这个数
        for k=1:SearchAgents_no% 计算每个搜索代理的适应度，并存入f_f数组
            tempfobj(k)=fobj(Positions(k,:));
            f_f(k)=tempfobj(k);%f_f(1)是2.1672e+19
        end    
        %2.搜索代理适应度排序:sort函数
        f_f=sort(f_f);%对数组中的元素升序排列
        %3.求出较小段的平均适应度和较大段的平均适应度
        middle=round(SearchAgents_no/2);%取排序后的适应度的中间位置坐标
        temp1=f_f(1:middle);%暂取小适应度部分数组
        low_aver=mean(temp1);
        temp2=f_f(middle+1:SearchAgents_no);
        high_aver=mean(temp2);
        %4.计算原始权值：（1）计算适应度最小的位置坐标：low_value=f_f(1)（2）计算适应度最大的位置坐标high_value=f_f(SearchAgents_no)
        low_value=f_f(1);
        high_value=f_f(SearchAgents_no);
        best_pos=[];
        worse_pos=[];
        for n=1:SearchAgents_no% 有可能存在适应度一样的位置
            if tempfobj(n)==low_value
                best_pos=Positions(n,:);
            end
            if tempfobj(n)==high_value
                worse_pos=Positions(n,:);
            end
        end
       %当仅仅改了自适应权重时，两个参数均等于0.048时挺好的，0.002时也挺好的。0.05时结果是0.023359，0.052时0.013213，0.053时0.018557
       %0.054时0.025808
        reduce_pos=worse_pos-best_pos;
        sm=sum(reduce_pos.*reduce_pos,2);
        sk=power(sm,0.5);
        d1=0.005;
        d2=0.005;
        w_original=sk.*d1+d2.*(ub-lb)/t;
        %5.判断当前个体的适应度所处范围后，计算当前个体的扰动权值
        if f_f(i)<=low_aver
            w=w_original*1.1;
        end
        if f_f(i)>low_aver&&f_f(i)<high_aver
            w=w_original*1;
        end
        if f_f(i)>=high_aver
            w=w_original.*0.8;
        end
        %针对改进2，新加入的——end
        for j=1:size(Positions,2)
            if p<0.5   
                if abs(A)>=1
                    rand_leader_index = floor(SearchAgents_no*rand()+1);
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand=abs(C*X_rand(j)-Positions(i,j)); % Eq. (2.7)
                    Positions(i,j)=X_rand(j)-A*D_X_rand;      % Eq. (2.8)
                elseif abs(A)<1
                    D_Leader=abs(C*Leader_pos(j)-Positions(i,j)); % Eq. (2.1)
                    Positions(i,j)=w*Leader_pos(j)-A*D_Leader;      % Eq. (2.2)
                end
            elseif p>=0.5
                distance2Leader=abs(Leader_pos(j)-Positions(i,j));
                % Eq. (2.5)
                Positions(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+w*Leader_pos(j);%exp(1)就是常数e= 2.7183   
            end
        end
    end
    t=t+1;
    Convergence_curve(t)=Leader_score;
    [t Leader_score]
end



