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
    for i=1:size(Positions,1)%1����>30
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;%��һ������
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        % Calculate objective function for each search agent
        fitness=fobj(Positions(i,:));%�������ǰ����ֵ
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
        %��ԸĽ�2������ġ���start
        f_f=zeros(1,SearchAgents_no);% ����һ��1*SearchAgents_no��Ԫ�ؾ�Ϊ0�ľ���
        % ע�⣺ ��ӡ����ʱ1.0e+05 *��ʾ�����������Ҫ���������
        for k=1:SearchAgents_no% ����ÿ�������������Ӧ�ȣ�������f_f����
            tempfobj(k)=fobj(Positions(k,:));
            f_f(k)=tempfobj(k);%f_f(1)��2.1672e+19
        end    
        %2.����������Ӧ������:sort����
        f_f=sort(f_f);%�������е�Ԫ����������
        %3.�����С�ε�ƽ����Ӧ�Ⱥͽϴ�ε�ƽ����Ӧ��
        middle=round(SearchAgents_no/2);%ȡ��������Ӧ�ȵ��м�λ������
        temp1=f_f(1:middle);%��ȡС��Ӧ�Ȳ�������
        low_aver=mean(temp1);
        temp2=f_f(middle+1:SearchAgents_no);
        high_aver=mean(temp2);
        %4.����ԭʼȨֵ����1��������Ӧ����С��λ�����꣺low_value=f_f(1)��2��������Ӧ������λ������high_value=f_f(SearchAgents_no)
        low_value=f_f(1);
        high_value=f_f(SearchAgents_no);
        best_pos=[];
        worse_pos=[];
        for n=1:SearchAgents_no% �п��ܴ�����Ӧ��һ����λ��
            if tempfobj(n)==low_value
                best_pos=Positions(n,:);
            end
            if tempfobj(n)==high_value
                worse_pos=Positions(n,:);
            end
        end
       %��������������ӦȨ��ʱ����������������0.048ʱͦ�õģ�0.002ʱҲͦ�õġ�0.05ʱ�����0.023359��0.052ʱ0.013213��0.053ʱ0.018557
       %0.054ʱ0.025808
        reduce_pos=worse_pos-best_pos;
        sm=sum(reduce_pos.*reduce_pos,2);
        sk=power(sm,0.5);
        d1=0.005;
        d2=0.005;
        w_original=sk.*d1+d2.*(ub-lb)/t;
        %5.�жϵ�ǰ�������Ӧ��������Χ�󣬼��㵱ǰ������Ŷ�Ȩֵ
        if f_f(i)<=low_aver
            w=w_original*1.1;
        end
        if f_f(i)>low_aver&&f_f(i)<high_aver
            w=w_original*1;
        end
        if f_f(i)>=high_aver
            w=w_original.*0.8;
        end
        %��ԸĽ�2���¼���ġ���end
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
                Positions(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+w*Leader_pos(j);%exp(1)���ǳ���e= 2.7183   
            end
        end
    end
    t=t+1;
    Convergence_curve(t)=Leader_score;
    [t Leader_score]
end



