M = csvread('data.csv');
for i=1:24
    m(i)=mean(M(:,i));
    sigma(i)=std(M(:,i));
    minim(i)=min(M(:,i));
    maxim(i)=max(M(:,i));
end

normM=zeros(43,24);
for i=1:24
    for k=1:43
        normM(k,i)=(M(k,i)-m(i))/sigma(i);
    end
end

for i=1:24
 %   plot(normM(:,i))
 %   hold on
end

for i=1:24
    f(i)=kstest(normM(:,i), 'Alpha', 0.005);
end

% Experiment 1. 
% Empirical probability nu(U) 
% U_1
gam_vector=[0:0.5:4];
itermax=1000;
sum1=0; 
for iter=1:itermax
    for i=1:24
        x(i)=normrnd(m(i),sigma(i)); % generate sample point
    end
    test1=1; % by default point x belongs to U_1
    for i=1:24
        if x(i)>maxim(i) || x(i)<minim(i)
            test1=0;
        end
    end
    sum1=sum1+test1;
end
nu1=sum1/itermax;
plot(gam_vector, nu1*ones(1,9))


%U_2
j=1; % counter
for gamma=0:0.5:4
    sum=0;
    for iter=1:itermax
        for i=1:24
            x(i)=normrnd(m(i),sigma(i)); % generate sample point
        end
     test=1;
     for i=1:24
        if x(i)>gamma*sigma(i)+m(i) || x(i)<-gamma*sigma(i)+m(i)
            test=0;
        end
     end
     sum=sum+test;
    end
    nu2(j)=sum/itermax;
    j=j+1;
end
 plot(gam_vector,nu2)

%U_3
j=1; % counter
for gamma=0:0.5:4
    sum=0;
    for iter=1:itermax
     t=0; % summation variable
        for i=1:24
            x(i)=normrnd(m(i),sigma(i)); % generate sample point
            t=t+abs((x(i)-m(i))/sigma(i));
        end
     test=1;
     for i=1:24
        if t>gamma*sqrt(24)
            test=0;
        end
     end
     sum=sum+test;
    end
    nu3(j)=sum/itermax;
    j=j+1;
end
 plot(gam_vector,nu3)
 
 
 plot(gam_vector, nu1*ones(1,9),'-.b', 'Linewidth',4)
 hold on 
 plot(gam_vector,nu2,'--','Linewidth',4)
 plot(gam_vector,nu3,'Linewidth',4)
 grid on
tx=xlabel({' ', 'Robustness level  \Gamma'}, 'FontSize', 40)
ty=ylabel({'Empirical probability \nu(U)',' '}, 'FontSize', 40)
set(gca,'FontSize',35)
h=legend({'Set U_1', 'Set U_2(\Gamma)','Set U_3(\Gamma)'},'Location', 'west','FontSize',40)

 



