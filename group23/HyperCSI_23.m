%======================================================================
%  Input
%  X is M-by-L data matrix, where M is the number of spectral bands and L is the number of pixels.
%  N is the number of endmembers.
%----------------------------------------------------------------------
%  Output
%  A_est is M-by-N mixing matrix whose columns are estimated endmember signatures.
%  S_est is N-by-L source matrix whose rows are estimated abundance maps.
%  time is the computation time (in secs).
%========================================================================

function [A_est, S_est, time] = HyperCSI_23(X,N)
t0 = clock;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% you need to implement HyperCSI here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dimension
M=224;
L=10000;
%step1
%calculate mean d
d=zeros(M,1);%ミM*10痻皚
for i=1:M
    d(i,1)=sum(X(i,:))/L;%盢X[n]羆癬ㄓキА
end

%calculate mean-remove data matrix
U=zeros(M,L);%ミM*L0痻皚
for i=1:M
    for j=1:L
        U(i,j)=X(i,j)-d(i,1);
    end
end
U_t=transpose(U); %L*M
U_x=U*U_t; %M*M

%U_x*V=V*D
[V D]=eig(U_x);

%calculate C
C=V(:,M-N+2:M);%程ㄢcolumn局Τ程癟璝ㄤcolumn穦癟ぃìτ礚猭临
C_t=transpose(C);

%calculate x_tilde
x_tilde=C_t*U; %(N-1)*L

%end of step1

%step2
%find lifted DR data
x_lifted=ones(N,L);
for j=1:(N-1) %add one row with all one
    for k=1:L
        x_lifted(j,k)=x_tilde(j,k);
    end
end

%SPA(maximization and projection)
alpha_bar=zeros(N);
x_2norm=zeros(1,L);
x_temp=x_lifted;%use for iteration
for i=1:N  
    %maximization
    %find 2-norm of x_temp of each column
    for j=1:L
        x_2norm(1,j)=norm(x_temp(:,j));
    end
    [value,index]=max(x_2norm);%return the maximum of x_2norm with value and column index
    alpha_bar(:,i)=x_lifted(:,index);
    alpha_bar_t=transpose(alpha_bar);% transpose of a
    %end of maximization
    %projection
    P=eye(N)-((alpha_bar(:,i)*alpha_bar_t(i,:))/(value^2));%projection matrix
    x_temp=P*x_temp;%iteration
end

%remove N-th row(all ones) of alpha_bar and obtain alpha_tilde
alpha_tilde=alpha_bar(1:N-1,:);

%end of SPA
%end of step2

%step3
%find b_tilde
b_tilde=find_b_tilde(alpha_tilde,N);
%end of finding b_tilde

%calculate the radius of Euclidean norm ball
r_temp=[];%r_temp saves the distance from two distinct points of alpha_tilde
for i=1:(N-1)
    for j=1:(N-i)
        r_temp=[r_temp ;norm((alpha_tilde(:,i)-alpha_tilde(:,i+j)))];
    end
end
r=min(r_temp)/2;%r is the radius of Euclidean norm ball
%calculate the domain {x_tilde&B}
%–Rcolumn计ぃ妓ぃ俱Θ蝴痻皚
R1=[];%domain of Euclidean norm ball with center alpha_tilde(:,1)
R2=[];%domain of Euclidean norm ball with center alpha_tilde(:,2)
R3=[];%domain of Euclidean norm ball with center alpha_tilde(:,3)
for j=1:L
    %find the x_tilde in Euclidean norm ball with center alpha_tilde(:,1)
    if(norm(x_tilde(:,j)-alpha_tilde(:,1))<r)
        R1=[R1  x_tilde(:,j)];
    end
    %find the x_tilde in Euclidean norm ball with center alpha_tilde(:,2)
    if(norm(x_tilde(:,j)-alpha_tilde(:,2))<r)
        R2=[R2  x_tilde(:,j)];
    end
    %find the x_tilde in Euclidean norm ball with center alpha_tilde(:,3)
    if(norm(x_tilde(:,j)-alpha_tilde(:,3))<r)
        R3=[R3  x_tilde(:,j)];
    end
end

%end of step3


%step4
p_active1=[]; %active pixel for i=1
p_active2=[]; %active pixel for i=2
p_active3=[]; %active pixel for i=3
for i=1:N
    if(i==1) %find the active pixel in R2 and R3
        [~,index]=max(transpose(b_tilde(:,i))*R2);
        p_active1= [p_active1 R2(:,index)];
        [~,index]=max(transpose(b_tilde(:,i))*R3);
        p_active1=[p_active1 R3(:,index)];  
    end
    if(i==2) %find the active pixel in R1 and R3
        [~,index]=max(transpose(b_tilde(:,i))*R1);
        p_active2= [p_active2 R1(:,index)];
        [~,index]=max(transpose(b_tilde(:,i))*R3);
        p_active2=[p_active2 R3(:,index)];  
    end
    if(i==3) %find the active pixel in R1 and R2
        [~,index]=max(transpose(b_tilde(:,i))*R1);
        p_active3=[p_active3 R1(:,index)]; 
        [~,index]=max(transpose(b_tilde(:,i))*R2);
        p_active3=[p_active3 R2(:,index)];  
    end
end

%find b_hat
b_hat(:,1)=find_b_hat([p_active1 alpha_tilde(:,1)],N);
b_hat(:,2)=find_b_hat([p_active2 alpha_tilde(:,2)],N);
b_hat(:,3)=find_b_hat([p_active3 alpha_tilde(:,3)],N);

%end of finding b_hat

%find h_hat
h_hat=[];
for i=1:N
    h_hat=[h_hat;max(transpose(b_hat(:,i))*x_tilde)];
end
%end of finding h_hat
%end of step4

%step5
%since the data set doesn't have noise,we set c=1 and =0.9
c=1;
%end of step5

%step6
%calulate alpha_hat
alpha_hat=[];
b_hat_temp=[];
h_hat_temp=[];
for i=1:N
   for j=1:N
       if(j~=i)%remove i-th data
           b_hat_temp=[b_hat_temp b_hat(:,j)];
           h_hat_temp=[h_hat_temp;h_hat(j,:)];
       end   
   end
   alpha_hat=[alpha_hat inv(transpose(b_hat_temp))*(h_hat_temp/c)];
   b_hat_temp=[];%reset for next iteration
   h_hat_temp=[];
end
%end of calulating alpha_hat

%calculate a_hat
a_hat=C*alpha_hat+d;
%end of calculating a_hat
%end of step6

%step7
%calculate the abundace vector
s=zeros(N,L);
for i=1:N
    den=h_hat(i,:)-transpose(b_hat(:,i))*alpha_hat(:,i);%denominator
    for j=1:L
        num=h_hat(i,:)-transpose(b_hat(:,i))*x_tilde(:,j);% numerator
        s(i,j)=num/den;
        if(s(i,j)<0)%s must be non-negative
            s(i,j)=0;
        end
    end
end
%end of step7

%output
A_est=a_hat;
S_est=s;
%end of output

time = etime(clock,t0);
function [b] = find_b_tilde(a,N)
for i=1:N
    j=mod((i+1),N)+1;%set j!=i and avoid "j" being 0
    a_remove=[a(:,i) a(:,j)];%the data we want to remove
    Q=setdiff(a,a_remove,'stable');%remove i,j column and not being sorted
    P=Q-a(:,j)*transpose(ones(N-2,1));
    P_t=transpose(P);
    b(:,i)=(eye(N-1)-P*inv(P_t*P)*P_t)*(a(:,j)-a(:,i));
end
return
function [b] = find_b_hat(a,N) 
    %a=[p_active a_tilde(:,i)],we remove a_tilde(:,i) and second column of p_active
    a_remove=[a(:,2) a(:,3)];%the data we want to remove
    Q=setdiff(a,a_remove,'stable');%remove 2,3 column and not being sorted
    P=Q-a(:,2)*transpose(ones(N-2));
    P_t=transpose(P);
    b=(eye(N-1)-P*inv(P_t*P)*P_t)*(a(:,2)-a(:,3));
return
