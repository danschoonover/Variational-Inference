%function [U,Cov,pi]=Varyashonal(X)
clear all;
load 'X.mat'
X=data;


K=3;
[N,D]=size(X);
P=zeros(N,K);
S=cell(K,1);
v=1+rand(K,1);
v0=D;
a0=.01*ones(K,1);

m=ones(K,D);    %initialize with rand*range(x,y)


mold=rand(size(m));
m0=5*ones(1,D); %cell



B=.1*ones(K,1);
a=(rand)*ones(K,1);
count=0;
B0=.01;
Euv=zeros(K,1);
LonA=zeros(K,1);
LonPi=zeros(K,1);

xhat=cell(K,1);
R=(1/K)*ones(N,K);
Rold=R;
W0=eye(D);

W=cell(K,1);
for k=1:K,
    W{k}=eye(D);
end

Nkold=0;

while(1)
    for n=1:N,
        for k=1:K                        %E step
            
            Euv(k)=(D/B(k))+v(k)*(X(n,:)-m(k,:))*W{k}*(X(n,:)-m(k,:))';  %problem here!
            
            sumofall=0;
            for i =1:D,
                sumofall=sumofall+psi((v(k)+1-i)/2);
            end
            LonA(k)=sumofall+D*log(2)+log(det(W{k}));
            
            aHat=sum(a);
            LonPi(k)=psi(a(k))-psi(aHat);
            
            P(n,k)=exp(LonPi(k)+.5*LonA(k)+-.5*Euv(k)+-D*log(2*pi)/2);
            
        end
        
    end
    
    for n=1:N,
        for k=1:K
            R(n,k)=P(n,k)/sum(P(n,:));      %E step continued
        end
    end
    
    for k=1:K,
        Nk(k)=sum(R(:,k));
    end
    for k=1:K,
        xhat{k}=(R(:,k)'*X)./Nk(k);
    end
    
    for k=1:K,
        eX=repmat(R(:,k),1,D).*(X-repmat(xhat{k},N,1));
        S{k}=eX'*(X-repmat(xhat{k},N,1))./Nk(k);
    end
    
    for k=1:K,                                             %M step
        a(k)=a0(k)+Nk(k);
        B(k)=B0+Nk(k);
        m(k,:)=(B0*m0+Nk(k)*xhat{k})./B(k);
        
        bigone=((xhat{k}-m0)'*(xhat{k}-m0))*(B0*Nk(k)/(B0+Nk(k)));               %last error is here!!
        W{k}=inv(inv(W0)+(Nk(k)*S{k})+bigone);
        v(k)=v0+Nk(k)+1;
                
    end
    
    if all(all(Nk==Nkold))% all(all(R==Rold)),%count>80%(R==Rold)%end condition
        break
    else
        clc
        count=count+1
        Nkold=Nk;
        
        for k=1:K,
            temp(k)=randg(a(k));
        end
        temp=temp./sum(temp)
        m
    end
    
end


for k=1:K,
    Cov{k}=wishrnd(W{k},v(k));
end

for k=1:K,
    Mu{k}=m(k,:);
end

for k=1:K,
    temp(k)=randg(a(k));
end
temp=temp./sum(temp)

obj = gmdistribution(Mu{1},S{1},temp(1));
ezcontour(@(x,y)pdf(obj,[x y]),[-2 2],[-2 2]);

obj = gmdistribution(Mu{3},S{3},temp(3));
ezcontour(@(x,y)pdf(obj,[x y]),[-2 2],[-2 2]);
