function [opt,data,status]=blockpop_uncons_first(poly,n,pd,solver,varargin)
d=pd/2;
x = sym('x',[1 n]);
[coe, terms] = coeffs(poly,x);
coe=double(coe);
lt=length(terms);
supp=zeros(lt,n);
for i=1:lt
    for j=1:n
        supp(i,j)=feval(symengine,'degree',terms(i),x(j));
    end
end
lsize=size(supp);
lsupp=lsize(1);
basis=get_basis(n,d);
lb=nchoosek(n+d,d);

% compute a monomial basis by the Newton polytope method
if nargin==4||varargin{2}==1
if sum(supp(lsupp,:))~=0
   supp=[supp;zeros(1,n)];
   coe=[coe 0];
   lsupp=lsupp+1;
end
A0=[-1/2*supp ones(lsupp,1)];
t=1;
indexb=1:lb;
temp=unique(supp,'rows');
while t<=lb
      i=indexb(t);
      if bfind(supp,lsupp,2*basis(i,:),n)~=0
         t=t+1;
      else  
         A=[A0;[basis(i,:),-1]];
         [vx,~,exitflag]=linprog([basis(i,:),-1],A,zeros(1,lsupp+1));
         if exitflag==1
            t=t+1;
         else
            lb=lb-1;
            indexb(t)=[];
            r=t;
            while lb>=r
                  j=indexb(r);
                  if [basis(j,:),-1]*vx<=-0.001&&bfind(supp,lsupp,2*basis(i,:),n)==0
                     lb=lb-1;
                     indexb(r)=[];
                  else
                      r=r+1;
                  end
            end
         end
      end
end
basis=basis(indexb,:);
end

osupp=odd_supp(n,supp);
osupp=unique(osupp,'rows');
losize=size(osupp);
lo=losize(1);
edges=[];    
for i = 1:lb
    for j = i:lb
        bi=basis(i,:)+basis(j,:);
        if length(bi(~mod(bi,2)))==n||bfind(osupp,lo,bi,n)~=0
           edges=[edges;[i j]];
        end
    end
end
G=graph(edges(:,1),edges(:,2));
blocks=conncomp(G,'OutputForm','cell');
cl=length(blocks);
if cl==1
   opt=0;
   data=0;
   status=0;
   disp('Unblockable');
   return
end
blocksize=zeros(1,cl);
for i=1:cl
    blocksize(i)=length(blocks{i});
end
sizes=[unique(blocksize);hist(blocksize,unique(blocksize))];
disp('blocksizes=');
disp(sizes);
[supp1,opt,status]=blockupop(n,supp,coe,basis,blocks,cl,blocksize,solver);
data.supp=supp;
data.basis=basis;
data.coe=coe;
data.supp1=supp1;
data.sizes=sizes;
end
