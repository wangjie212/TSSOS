function [opt,data,status]=blockpop_uncons_higher(n,data,solver)
basis=data.basis;
sizes=data.sizes;
supp=data.supp;
coe=data.coe;
supp1=data.supp1;

%obtain the blocks
osupp=odd_supp(n,supp1);
losize=size(osupp);
lo=losize(1);
lb=length(basis);
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
blocksize=zeros(1,cl);
for i=1:cl
    blocksize(i)=length(blocks{i});
end
nsizes=[unique(blocksize);hist(blocksize,unique(blocksize))];
s1=size(sizes);
s2=size(nsizes);             
if s1(2)~=s2(2)||~isempty(setdiff(nsizes,sizes,'rows'))
   data.sizes=nsizes;
   disp('blocksizes=');
   disp(nsizes);
else
   opt=0;
   data=0;
   status=0;
   disp('No higher blocking hierarchy'); 
   return
end

[supp1,opt,status]=blockupop(n,supp,coe,basis,blocks,cl,blocksize,solver);
data.supp1=supp1;
end