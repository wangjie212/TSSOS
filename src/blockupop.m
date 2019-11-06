function [supp1,opt,status]=blockupop(n,supp,coe,basis,blocks,cl,blocksize,solver)
lsize=size(supp);
lsupp=lsize(1);
supp1=[];
for i=1:cl
    for j=1:blocksize(i)
        for r=j:blocksize(i)
            supp1=[supp1;basis(blocks{i}(j),:)+basis(blocks{i}(r),:)];
        end
    end
end
supp1=unique(supp1,'rows');
lsupp1=length(supp1);
y=sdpvar(1,lsupp1);
Locb=bfind(supp1,lsupp1,supp(1,:),n);
obj=coe(1)*y(Locb);
for i=2:lsupp
    Locb=bfind(supp1,lsupp1,supp(i,:),n);
    obj=obj+coe(i)*y(Locb);
end
mom=cell(1,cl);
for i=1:cl  
    mom{i}=sdpvar(blocksize(i));
    for j=1:blocksize(i)
        for r=j:blocksize(i)
            bi=basis(blocks{i}(j),:)+basis(blocks{i}(r),:);
            Locb=bfind(supp1,lsupp1,bi,n);
            mom{i}(j,r)=y(Locb);
            mom{i}(r,j)=y(Locb);
        end
    end
end
momCons=[];
for i=1:cl
    momCons=[momCons,mom{i}>=0];      
end
ops=sdpsettings('solver',solver);
sol=optimize([momCons,y(1)==1],obj,ops);
status=sol.problem;
opt=value(obj);
if status==0
    disp('Optimal=');
    disp(opt);
else
    disp('Failed');
end
end