function [fsupp,gsupp,opt,status]=blockcpop(n,m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,solver)
    fsupp=[];
    for i=1:fcl
        for j=1:fblocksize(i)
            for r=j:fblocksize(i)
                bi=fbasis(fblocks{i}(j),:)+fbasis(fblocks{i}(r),:);
                fsupp=[fsupp;bi];
            end
        end
    end
    gsupp=cell(1,m);
    supp1=fsupp;
    for k=1:m
        gsupp{k}=[];
        for i=1:gcl(k)
            for j=1:gblocksize{k}(i)
                for r=j:gblocksize{k}(i)
                    for s=1:lt(k+1)
                        bi=ssupp{k+1}(s,:)+gbasis{k}(gblocks{k}{i}(j),:)+gbasis{k}(gblocks{k}{i}(r),:);
                        gsupp{k}=[gsupp{k};bi];
                    end
                end
            end
        end
        supp1=[supp1;gsupp{k}];
    end
    
supp1=unique(supp1,'rows');
lsupp1=length(supp1);
y=sdpvar(1,lsupp1);
Locb=bfind(supp1,lsupp1,ssupp{1}(1,:),n);
obj=coe{1}(1)*y(Locb);
for i=2:lt(1)
    Locb=bfind(supp1,lsupp1,ssupp{1}(i,:),n);
    obj=obj+coe{1}(i)*y(Locb);
end
fmom=cell(1,fcl);
for i=1:fcl  
    fmom{i}=sdpvar(fblocksize(i));
    for j=1:fblocksize(i)
        for r=j:fblocksize(i)      
            bi=fbasis(fblocks{i}(j),:)+fbasis(fblocks{i}(r),:);
            Locb=bfind(supp1,lsupp1,bi,n);
            fmom{i}(j,r)=y(Locb);
            fmom{i}(r,j)=y(Locb);
        end
    end
end
gmom=cell(1,m);
for k=1:m
    gmom{k}=cell(1,gcl(k));
    for i=1:gcl(k)  
        gmom{k}{i}=sdpvar(gblocksize{k}(i));
        for j=1:gblocksize{k}(i)
            for r=j:gblocksize{k}(i)
                Locb=zeros(1,lt(k+1));
                for s=1:lt(k+1)
                    bi=ssupp{k+1}(s,:)+gbasis{k}(gblocks{k}{i}(j),:)+gbasis{k}(gblocks{k}{i}(r),:);
                    Locb(s)=bfind(supp1,lsupp1,bi,n);
                end
                gmom{k}{i}(j,r)=coe{k+1}*y(Locb)';
                gmom{k}{i}(r,j)=gmom{k}{i}(j,r);
            end
        end
    end
end
momCons=[];
for i=1:fcl
    momCons=[momCons,fmom{i}>=0];      
end
for k=1:m
    for i=1:gcl(k)
        momCons=[momCons,gmom{k}{i}>=0];   
    end
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