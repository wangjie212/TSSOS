function [opt,data,status]=blockpop_cons_first(f,g,n,m,d,dg,solver)
x = sym('x',[1 n]);
pop=[f,g];
coe=cell(1,m+1);
terms=cell(1,m+1);
ssupp=cell(1,m+1);
supp=[];
lt=zeros(1,m+1);
for k=1:m+1
    [coe{k}, terms{k}] = coeffs(pop(k),x);
    lt(k)=length(terms{k});
    ssupp{k}=zeros(lt(k),n);
    coe{k}=double(coe{k});
    for i=1:lt(k)
        for j=1:n
            ssupp{k}(i,j)=feval(symengine,'degree',terms{k}(i),x(j));
        end
    end
    supp=[supp;ssupp{k}];
end
supp=unique(supp,'rows');

fbasis=get_basis(n,d);
flb=nchoosek(n+d,d);
glb=zeros(1,m);
gbasis=cell(1,m);
for k=1:m
    glb(k)=nchoosek(n+d-ceil(dg(k)/2),d-ceil(dg(k)/2));
    gbasis{k}=get_basis(n,d-ceil(dg(k)/2));
end

%obtain the blocks
osupp=odd_supp(n,supp);
losize=size(osupp);
lo=losize(1);
fedges=[];
for i = 1:flb
    for j = i:flb
        bi=fbasis(i,:)+fbasis(j,:);
        if length(bi(~mod(bi,2)))==n||bfind(osupp,lo,bi,n)~=0
           fedges=[fedges;[i j]];
        end
    end
end
fG=graph(fedges(:,1),fedges(:,2));
fblocks=conncomp(fG,'OutputForm','cell');
fcl=length(fblocks);
if fcl==1
   opt=0;
   data=0;
   status=0;
   disp('Unblockable');
   return
end
fblocksize=zeros(1,fcl);
for i=1:fcl
    fblocksize(i)=length(fblocks{i});
end
fsizes=[unique(fblocksize);hist(fblocksize,unique(fblocksize))];
disp('fblocksizes=');
disp(fsizes);
disp('gblocksizes=');
gblocks=cell(1,m); 
gcl=zeros(1,m);
gblocksize=cell(1,m);
for k=1:m
    gedges=[];  
    for i = 1:glb(k)
        for j = i:glb(k)
            r=1;
            while r<=lt(k+1)
                  bi=ssupp{k+1}(r,:)+gbasis{k}(i,:)+gbasis{k}(j,:);
                  if length(bi(~mod(bi,2)))==n||bfind(osupp,lo,bi,n)~=0
                     break
                  else
                     r=r+1;
                  end
            end
            if r<=lt(k+1)
               gedges=[gedges;[i j]];
            end
        end
    end
    gG=graph(gedges(:,1),gedges(:,2));
    gblocks{k}=conncomp(gG,'OutputForm','cell');
    gcl(k)=length(gblocks{k});
    gblocksize{k}=zeros(1,gcl(k));
    for i=1:gcl(k)
        gblocksize{k}(i)=length(gblocks{k}{i});
    end
    gsizes=[unique(gblocksize{k});hist(gblocksize{k},unique(gblocksize{k}))];
    disp(gsizes);
end
[fsupp,gsupp,opt,status]=blockcpop(n,m,ssupp,coe,lt,fbasis,gbasis,fblocks,fcl,fblocksize,gblocks,gcl,gblocksize,solver);
data.fbasis=fbasis;
data.gbasis=gbasis;
data.lt=lt;
data.coe=coe;
data.ssupp=ssupp;
data.fsizes=fsizes;
data.fsupp=fsupp;
data.gsupp=gsupp;
end