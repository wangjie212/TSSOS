function [opt,data,status]=blockpop_cons_higher(n,m,data,solver)
fbasis=data.fbasis;
gbasis=data.gbasis;
fsizes=data.fsizes;
fsupp=data.fsupp;
gsupp=data.gsupp;
ssupp=data.ssupp;
coe=data.coe;
lt=data.lt;

%obtain the blocks
ofsupp=odd_supp(n,fsupp);
ofsupp=unique(ofsupp,'rows');
lofsize=size(ofsupp);
lfo=lofsize(1);
flb=length(fbasis);
fedges=[];
for i = 1:flb
    for j = i:flb
        bi=fbasis(i,:)+fbasis(j,:);
        if length(bi(~mod(bi,2)))==n||bfind(ofsupp,lfo,bi,n)~=0
           fedges=[fedges;[i j]];
        end
    end
end
fG=graph(fedges(:,1),fedges(:,2));
fblocks=conncomp(fG,'OutputForm','cell');
fcl=length(fblocks);
fblocksize=zeros(1,fcl);
for i=1:fcl
    fblocksize(i)=length(fblocks{i});
end
nfsizes=[unique(fblocksize);hist(fblocksize,unique(fblocksize))];
s1=size(fsizes);
s2=size(nfsizes);             
if s1(2)~=s2(2)||~isempty(setdiff(nfsizes,fsizes,'rows'))
   data.fsizes=nfsizes;
   disp('fblocksizes=');
   disp(nfsizes);
else
   opt=0;
   data=0;
   status=0;
   disp('No higher blocking hierarchy'); 
   return
end
disp('gblocksizes=');
glb=zeros(1,m);
gblocks=cell(1,m); 
gcl=zeros(1,m);
gblocksize=cell(1,m);
for k=1:m
    glb(k)=length(gbasis{k});
    ggsupp=[fsupp;gsupp{k}];
    ogsupp=odd_supp(n,ggsupp);
    ogsupp=unique(ogsupp,'rows');
    logsize=size(ogsupp);
    lgo=logsize(1);
    gedges=[];
       for i = 1:glb(k)
           for j = i:glb(k)
               r=1;
               while r<=lt(k+1)
                     bi=ssupp{k+1}(r,:)+gbasis{k}(i,:)+gbasis{k}(j,:);
                     if length(bi(~mod(bi,2)))==n||bfind(ogsupp,lgo,bi,n)~=0
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
data.fsupp=fsupp;
data.gsupp=gsupp;
end