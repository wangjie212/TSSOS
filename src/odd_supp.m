function osupp=odd_supp(n,supp)
i=1;
lsize=size(supp);
lo=lsize(1);
indexb=1:lo;
while lo>=i
      bi=supp(indexb(i),:);
      if length(bi(~mod(bi,2)))==n
         indexb(i)=[];
         lo=lo-1;
      else
          i=i+1;
      end
end
osupp=supp(indexb,:);
end