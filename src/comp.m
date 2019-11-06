function sign=comp(a,b,n)
    i=1;
    while i<=n
          if a(i)<b(i)
             sign=-1;
             return
          elseif a(i)>b(i)
             sign=1;
             return
          else
             i=i+1;
          end
    end
    if i==n+1
       sign=0;
       return
    end
end