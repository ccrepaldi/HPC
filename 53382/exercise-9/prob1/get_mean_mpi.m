Data = load('prob1.out');

Data=sort(Data);
i=0;
for n=0:17:374
    i=i+1;
    avg(i,2)= mean(Data(n+1:n+17,2));
    avg(i,3) = std(Data(n+1:n+17,2))/sqrt(17);
end

avg(:,1)=2:24;