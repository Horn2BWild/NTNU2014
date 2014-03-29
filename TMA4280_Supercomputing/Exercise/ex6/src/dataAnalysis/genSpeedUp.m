function [  ] = genSpeedUp( )
%GENSPEEDUP Summary of this function goes here
%   Detailed explanation goes here

baseline = dlmread('proc1Poisson.txt');
bsize = baseline(:,3);
btime = baseline(:,4);

vary = dlmread('varyPoisson.txt');
vsize = vary(:,3);
vtime = vary(:,4);
vproc = vary(:,1);

[uniqueSizes,~,sInd] = unique(vsize);
[uniqueProc,~,pInd] = unique(vproc);
datalen = length(uniqueProc);
data = zeros(datalen, length(uniqueSizes));

for i = 1:uniqueProc
    loc = 1;
    for j = 1:length(vtime)
        if (pInd(j) == i)
            data(i,loc) = vtime(j);
            loc = loc + 1;
        end
    end
end

speedup = data';

figure, hold, plot ([1;uniqueProc],[1;uniqueProc]);
for i = 1:length(uniqueSizes)
    plot([1;uniqueProc], [1,btime(i)./speedup(i,:)]);
end


figure, hold, plot(uniqueSizes, btime);
for i = 1:uniqueProc
    plot(uniqueSizes, data(i,:));
end

end

