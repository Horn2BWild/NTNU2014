function [  ] = genSpeedUp( )
%GENSPEEDUP Summary of this function goes here
%   Detailed explanation goes here

baseline = dlmread('proc1Poisson.txt');
bsize = baseline(:,3);
btime = baseline(:,4);

vary = dlmread('varying2.txt');
vsize = vary(:,3);
vtime = vary(:,4);
vproc = vary(:,1);

[uniqueSizes,~,sInd] = unique(vsize);
[uniqueProc,~,pInd] = unique(vproc);
datalen = length(uniqueProc);
data = zeros(datalen, length(uniqueSizes));
datac = zeros(datalen, length(uniqueSizes));

for i = 1:length(uniqueProc)
    loc = 1;
    for k = 1:length(vsize)
        for j = 1:length(vtime)
            if (pInd(j) == i && sInd(j) == loc)
                datac(i,loc) = datac(i,loc) + 1;
                data(i,loc) = (data(i,loc) * (datac(i,loc)-1) + vtime(j))/datac(i,loc);
            end
        end
        loc = loc + 1;
    end
end

speedup = data';

figure, hold, plot ([1;uniqueProc],[1;uniqueProc]);
% for i = 1:length(uniqueSizes)
i= 1;
plot([1;uniqueProc], [1,btime(1)./speedup(i,:)]);
i= 2;
plot([1;uniqueProc], [1,btime(5)./speedup(i,:)]);
i= 3;
plot([1;uniqueProc], [1,btime(9)./speedup(i,:)]);
i= 4;
plot([1;uniqueProc], [1,btime(13)./speedup(i,:)]);
% end


% figure, hold, plot(uniqueSizes, btime);
% for i = 1:uniqueProc
%     plot(uniqueSizes, data(i,:));
% end

end

