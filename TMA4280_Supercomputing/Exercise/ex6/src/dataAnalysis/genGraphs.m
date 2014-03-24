function [] = genGraphs()
%GENGRAPHCS Summary of this function goes here
%   Detailed explanation goes here

data = dlmread('data.txt');

processes = data(:,1);
threads = data(:,2);
size = data(:,3);
time = data(:,4);
error = data(:,5);

[uniqueSizes, avgError] = averageErrorPerSize(size,error);
figure, loglog(uniqueSizes, avgError, uniqueSizes, 1./(uniqueSizes).^2);

[uniqueSizes, minTime] = minTimes(size, time);
figure, plot(uniqueSizes, minTime);

[uniqueSizes, optProc, optThrd] = optimalParallel(size, processes, threads, time);
figure, plot(uniqueSizes, optProc, uniqueSizes, optThrd);
end

function [uniqueSizes, avgErrors] = averageErrorPerSize(size, error)
[uniqueSizes,~,ind] = unique(size);
avgErrors = accumarray(ind,error)./accumarray(ind,1);
end

function [uniqueSizes,minTime] = minTimes(size, time)
[uniqueSizes,~,ind] = unique(size);
minTime = accumarray(ind, time, [] , @min);
end

function [uniqueSizes, optProc, optThrd] = optimalParallel(sizes, processes, threads, time)
[uniqueSizes,~,ind] = unique(sizes);
optProc = uniqueSizes;
optThrd = uniqueSizes;
minTime = accumarray(ind, time, [] , @min);
j = 1;
for i = 1:size(uniqueSizes)
    currentSize = uniqueSizes(i);
    while (j <= length(sizes)) && (currentSize == sizes(j))
        if time(j) == minTime(i)
            optProc(i) = processes(j);
            optThrd(i) = threads(j);
        end
        j = j+1;
    end
end

end

