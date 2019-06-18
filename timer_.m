function [record, total] = timer_ (n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   use 'timer_(n)' to record the time used for type-n operation
%   the same operation can be recorded many times
%

    persistent record_ total_

    if n == 0
    % reset timer

        tic
    elseif n == -1
    % reset record_ total_

        record_ = [];
        total_ = tic;
    elseif n == -2
    % retrive record total

        record = record_;
        total = total_;
    elseif n == -3
    % end of 'total_' timer

        total_ = toc(total_);
    elseif n == -9
    % plot current record
    
        figure;
        hold on;
        index = unique(record_(:, 1));
        avgArray = zeros(length(index), 1);
        sumArray = zeros(length(index), 1);
        maxTime = max(record_(:, 2));
        axis([-0.5 length(index)+0.5 0 maxTime*1.2]);
        for i = 1:length(index)
        
            tempArray = record_(record_(:, 1) == index(i), 2);
            avgArray(i) = mean(tempArray);
            sumArray(i) = sum(tempArray);
            scatter(ones(length(tempArray), 1)*i, tempArray, 'filled', 'd');
        end
        plot(1:length(index), avgArray');
        plot(1:length(index), sumArray'./max(sumArray).*maxTime);
        xlabel('operation index');
        ylabel('second');
        % (recorded time) / (total time) percentage
        sumTime = sum(record_(:, 2));
        ratio = sumTime / total_;
        plot([0 0], [0 maxTime], ':');
        scatter([0 0 0], [0 maxTime maxTime*ratio], 'filled', 'o');
        text(-0.32, maxTime, num2str(total_));
        text(0.05, maxTime*ratio, num2str(sumTime));
    else
  
        record_ = [record_; n, toc];
        tic
    end
end
