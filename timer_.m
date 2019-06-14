function [record] = timer_ (n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   use 'timer_(n)' to record the time used for type-n operation
%   the same operation can be recorded many times
%

    persistent record_

    if n == 0
    % reset timer
        
        tic
    elseif n == -1
    % reset record_
        
        record_ = [];
    elseif n == -2
    % retrive record
    
        record = record_;
    elseif n == -9
    % plot current record
    
        figure;
        hold on;
        index = unique(record_(:, 1));
        avgArray = zeros(length(index), 1);
        axis([0.5 length(index)+0.5 0 max(record_(:, 2)*1.2)]);
        for i = 1:length(index)
        
            tempArray = record_(record_(:, 1) == index(i), 2);
            avgArray(i) = mean(tempArray);
            scatter(ones(length(tempArray), 1)*i, tempArray, 'filled', 'd');
        end
        plot(1:length(index), avgArray');
        xlabel('operation index');
        ylabel('second');
    else
        
        record_ = [record_; n, toc];
        tic
    end
end
