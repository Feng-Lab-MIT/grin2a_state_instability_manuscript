%% this script is for selecting valid blocks

function valid=testBlock(block)

if length(block.RewardedTimes)<10
    valid='block too short (<10)';
else
    if length(block.RewardedTimes)/length(block.InitiationTimes)<0.5
        valid='too many aborted trials';
    else
        if length(block.LRpress)<99
            if length(unique(block.LRpress(end-5:end)))>1
                valid='noLRcommitment';
            else if unique(block.LRpress(end-5:end))==6
        valid=1;
                else
                    valid='need to check';
                end
            end
        else 
            valid='no LR commitment within 100trials';
        end
    end
end

end