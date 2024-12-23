

load('hmm_2_3_4stateresult.mat', 'Outpurparams4', 'NLLout4');


sortedemiss4=zeros(100,5);
for i=1:100
sortedemiss4(i,:)=[sort([Outpurparams4(i,13:16)]),NLLout4(i)];


end


sortrowsortedemiss4=sortrows(sortedemiss);

%%


load('hmm_2_3_4stateresult.mat', 'Outpurparams', 'NLLout');

%%
sortedemiss2=zeros(100,3);
for i=1:100
sortedemiss2(i,:)=[sort([Outpurparams(i,3:4)]),NLLout(i)];


end


sortrowsortedemiss2=sortrows(sortedemiss2);