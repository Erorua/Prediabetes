function Si_rate_change_male


% ============================================================================================
% Description
% ============================================================================================


%%% Author: Aurore Woller

%%% Date: December 2022

%%% Uni: Weizmann institute of Science

%%% Description: construct histogram of rate of change of insulin resistance 

%%% For male cohort

%%% Bootstrapping to calculate variability of the mode

%%% input: Si rate of change file : "Male_cohort_Si_rate_change.txt"

%%% output: timescale of rate of change of insulin resistance in days



% ============================================================================================
% Histogram
% ============================================================================================


%%%% Si rate of change file

m=load('./Male_cohort_Si_rate_change.txt');


mmmm=[];


% Remove zero values

for i=1:length(m(:,1))
    
    if m(i)==0
        
        
    else
        
     mmmm=[mmmm;m(i,1)];   
     
    end
    
    
end








mmmm

m=log10(mmmm);



%%%%% Kernel density estimation
% 
 [f,xj] = ksdensity(m); 


f'; %proba/occurence

xj';% rate of change of SI

[mx,ix]=max(f') % max of proba

trapz(xj, f) % normalization check

mode_ksdensity=10^(xj(ix)) % mode

tau=1./mode_ksdensity % caracteristic timescale





% %%% Histogram for Delta SI (rate of change of insulin sensitivity)
% 
%  figure(5)
% % 
% plot(xj,f,'LineWidth',4,'Color',[0 0 1]);
% 
% hold on;
% 
% 
% xlabel('Insulin sensitivity rate of change, |\alpha| (day-1, log10)')
% 
% ylabel('Probability')





 
%%% Histogram for Delta T_IR (rate of change of insulin resistance)
 
figure(4)


tauhist=log10(1./(mmmm));



nbins=20;

[n]=histogram(tauhist,nbins,'FaceColor',[0 0 1]);

xlabel('Insulin sensitivity rate of change, |\alpha| (day-1, log10)')

ylabel('Probability')



set(gca,'FontName','Arial','FontSize',20);


n


values=n.Values; % proba occurence

edges=(n.BinEdges); % bin edges

widtth=n. BinWidth; % bin width

bar(edges(1:end-1)+widtth*0.5,values/sum(values.*widtth))

hold on;

% 
xlabel('Insulin resistance timescale, T_{IR} (day, log10)')

ylabel('Probability')



set(gca,'FontName','Arial','FontSize',20);





hold on;
















%%%%% Bootstrapping of modes to find the variability of mode


taub=[];

modeb=[];

for i=1:5000
    
    
    
n=length(mmmm(:,1));

mmmmb=mmmm(randi(n,n,1),:);

mb=log10(mmmmb);

% % Kernel density estimation

[f,xj] = ksdensity(mb); 


f'; %proba/occurence

xj';% rate of change of SI

[mx,ix]=max(f'); % max of proba

trapz(xj, f); % normalization check

mode_ksdensity=10^(xj(ix)); % mode

tau=1./mode_ksdensity; % caracteristic timescale


taub=[taub;tau];

modeb=[modeb;mode_ksdensity];
    
    
    
end



taub;

taubm=mean(taub)

taubsd=std(taub)

 Store=[taubm taubsd]; % mean and standard deviation of the modes obtained by bootstrapping (expressed in days)



% save Delta_Si_bootstraping Store -ASCII 




