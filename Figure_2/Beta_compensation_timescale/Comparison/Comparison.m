function Timescales_comparison


% ============================================================================================
% Description
% ============================================================================================

%%% Author: Aurore Woller

%%% Date: December 2022

%%% Uni: Weizmann institute of Science

%%% Description: T_IR and T_B timescales comparison (in days)

%%% with T_IR : insulin resistance timescale

%%% and with T_B: beta cell compensation timescale

%%% Seasonal input file: "" from seasonal T_B estimation

%%% Postpartum input file: "Postpartum_Fit.txt" from postpartum T_B estima.

%%% Insulin sens. input file:"All_cohort_Si_rate_change.txt" from SI estim






% ============================================================================================
% T_B from seasonal Clalit data
% ============================================================================================


% Mean between Female and Male delay




phif=0.5*(44.7+49.8); % mean delay in days (see Tendler et al. PNAS 2021))

phifl=0.5*(44.7+49.8)-0.54; % (see Tendler et al. PNAS 2021)

phifh=0.5*(44.7+49.8)+0.54; % (see Tendler et al. PNAS 2021)
















figure(2)

hold on;

bar(1,phif,'FaceColor',[0 0 1]);

hold on;

%errorbar(1,taumeanf,taumeanf-tauvarlf,tauvarhf-taumeanf,'LineWidth',2,'Color',[0 0 0])




errorbar(1,phif,phif-phifl,phifh-phif,'LineWidth',2,'Color',[0 0 0])


hold on;





hold on;




ylabel('Timescale T (days)','FontName','Arial')

hold on;

set(gca,'FontName','Arial','FontSize',20);

hold on;


% ============================================================================================
% T_B from postpartum Clalit data
% ============================================================================================


%%% file with matrix [p(1) p(2) p(3) from postpartum estimation


Malle=load('./Postpartum_Fit.txt'); % T_B=tau=1./(3*median(p(2))] (in days)



G0me=median(Malle(:,2));


G0q25 = quantile(Malle(:,2),0.25);


G0q75 = quantile(Malle(:,2),0.75);



taume=1./(G0me);

tau25=1./(G0q25);

tau75=1./(G0q75);

taume-tau25;

tau75-taume;

bar(4,taume,'FaceColor',[0.125490196078431 0.47843137254902 0.250980392156863]);

hold on;


errorbar(4,taume,taume-tau25,tau75-taume,'MarkerSize',10,'LineWidth',2,'Color',[0 0 0])

hold on;



% % ============================================================================================
% % Insulin resistance timescale T_IR from Eran Segal entire cohort
% % ============================================================================================

%%% file with insulin sensitivity rate of change from Eran Segal cohort


m=load('./All_cohort_Si_rate_change.txt');


mmmm=[];

for i=1:length(m(:,1))
    
    if m(i)==0
        
        
    else
        
     mmmm=[mmmm;m(i,1)];   
     
    end
    
    
end








mmmm;

m=log10(mmmm);





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

taubm=mean(taub);

taubsd=std(taub);

 Store=[taubm taubsd];

bar(7,taubm,'FaceColor',[1 0 0]);

hold on;

   
errorbar(7,taubm,taubsd,taubsd,'LineWidth',2,'MarkerSize',10,'Color',[0 0 0])


set(gca,'FontName','Arial','FontSize',20,'YMinorTick','on','YScale','linear');


 set(gca,'FontName','Arial','FontSize',20,'XTick',[1 4 7],'XTickLabel', {'T_{comp}_{ season}','T_{comp}_{ postpartum}','T_{IR}'});

 
 
 xlim([0 8])
 
  
 ylim([0 450])
 
 pbaspect([1 1 1])





