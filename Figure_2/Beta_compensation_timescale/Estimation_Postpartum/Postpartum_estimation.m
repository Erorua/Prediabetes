function Beta_cell_timescale_estimation_from_postpartum_data

% ============================================================================================
% Description
% ============================================================================================

%%% Author: Aurore Woller

%%% Date: December 2022

%%% Uni: Weizmann institute of Science

%%% Description: Beta cell timescale estimation T_B (in days)

%%% From postpartum Clalit data

%%% Fits the function G(t)=p(3)*(1-p(1)*exp(-p(2)*t) to the data;

%%% with T_B=1/(3*p(2))

%%% Bootstrapping method

%%% input file: glucose and insulin Clalit postpartum data

%%% output file : "Postpartum_Fit.txt"

%%% output file contains parameters p(1), p(2) and p(3) from bootstrapping

%%% T_B=tau=1./(3*median(p(2))] (in days)









% ============================================================================================
% Load and process Clalit data 
% ============================================================================================

Gc=readtable('./GLUCOSE_BLOOD.csv');% table with pregancy and postpartum G 


Ic=readtable('./INSULIN_2w.csv');% table with pregancy and postpartum G 


Aa = table2array(Gc(:,1)); %times of measurements

Times=categorical(Aa(:,1));% times without bracket

Data=[];

time= -59:2:79; % time vector

time=time'*7; % time in days

% Maxtrix with [Time (days)  Glucose (mmol/L)  std (mmol/L)  people #  Insulin (pmol/L)  std (pmol/L)  people #]

Data=[Data;time 0.0555*table2array(Gc(:,3)) 0.0555*table2array(Gc(:,5)) table2array(Gc(:,2)) 6*table2array(Ic(:,3)) 6*table2array(Ic(:,5)) table2array(Ic(:,2))] % 0.0555 =conversion factor to mmol/L



% To plot data  from t > 2 weeks after pregancy (because SI=+- cste then)

kpost=find(Data(:,1) > 21); % because the two first data points are t= 7 days and t =21 days

tpost=Data(kpost:end,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start of the bootstrapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Malle=[];

for ibs=1:5000




% ============================================================================================
% Fit of mu
% ============================================================================================


 nbst=length(Data(kpost:end,1));


 ydata=[];

 
 % Draw random # in normal distribution to get postpartum time series of glucose 
 
for its=1:nbst
    
    moye=Data(kpost(its),2);
    
    standt=Data(kpost(its),3)./sqrt(Data(kpost(its),4));
    
   

    ran =standt.*randn + moye;

    
    ydata=[ydata;ran];






end


xdata=Data(kpost:end,1)-Data(kpost(1),1); % time 
 
ydata=ydata;





% initial  parameter guess


 p0 = [1 0.01 4.87];
  
  
% lower limit of parameters  

 
lb = [0.01 0.001 4.8];



% higher limit of parameters
 


ub = [10 0.1 4.9];




options = optimoptions('lsqcurvefit','Display','iter');
 
 
options=optimset('disp','iter','LargeScale','off','TolFun',.001,'MaxIter',1000,'MaxFunEvals',1000);
 
 
[p,resnorm,~,exitflag,output] = lsqcurvefit(@KineticFit,p0,xdata,ydata,lb,ub,options)


mufit=p(1);

G0fit=p(2);

Gend=p(3);


Malle=[Malle;mufit G0fit Gend]; 



end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of the bootstrapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Malle

%%% Output with parameters obtained from bootstrapping

save Postpartum_Fit.txt Malle -ASCII %%% output with [p(1) p(2) p(3)]


mum=mean(Malle(:,1));


musd=std(Malle(:,1));

mume=median(Malle(:,1));


mumad=mad(Malle(:,1));



muq25 = quantile(Malle(:,1),0.25);


muq75 = quantile(Malle(:,1),0.75);


mup5=prctile(Malle(:,1),5);


mup95=prctile(Malle(:,1),95);


G0m=mean(Malle(:,2));


G0sd=std(Malle(:,2));

G0me=median(Malle(:,2));


G0q25 = quantile(Malle(:,2),0.25);


G0q75 = quantile(Malle(:,2),0.75);


G0p5=prctile(Malle(:,2),5);


G0p95=prctile(Malle(:,2),95);


Gendm=mean(Malle(:,3));


Gendsd=std(Malle(:,3));

Gendme=median(Malle(:,3));

Gendq25 = quantile(Malle(:,3),0.25);


Gendq75 = quantile(Malle(:,3),0.75);


tau=[1./(3*Malle(:,2))];

taumed=median(tau) % median of beta  compensation timescale T_B

taumad=mad(tau) % median abs devia. of beta  compensation timescale T_B


% ============================================================================================
% Plot
% ============================================================================================
figure(1)


errorbar(Data(kpost:end,1),Data(kpost:end,2),Data(kpost:end,3)./sqrt(Data(kpost:end,4)),'Marker','o','LineWidth',2,'Color',[0 0 1])


hold on;

pmed=[mume G0me Gendme];


plot(xdata+Data(kpost(1),1),KineticFit(pmed,xdata),'LineWidth',2,'Color',[0.65 0.65 0.65])

hold on;   

p25=[muq25 G0q25 Gendq25];



p75=[muq75 G0q75 Gendq75];


hold on; 


hold on

x1=(xdata+Data(kpost(1),1))';  
x2=(xdata+Data(kpost(1),1))'; 
y1=(KineticFit(pmed,xdata))';                      %#create first curve
y2=(KineticFit(p25,xdata))';                          %#create second curve
X=[x1,fliplr(x2)];                %#create continuous x value array for plotting
Y=[y1,fliplr(y2)];              %#create y values for out and then back

patch(X,Y,'y'); 

hold on;
% 
x1=(xdata+Data(kpost(1),1))';  
x2=(xdata+Data(kpost(1),1))';  
y1=(KineticFit(p75,xdata))';                      %#create first curve
y2=(KineticFit(pmed,xdata))';                          %#create second curve
X=[x1,fliplr(x2)];                %#create continuous x value array for plotting
Y=[y1,fliplr(y2)];              %#create y values for out and then back

patch(X,Y,'y'); 

hold on;

 




xlabel('Time (days)','FontName','Arial')

ylabel('Glucose (mmol/L)','FontName','Arial')

hold on;

set(gca,'FontName','Arial','FontSize',20);

hold on;



pbaspect([1 1 1])

legend('Postpartum Data','Fit')





% ============================================================================================
% Kinetic
% ============================================================================================


function vv = KineticFit(p,xdata)



vv=p(3)*(1-p(1)*exp(-p(2)*xdata));





