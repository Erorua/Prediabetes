function  Si_Rate_of_change_all_cohort



% ============================================================================================
% Description
% ============================================================================================



%%% Author: Aurore Woller & Yuval Tamir

%%% Date: December 2022

%%% Uni: Weizmann institute of Science

%%% Description: compute rate of change of insulin sensitivity 

%%% Choose if compute for upper/lower quartile of bmi/age of patients:
%%% quartile = 1 for lower quartile, quartile = 4 for upper.
%%% age = 1 for quartile on age values, age = 0 for bmi values.

%%% Compute as abs(Delta(log(Si)/Delta(t))

%%% Data from Segal Cohort

%%% output: file "All_cohort_Si_rate_change.txt"







% ============================================================================================
% Data
% ============================================================================================


%%% Eran Segal human data


Datap=readtable('./Prediabetic_Infos_Header.txt');


Datad=readtable('./Diabetic_Infos_Header.txt');


Data=[table2array(Datap);table2array(Datad)]; % matrix [Reg Date G I gender Age BMI]
%%
%User input
age = 0;%0 for quratile of bmi values and 1 for age values
quartile = 1;% 1 for lower quartile, 4 for upper.
%%
%Get quartiles out of data



temp = diff(Data(:,1));%Find where next patient appears in the list
first_patient_indx=[1;find(temp)+1];
initial_BMIs=Data(first_patient_indx,7);
initial_ages=Data(first_patient_indx,6);

qrtl1_bmi_thresh=prctile(initial_BMIs,25);
qrtl4_bmi_thresh=prctile(initial_BMIs,75);
qrtl1_bmi_full=initial_BMIs<qrtl1_bmi_thresh;
qrtl4_bmi_full=initial_BMIs>qrtl4_bmi_thresh;

qrtl1_age_thresh=prctile(initial_ages,25);
qrtl4_age_thresh=prctile(initial_ages,75);
qrtl1_age_full=initial_ages<qrtl1_age_thresh;
qrtl4_age_full=initial_ages>qrtl4_age_thresh;



[x,a,c]=unique(Data(:,1));

temp=c(first_patient_indx(qrtl1_bmi_full));
[indx_qrtl1_bmi, ii]=ismember(c,temp);
temp=c(first_patient_indx(qrtl4_bmi_full));
[indx_qrtl4_bmi, ii]=ismember(c,temp);
temp=c(first_patient_indx(qrtl1_age_full));
[indx_qrtl1_age, ii]=ismember(c,temp);
temp=c(first_patient_indx(qrtl4_age_full));
[indx_qrtl4_age, ii]=ismember(c,temp);

if age == 0 && quartile ==1
    Data=Data(indx_qrtl1_bmi,:);
elseif age == 0 && quartile ==4
    Data=Data(indx_qrtl4_bmi,:);
elseif age == 1 && quartile ==1
   Data=Data(indx_qrtl1_age,:);
else      
    Data=Data(indx_qrtl4_age,:);
end
%%
%%% Insulin resistance

 Homa_ir=Data(:,4).*Data(:,3)/22.5;
 
 
 Sim=1./Homa_ir;
 
 Storem=[Sim Data(:,3) Data(:,1) Data(:,2) Data(:,4)]% matrix with [Si G Reg Date I]
 
 
 
 
% ============================================================================================
% Calculations
% ============================================================================================
 
 
 
%%%% Compute rate of change of Insulin sensitibity
 
 Hero=[];

 Same=[];
 
 Stock=[];
 
 Fridge=[];
 
 Hero=[];
 
 rrr=1;
 z=1;
 
 zj=1;
 
 format longG

     
Same=[Same;Storem(:,1) Storem(:,2) Storem(:,4) Storem(:,3)];
     
 aaa=1;

 for i=1:length(Storem(:,1))
     
     
       i;
       
    if aaa ==1
       
     Same=[Storem(i,1) Storem(i,2) Storem(i,4) Storem(i,3)];
     
     
    else
        
        aaa;
        
    Same=Same;
        
        
    end
    
    
if i == length(Storem(:,1))
    
    
    Storem=[Storem;0 0 0 0 0 ];
    
    
    
end
    
    
    
    
    
    
if Storem(i+1,3) == Storem(i,3)
         
         aaa=aaa+1;
         
         
         
         
         Same=[Same;Storem(i+1,1) Storem(i+1,2) Storem(i+1,4) Storem(i+1,3)];
         
         
         
         
     else 
         
         
         aaa=1;
         
       
         
         Same;
         
         Sami=Same;
         
         ll=length(Same(:,1));
         
         Samee=[];
         
         if ll >1
             
         bbb=1;
             
         for jj=1:ll-1
             
        jj;
             
             
         xxx=log(Same(jj+1,1))-log(Same(jj,1));  
         
         yyy=(Same(jj+1,3)-Same(jj,3));
           
         
         Samee=[Samee;yyy xxx];
         
         zzzzzx=abs(xxx/yyy);
         
          Hero=[Hero;zzzzzx]; %%%% rate of change of SI: abs(Delta(log(Si)/Delta(t))
             
         end
         
         else
           Samee=[Same(:,1) Same(:,3)];
             
         end
         
        Samee;
        
             

         

         
         
         
         

          
        
         
         Same=[];

 

 
         
     end
     
     
 end
 
 
    
     Same=[Storem(:,1) Storem(:,2) Storem(:,4) Storem(:,3)];
     
     
     length(Same);
     
     
Hero
 
% 
 save All_cohort_Si_rate_change.txt Hero -ASCII 
