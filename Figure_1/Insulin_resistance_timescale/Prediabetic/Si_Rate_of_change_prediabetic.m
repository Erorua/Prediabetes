function  delta_SI


% ============================================================================================
% Description
% ============================================================================================

%%% Author: Aurore Woller

%%% Date: December 2022

%%% Uni: Weizmann institute of Science

%%% Description: compute rate of change of  insulin sensitivity  

%%% Compute as abs(Delta(log(Si)/Delta(t))

%%% Data from Segal Cohort - Prediabetic people

%%% output: file "Prediabetic_cohort_Si_rate_change.txt"



% ============================================================================================
% Data
% ============================================================================================


%%% Eran Segal human data


Datap=readtable('./Prediabetic_Infos_Header.txt');



Data=[table2array(Datap)]; % matrix [Reg Date G I gender Age BMI]


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
 save Prediabetic_cohort_Si_rate_change.txt Hero -ASCII 
