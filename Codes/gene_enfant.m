function [enfant,param]=gene_enfant(population,kp_k0,k0, nombre_pixels_conducteurs, prob_mutation, prob_croisement, meilleurs,indice)

%Sélection aléatoire des parents parmis les meilleurs
    parent_1=population(:,:,indice(ceil(rand*meilleurs)));
    parent_2=population(:,:,indice(ceil(rand*meilleurs)));
        
    test=rand;
    var1=prob_croisement*rand;
    var2=prob_mutation*rand;
    param=[var1 var2];
    
    if test<0.5; 
    %Croisement vertical et correction
    enfant = croisement_vert(parent_1,parent_2,var1);    
    enfant = correction_mutation(enfant, kp_k0,k0, nombre_pixels_conducteurs,var2);
    else  
    %Croisement horizontal et correction
    enfant = croisement_hori(parent_1,parent_2,var1); 
    enfant = correction_mutation(enfant, kp_k0,k0, nombre_pixels_conducteurs,var2);
    end;