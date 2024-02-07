%*****************Automate cellulaire*********************INPG/BOICHOT/2008
function [distance,somme_entropie, entropie, border_variance,variance, moyenne_temp,t_max_sortie,temp,grad, variance_grad]=finite_temp_direct_sparse(cond_haute,cond_basse,temp_puits,pas_x,p_vol,condi_limites)

%****Pré-allocation de la taille des matrices utilisées dans les boucles***
[hauteur,largeur,~]=size(condi_limites);
temp=ones(hauteur,largeur).*temp_puits;
condu_tab=zeros(hauteur,largeur,5);
    
%****création de la matrice des tables de conductances entre noeuds********    
        %table équivalente provisoire
        condi_temp=condi_limites;
        
        for k=1:1:(hauteur)
        for l=1:1:(largeur)
        if condi_temp(k,l)==-2; condi_temp(k,l)=1e-9;end
        if condi_temp(k,l)==-3; condi_temp(k,l)=cond_haute;end
        end
        end
        
        for k=2:1:(hauteur-1)
        for l=2:1:(largeur-1)

            %calcul des conductances entre deux noeuds
            condu_tab(k, l, 2) = pas_x / ((pas_x / 2) / condi_temp(k-1,l) + (pas_x / 2) / condi_temp(k, l));
            condu_tab(k, l, 3) = pas_x / ((pas_x / 2) / condi_temp(k,l+1) + (pas_x / 2) / condi_temp(k, l));
            condu_tab(k, l, 4) = pas_x / ((pas_x / 2) / condi_temp(k+1,l) + (pas_x / 2) / condi_temp(k, l));
            condu_tab(k, l, 5) = pas_x / ((pas_x / 2) / condi_temp(k,l-1) + (pas_x / 2) / condi_temp(k, l));
        
        end
        end


%******************************************Résolution équation chaleur%différences finies et méthode directe**

adresse = zeros(hauteur, largeur);
k=1;
%création de la matrice d'adressage
        for i=1:1:hauteur
        for j=1:1:largeur

            adresse(i,j)=k;
            k=k+1;
            
        end
        end
        
%création des matrices vides

B=zeros((hauteur)*(largeur),1);
k=0;

 %construction du système linéaire pour la méthode directe
    for i=1:1:(hauteur)
    for j=1:1:(largeur)
    
        ind_1=adresse(i,j);
        
        if not (condi_limites(i,j)==-2||condi_limites(i,j)==-3)
        ind_2=adresse(i-1,j);
        ind_3=adresse(i,j+1);
        ind_4=adresse(i+1,j);
        ind_5=adresse(i,j-1);
        somme_cond = condu_tab(i, j, 2)+condu_tab(i, j, 3)+condu_tab(i, j, 4)+condu_tab(i, j, 5);
        end
        
        fin=0;
        %On balaye les 5 cellules et leurs pondérations
        %Cellule centrale (1)
        if condi_limites(i,j)==-2||condi_limites(i,j)==-3
            B(ind_1)=temp(i,j);
            k=k+1;
            ligne(k)=ind_1;
            colonne(k)=ind_1;
            valeur(k)=1;
            fin=1;
        else
            k=k+1;
            ligne(k)=ind_1;
            colonne(k)=ind_1;
            valeur(k)=-1*somme_cond;
            if condi_limites(i,j)==cond_basse; B(ind_1)=B(ind_1)-p_vol*pas_x^2;end
        end
        
        
        if not (fin==1)
        %Cellule N (2)
        if condi_limites(i-1,j)==-2||condi_limites(i-1,j)==-3
        B(ind_1)=B(ind_1)-temp(i-1,j)*condu_tab(i, j, 2);
        else
            k=k+1;
            ligne(k)=ind_1;
            colonne(k)=ind_2;
            valeur(k)=1*condu_tab(i, j, 2);
        end
        
        %Cellule E (3)
        if condi_limites(i,j+1)==-2||condi_limites(i,j+1)==-3
        B(ind_1)=B(ind_1)-temp(i,j+1)*condu_tab(i, j, 3);
        else
            k=k+1;
            ligne(k)=ind_1;
            colonne(k)=ind_3;
            valeur(k)=1*condu_tab(i, j, 3);  
        end
        
        %Cellule N (4)
        if condi_limites(i+1,j)==-2||condi_limites(i+1,j)==-3
        B(ind_1)=B(ind_1)-temp(i+1,j)*condu_tab(i, j, 4);
        else
            k=k+1;
            ligne(k)=ind_1;
            colonne(k)=ind_4;
            valeur(k)=1*condu_tab(i, j, 4);
        end
        
        %Cellule O (5)
        if condi_limites(i,j-1)==-2||condi_limites(i,j-1)==-3
        B(ind_1)=B(ind_1)-temp(i,j-1)*condu_tab(i, j, 5);
        else
            k=k+1;
            ligne(k)=ind_1;
            colonne(k)=ind_5;
            valeur(k)=1*condu_tab(i, j, 5);
                    
        end
        end

    end
    end    

A=sparse(ligne,colonne,valeur);

T=A\B;

%création de la matrice des températures
        for i=1:1:hauteur
        for j=1:1:largeur

            temp(i,j)=T(adresse(i,j));
           
        end
        end


t_max_sortie = max(max(temp));
variance=std(std(temp))^2;
[gradX, gradY]=gradient(temp(2:end-1,2:end-1));
grad=(gradX.^2+gradY.^2).^0.5;
variance_grad=std(std(grad))^2;
moyenne_temp=mean(mean(temp));
border_variance=std([temp(2:end-1,2)' temp(2,2:end-1)]);

%*****Calcul de la table d'entropie******
somme_entropie=0;

entropie=zeros(hauteur,largeur);

for k=1:1:hauteur
for l=1:1:largeur

entropie(k,l)=0;
    
if not(condi_limites(k,l)==-2||condi_limites(k,l)==-3) 

    flux1=condu_tab(k,l,2)*(temp(k-1,l)-temp(k,l));
    flux2=condu_tab(k,l,3)*(temp(k,l+1)-temp(k,l));
    flux3=condu_tab(k,l,4)*(temp(k+1,l)-temp(k,l));
    flux4=condu_tab(k,l,5)*(temp(k,l-1)-temp(k,l));
    
    flux5=0;
    if condi_limites(k,l)==cond_basse; flux5=p_vol*pas_x*pas_x; end
    if not(condi_limites(k-1,l)==-2)
        entropie(k,l)=entropie(k,l)+(abs(flux1/temp(k,l)-abs(flux1/temp(k-1,l))))*sign(temp(k,l)-temp(k-1,l));
    end
    if not(condi_limites(k,l+1)==-2)
        entropie(k,l)=entropie(k,l)+(abs(flux2/temp(k,l)-abs(flux2/temp(k,l+1))))*sign(temp(k,l)-temp(k,l+1));
    end
    if not(condi_limites(k+1,l)==-2)
        entropie(k,l)=entropie(k,l)+(abs(flux3/temp(k,l)-abs(flux3/temp(k+1,l))))*sign(temp(k,l)-temp(k+1,l));
    end
    if not(condi_limites(k,l-1)==-2)
        entropie(k,l)=entropie(k,l)+(abs(flux4/temp(k,l)-abs(flux4/temp(k,l-1))))*sign(temp(k,l)-temp(k,l-1));
    end
    
    entropie(k,l)=entropie(k,l)+flux5/temp(k,l);
    somme_entropie=somme_entropie+entropie(k,l);

end

end
end

[row,col] = find(temp==t_max_sortie);
row_ref=hauteur;
col_ref=largeur;
distance=(row-row_ref)^2+(col-col_ref)^2;




