function enfant=croisement_hori(parent1,parent2,prob_croisement)

%récupère les dimensions de l'image
[hauteur,largeur,profondeur]=size(parent1);

test=1;

for j=1:1:largeur
    for i=1:1:hauteur
    
        if rand<prob_croisement; test=test*-1;end;
        
        if test==1; enfant(i,j,:)=parent1(i,j,:);end;
        if test==-1;enfant(i,j,:)=parent2(i,j,:);end;
        
    end;
end;