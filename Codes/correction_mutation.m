function enfant2=correction_mutation(enfant,kp_k0,k0,nombre_pixels_conducteurs, p_mutation);

[hauteur,largeur,profondeur]=size(enfant);
pixels_cond=0;

for j=2:1:largeur-1
    for i=2:1:hauteur-1
        if enfant(i,j)==k0*kp_k0; pixels_cond=pixels_cond+1;end;
    end;
end;


%****************calcul de la liste des zones où la correction est possible**************
n=0;
m=0;
for j=2:1:largeur-1
    for i=2:1:hauteur-1
        
        vert=0;
        blanc=0;
        
        if enfant(i+1,j,1)==k0; blanc=blanc+1;end;
        if enfant(i-1,j,1)==k0; blanc=blanc+1;end;
        if enfant(i,j+1,1)==k0; blanc=blanc+1;end;
        if enfant(i,j-1,1)==k0; blanc=blanc+1;end;
        if enfant(i+1,j+1,1)==k0; blanc=blanc+1;end;
        if enfant(i-1,j-1,1)==k0; blanc=blanc+1;end;
        if enfant(i+1,j-1,1)==k0; blanc=blanc+1;end;
        if enfant(i-1,j+1,1)==k0; blanc=blanc+1;end;
        
        if enfant(i+1,j,1)==k0*kp_k0; vert=vert+1;end;
        if enfant(i-1,j,1)==k0*kp_k0; vert=vert+1;end;
        if enfant(i,j+1,1)==k0*kp_k0; vert=vert+1;end;
        if enfant(i,j-1,1)==k0*kp_k0; vert=vert+1;end;
        if enfant(i+1,j+1,1)==k0*kp_k0; vert=vert+1;end;
        if enfant(i-1,j-1,1)==k0*kp_k0; vert=vert+1;end;
        if enfant(i+1,j-1,1)==k0*kp_k0; vert=vert+1;end;
        if enfant(i-1,j+1,1)==k0*kp_k0; vert=vert+1;end;      
        
        %Futures positions corrigeables pour les pixels conducteurs->non
        %conducteurs
        if enfant(i,j,1)==k0 && not(vert==0);
        n=n+1;    
        mut_pos_blanc(n,1)=i;
        mut_pos_blanc(n,2)=j;
        mut_pos_blanc(n,3)=1;
        end;
        
        %Futures positions corrigeables pour les pixels non
        %conducteurs->conducteurs
        if enfant(i,j,1)==k0*kp_k0 && not(blanc==0);
        m=m+1;    
        mut_pos_vert(m,1)=i;
        mut_pos_vert(m,2)=j;
        mut_pos_vert(m,3)=1;
        end;        

    end;
end;

difference=pixels_cond-nombre_pixels_conducteurs;

%correction pixels conducteur vers non conducteur
if difference<0;
o=0;
while o<abs(difference);
    pos=ceil(rand*n);
    if mut_pos_blanc(pos,3)==1;
        i=mut_pos_blanc(pos,1);
        j=mut_pos_blanc(pos,2);
        enfant(i,j)=k0*kp_k0;
        mut_pos_blanc(pos,3)=0;
        o=o+1;    
    end;
end;
end;
%correction pixels verts vers blancs
if difference>0;
o=0;
while o<abs(difference);
    pos=ceil(rand*m);
    if mut_pos_vert(pos,3)==1;
        i=mut_pos_vert(pos,1);
        j=mut_pos_vert(pos,2);
        enfant(i,j)=k0;
        mut_pos_vert(pos,3)=0;
        o=o+1;    
    end;
end;
end;

%Mutation sur les positions restantes après correction
m_restant=sum(mut_pos_blanc(:,3));
n_restant=sum(mut_pos_vert(:,3));
nombre_mutations=max(round(p_mutation*min(m_restant, n_restant)),0);
%mutation pixels blancs vers verts
o=0;
while o<nombre_mutations;
    pos=ceil(rand*n);
    if mut_pos_blanc(pos,3)==1;
        i=mut_pos_blanc(pos,1);
        j=mut_pos_blanc(pos,2);
        enfant(i,j)=k0*kp_k0;
        mut_pos_blanc(pos,3)=0;
        o=o+1;    
    end;
end;

%mutation pixels verts vers blancs
o=0;
while o<nombre_mutations;
    pos=ceil(rand*m);
    if mut_pos_vert(pos,3)==1;
        i=mut_pos_vert(pos,1);
        j=mut_pos_vert(pos,2);
        enfant(i,j)=k0;
        mut_pos_vert(pos,3)=0;
        o=o+1;    
    end;
end;


enfant2=enfant;