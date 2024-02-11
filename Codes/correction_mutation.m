%https://doi.org/10.1016/j.ijthermalsci.2016.05.015
%https://github.com/Raphael-Boichot/A-genetic-algorithm-for-topology-optimization-of-area-to-point-heat-conduction-problem
function child2=correction_mutation(child,kp_k0,k0,number_conductive_cells, p_mutation)
[height,width,~]=size(child);
pixels_cond=0;

%Mutation, here we just add conductive matter randomly around existing conductive matter
number_of_mutations=ceil(p_mutation*number_conductive_cells);
cell_counter=0;
while not(number_of_mutations==cell_counter)
    row=ceil(rand*height);
    column=ceil(rand*width);
    if (child(row,column)==k0)&&...
            ((child(row-1,column)==k0*kp_k0)||...
            (child(row+1,column)==k0*kp_k0)||...
            (child(row,column-1)==k0*kp_k0)||...
            (child(row,column+1)==k0*kp_k0)||...
            (child(row+1,column+1)==k0*kp_k0)||...
            (child(row+1,column-1)==k0*kp_k0)||...
            (child(row-1,column-1)==k0*kp_k0)||...
            (child(row-1,column+1)==k0*kp_k0))
        child(row,column)=k0*kp_k0;
        cell_counter=cell_counter+1;
    end
end

%here we verify the actual number of conductive cells as crossover was not conservative anyway
for j=2:1:width-1
    for i=2:1:height-1
        if child(i,j)==k0*kp_k0; pixels_cond=pixels_cond+1;end
    end
end

%random correction to keep the number of conductive cells constant with time that acts as mutation too
if pixels_cond<number_conductive_cells
    while not(pixels_cond==number_conductive_cells)
        row=ceil(rand*height);
        column=ceil(rand*width);
        if (child(row,column)==k0)&&...
                ((child(row-1,column)==k0*kp_k0)||...
                (child(row+1,column)==k0*kp_k0)||...
                (child(row,column-1)==k0*kp_k0)||...
                (child(row,column+1)==k0*kp_k0)||...
                (child(row+1,column+1)==k0*kp_k0)||...
                (child(row+1,column-1)==k0*kp_k0)||...
                (child(row-1,column-1)==k0*kp_k0)||...
                (child(row-1,column+1)==k0*kp_k0))
            child(row,column)=k0*kp_k0;
            pixels_cond=pixels_cond+1;
        end
    end
end

if pixels_cond>number_conductive_cells
    while not(pixels_cond==number_conductive_cells)
        row=ceil(rand*height);
        column=ceil(rand*width);
        if (child(row,column)==k0*kp_k0)&&...
                ((child(row-1,column)==k0)||...
                (child(row+1,column)==k0)||...
                (child(row,column-1)==k0)||...
                (child(row,column+1)==k0)||...
                (child(row+1,column+1)==k0)||...
                (child(row+1,column-1)==k0)||...
                (child(row-1,column-1)==k0)||...
                (child(row-1,column+1)==k0))
            child(row,column)=k0;
            pixels_cond=pixels_cond-1;
        end
    end
end

% old code from 2015, very inefficient
% %****************calcul de la liste des zones où la correction est possible**************
% n=0;
% m=0;
% for j=2:1:width-1
%     for i=2:1:height-1
%
%         vert=0;
%         blanc=0;
%
%         if child(i+1,j,1)==k0; blanc=blanc+1;end;
%         if child(i-1,j,1)==k0; blanc=blanc+1;end;
%         if child(i,j+1,1)==k0; blanc=blanc+1;end;
%         if child(i,j-1,1)==k0; blanc=blanc+1;end;
%         if child(i+1,j+1,1)==k0; blanc=blanc+1;end;
%         if child(i-1,j-1,1)==k0; blanc=blanc+1;end;
%         if child(i+1,j-1,1)==k0; blanc=blanc+1;end;
%         if child(i-1,j+1,1)==k0; blanc=blanc+1;end;
%
%         if child(i+1,j,1)==k0*kp_k0; vert=vert+1;end;
%         if child(i-1,j,1)==k0*kp_k0; vert=vert+1;end;
%         if child(i,j+1,1)==k0*kp_k0; vert=vert+1;end;
%         if child(i,j-1,1)==k0*kp_k0; vert=vert+1;end;
%         if child(i+1,j+1,1)==k0*kp_k0; vert=vert+1;end;
%         if child(i-1,j-1,1)==k0*kp_k0; vert=vert+1;end;
%         if child(i+1,j-1,1)==k0*kp_k0; vert=vert+1;end;
%         if child(i-1,j+1,1)==k0*kp_k0; vert=vert+1;end;
%
%         %Futures positions corrigeables pour les pixels conducteurs->non
%         %conducteurs
%         if child(i,j,1)==k0 && not(vert==0);
%         n=n+1;
%         mut_pos_blanc(n,1)=i;
%         mut_pos_blanc(n,2)=j;
%         mut_pos_blanc(n,3)=1;
%         end;
%
%         %Futures positions corrigeables pour les pixels non
%         %conducteurs->conducteurs
%         if child(i,j,1)==k0*kp_k0 && not(blanc==0);
%         m=m+1;
%         mut_pos_vert(m,1)=i;
%         mut_pos_vert(m,2)=j;
%         mut_pos_vert(m,3)=1;
%         end;
%
%     end;
% end;
%
% difference=pixels_cond-number_conductive_cells;
%
% %correction pixels conducteur vers non conducteur
% if difference<0;
% o=0;
% while o<abs(difference);
%     pos=ceil(rand*n);
%     if mut_pos_blanc(pos,3)==1;
%         i=mut_pos_blanc(pos,1);
%         j=mut_pos_blanc(pos,2);
%         child(i,j)=k0*kp_k0;
%         mut_pos_blanc(pos,3)=0;
%         o=o+1;
%     end;
% end;
% end;
% %correction pixels verts vers blancs
% if difference>0;
% o=0;
% while o<abs(difference);
%     pos=ceil(rand*m);
%     if mut_pos_vert(pos,3)==1;
%         i=mut_pos_vert(pos,1);
%         j=mut_pos_vert(pos,2);
%         child(i,j)=k0;
%         mut_pos_vert(pos,3)=0;
%         o=o+1;
%     end;
% end;
% end;
%
% %Mutation sur les positions restantes après correction
% m_restant=sum(mut_pos_blanc(:,3));
% n_restant=sum(mut_pos_vert(:,3));
% nombre_mutations=max(round(p_mutation*min(m_restant, n_restant)),0);
% %mutation pixels blancs vers verts
% o=0;
% while o<nombre_mutations;
%     pos=ceil(rand*n);
%     if mut_pos_blanc(pos,3)==1;
%         i=mut_pos_blanc(pos,1);
%         j=mut_pos_blanc(pos,2);
%         child(i,j)=k0*kp_k0;
%         mut_pos_blanc(pos,3)=0;
%         o=o+1;
%     end;
% end;
%
% %mutation pixels verts vers blancs
% o=0;
% while o<nombre_mutations;
%     pos=ceil(rand*m);
%     if mut_pos_vert(pos,3)==1;
%         i=mut_pos_vert(pos,1);
%         j=mut_pos_vert(pos,2);
%         child(i,j)=k0;
%         mut_pos_vert(pos,3)=0;
%         o=o+1;
%     end;
% end;


child2=child;