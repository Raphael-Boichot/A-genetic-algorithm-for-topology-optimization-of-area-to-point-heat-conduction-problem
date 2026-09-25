%https://doi.org/10.1016/j.ijthermalsci.2016.05.015
%https://github.com/Raphael-Boichot/A-genetic-algorithm-for-topology-optimization-of-area-to-point-heat-conduction-problem
function [child,param]=generate_child(population,kp_k0,k0, number_conductive_cells, prob_mutation, prob_crossover, population_best,indice)

%We select randomly two good parents
parent_1=population(:,:,indice(ceil(rand*population_best)));
parent_2=population(:,:,indice(ceil(rand*population_best)));

test=rand;
prob_crossover_local=prob_crossover*rand;
prob_mutation_local=prob_mutation*rand;
param=[prob_crossover_local prob_mutation_local];
[height,width,~]=size(parent_1);
direction=1;

if test<0.5
    %vertical crossover and mutation
    for i=1:1:height
        for j=1:1:width
            if rand<prob_crossover_local; direction=direction*-1;end
            if direction==1; child(i,j,:)=parent_1(i,j,:);end
            if direction==-1;child(i,j,:)=parent_2(i,j,:);end
        end
    end
    child = correction_mutation(child, kp_k0,k0, number_conductive_cells,prob_mutation_local);
else

    %horizontal crossover and mutation
    for j=1:1:width
        for i=1:1:height
            if rand<prob_crossover_local; direction=direction*-1;end
            if direction==1; child(i,j,:)=parent_1(i,j,:);end
            if direction==-1;child(i,j,:)=parent_2(i,j,:);end
        end
    end
    child = correction_mutation(child, kp_k0,k0, number_conductive_cells,prob_mutation_local);
end