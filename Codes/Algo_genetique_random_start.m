%https://doi.org/10.1016/j.ijthermalsci.2016.05.015
clc;
clear;
close all;
%---------Conditions for thermal science-----------------------------------
k0=1;
kp_k0=10;
p=1e6;
pas_x=0.001;
T_ref=298;
taux_remplissage=0.3;
imname='50x100.bmp';
%--------------------------------------------------------------------------

%--------Hyper parameters for genetic algorithm----------------------------
taille_pop=1000;
meilleurs=200;
nb_generations=10000;
prob_croisement=0.2;
probabilite_mutation_maximale=0.1;
%--------------------------------------------------------------------------

rng('shuffle', 'twister')
mkdir('Figure');
mkdir('Best_topology');
mkdir('Average_topology');


T_comp=0;
table=zeros(taille_pop,2);
couleurs=imread(imname);
[hauteur,largeur,profondeur]=size(couleurs);
conditi_limites_ini=zeros(hauteur,largeur);
pixels_blancs=0;
for k = 1:1:hauteur
    for l = 1:1:largeur
        
        rouge = couleurs(k,l,1);
        vert = couleurs(k,l,2);
        bleu = couleurs(k,l,3);
        
        if (rouge == 127) && (vert == 127) && (bleu == 127)
            conditi_limites_ini(k,l)=-2;
        end
        
        if (rouge == 0) && (vert == 0) && (bleu == 255)
            conditi_limites_ini(k,l)=-3;
        end
        
        if (rouge == 255) && (vert == 255) && (bleu == 255)
            conditi_limites_ini(k,l)=k0;
            pixels_blancs=pixels_blancs+1;
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial population creation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

population=zeros(hauteur, largeur,taille_pop);
nombre_pixels_conducteurs=ceil(pixels_blancs*taux_remplissage);
disp('Creating the initial population frow scratch...');
tic
parfor i=1:1:taille_pop
    population(:,:,i)=init_image(conditi_limites_ini,nombre_pixels_conducteurs, k0, kp_k0);
end
topology_history=zeros(hauteur, largeur, nb_generations);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fitness evaluation and sorting by fitness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Position',[100 100 800 800]);
for g=1:1:nb_generations
    tic
    
    %Mutation rate is decreased with epoch following an empirical law that
    %works well on this problem
    prob_mutation=probabilite_mutation_maximale*exp(-5.8*g/nb_generations);
    
    %Fitness calculation
    disp('Calculating fitness for each individual...');
    % Variables in this order : variance, moyenne, temperature_max, map
    % temperatures, map gradients, variance gradients
    parfor i=1:1:taille_pop
        [~,~, ~, variance_border,~, ~,t_max, temp, ~,var_grad]=finite_temp_direct_sparse(k0*kp_k0,k0,T_ref,pas_x,p,population(:,:,i));
        fitness(i,g)=t_max;
    end
    
    %Best topologies are kept for next step
    temp_temp=[(1:1:taille_pop)',fitness(:,g)];
    pop_classe=sortrows(temp_temp, 2);
    indice=pop_classe(1:meilleurs,1);
    opti_p_crois=table(indice(1),1);
    opti_p_mut=table(indice(1),2);
    
    if (g>1)
        if (fitness(1,g)==fitness(1,g-1))
            disp('No best configuration detected')
        else
            disp('Best configuration detected')
        end
    end
    
    %elitism with two times the best topology, you'll see why after
    nouvelle_population(:,:,1)=population(:,:,indice(1));
    nouvelle_population(:,:,2)=population(:,:,indice(1));
    topology_history(:,:,g)=population(:,:,indice(1));
    %For keeping trace of the hyper parameters leading to success
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %New population of childrens
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp('Applying the mutation/crossover algorithm...');
    parfor m=3:1:(taille_pop)
        [nouvelle_population(:,:,m),table(m,:)]=gene_enfant(population,kp_k0,k0, nombre_pixels_conducteurs, prob_mutation, prob_croisement, meilleurs,indice);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get the best topology
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    meilleure_image=nouvelle_population(:,:,1);
    best_image=zeros(hauteur,largeur,3);
    image_moyenne=zeros(hauteur,largeur,3);
    
    somme_controle=0;
    for k = 1:1:hauteur
        for l = 1:1:largeur
            
            if meilleure_image(k,l)==k0
                best_image(k,l,1)=255;
                best_image(k,l,2)=255;
                best_image(k,l,3)=255;
                image_moyenne(k,l,1)=255-((mean(nouvelle_population(k,l,:))-k0)/(k0*kp_k0-k0))*255;
                image_moyenne(k,l,2)=255-((mean(nouvelle_population(k,l,:))-k0)/(k0*kp_k0-k0))*255;
                image_moyenne(k,l,3)=255-((mean(nouvelle_population(k,l,:))-k0)/(k0*kp_k0-k0))*255;
            end
            
            if meilleure_image(k,l)==k0*kp_k0
                best_image(k,l,1)=0;
                best_image(k,l,2)=0;
                best_image(k,l,3)=0;
                image_moyenne(k,l,1)=255-((mean(nouvelle_population(k,l,:))-k0)/(k0*kp_k0-k0))*255;
                image_moyenne(k,l,2)=255-((mean(nouvelle_population(k,l,:))-k0)/(k0*kp_k0-k0))*255;
                image_moyenne(k,l,3)=255-((mean(nouvelle_population(k,l,:))-k0)/(k0*kp_k0-k0))*255;
                somme_controle=somme_controle+1;
            end
            
            if meilleure_image(k,l)==-2
                best_image(k,l,1)=127;
                best_image(k,l,2)=127;
                best_image(k,l,3)=127;
                image_moyenne(k,l,1)=127;
                image_moyenne(k,l,2)=127;
                image_moyenne(k,l,3)=127;
            end
            
            if meilleure_image(k,l)==-3
                best_image(k,l,1)=0;
                best_image(k,l,2)=0;
                best_image(k,l,3)=255;
                image_moyenne(k,l,1)=0;
                image_moyenne(k,l,2)=0;
                image_moyenne(k,l,3)=255;
                
            end
            
        end
    end
    
    best_image=uint8(best_image);
    image_moyenne=uint8(image_moyenne);
    
    miroir_best=fliplr(best_image(1:hauteur,1:largeur-1,:));
    miroir_best2=fliplr(miroir_best);
    miroir_mean=fliplr(image_moyenne(1:hauteur,1:largeur-1,:));
    miroir_mean2=fliplr(miroir_mean);
    
    if g==1
        imwrite([miroir_best2,miroir_best],['Best_topology\Best_topology_',num2str(g,'%06.f'),'.png']);
    end
    if (g>1)
        if not(sum(abs(sum(topology_history(:,:,g)-topology_history(:,:,g-1)))))==0
            imwrite([miroir_best2,miroir_best],['Best_topology\Best_topology_',num2str(g,'%06.f'),'.png']);
        end
    end
    
    imwrite([miroir_best2,miroir_best],'Best_topology.png');
    imwrite([miroir_mean2,miroir_mean],'Average_topology.png');
    imwrite([miroir_mean2,miroir_mean],['Average_topology\Average_topology_',num2str(g,'%06.f'),'.png']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %output to console
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %clc;
    disp(['---------Epoch: ',num2str(g),'---------']);
    disp(['Sum of conductive cells: ', num2str(nombre_pixels_conducteurs), ' must be equal to : ', num2str(somme_controle)]);
    disp(['Current maximal mutation probability: ', num2str(prob_mutation)]);
    disp(['Last successfull - crossover rate: ', num2str(opti_p_crois), ' / mutation rate : ', num2str(opti_p_mut)]);
    disp(['Best fitness: ', num2str(fitness(1,g))]);
    save Etat_courant.mat
    toc
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %output to plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if g==2
        norme_iteration_1=fitness(1,1)-fitness(1,2);
    end
    
    if g>1
        residus(1,g-1) = (fitness(1,g-1)-fitness(1,g))/norme_iteration_1;
        subplot(2,4,1);
        plot(log10(residus), '.r');
        title('Residuals');
        xlabel('Generation');
        ylabel('log10 value');
    end
    
    subplot(2,4,2);
    imagesc(best_image);
    title('Best topology');
    
    subplot(2,4,3);
    imagesc(image_moyenne);
    title('Average topology');
    
    subplot(2,4,4);
    [distance,somme_entropie, entropie, border_variance, variance, moyenne,t_max,temp_affichage, grad, variance_grad]=finite_temp_direct_sparse(k0*kp_k0,k0,T_ref,pas_x,p,population(:,:,1));
    imagesc(temp_affichage);
    title('Objective function');
    
    P_1(g,1)=opti_p_crois;
    P_1(g,2)=opti_p_mut;
    P_1(g,3)=prob_mutation;
    
    subplot(2,4,5);
    plot(P_1(:,1), '.b');
    title('Mutation rate');
    xlabel('Generation');
    ylabel('Value');
    
    subplot(2,4,6);
    plot(log10(P_1(:,2)), '.b')
    hold on
    plot(log10(P_1(:,3)), '.r');
    hold off
    title('log10 mutation rate');
    xlabel('Generation');
    
    subplot(2,4,7);
    imagesc(log10(entropie(2:end-1,2:end-1)));
    title('Log10 Entropy');
    
    subplot(2,4,8);
    imagesc(grad);
    title('Thermal gradients');
    colormap jet
    
    pause(0.01);
    population=nouvelle_population;
    
    saveas(gcf,['Figure\Figure_',num2str(g,'%06.f')],'png');
    saveas(gcf,'Figure.png');
end

