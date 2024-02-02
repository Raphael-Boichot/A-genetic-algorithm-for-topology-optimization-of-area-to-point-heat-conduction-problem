load Etat_courant.mat
prob_mutation_max=0.1;
m=g;
close all
figure('Position',[100 100 800 800]);
for g=m:1:nb_generations
    tic
    
    %Mutation rate is decreased with epoch following an empirical law that
    %works well on this problem
    prob_mutation=prob_mutation_max*exp(-5.8*g/nb_generations);
    
    %Fitness calculation
    disp('Calculating fitness for each individual...');
    % Variables in this order : variance, moyenne, temperature_max, map
    % temperatures, map gradients, variance gradients
    parfor i=1:1:population_size
        [~,~, ~, variance_border,~, ~,t_max, temp, ~,var_grad]=finite_temp_direct_sparse(k0*kp_k0,k0,T_ref,step_x,p,population(:,:,i));
        fitness(i,g)=t_max;
    end
    
    %Best topologies are kept for next step
    temp_temp=[(1:1:population_size)',fitness(:,g)];
    pop_classe=sortrows(temp_temp, 2);
    indice=pop_classe(1:population_best,1);
    opti_p_crois=table(indice(1),1);
    opti_p_mut=table(indice(1),2);

    %elitism
    new_population(:,:,1)=population(:,:,indice(1));
    new_population(:,:,2)=population(:,:,indice(1));
    
    if (g>1)
        if (fitness(1,g)==fitness(1,g-1))
             disp('---------No better configuration detected---------')
            %disp('Applying the ESO algorithm')
            %[ESO_shape,growth,etching] = fun_ESO_algorithm(population(:,:,indice(1)),k0*kp_k0,k0,T_ref,step_x,p);
            %new_population(:,:,1)=ESO_shape;
        else
            disp('------------Better configuration detected----------')
        end
    end
    
    topology_history(:,:,g)=new_population(:,:,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %New population of childrens
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp('Applying the mutation/crossover algorithm...');
    parfor m=3:1:(population_size)
        [new_population(:,:,m),table(m,:)]=gene_enfant(population,kp_k0,k0, conductive_pixels, prob_mutation, prob_croisement, population_best,indice);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get the best topology
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    best_topology=new_population(:,:,1);
    best_image=zeros(height,width,3);
    mean_topology=zeros(height,width,3);
    
    checksum=0;
    for k = 1:1:height
        for l = 1:1:width
            
            if best_topology(k,l)==k0
                best_image(k,l,1)=255;
                best_image(k,l,2)=255;
                best_image(k,l,3)=255;
                mean_topology(k,l,1)=255-((mean(new_population(k,l,:))-k0)/(k0*kp_k0-k0))*255;
                mean_topology(k,l,2)=255-((mean(new_population(k,l,:))-k0)/(k0*kp_k0-k0))*255;
                mean_topology(k,l,3)=255-((mean(new_population(k,l,:))-k0)/(k0*kp_k0-k0))*255;
            end
            
            if best_topology(k,l)==k0*kp_k0
                best_image(k,l,1)=0;
                best_image(k,l,2)=0;
                best_image(k,l,3)=0;
                mean_topology(k,l,1)=255-((mean(new_population(k,l,:))-k0)/(k0*kp_k0-k0))*255;
                mean_topology(k,l,2)=255-((mean(new_population(k,l,:))-k0)/(k0*kp_k0-k0))*255;
                mean_topology(k,l,3)=255-((mean(new_population(k,l,:))-k0)/(k0*kp_k0-k0))*255;
                checksum=checksum+1;
            end
            
            if best_topology(k,l)==-2
                best_image(k,l,1)=127;
                best_image(k,l,2)=127;
                best_image(k,l,3)=127;
                mean_topology(k,l,1)=127;
                mean_topology(k,l,2)=127;
                mean_topology(k,l,3)=127;
            end
            
            if best_topology(k,l)==-3
                best_image(k,l,1)=0;
                best_image(k,l,2)=0;
                best_image(k,l,3)=255;
                mean_topology(k,l,1)=0;
                mean_topology(k,l,2)=0;
                mean_topology(k,l,3)=255;
                
            end
            
        end
    end
    
    best_image=uint8(best_image);
    mean_topology=uint8(mean_topology);
    
    miroir_best=fliplr(best_image(1:height,1:width-1,:));
    miroir_best2=fliplr(miroir_best);
    miroir_mean=fliplr(mean_topology(1:height,1:width-1,:));
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
    disp(['Sum of conductive cells: ', num2str(conductive_pixels), ' must be equal to : ', num2str(checksum)]);
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
        residuals(1,g-1) = (fitness(1,g-1)-fitness(1,g))/norme_iteration_1;
        subplot(2,4,1);
        plot(log10(residuals), '.r');
        title('Residuals');
        xlabel('Generation');
        ylabel('log10 value');
    end
    
    subplot(2,4,2);
    imagesc(best_image);
    title('Best topology');
    
    subplot(2,4,3);
    imagesc(mean_topology);
    title('Average topology');
    
    subplot(2,4,4);
    [distance,somme_entropie, entropie, border_variance, variance, moyenne,t_max,temp_affichage, grad, variance_grad]=finite_temp_direct_sparse(k0*kp_k0,k0,T_ref,step_x,p,population(:,:,1));
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
    population=new_population;
    
    saveas(gcf,['Figure\Figure_',num2str(g,'%06.f')],'png');
    saveas(gcf,'Figure.png');
end

