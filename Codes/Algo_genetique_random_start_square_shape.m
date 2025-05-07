%https://doi.org/10.1016/j.ijthermalsci.2016.05.015
%https://github.com/Raphael-Boichot/A-genetic-algorithm-for-topology-optimization-of-area-to-point-heat-conduction-problem
clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%User parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------Conditions for thermal science-----------------------------------
k0=1;                           %conductivity of the heating matter
kp_k0=10;                       %conductivity of the draining material
p=1e6;                          %surface of volume power density
delta_x=0.001;                  %size of x/y square cells
Heat_sink_temperature=298;      %self explanatory
filling_ratio=0.3;              %ratio of conductive matter on the surface
starting_image='50x100_spicy.bmp';
%--------Hyper parameters for genetic algorithm----------------------------
%hyper parameters have been optimized with blood, sweat, and tears
%Believe me, they are efficient for tackling this problem
population_size=1000;           %size of the topology dataset
population_best=200;            %rank of the last topology allowed to survive
split_crossover=5;              %Image will be splitted split_crossover times in average during crossover (constant with epoch)
prob_mutation_max=0.05;         %mutation probability at each cell (decrases with epoch, see code)
convergence_criterion=500;      %after n steps without better configuration, code stops
%--------------------------------------------------------------------------

rng('shuffle', 'twister')
mkdir('Figure');
mkdir('Best_topology');
mkdir('Average_topology');
figure('Position',[100 100 1300 1000]);

T_comp=0;
table=zeros(population_size,2);
pixels=imread(starting_image);
[height,width,layers]=size(pixels);
prob_crossover=split_crossover/(height*width);
Initial_boundary_limits=zeros(height,width);
non_conductive_pixels=0;
yet_conductive_pixels=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Image translation into boundary limits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:1:height
    for l = 1:1:width
        red = pixels(k,l,1);
        green = pixels(k,l,2);
        blue = pixels(k,l,3);
        if (red == 127) && (green == 127) && (blue == 127)
            Initial_boundary_limits(k,l)=-2;
        end
        if (red == 0) && (green == 0) && (blue == 255)
            Initial_boundary_limits(k,l)=-3;
        end
        if (red == 0) && (green == 0) && (blue == 0)
            Initial_boundary_limits(k,l)=k0*kp_k0;
            yet_conductive_pixels=yet_conductive_pixels+1;
        end
        if (red == 255) && (green == 255) && (blue == 255)
            Initial_boundary_limits(k,l)=k0;
            non_conductive_pixels=non_conductive_pixels+1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial population creation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
population=zeros(height, width,population_size);
conductive_pixels=ceil((non_conductive_pixels+yet_conductive_pixels)*filling_ratio);
checksum=conductive_pixels;
%topology_history=zeros(height, width, nb_generations);

%Check is a preceding session was crashed and reload it
if isfile('Current_state.mat')
    disp('Reloading last session and continuing from last known configuration...');
    disp('To start from scratch erase manually the *.mat file in directory');
    load Current_state.mat
    m=g;
else
    disp('Creating the initial population frow scratch and cleaning directories...');
    delete('Figure/*.png');
    delete('Best_topology/*.png');
    delete('Average_topology/*.png');
    parfor i=1:population_size
        population(:,:,i)=init_image(Initial_boundary_limits,conductive_pixels, yet_conductive_pixels, k0, kp_k0);
    end
    m=1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=m;
convergence_counter=0;
while convergence_counter<convergence_criterion

    %Mutation rate is decreased with epoch following an empirical law that works well on this problem
    prob_mutation=prob_mutation_max*exp(-5*g/10000);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %output to console
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    disp(' ');
    disp(['---------Epoch: ',num2str(g),'---------']);
    disp(['Checksum: ',num2str(checksum-conductive_pixels),' (must be 0)']);
    disp(['Current maximal mutation probability: ', num2str(prob_mutation)]);
    disp('Calculating fitness for each individual...');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %evaluate the fitness / choose objective function
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    parfor i=1:population_size
        % Variables output in this order :
        % 1. Distance of the hotest cell to the heat sink (scalar)
        % 2. Sum of cell entropy (scalar)
        % 3. Entropy map (matrix)
        % 4. Variance of temperatures accross the 1D adabatic borders (scalar)
        % 5. Variance of temperatures accross the 2D domain (scalar)
        % 6. Mean temperature (scalar)
        % 7. Maximal temperature accross the 2D domain (scalar)
        % 8. Map of temperatures (matrix)
        % 9. map of thermal gradients (matrix)
        % 10. Variance of gradients across the 2D domain (scalar)
        [~,~,~,~,~,~,t_max,~,~,~]=finite_temp_direct_sparse(k0*kp_k0,k0,Heat_sink_temperature,delta_x,p,population(:,:,i));
        fitness(i,g)=t_max;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get the best topology
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp_temp=[(1:1:population_size)',fitness(:,g)];
    pop_classe=sortrows(temp_temp, 2);
    indice=pop_classe(1:population_best,1);
    opti_p_crois=table(indice(1),1);
    opti_p_mut=table(indice(1),2);

    %Best topologies are kept for next step because eugenics
    new_population(:,:,1)=population(:,:,indice(1));
    new_population(:,:,2)=population(:,:,indice(1));%this one can be used as target for another algorithm

    if (g>1)
        if (fitness(1,g)==fitness(1,g-1))
            disp('---------No better configuration detected---------')
            convergence_counter=convergence_counter+1;
        else
            disp('>>>>>>>>>Better configuration detected<<<<<<<<<<<<')
            convergence_counter=0;
        end
    end

    topology_history(:,:,g)=new_population(:,:,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %New population of childrens
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Applying the mutation/crossover algorithm...');
    parfor m=3:(population_size)
        [new_population(:,:,m),table(m,:)]=generate_child(population,kp_k0,k0, conductive_pixels, prob_mutation, prob_crossover, population_best,indice);
    end

    best_topology=new_population(:,:,1);
    best_image=zeros(height,width,3);
    mean_topology=zeros(height,width,3);
    checksum=0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Make fancy display
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    miroir_best=fliplr(best_image(1:height-1,1:width-1,:));
    miroir_best2=best_image(1:height-1,1:width-1,:);
    miroir_mean=fliplr(mean_topology(1:height-1,1:width-1,:));
    miroir_mean2=mean_topology(1:height-1,1:width-1,:);

    if g==1
        imwrite([miroir_best2,miroir_best],['Best_topology\Best_topology_',num2str(g,'%06.f'),'.png']);
    end
    if (g>1)
        if not(sum(abs(sum(topology_history(:,:,g)-topology_history(:,:,g-1)))))==0
            imwrite([miroir_best2,miroir_best],['Best_topology\Best_topology_',num2str(g,'%06.f'),'.png']);
        end
    end

    imwrite(best_image,'Best_topology.png');
    imwrite(mean_topology,'Average_topology.png');
    imwrite([miroir_mean2,miroir_mean],['Average_topology\Average_topology_',num2str(g,'%06.f'),'.png']);

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
    [~,~,entropy_map,~,~,~,t_max,map_temperature,grad,~]=finite_temp_direct_sparse(k0*kp_k0,k0,Heat_sink_temperature,delta_x,p,population(:,:,1));
    imagesc(map_temperature);
    title('Objective function');

    disp(['Maximal temperature: ',num2str(t_max)])
    disp(['Epochs without change: ',num2str(convergence_counter),'/',num2str(convergence_criterion), ' to convergence'])
    P_1(g,1)=opti_p_crois;
    P_1(g,2)=opti_p_mut;
    P_1(g,3)=prob_mutation;

    subplot(2,4,5);
    plot(P_1(:,1), '.b');
    title('Crossover rate');
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
    imagesc(log10(entropy_map(2:end-1,2:end-1)));
    title('Log10 Entropy');

    subplot(2,4,8);
    imagesc(grad);
    title('Thermal gradients');
    colormap jet

    pause(0.01);
    population=new_population;

    saveas(gcf,['Figure\Figure_',num2str(g,'%06.f')],'png');
    saveas(gcf,'Figure.png');

    if (rem(g,10)==0)||(g==1)
        disp('Saving current memory state...');
        save Current_state.mat
    end
    g=g+1;
    toc
end
disp('Shape convergence !')
