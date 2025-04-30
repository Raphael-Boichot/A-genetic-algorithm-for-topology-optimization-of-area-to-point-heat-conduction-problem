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
pixels=imread('Best_topology.png');
%--------------------------------------------------------------------------

[height, width, deepness]=size(pixels);
Boundary_conditions=zeros(height, width);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Image translation into boundary limits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:1:height
    for l = 1:1:width
        red = pixels(k,l,1);
        green = pixels(k,l,2);
        blue = pixels(k,l,3);
        if (red == 127) && (green == 127) && (blue == 127)
            Boundary_conditions(k,l)=-2;
        end
        if (red == 0) && (green == 0) && (blue == 255)
            Boundary_conditions(k,l)=-3;
        end
        if (red == 255) && (green == 255) && (blue == 255)
            Boundary_conditions(k,l)=k0;
        end
        if (red == 0) && (green == 0) && (blue == 0)
            Boundary_conditions(k,l)=kp_k0;
        end
    end
end

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
[distance,sum_of_entropy, entropy, border_variance,variance, mean_temperature,maximal_temperature,temp,grad,variance_grad]=finite_temp_direct_sparse(k0*kp_k0,k0,Heat_sink_temperature,delta_x,p,Boundary_conditions);

map_to_colorize=log10(entropy+1e-3);
rgbImage = grayscale_to_colormap(map_to_colorize, jet(256));
imshow(rgbImage);
imwrite(rgbImage,'Entropy.png');
drawnow
pause(1)

map_to_colorize=temp;
rgbImage = grayscale_to_colormap(map_to_colorize, jet(256));
imshow(rgbImage);
imwrite(rgbImage,'Temperature.png');
drawnow
pause(1)

map_to_colorize=grad;
rgbImage = grayscale_to_colormap(map_to_colorize, jet(256));
imshow(rgbImage);
imwrite(rgbImage,'Thermal_gradients.png');
drawnow
pause(1)