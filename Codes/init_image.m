%https://doi.org/10.1016/j.ijthermalsci.2016.05.015
%https://github.com/Raphael-Boichot/A-genetic-algorithm-for-topology-optimization-of-area-to-point-heat-conduction-problem
function boundary_conditions_ini=init_image(boundary_conditions_ini,number_conductive_cells, yet_conductive_pixels, k0, kp_k0)

%récupère les dimensions de l'image
[height,width,~]=size(boundary_conditions_ini);
k=0;

if number_conductive_cells>yet_conductive_pixels
    while k<(number_conductive_cells-yet_conductive_pixels)
        h=ceil(rand*height);
        l=ceil(rand*width);
        if boundary_conditions_ini(h,l)==k0
            boundary_conditions_ini(h,l)=k0*kp_k0;
            k=k+1;
        end
    end
end

if number_conductive_cells<yet_conductive_pixels
    while k<(yet_conductive_pixels-number_conductive_cells)
        h=ceil(rand*height);
        l=ceil(rand*width);
        if boundary_conditions_ini(h,l)==k0*kp_k0
            boundary_conditions_ini(h,l)=k0;
            k=k+1;
        end
    end
end