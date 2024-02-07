function boundary_conditions_ini=init_image(boundary_conditions_ini,number_conductive_cells, k0, kp_k0)

%récupère les dimensions de l'image
[height,width,~]=size(boundary_conditions_ini);
k=0;

while k<number_conductive_cells
    h=ceil(rand*height);
    l=ceil(rand*width);
    if boundary_conditions_ini(h,l)==k0 
    boundary_conditions_ini(h,l)=k0*kp_k0; 
    k=k+1;
    end
        
end