function condition_limites_ini=init_image(condition_limites_ini,nombre_pixels_conducteurs, k0, kp_k0);

%récupère les dimensions de l'image
[hauteur,largeur,profondeur]=size(condition_limites_ini);
k=0;

while k<nombre_pixels_conducteurs;
    h=ceil(rand*hauteur);
    l=ceil(rand*largeur);
    
    if condition_limites_ini(h,l)==k0; 
    condition_limites_ini(h,l)=k0*kp_k0; 
    k=k+1;
    end;
        
end;