%*****************Automate cellulaire*********************INPG/BOICHOT/2008
function [distance,sum_of_entropy, entropy, border_variance,variance, mean_temperature,maximal_temperature,temp,grad,variance_grad]=finite_temp_direct_sparse(kp_k0,k0,heat_sink_temperature,pas_x,p_vol,boundary_conditions)

[height,width,~]=size(boundary_conditions);
temp=ones(height,width).*heat_sink_temperature;
conductivity_table=zeros(height,width,5);
temporary_boundaries=boundary_conditions;

for k=1:1:(height)
    for l=1:1:(width)
        if temporary_boundaries(k,l)==-2; temporary_boundaries(k,l)=1e-9;end
        if temporary_boundaries(k,l)==-3; temporary_boundaries(k,l)=kp_k0;end
    end
end

for k=2:1:(height-1)
    for l=2:1:(width-1)
        %thermal conductance table
        conductivity_table(k, l, 2) = pas_x / ((pas_x / 2) / temporary_boundaries(k-1,l) + (pas_x / 2) / temporary_boundaries(k, l));
        conductivity_table(k, l, 3) = pas_x / ((pas_x / 2) / temporary_boundaries(k,l+1) + (pas_x / 2) / temporary_boundaries(k, l));
        conductivity_table(k, l, 4) = pas_x / ((pas_x / 2) / temporary_boundaries(k+1,l) + (pas_x / 2) / temporary_boundaries(k, l));
        conductivity_table(k, l, 5) = pas_x / ((pas_x / 2) / temporary_boundaries(k,l-1) + (pas_x / 2) / temporary_boundaries(k, l));
    end
end

address = zeros(height, width);
k=1;
%address matrix
for i=1:1:height
    for j=1:1:width
        address(i,j)=k;
        k=k+1;
    end
end

%Initialization
B=zeros((height)*(width),1);
k=0;

%direct method for solving heat equation by finite difference
for i=1:1:(height)
    for j=1:1:(width)
        ind_1=address(i,j);
        if not (boundary_conditions(i,j)==-2||boundary_conditions(i,j)==-3)
            ind_2=address(i-1,j);
            ind_3=address(i,j+1);
            ind_4=address(i+1,j);
            ind_5=address(i,j-1);
            somme_cond = conductivity_table(i, j, 2)+conductivity_table(i, j, 3)+conductivity_table(i, j, 4)+conductivity_table(i, j, 5);
        end
        fin=0;

        %Central cell (1)
        if boundary_conditions(i,j)==-2||boundary_conditions(i,j)==-3
            B(ind_1)=temp(i,j);
            k=k+1;
            row(k)=ind_1;
            column(k)=ind_1;
            value(k)=1;
            fin=1;
        else
            k=k+1;
            row(k)=ind_1;
            column(k)=ind_1;
            value(k)=-1*somme_cond;
            if boundary_conditions(i,j)==k0; B(ind_1)=B(ind_1)-p_vol*pas_x^2;end
        end
        if not (fin==1)
            %Cell North (2)
            if boundary_conditions(i-1,j)==-2||boundary_conditions(i-1,j)==-3
                B(ind_1)=B(ind_1)-temp(i-1,j)*conductivity_table(i, j, 2);
            else
                k=k+1;
                row(k)=ind_1;
                column(k)=ind_2;
                value(k)=1*conductivity_table(i, j, 2);
            end
            %Cell East (3)
            if boundary_conditions(i,j+1)==-2||boundary_conditions(i,j+1)==-3
                B(ind_1)=B(ind_1)-temp(i,j+1)*conductivity_table(i, j, 3);
            else
                k=k+1;
                row(k)=ind_1;
                column(k)=ind_3;
                value(k)=1*conductivity_table(i, j, 3);
            end
            %Cell North (4)
            if boundary_conditions(i+1,j)==-2||boundary_conditions(i+1,j)==-3
                B(ind_1)=B(ind_1)-temp(i+1,j)*conductivity_table(i, j, 4);
            else
                k=k+1;
                row(k)=ind_1;
                column(k)=ind_4;
                value(k)=1*conductivity_table(i, j, 4);
            end
            %Cell West (5)
            if boundary_conditions(i,j-1)==-2||boundary_conditions(i,j-1)==-3
                B(ind_1)=B(ind_1)-temp(i,j-1)*conductivity_table(i, j, 5);
            else
                k=k+1;
                row(k)=ind_1;
                column(k)=ind_5;
                value(k)=1*conductivity_table(i, j, 5);
            end
        end
    end
end

A=sparse(row,column,value);
T=A\B;

%Temperature map formating
for i=1:1:height
    for j=1:1:width
        temp(i,j)=T(address(i,j));
    end
end

maximal_temperature = max(max(temp));
variance=std(std(temp))^2;
[gradX, gradY]=gradient(temp(2:end-1,2:end-1));
grad=(gradX.^2+gradY.^2).^0.5;
variance_grad=std(std(grad))^2;
mean_temperature=mean(mean(temp));
border_variance=std([temp(2:end-1,2)' temp(2,2:end-1)]);

%Entropy table calculation
sum_of_entropy=0;
entropy=zeros(height,width);

for k=1:1:height
    for l=1:1:width
        entropy(k,l)=0;
        if not(boundary_conditions(k,l)==-2||boundary_conditions(k,l)==-3)
            flux1=conductivity_table(k,l,2)*(temp(k-1,l)-temp(k,l));
            flux2=conductivity_table(k,l,3)*(temp(k,l+1)-temp(k,l));
            flux3=conductivity_table(k,l,4)*(temp(k+1,l)-temp(k,l));
            flux4=conductivity_table(k,l,5)*(temp(k,l-1)-temp(k,l));
            flux5=0;
            if boundary_conditions(k,l)==k0; flux5=p_vol*pas_x*pas_x; end
            if not(boundary_conditions(k-1,l)==-2)
                entropy(k,l)=entropy(k,l)+(abs(flux1/temp(k,l)-abs(flux1/temp(k-1,l))))*sign(temp(k,l)-temp(k-1,l));
            end
            if not(boundary_conditions(k,l+1)==-2)
                entropy(k,l)=entropy(k,l)+(abs(flux2/temp(k,l)-abs(flux2/temp(k,l+1))))*sign(temp(k,l)-temp(k,l+1));
            end
            if not(boundary_conditions(k+1,l)==-2)
                entropy(k,l)=entropy(k,l)+(abs(flux3/temp(k,l)-abs(flux3/temp(k+1,l))))*sign(temp(k,l)-temp(k+1,l));
            end
            if not(boundary_conditions(k,l-1)==-2)
                entropy(k,l)=entropy(k,l)+(abs(flux4/temp(k,l)-abs(flux4/temp(k,l-1))))*sign(temp(k,l)-temp(k,l-1));
            end
            entropy(k,l)=entropy(k,l)+flux5/temp(k,l);
            sum_of_entropy=sum_of_entropy+entropy(k,l);
        end
    end
end

[row,col] = find(temp==maximal_temperature);
row_ref=height;
col_ref=width;
distance=(row-row_ref)^2+(col-col_ref)^2;




