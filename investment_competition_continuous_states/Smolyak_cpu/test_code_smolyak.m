clear
tic
% 1. Smolyak anisotropic method for 2 dimensions;
% -----------------------------------------------
d=10;
vector_mus_dimensions = 2*ones(1,d); % Introduce the level of approximation in every dimension from 1 to 10; see Section 4 of JMMV (2014)

d=length(vector_mus_dimensions); % Number of dimensions
mu_max  = max(vector_mus_dimensions); % Compute the maximum level of approximation across all dimensions

Smolyak_elem_iso = Smolyak_Elem_Isotrop(d,mu_max);
    % Construct the matrix of indices of multidimesional Smolyak elements (grid points and polynomial basis functions) that satisfy the usual
    % isotropic Smolyak rule for the approximation level equal to "mu_max"
Smol_elem_ani = Smolyak_Elem_Anisotrop(Smolyak_elem_iso,vector_mus_dimensions);
    % Select from the matrix of isotropic indices "Smol elem_iso" a subset of indices that correspond to the given anisotropic "vector_mus_dimensions"
Smol_grid_ani = Smolyak_Grid(d,mu_max,Smol_elem_ani);
    % Construct the Smolyak grid for the given subindices of anisotropic Smolyak elements "Smol_elem_ani"
Smol_polynom = Smolyak_Polynomial(Smol_grid_ani,d,mu_max,Smol_elem_ani);
    % Matrix of the polynomial basis functions evaluated in the grid points
toc
if d==2
    figure(6)
    scatter(Smol_grid_ani(:,1),Smol_grid_ani(:,2),'filled'),title('Smolyak grid') % Enable it to draw the Smolyak grid
end
% 2. Interpolation of a function y=2*x1 .*exp(-4*x1.^2-16*x2.^2);
% ---------------------------------------------------------------
y_Smolyak = 2*Smol_grid_ani(:,1) .* exp(-4*Smol_grid_ani(:,1).^2 - 16*Smol_grid_ani(:,2).^2);
    % Evaluate the function on the Smolyak grid

b = Smol_polynom\y_Smolyak;
    % Compute the coefficients of Smolyak interpolating polynomial

size(b)
    % Take a look to see how many coefficients this involves
toc

% 3. Compare the true and interpolated functions on a dense grid
% ---------------------------------------------------------------
[x1,x2] = meshgrid(-1:0.05:1, -1:0.1:1);
    % Create a uniformly spaced grid of 41 points on [-1,1]x[-1,1]

y_true = 2*x1.* exp(-4*x1.^2 - 16*x2.^2);
    % Evaluate the true function on the grid

figure(7)
subplot(1,2,1), surf(x1,x2,y_true),title('True function')
    % Plot the true function

for j=1:size(x1(1,:),2)
    y_fitted(:,j) = Smolyak_Polynomial([x1(:,j),x2(:,j)],d,mu_max,Smol_elem_ani)*b;
    % Evaluate Smolyak interpolating polynomial on the grid
end

subplot(1,2,2), surf(x1,x2,y_fitted),title('Smolyak interpolation')
    % Plot Smolyak interpolation

%
% Note: the commands being called here to do the Smolyak anisotropic grids
% based on chebyshev polynomials are not limited to the two dimensions that
% we use here. They work for any number 'd' of dimensions.
%