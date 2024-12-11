function [basis_reshaped,basis_single_reshaped,...
    basis_diff_reshaped,basis_diff_single_reshaped,...
    basis_diff2_reshaped,basis_diff2_single_reshaped]=...
    base_func(state,state_min,state_max,Smol_elem,mu_max,d,ind)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input:
%%% state: n_pts*n_state*n_node

%% Output:
%%% basis_reshaped: n_pts*n_coef*n_node
%%% basis_single_reshaped: n_pts*n_coef*n_node*n_state
%%% basis_diff_reshaped: n_pts*n_coef*n_node*n_state
%%% basis_diff_single_reshaped: n_pts*n_coef*n_node*n_state
%%% basis_diff2_reshaped: n_pts*n_coef*n_node*n_state
%%% basis_diff2_single_reshaped: n_pts*n_coef*n_node*n_state

%% overview:
% Step1. Reshape state variables
% Step2. Construct Smolyak Polynomial
% Step3. Reshape basis functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n_pts,n_state,n_node]=size(state);

%% Reshape state: n_pts*n_node*n_state => (n_pts*n_node)*n_state
state_reshaped=reshape(permute(state,[1,3,2]),n_pts*n_node,n_state);

%% Construct Smolyak polynomial
Smol_point=Rescale(state_reshaped,-1,1,state_min,state_max);

if nargout>=5
	[basis,basis_single,basis_diff,basis_diff_single,...
        basis_diff2,basis_diff2_single]=...
        Smolyak_Polynomial(Smol_point,d,mu_max,Smol_elem,ind);

    %%  Adjust rescale size for first derivatives
	% Let f and Psi denotes the rescale function and Smolyak function.
	% Note that g''(x)=Psi''(f(x))(f'(x))^2.
	% f'(x) corresponds to the following rescale_rate.

	Rescale_rate=(state_max-state_min)/2; %1*n_state
	basis_diff2=basis_diff2./reshape(Rescale_rate.^2,1,1,n_state);
	basis_diff2_single=basis_diff2_single./reshape(Rescale_rate.^2,1,1,n_state);
	basis_diff=basis_diff./reshape(Rescale_rate,1,1,n_state);
	basis_diff_single=basis_diff_single./reshape(Rescale_rate,1,1,n_state);

elseif nargout>=3
	[basis,basis_single,basis_diff,basis_diff_single]=...
        Smolyak_Polynomial(Smol_point,d,mu_max,Smol_elem,ind);

    %%  Adjust rescale size for first derivatives
	% Let f and Psi denotes the rescale function and Smolyak function.
	% Note that g'(x)=Psi'(f(x))(f'(x)).
	% f'(x) corresponds to the following rescale_rate.

	Rescale_rate=(state_max-state_min)/2; %1*n_state
	basis_diff=basis_diff./reshape(Rescale_rate,1,1,n_state);
	basis_diff_single=basis_diff_single./reshape(Rescale_rate,1,1,n_state);

else
	[basis,basis_single]=Smolyak_Polynomial(Smol_point,d,mu_max,Smol_elem,ind);
end

n_coef=size(basis,2);


%% Reshape basis: (n_pts*n_node)*n_coef =>  n_pts*n_coef*n_node
basis_reshaped=permute(reshape(basis,n_pts,n_node,n_coef),[1,3,2]);

%% Reshape basis_single: (n_pts*n_node)*n_coef*n_state =>  n_pts*n_coef*n_node*n_state
basis_single_reshaped=permute(reshape(basis_single,n_pts,n_node,n_coef,n_state),[1,3,2,4]);

if nargout>=3
	%% Reshape basis_diff: n_pts*n_node*n_coef*n_state => n_pts*n_coef*n_node*n_state
	basis_diff_reshaped=permute(reshape(basis_diff,n_pts,n_node,n_coef,n_state),[1,3,2,4]);

	%% Reshape basis_diff_single: n_pts*n_node*n_coef*n_state => n_pts*n_coef*n_node*n_state
	basis_diff_single_reshaped=permute(reshape(basis_diff_single,n_pts,n_node,n_coef,n_state),[1,3,2,4]);
end


if nargout>=5
	%% Reshape basis_diff: n_pts*n_node*n_coef*n_state => n_pts*n_coef*n_node*n_state
	basis_diff2_reshaped=permute(reshape(basis_diff2,n_pts,n_node,n_coef,n_state),[1,3,2,4]);

	%% Reshape basis_diff_single: n_pts*n_node*n_coef*n_state => n_pts*n_coef*n_node*n_state
	basis_diff2_single_reshaped=permute(reshape(basis_diff2_single,n_pts,n_node,n_coef,n_state),[1,3,2,4]);
end

return
