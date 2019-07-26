function [v_opt, d_opt] = ICQM(mtx_M, v_v, d_cc)
%% ICQM Integer Convex Quadratic Minimizer 
% Find a solution to the following problem:
% 
% minimize:      x'*mtx_M*x + 2*v_v'*x + d_cc
% subject to:    x with integer components
%
% This function implements a translation of the Schnorr-Euchner algorithm 
% for ordinary integer least squares to the described integer convex 
% quadratic minimzation problem. 
% 
% Taken from Algorithm 2.1.1 `SEARCH' in:
% Borno, Mazen Al. "Reduction in solving some integer least squares problems." arXiv preprint arXiv:1101.0382 (2011).
%
% Inputs:
%    mtx_M = K*K real positive-semidefinite matrix
%    v_v = K*1 real vector
%    d_cc = real scalar
%
% Outputs: 
%    v_opt = K*1 integer vector, optimal point 
%    d_opt = real scalar, optimal value
%

%% Setup
	n_dim = size(mtx_M,1); 

	% Stated problem is equivalent to instead minimizing:
	% ||mtx_R*v_x - v_y||^2 + d_cc2
	mtx_R = cholcov(mtx_M);
	v_y = -pinv(mtx_R')*v_v;
	%d_cc2 = d_cc - v_y'*v_y;
	v_diag_R = diag(mtx_R);

	% Auxillary functions
	% Equation (2.4)
	fn_c = @(v_z, n_k) ( v_y(n_k) - ...
		sum(mtx_R(n_k,(n_k+1):end)*v_z((n_k+1):end,1)) )/mtx_R(n_k,n_k);
	% beta^2 formula
	fn_beta2_partial = @(v_z, v_c, n_k) ...
		sum( ( v_diag_R((n_k+1):end,1).* ...
		(v_z((n_k+1):end,1)-v_c((n_k+1):end,1)) ).^2 );
	
	% Globals
	v_z = zeros(n_dim,1);
	v_hatz = zeros(n_dim,1);
	v_c = zeros(n_dim,1);
	v_Delta = zeros(n_dim,1);

%% Schnorr-Euchner
	% Initialization
	n_k = n_dim;
	d_beta2 = inf;

	STEP = 2;
	while(STEP > 0)
	if(STEP == 2)
		v_c(n_k) = fn_c(v_z, n_k);
		v_z(n_k) = round(v_c(n_k)); 
		v_Delta(n_k) = sign(v_c(n_k)-v_z(n_k));

		STEP = 3;
		continue; 
	elseif(STEP == 3) % Main step
		d_A = mtx_R(n_k,n_k)^2*(v_z(n_k)-v_c(n_k))^2;
		d_B = d_beta2 - fn_beta2_partial(v_z, v_c, n_k);
		if(d_A > d_B)
			STEP = 4;
			continue;
		elseif( n_k > 1 )
			n_k = n_k-1;

			STEP = 2; 
			continue;
		else % case n_k = 1;
			STEP = 5;
			continue;
		end
	elseif(STEP == 4) % Invalid point
		if(n_k == n_dim)
			STEP = -1;
			continue;
		else
			n_k = n_k+1;

			STEP = 6;
			continue;
		end
	elseif(STEP == 5) % Found valid point
		v_hatz = v_z;
		d_beta2 = fn_beta2_partial(v_hatz, v_c, n_dim);
		n_k = n_k+1;

		STEP = 6;
		continue;
	elseif(STEP == 6) % Enumeration at level n_k
		v_z(n_k) = v_z(n_k) + v_Delta(n_k);
		v_Delta(n_k) = -v_Delta(n_k) - sign(v_Delta(n_k));

		STEP = 3; 
		continue;
	end
	end

	v_opt = v_hatz;
	d_opt = v_opt'*mtx_M*v_opt + 2*v_v'*v_opt + d_cc;
	return;
end
