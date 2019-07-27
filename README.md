# icqm
Integer Convex Quadratic Minimization

    %% icqm Integer Convex Quadratic Minimizer                                                       
    % Find a solution to the following problem:                                                      
    %                                                                                                
    % minimize:      x'*mtx_M*x + 2*v_v'*x + d_cc                                                    
    % subject to:    x with integer components                                                       
    %                                                                                                
    % This function implements a translation of the Schnorr-Euchner algorithm                        
    % for ordinary integer least squares to the described integer convex                             
    % quadratic minimzation problem.                                                                 
    %                                                                                                
    % Taken from Algorithm 5, Figure 2 in:
    % GHASEMMEHDI AND AGRELL: FASTER RECURSIONS IN SPHERE DECODING.
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
                                            
