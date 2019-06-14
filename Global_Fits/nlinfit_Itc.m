function [ p, err ] = nlinfit_Itc( param_Tot, fixed, TitrPar, I0, x, T_rm, dQm_exp )
% modified version of nlinfit to allow for fixed parameters
% following, but does not use varargin
% learn-one-thing-a-day.blogspot.com/2012/11/
% how-to-hold-parameters-in-nlinfit-matlab
% NB in order to deal with imaginary numbers during fitting 
% procedure, i have taken real parts of solved parameters, 
% residuals, and jacobians from nlinfit for purposes of 
% estimating errors.  imag numbers can easily arise due to 
% a square root in function fcn

% copy ALL (both free and fixed) parameters in beta0 into p
p = param_Tot;

% assign zero uncertainty to fixed paraemters
err = zeros( size(param_Tot) );

% extract free parameters, i.e., those NOT (~) fixed
param_Free = param_Tot( ~fixed );

% NB only free parameters are passed to nlinfit
[p1, r1, JJ1] = nlinfit( x, dQm_exp, @nested_fcn, param_Free );

% get error bars for free parameters
%delta = nlparci( p1Real, r1Real, JJ1Real ); % wgn-28Sept2013
delta = nlparci( p1, r1, JJ1 );

err1 = 0.5 * ( delta(:,2) - delta(:,1) ); 

% update free parameters, i.e., those NOT (~) fixed
p(~fixed) = p1; 

% update according error bars
err(~fixed) = err1;

	% Nested function using only parameters as inputs
	% while inheriting addl parameters from outer function

	function yy = nested_fcn( param_Passed, x ) 
        
        param_Nest(  fixed ) = param_Tot( fixed );
        param_Nest( ~fixed ) = param_Passed;
        
        yy = Fcn_Calc_dQm_N( param_Nest, TitrPar, I0, x, T_rm );

    end

end

