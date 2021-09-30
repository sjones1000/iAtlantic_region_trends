function [specvol_CT_exact, alpha_CT_exact, beta_CT_exact] = gsw_specvol_alpha_beta_CT_exact(SA,CT,p)

% gsw_specvol_alpha_beta_CT_exact        specific volume, thermal expansion 
%                                          & saline contraction coefficient
%==========================================================================
% 
% USAGE:  
%  [specvol_CT_exact, alpha_CT_exact, beta_CT_exact] = ...
%                                  gsw_specvol_alpha_beta_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates specific volume, the appropiate thermal expansion coefficient
%  and the appropriate saline contraction coefficient of seawater from 
%  Absolute Salinity and Conservative Temperature.  
%
%  Note that this function uses the full Gibbs function.  There is an 
%  alternative to calling this function, namely 
%  gsw_specvol_alpha_beta(SA,CT,p), which uses the computationally
%  efficient 75-term expression for density in terms of SA, CT and p
%  (Roquet et al., 2015).  
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  specvol_CT_exact =  specific volume                             [ m/kg ]
%  alpha_CT_exact   =  thermal expansion coefficient                [ 1/K ]
%                      with respect to Conservative Temperature
%  beta_CT_exact    =  saline (i.e. haline) contraction            [ kg/g ]
%                      coefficient at constant Conservative Temperature
%
% AUTHOR: 
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05 (17th January, 2015)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix A.20 and appendix K of this TEOS-10 Manual. 
%
%  Roquet, F., G. Madec, T.J. McDougall and P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_specvol_alpha_beta_CT_exact:  Requires three inputs')
end %if

if ~(nargout == 3)
   error('gsw_specvol_alpha_beta_CT_exact:  Requires three outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_specvol_alpha_beta_CT_exact: SA and CT must have same dimensions')
end

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,ms), :);              % copy down each column.
elseif (ms == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,ns));               % copy across each row.
elseif (ns == mp) & (np == 1)          % p is a transposed row vector,
    p = p.';                              % transposed then
    p = p(ones(1,ms), :);                % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_specvol_alpha_beta_CT_exact: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    CT = CT.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

t = gsw_t_from_CT(SA,CT,p);
specvol_CT_exact = gsw_specvol_t_exact(SA,t,p);
alpha_CT_exact = gsw_alpha_wrt_CT_t_exact(SA,t,p);
beta_CT_exact = gsw_beta_const_CT_t_exact(SA,t,p);

if transposed
    specvol_CT_exact = specvol_CT_exact.';
    alpha_CT_exact = alpha_CT_exact.';
    beta_CT_exact = beta_CT_exact.';
end

end
