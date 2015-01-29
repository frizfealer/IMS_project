
% harmonic.m by David Terr, Raytheon, 9-10-04

% Given a complex number z, estimate H_z, the zth Harmonic number. 
% For z a positive integer, this is the sum 1 + 1/2 + 1/3 + ... + 1/z.

% Warning: This may not be accurate if |z| is small.

% Reference: MathWorld, Harmonic Number entry

function h = user_harmonic(z)

if ( z == 1 ) 
    h = 1;
else
%     h = log(z) + 0.577215664901532 + 1/(2*z) - 1/(12*z^2) + 1/(120*z^4) - 1/(252*z^6);
    h = psi(z+1) + 0.577215664901532;
end