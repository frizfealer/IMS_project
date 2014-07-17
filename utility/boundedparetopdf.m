function [ val ] = boundedparetopdf( L, H, alpha, x )
%boundedparetopdf from wikipedia 
%http://en.wikipedia.org/wiki/Pareto_distribution#Generalized_Pareto_distributions
val = (alpha * L^alpha).*x.^(-alpha-1) / (1-(L/H)^alpha);

end

