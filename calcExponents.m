function lambda = calcExponents(pol, tau)
    pos_ind = real(pol)>0 & abs(imag(pol))<1e-8;
    neg_ind = real(pol)<0 & abs(imag(pol))<1e-8;

    [ccPol, ~] = sort(log(pol(~pos_ind&~neg_ind)));

    negPol = [real(log(pol(neg_ind))) + 1i*pi; real(log(pol(neg_ind))) - 1i*pi;];
    [sorted_negPol, ~] = sort(negPol);
    
    logpol = [log(pol(pos_ind)); ccPol; sorted_negPol];
    lambda = logpol/tau;
end