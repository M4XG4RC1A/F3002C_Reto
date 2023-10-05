JSA = 0;
for j = 1: lenght(omgp)
    KP2 = k_lambda(2.*pi.*3e14./(OMS+OMI-omgp(j)),coeff);
    DK = KP1+KP2-KS-KI-NL;

    JSA = JSA+ dwp.*(exp(-((omgp(j)-omgp0).^2)./(sigma.^2))
    .*exp(-(OMS+OMI-omgp(j)-omgp0).^2./(sigma.^2))).*sinc( L.*(KP1(j)+KP2-KS-KI-NL)./2)
    .*exp(1i.*L.*(KP1(j)+KP2-KS-KI-NL)./2)
    