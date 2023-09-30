function [lams,lami] = lsli(lamp,NL,as,coeff)

options = optimset('Display', 'on', 'TolFin', 1e-20);
optnew = optimset(options,'TolX', 1e-20);

lams = fzero(@(lams)DK(lamp,lams,coeff,NL),as,optnew);
lami = (2*pi*3e14)./((2*pi*3e14)./lamp + (2*pi*3e14)./lamp - (2*pi*3e14)./lams);

end




