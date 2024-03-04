function out=gfun_cubicquintic(u,alpha3,beta3,alpha4,beta4)
  us = u.*u;
  cu = conj(u);
  ucub = us.*cu;

  out = (alpha3+1i*beta3)*ucub+(alpha4+1i*beta4)*(ucub.*(u.*cu));

end
