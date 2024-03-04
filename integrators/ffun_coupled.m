function out=ffun_coupled(u,M,alpha3,beta3,alpha4,beta4,alpha5,ft,ift)

  mod2{1} = conj(u{1}).*u{1};
  mod4{1} = mod2{1}.*mod2{1};

  mod2{2} = conj(u{2}).*u{2};
  mod4{2} = mod2{2}.*mod2{2};

  out{1} = ift(M{1}.*ft(u{1}))+...
  ((alpha3+1i*beta3)*mod2{1}+(alpha4+1i*beta4)*mod4{1}+alpha5*mod2{2}).*u{1};
  out{2} = ift(M{2}.*ft(u{2}))+...
  ((alpha3+1i*beta3)*mod2{2}+(alpha4+1i*beta4)*mod4{2}+alpha5*mod2{1}).*u{2};
