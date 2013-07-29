function out = bringBZ(k);

  %# check if k in [0,1);
  %# otherwise bring into BZ
  global tlc;

  for n = 1:3;

    if(k(n) < -tlc);
      k(n) += 1;
    endif;

    if( k(n) >= 1-tlc);
      k(n) -= 1;
    endif;
    
  endfor

  out = k;

endfunction;