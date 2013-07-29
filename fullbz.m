function [partner_irr,ks] = fullbz(syms,kirr);
  
  global tlc;
  Nkirr = size(kirr)(1);
  Nsyms = size(syms)(1);
 
  lenk = 0;
  ks   = [];
  partner_irr = [];
  
  for ii = 1:Nkirr;
    k = kirr(ii,:);
    for jj = 1:Nsyms;
      
      S = reshape(syms(jj,:),3,3);
      kt = transpose( S*k' );
      kt(kt < tlc) += 1;
      
      if(lenk == 0); 
	ks = [ks; kt];
	partner_irr = [partner_irr; ii];
	lenk += 1;
      else;
	diff = min(sqrt(sum( (ks - ones(lenk,1)*kt).^2,2)));
	if( diff > tlc)
	  ks = [ks; kt];
	  partner_irr = [partner_irr; ii];
	  lenk += 1;
	endif;
      endif;
      
    endfor;
  endfor;
  
  for ii = 1:size(ks)(1);
    ks(ii,:) = bringBZ(ks(ii,:));
  endfor;

endfunction;