function out = getxi(nq,bandv,kcut,S,ws);

  global Gmat;
  global Nkx;
  global vF;
  global EFdirac;
  global EFmid;
  global energies;
  global syms;
  global kirr;
  
  %# run fullbz: gives mapping between wedge and full zone
  [partner_irr,ks]=fullbz(syms,kirr);
  Nk = size(ks)(1);

  %# derived parameters:
  %#--------------------
  EF = EFdirac + EFmid;
  ryd = 27.211/2;
  q = kirr(nq,:); 
  printf("getxi: doing q=(%f %f) \n",q(1),q(2));
  d = 1/Nkx;
  tlc =0.001;
  lenS = prod(S);
  bandc = 5; %# which empty state
  if(bandv == 5);
    brd_fac = 1;
  else;
    brd_fac = 4;
  endif;
  eta = brd_fac*norm(Gmat(1,:))/Nkx/S(1)*vF*27.21; %# in eV

  ms=[0:prod(S)-1];
  ms=ms.';
  m1=rem(ms,S(1));
  m2=rem(floor(ms/S(1)),S(2));
  t = m1/S(1);
  u = m2/S(2);
  
  indxv = [];
  avgm  = [];
  
  %# which kpoint am I doing - specified by indx1
  for indx1 = 1:Nk;
    k1 = ks(indx1,:);
    dk1 = k1 - [2/3 1/3 0];
    dk2 = k1 - [1/3 2/3 0];
    dis1 = norm(dk1*Gmat);
    dis2 = norm(dk2*Gmat);

    if( dis1 < kcut || dis2 < kcut);
      
      #printf("found canidate k1=%f %f %d\n",k1(1),k1(2),indx1);
      k2 = bringBZ( k1 + [d 0 0] );
      k3 = bringBZ( k1 + [d d 0] );
      k4 = bringBZ( k1 + [0 d 0] );
      
      e1 = energies(partner_irr(indx1),bandc);
      [dis2,ind2] = min( sqrt( sum( (ks-ones(Nk,1)*k2).^2,2) ) ); 
      if( dis2(1) > tlc);
	printf("could not find k2 in list of kpoints! \n");
      else;
	indx2 = ind2(1);
	e2 = energies(partner_irr(indx2),bandc);
      endif;
      
      [dis3,ind3] = min( sqrt( sum( (ks-ones(Nk,1)*k3).^2,2) ) );
      if( dis3(1) > tlc);
	printf("could not find k3 in list of kpoints! \n");
      else;
	indx3 = ind3(1);
	e3 = energies(partner_irr(indx3),bandc);
      endif;
      
      [dis4,ind4] = min( sqrt( sum( (ks-ones(Nk,1)*k4).^2,2) ) );
      if( dis4(1) > tlc);
	printf("could not find k4 in list of kpoints! \n");
      else;
	indx4 = ind4(1);
	e4 = energies(partner_irr(indx4),bandc);
      endif;
      
      kq1 = bringBZ( k1 + q );
      kq2 = bringBZ( k2 + q );
      kq3 = bringBZ( k3 + q );      
      kq4 = bringBZ( k4 + q );
      
      [dis1,ind1] = min( sqrt( sum( (ks-ones(Nk,1)*kq1).^2,2) ) ); 
      if( dis1(1) > tlc);
	printf("could not find kq1=%f %f in list of kpoints! \n",kq1(1),kq1(2));
      else;
	indx1q = ind1(1);
	e1q = energies(partner_irr(indx1q),bandv);
      endif;

      [dis2,ind2] = min( sqrt( sum( (ks-ones(Nk,1)*kq2).^2,2) ) ); 
      if( dis2(1) > tlc);
	printf("could not find kq2=%f %f in list of kpoints! \n",kq2(1),kg2(2));
      else;
	indx2q = ind2(1);
	e2q = energies(partner_irr(indx2q),bandv);
      endif;
      
      [dis3,ind3] = min( sqrt( sum( (ks-ones(Nk,1)*kq3).^2,2) ) );
      if( dis3(1) > tlc);
	printf("could not find kq3=%f %f in list of kpoints! \n",kq3(1),kq3(2));
      else;
	indx3q = ind3(1);
	e3q = energies(partner_irr(indx3q),bandv);
      endif;
      
      [dis4,ind4] = min( sqrt( sum( (ks-ones(Nk,1)*kq4).^2,2) ) );
      if( dis4(1) > tlc);
	printf("could not find kq4=%f %f in list of kpoints! \n",kq4(1),kq4(2));
      else;
	indx4q = ind4(1);
	e4q = energies(partner_irr(indx4q),bandv);
      endif;
      
      
      #%Eks is empty, Ekqs is occ
      Eks = (1-t).*(1-u)*e1  + t.*(1-u)*e2  + t.*u*e3  + (1-t).*u*e4;
      Ekqs= (1-t).*(1-u)*e1q + t.*(1-u)*e2q + t.*u*e3q + (1-t).*u*e4q;
      dEs = Ekqs - Eks;
      occs = (Eks > EF).*(Ekqs < EF);

      if( sum(occs) > tlc);
	
	#printf("getxi: found transition! \n")
	avgv = zeros(size(ws));
	for n = 1:length(ws);      
	  avgv(n)  = sum( occs./(dEs - ws(n) - I*eta) );
	  avgv(n) += sum( occs./(dEs + ws(n) + I*eta) );
	endfor;
	avgv *= -ryd*0.5/lenS;
	indxv = [indxv; indx1];
	avgm  = [avgm avgv];
      endif;
      
    endif; %# kcut-if
  endfor; %# indx1-loop
  
  %# write output:
  %#--------------
  Ntrans = length(indxv);
  avg_r  = reshape(avgm,Ntrans*length(ws),1);
  
  avg_r2 = [real(avg_r) imag(avg_r)];
  f1name = strcat("avg",int2str(bandv),"_file");
  fid = fopen(f1name,'w');
  fprintf(fid,"%20.12f %20.12f \n",avg_r2');
  fclose(fid);
  
  %# first entry in indxv is number of transitions
  %# other entries are their indices
  indxv = [Ntrans; indxv];
  f2name = strcat("indx",int2str(bandv),"_file");
  fid = fopen(f2name,'w');
  fprintf(fid,"%d \n",indxv');
  fclose(fid);

endfunction;