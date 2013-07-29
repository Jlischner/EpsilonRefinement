more off;
addpath("/auto/jlischner/Dropbox/RefinementCode/EpsRefinement");
input;
global tlc = 1/Nkx/2;
[partner_irr,fkull] = fullbz(syms,kirr);
Nkirr = size(kirr)(1);

ws = [];
tmpfreq = init_frequency;
ifreqcounter = 0;
freqstep = delta_frequency;
while ( tmpfreq <= frequency_high_cutoff);

  ws = [ws tmpfreq];
  ifreqcounter += 1;
  if(tmpfreq < frequency_low_cutoff);
    tmpfreq += delta_frequency;
  else;
    freqstep += delta_frequency_step;
    tmpfreq  += freqstep;
  endif;
endwhile;

alpha = sqrt(2*EFdirac/27.21/kappa); %# plasmon dispersion:alpha*sqrt(q)
global vF = 0.38859; %# corresponding to 0.85*10^8 cm/s
qtgc = alpha^2/4/vF^2; %# tangent criterion q in bohr
qcut_tgc = qtgc * 1.4; %# maximum q included in refinement
qpad = norm(Gmat(1,:))/Nkx; %# small padding
kF = EFdirac/27.21/vF; %# fermi wave vector in bohr
kcut = kF + qcut_tgc + qpad; %# distance from kdirac being searched for transitions

%# find relevant points in irrBZ
for ii = 2:Nkirr;

  printf("doing %d out of %d \n",ii,Nkirr);
  normq = norm(kirr(ii,:)*Gmat);
  if( normq < qcut_tgc);
    
    printf("found small qpoint! \n");
    printf("%d %f %f \n",ii,kirr(ii,1),kirr(ii,2));
    getxi(ii,4,kcut,S,ws);
    getxi(ii,5,kcut,S,ws);
    system(strcat("mkdir q",int2str(ii) ));
    system(strcat("mv avg4_file q",int2str(ii) ));
    system(strcat("mv avg5_file q",int2str(ii) ));
    system(strcat("mv indx4_file q",int2str(ii) ));
    system(strcat("mv indx5_file q",int2str(ii) ));
  endif;
  
endfor;

printf("using kcut=%f \n",kcut);
save epsilon_output kcut EFmid EFdirac
more on;