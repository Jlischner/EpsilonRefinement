%# file containing symmetry matrices:
global syms =  load("Derek72grid/syms_para"); 

%# irreducible kgrid (wrt full group):
global kirr = load("Derek72grid/k_irr"); 

%# energies in irr. wedge; assumes dirac bands are 4th and 5th col.:
global energies = load("Derek72grid/en_para"); 

%# linear dimension of BGW kgrid:
global Nkx = 72; 

%# undoped fermi level in eV:
global EFmid = -0.775;

%# fermi level wrt to dirac point (=EFmid) in eV:
global EFdirac = 1.0 ; 

%# refined kgrid (new kgrid is Nkx*S(1) ):
S = [4 4]; 

%# freq. grid for BGW epsilon calculation:
init_frequency = 0.0
delta_frequency = 0.005
delta_frequency_step = 100
frequency_low_cutoff = 1.0 
frequency_high_cutoff = 100.0

%# metric tensor in Gspace (atomic units):
global Gmat = [0.78013  -1.35122   0.00000;
               0.78013   1.35122   0.00000;
               0.00000   0.00000   0.33069];

%# approx. background dielectric constant, used to calculate refinement radius:
kappa = 1.4;