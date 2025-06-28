#!/bin/bash
#SBATCH -N 1
#SBATCH -n 48                
#SBATCH -p normal
#SBATCH -J NbMoCO2




kp=10
ecut=80
natom=5
mat=NbMoCO2
ntyp=4
nbnd=30
pseudo=/home/pseudo/sg15



cat >scf.in <<EOF
&control
 calculation  ='scf',
 prefix       ='$mat',
 pseudo_dir   = "$pseudo"
 outdir       = './out',
 forc_conv_thr =  1.00000e-04
verbosity='high'
/
&system
    ibrav       = 0
    degauss     =  1.00000e-02
    ecutwfc     = $ecut
    !ecutrho     = $((12*ecut))
    nat         = $natom
    ntyp        = $ntyp
    nbnd        = $nbnd
    occupations = "smearing"
    smearing    = "gaussian"
    input_dft   = "PBE"
    nosym       = .TRUE.
    vdw_corr    = "DFT-D3"
/
&electrons
    conv_thr         =  1.00000e-10
    diago_david_ndim = 24
    diagonalization  = "david"
    electron_maxstep = 200
    mixing_beta      =  7.00000e-01
    mixing_mode      = "plain"
    mixing_ndim      = 24
    startingpot      = "atomic"
    startingwfc      = "atomic+random"

/
 

ATOMIC_SPECIES
Nb     92.90638  Nb_ONCV_PBE_sr.upf
Mo     95.94000  Mo_ONCV_PBE_sr.upf
C      12.01070  C_ONCV_PBE_sr.upf
O      15.99940  O_ONCV_PBE_sr.upf

CELL_PARAMETERS (angstrom)
   3.084938115   0.000000000   0.000000000
  -1.542469057   2.671634359   0.000000000
   0.000000000   0.000000000  25.000000000

ATOMIC_POSITIONS (crystal)
Nb            0.6666666666        0.3333333333        0.5487560083
Mo            0.3333333333        0.6666666666        0.4528740690
C            -0.0000000000       -0.0000000000        0.4996423876
O             0.3333333333        0.6666666666        0.5926987669
O             0.6666666666        0.3333333333        0.4083327682



K_POINTS {automatic}
$((2*kp)) $((2*kp)) 1 0 0 0

EOF
cat >nscf.in <<EOF
&control
 calculation  ='nscf',
 prefix       ='$mat',
 pseudo_dir   = "$pseudo"
 outdir       = './out',
 forc_conv_thr =  1.00000e-04
verbosity='high'
/
&system
    ibrav       = 0
    degauss     = 1.00000e-02
    ecutwfc     = $ecut
    !ecutrho     = $((12*ecut))
    nat         = $natom
    ntyp        = $ntyp
    nbnd        = $nbnd
    occupations = "fixed"
    smearing    = "gaussian"
    input_dft   = "PBE"
    nosym       = .TRUE.
    vdw_corr    = "DFT-D3"
/
&electrons
    conv_thr         =  1.00000e-10
    diago_david_ndim = 24
    diagonalization  = "david"
    electron_maxstep = 200
    mixing_beta      =  7.00000e-01
    mixing_mode      = "plain"
    mixing_ndim      = 24
    startingpot      = "atomic"
    startingwfc      = "atomic+random"

/
 

ATOMIC_SPECIES
Nb     92.90638  Nb_ONCV_PBE_sr.upf
Mo     95.94000  Mo_ONCV_PBE_sr.upf
C      12.01070  C_ONCV_PBE_sr.upf
O      15.99940  O_ONCV_PBE_sr.upf

CELL_PARAMETERS (angstrom)
   3.084938115   0.000000000   0.000000000
  -1.542469057   2.671634359   0.000000000
   0.000000000   0.000000000  25.000000000

ATOMIC_POSITIONS (crystal)
Nb            0.6666666666        0.3333333333        0.5487560083
Mo            0.3333333333        0.6666666666        0.4528740690
C            -0.0000000000       -0.0000000000        0.4996423876
O             0.3333333333        0.6666666666        0.5926987669
O             0.6666666666        0.3333333333        0.4083327682

K_POINTS {automatic}
$((3*kp)) $((3*kp)) 1 0 0 0
EOF



cat >bands.in <<EOF
&control
 calculation   ='bands',
 prefix        ='$mat',
 pseudo_dir    = "$pseudo"
 outdir        = './out',
 forc_conv_thr =  1.00000e-04
verbosity='high'
/
&system
 ibrav       = 0
 degauss     = 1.00000e-02
 nat         = $natom,
 ntyp        = $ntyp,
 ecutwfc     = $ecut,
 !ecutrho     = $((12*ecut)),
 nbnd        = $nbnd,
 occupations = "fixed"
 smearing    = "gaussian"
 nosym       = .TRUE.
    vdw_corr    = "DFT-D3"
/
&electrons
    conv_thr         =  1.00000e-10
    diago_david_ndim = 24
    diagonalization  = "david"
    electron_maxstep = 200
    mixing_beta      =  7.00000e-01
    mixing_mode      = "plain"
    mixing_ndim      = 24
    startingpot      = "atomic"
    startingwfc      = "atomic+random"
/

ATOMIC_SPECIES
Nb     92.90638  Nb_ONCV_PBE_sr.upf
Mo     95.94000  Mo_ONCV_PBE_sr.upf
C      12.01070  C_ONCV_PBE_sr.upf
O      15.99940  O_ONCV_PBE_sr.upf

CELL_PARAMETERS (angstrom)
   3.084938115   0.000000000   0.000000000
  -1.542469057   2.671634359   0.000000000
   0.000000000   0.000000000  25.000000000

ATOMIC_POSITIONS (crystal)
Nb            0.6666666666        0.3333333333        0.5487560083
Mo            0.3333333333        0.6666666666        0.4528740690
C            -0.0000000000       -0.0000000000        0.4996423876
O             0.3333333333        0.6666666666        0.5926987669
O             0.6666666666        0.3333333333        0.4083327682


K_POINTS {crystal_b}
5
   0.3333333333   0.3333333333   0.0000000000    100 !K 
   0.0000000000   0.0000000000   0.0000000000    100 !GAMMA                       
   0.5000000000   0.0000000000   0.0000000000    100 !M              
   0.3333333333   0.3333333333   0.0000000000    100 !K                           
   0.0000000000   0.0000000000   0.0000000000    100 !GAMMA          
 

EOF


cat >kpdos.in <<EOF
&PROJWFC
    prefix         = '$mat',
    outdir         = "./out"
    !filpdos        = '${mat}.k',
    filproj        = '${mat}proj.out'
    degauss        =  1.00000e-02
    deltae         =  1.00000e-02
    emax           =  5.00000e+01
    emin           = -5.00000e+01
    kresolveddos   = .true.
    ngauss=0
 /

EOF


cat >band.in <<EOF
&BANDS
    prefix         ="$mat"
    outdir         = "./out"
    filband        ="$mat.band"
 /
EOF


prun pw.x <scf.in > scf.out
prun pw.x <nscf.in > nscf.out
prun pw.x <bands.in > bands.out
prun bands.x <band.in > band.out
projwfc.x <kpdos.in > kpdos.out

