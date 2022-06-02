!=======================================================================
! A Neural Network Emulator of the Level-ice meltpond parametrisation
!
! This subroutine is an emulator of the Level-ice melt pond parametrisation
! The parametrisation which it emulates is described as:
!
! ----
!
! This meltpond parameterization was developed for use with the delta-
! Eddington radiation scheme, and only affects the radiation budget in
! the model.  That is, although the pond volume is tracked, that liquid
! water is not used elsewhere in the model for mass budgets or other
! physical processes.
!
! The input is the same as the original Level-ice parametrisation.
! The code builds on the work of the same module in the official Icepack 
! model, differing in how the variables modified by this parametrisation 
! (apnd, hpnd, ipnd, ffrac) are calculated.
!
! ----
!  For subroutines compute_ponds_lvl and selu:
!
! author: Simon Driscoll (University of Reading)
! contact: s.driscoll@pgr.reading.ac.uk
!
!=======================================================================

      module icepack_meltpond_lvl

      use icepack_kinds
      use icepack_parameters, only: c0, c1, c2, c10, p01, p5, puny
      use icepack_parameters, only: viscosity_dyn, rhoi, rhos, rhow, Timelt, Tffresh, Lfresh
      use icepack_parameters, only: gravit, depressT, rhofresh, kice, pndaspect, use_smliq_pnd
      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      implicit none

      private
      public :: compute_ponds_lvl

!=======================================================================

      contains

!=======================================================================

      subroutine compute_ponds_lvl(dt,     nilyr,        &
                                   ktherm,               &
                                   hi_min, dpscale,      &
                                   frzpnd,               &
                                   rfrac,  meltt, melts, &
                                   frain,  Tair,  fsurfn,&
                                   dhs,    ffrac,        &
                                   aicen,  vicen, vsnon, &
                                   qicen,  sicen,        &
                                   Tsfcn,  alvl,         &
                                   apnd,   hpnd,  ipnd,  &
                                   meltsliqn)

      integer (kind=int_kind), intent(in) :: &
         nilyr, &    ! number of ice layers
         ktherm      ! type of thermodynamics (0 0-layer, 1 BL99, 2 mushy)

      real (kind=dbl_kind), intent(in) :: &
         dt,       & ! time step (s)  
         hi_min,   & ! minimum ice thickness allowed for thermo (m)
         dpscale     ! alter e-folding time scale for flushing

      character (len=char_len), intent(in) :: &
         frzpnd      ! pond refreezing parameterization

      real (kind=dbl_kind), &
         intent(in) :: &
         Tsfcn, &    ! surface temperature (C)
         alvl,  &    ! fraction of level ice
         rfrac, &    ! water fraction retained for melt ponds
         meltt, &    ! top melt rate (m/s)
         melts, &    ! snow melt rate (m/s)
         frain, &    ! rainfall rate (kg/m2/s)
         Tair,  &    ! air temperature (K)
         fsurfn,&    ! atm-ice surface heat flux  (W/m2)
         aicen, &    ! ice area fraction
         vicen, &    ! ice volume (m)
         vsnon, &    ! snow volume (m)
         meltsliqn   ! liquid contribution to meltponds in dt (kg/m^2)

      real (kind=dbl_kind), &
         intent(inout) :: &
         apnd, hpnd, ipnd

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         qicen, &  ! ice layer enthalpy (J m-3)
         sicen     ! salinity (ppt)   

      real (kind=dbl_kind), &
         intent(in) :: &
         dhs       ! depth difference for snow on sea ice and pond ice

      real (kind=dbl_kind), &
         intent(out) :: &
         ffrac     ! fraction of fsurfn over pond used to melt ipond

      ! local temporary variables
      real (kind=dbl_kind), dimension (18) :: &
         x_t, &    ! x_t 
         layer_1   ! layer_1
      real (kind=dbl_kind), dimension (4) :: layer_2   ! layer_2

      REAL(KIND=dbl_kind) :: fortran_layer_1_consts(18), &
       fortran_layer_2_consts(4), fortran_scaler_const(18,2)

      REAL (kind=dbl_kind) :: fortran_layer_2_weights_1(4)
      REAL (kind=dbl_kind) :: fortran_layer_2_weights_2(4)
      REAL (kind=dbl_kind) :: fortran_layer_2_weights_3(4)
      REAL (kind=dbl_kind) :: fortran_layer_2_weights_4(4)
      REAL (kind=dbl_kind) :: fortran_layer_2_weights_5(4)
      REAL (kind=dbl_kind) :: fortran_layer_2_weights_6(4)
      REAL (kind=dbl_kind) :: fortran_layer_2_weights_7(4)
      REAL (kind=dbl_kind) :: fortran_layer_2_weights_8(4)
      REAL (kind=dbl_kind) :: fortran_layer_2_weights_9(4)
      REAL (kind=dbl_kind) :: fortran_layer_2_weights_10(4)
      REAL (kind=dbl_kind) :: fortran_layer_2_weights_11(4)
      REAL (kind=dbl_kind) :: fortran_layer_2_weights_12(4)
      REAL (kind=dbl_kind) :: fortran_layer_2_weights_13(4)
      REAL (kind=dbl_kind) :: fortran_layer_2_weights_14(4)
      REAL (kind=dbl_kind) :: fortran_layer_2_weights_15(4)
      REAL (kind=dbl_kind) :: fortran_layer_2_weights_16(4)
      REAL (kind=dbl_kind) :: fortran_layer_2_weights_17(4)
      REAL (kind=dbl_kind) :: fortran_layer_2_weights_18(4)

      REAL (kind=dbl_kind) :: fortran_layer_1_weights_1(18)
      REAL (kind=dbl_kind) :: fortran_layer_1_weights_2(18)
      REAL (kind=dbl_kind) :: fortran_layer_1_weights_3(18)
      REAL (kind=dbl_kind) :: fortran_layer_1_weights_4(18)
      REAL (kind=dbl_kind) :: fortran_layer_1_weights_5(18)
      REAL (kind=dbl_kind) :: fortran_layer_1_weights_6(18)
      REAL (kind=dbl_kind) :: fortran_layer_1_weights_7(18)
      REAL (kind=dbl_kind) :: fortran_layer_1_weights_8(18)
      REAL (kind=dbl_kind) :: fortran_layer_1_weights_9(18)
      REAL (kind=dbl_kind) :: fortran_layer_1_weights_10(18) 
      REAL (kind=dbl_kind) :: fortran_layer_1_weights_11(18)
      REAL (kind=dbl_kind) :: fortran_layer_1_weights_12(18)
      REAL (kind=dbl_kind) :: fortran_layer_1_weights_13(18) 
      REAL (kind=dbl_kind) :: fortran_layer_1_weights_14(18)
      REAL (kind=dbl_kind) :: fortran_layer_1_weights_15(18)
      REAL (kind=dbl_kind) :: fortran_layer_1_weights_16(18)
      REAL (kind=dbl_kind) :: fortran_layer_1_weights_17(18)
      REAL (kind=dbl_kind) :: fortran_layer_1_weights_18(18)

      INTEGER :: i, j, i_loop, j_loop
      REAL (kind=dbl_kind) :: selu_x, selu_x_dummy, x, tempmyreal
      REAL (kind=dbl_kind) :: fortran_scaler_stds(18)
      REAL (kind=dbl_kind) :: fortran_scaler_means(18)
      REAL (kind=dbl_kind) :: y_fortran_scaler_means(4)
      REAL (kind=dbl_kind) :: y_fortran_scaler_stds(4)
      !--
      !For every variable in we output its value at the start of the lvl melt pond parametrisation
      !--
      Print *, "++++++++++++++++++++"
      Print *, "Inside compute_ponds_lvl SubRoutine (1)."
      Print *, "++++++++++++++++++++"
      Print *, "nilyr, is:                    ", nilyr
      Print *, "ktherm, is:                   ", ktherm
      Print *, "dt, is:                       ", dt
      Print *, "hi_min, is:                   ", hi_min
      Print *, "dpscale, is:                  ", dpscale
      Print *, "frzpnd, is:                   ", frzpnd
      Print *, "Tsfcn, is:                    ", Tsfcn
      Print *, "alvl, is:                     ", alvl
      Print *, "rfrac, is:                    ", rfrac
      Print *, "meltt, is:                    ", meltt
      Print *, "melts, is:                    ", melts
      Print *, "frain, is:                    ", frain
      Print *, "Tair, is:                     ", Tair
      Print *, "fsurfn, is:                   ", fsurfn
      Print *, "aicen, is:                    ", aicen
      Print *, "vicen, is:                    ", vicen
      Print *, "vsnon, is:                    ", vsnon
      Print *, "meltsliqn, is:                ", meltsliqn
      Print *, "qicen, is:                    ", qicen
      Print *, "sicen, is:                    ", sicen
      Print *, "dhs, is:                      ", dhs
      Print *, "apnd (in), is:                ", apnd
      Print *, "hpnd (in), is:                ", hpnd
      Print *, "ipnd (in), is:                ", ipnd
      Print *, "ffrac (in), is:               ", ffrac
      Print *, "++++++++++++++++++++"
      Print *, "Inside compute_ponds_lvl SubRoutine (2)."
      Print *, "++++++++++++++++++++"

      !--
      !End of Outputting Values at the Start of the Melt Pond Parametrization
      !--
 
!      fortran_scaler_means=(/3600.,7.,0.01,0.2646114,&
!      9.71423986e-05,1.76247313e-05,6.66644749e-06,&
!      265.00853045,-2.07534244,-0.00170162,0.13483694,&
!      0.15309884,0.01481411,-4.9671632,0.30401578,&
!      0.04058217,0.00421382,0.00105196/)
!
!      fortran_scaler_stds=(/0.,0.,0.,0.23762826,&
!      0.00087106,0.00017715,3.1527981e-05,10.90969017,&
!      31.44275292,0.02271638,0.27956266,0.33244626,&
!      0.03937458,7.3650694,0.45935757,0.16684022,&
!      0.0177957,0.0062573/)

      fortran_scaler_means=(/3600.0, 7.0, 0.00999999999999987, 0.2646113956594903,&
9.714239864602315e-05, 1.762473127220455e-05,&
6.666447488584481e-06,&
265.00853044520557,&
-2.0753424353632535,&
-0.0017016244352081136,&
0.13483693606998826,&
0.15309884061482493,&
0.014814112832953533,&
-4.967163197005169,&
0.3040157849498504,&
0.040582169756724226,&
0.004213819393047367,&
0.0010519599091466264/)

fortran_scaler_stds=(/0.0,&
0.0,&
0.0,&
0.23762826100983728,&
0.0008710638655181999,&
0.0001771479373808342,&
3.152798099030236e-05,&
10.909690167632933,&
31.442752918818037,&
0.022716383865801285,&
0.2795626600115734,&
0.3324462560099669,&
0.03937457942062025,&
7.365069400018307,&
0.4593575671298159,&
0.1668402222713476,&
0.017795698641136903,&
0.006257295176157965/)
     
y_fortran_scaler_means=(/0.040593142607650026,&
0.004266672560698676,&
0.0009741549212324771,&
0.0029534839555229136/)

y_fortran_scaler_stds=(/0.16688955556357518,&
0.01791779601435871,&
0.0058788796742160815,&
0.0542442817958668/)

fortran_layer_1_weights_1 =         (/ 8.84567201e-02,  2.80413628e-01, -3.62941682e-01,&
           5.56962490e-02, -4.48166430e-02, -2.44742110e-01,&
           3.34501863e-01, -1.59155726e-02, -2.54723012e-01,&
          -3.01025271e-01,  1.65889740e-01,  2.81499624e-02,&
          -7.43213296e-03, -1.56211495e-01,  4.23392653e-03,&
          -3.90072167e-01,  1.76301420e-01,  1.43817663e-02/)
         
fortran_layer_1_weights_2 =         (/-1.03926212e-01, -3.98738593e-01, -3.18861157e-01,&
          -3.39277625e-01, -1.74841613e-01, -8.84881616e-02,&
          -1.73422679e-01, -1.97139278e-01, -2.61973619e-01,&
          -2.71975875e-01, -6.66919649e-02, -3.72957587e-02,&
          -4.00606513e-01,  3.73700559e-01, -3.07978332e-01,&
           3.47150087e-01, -4.32130694e-02, -6.93191290e-02/)
         
fortran_layer_1_weights_3 =         (/ 1.65355444e-01, -9.49888229e-02, -7.21848011e-03,&
           2.78173685e-01,  7.17760324e-02, -2.08396524e-01,&
           7.78534710e-02,  1.88621163e-01,  3.43812704e-02,&
           3.77780557e-01, -1.64494723e-01,  2.58259773e-01,&
          -4.06697005e-01,  3.42236876e-01,  8.85198116e-02,&
          -1.81261569e-01, -2.05504343e-01, -2.64911056e-01/)
         
fortran_layer_1_weights_4 =         (/ 1.67972684e-01, -2.32080072e-01,  4.10761327e-01,&
          -1.48775712e-01, -2.49420509e-01, -1.17248237e-01,&
          -1.57755554e-01, -3.36749554e-01, -1.87588900e-01,&
          -1.42073520e-02,  6.70602858e-01, -2.07156137e-01,&
          -2.05225036e-01, -3.36952120e-01,  2.03499511e-01,&
          -3.50860953e-02, -2.85749763e-01,  2.58120093e-02/)
         
fortran_layer_1_weights_5 =         (/ 2.08348107e+00, -8.47400784e-01, -8.75980198e-01,&
          -2.27244925e+00,  2.11672926e+00, -8.41158390e-01,&
           3.61004281e+00, -5.66395235e+00,  5.53503096e-01,&
          -1.81946754e-01,  4.72052187e-01,  6.81566596e-02,&
          -6.65201526e-03,  2.71286756e-01,  9.36920568e-02,&
           4.28004837e+00,  3.73928010e-01,  8.94674361e-01/)
         
fortran_layer_1_weights_6 =         (/ 6.06565654e-01, -8.21708739e-01, -3.69177982e-02,&
           1.12911880e+00,  6.11276925e-01, -8.95101070e-01,&
           1.01262882e-01,  8.02844405e-01, -1.67813748e-01,&
           3.61874282e-01,  1.39974177e+00,  1.65080763e-02,&
          -2.68039912e-01, -1.41877562e-01, -4.29339141e-01,&
          -4.81097214e-02, -1.24233596e-01,  1.68497890e-01/)
         
fortran_layer_1_weights_7 =         (/ 2.26911530e-01,  4.52979296e-01, -1.11035839e-01,&
           2.45157287e-01,  2.29914531e-01,  4.45190638e-01,&
          -6.34485781e-02, -1.24132168e-02, -9.57090333e-02,&
           1.42596876e-02,  6.00760728e-02, -9.36319877e-04,&
          -7.64716864e-02,  3.57236899e-02, -2.83878207e-01,&
           1.03447437e-01,  2.03727514e-01,  7.91217089e-02/)
         
fortran_layer_1_weights_8 =         (/-4.26621348e-01, -1.03619850e+00, -2.01276138e-01,&
          -3.90371978e-01, -4.35948014e-01, -1.01055086e+00,&
           2.67175138e-02, -2.98774093e-01,  1.25856578e+00,&
          -3.34533036e-01,  9.18932483e-02,  1.72262769e-02,&
           1.39116895e+00, -8.84691849e-02,  1.67916548e+00,&
          -1.72506452e-01,  4.26246375e-01, -2.63318270e-02/)
         
fortran_layer_1_weights_9 =         (/ 4.09091622e-01, -4.88320112e-01, -1.50129152e-02,&
          -7.87776947e+00,  4.13737237e-01, -5.06598651e-01,&
          -1.04859332e-02, -5.67750406e+00,  1.01361394e+00,&
          -2.82622242e+00, -2.27262452e-01, -1.10058384e-02,&
           1.32909334e+00, -3.26079759e-03,  1.79848802e+00,&
           3.12945426e-01,  1.94066372e-02,  3.24546546e-02/)
         
fortran_layer_1_weights_10 =         (/-3.48331556e-02, -8.22202936e-02,  3.94971657e+00,&
          -1.82308123e-01, -3.59739922e-02, -5.97311556e-02,&
          -1.02641983e+01,  9.29333878e+00, -4.57439721e-02,&
           3.47064042e+00,  8.79766226e-01, -1.05498405e-02,&
           1.19437706e+00,  6.96246386e-01,  8.23065341e-02,&
           9.01768780e+00,  4.48526049e+00, -4.11491729e-02/)
         
fortran_layer_1_weights_11 =         (/ 4.21420336e-02,  3.04841161e-01,  5.35141081e-02,&
          -7.89888203e-01,  4.65331972e-01,  1.94493219e-01,&
           9.80763137e-02, -3.01798254e-01, -2.36931935e-01,&
          -4.72382963e-01,  6.56127512e-01,  2.29932010e-01,&
           2.73187101e-01, -1.09646268e-01,  1.95009768e-01,&
           1.25181139e-01,  9.21006650e-02,  8.09260234e-02/)
         
fortran_layer_1_weights_12 =         (/ 9.34917778e-02, -8.40854943e-02, -8.38572532e-02,&
           5.21262825e-01,  8.64438489e-02, -8.12920779e-02,&
           1.62888438e-01,  5.07464468e-01,  3.18123907e-01,&
           5.56969225e-01,  8.53760302e-01,  1.00015122e-02,&
          -8.38373378e-02, -4.28444564e-01, -3.13030720e-01,&
          -9.66191962e-02,  4.67096176e-03, -2.84041334e-02/)
         
fortran_layer_1_weights_13 =         (/-9.09199491e-02,  1.43443542e-02,  4.82167929e-01,&
           1.21306944e+00, -7.83076063e-02,  1.59397162e-02,&
           7.10728467e-01,  7.83577979e-01, -5.94654989e+00,&
           3.60500038e-01,  9.06786993e-02, -2.89590983e-03,&
          -3.50112438e+00,  1.79933310e+00, -2.73456335e+00,&
          -1.60060436e-01, -3.04748917e+00,  5.14190644e-02/)
         
fortran_layer_1_weights_14 =         (/ 2.08862162e+00,  1.06297676e-02,  8.50548208e-01,&
           1.43406272e-01,  2.18611670e+00, -1.39223167e-03,&
          -5.31938672e-01,  1.69247553e-01,  3.01398015e+00,&
           9.33263898e-02, -3.34185690e-01, -2.38537788e-02,&
           3.04586148e+00,  3.44506562e-01,  3.56360698e+00,&
          -6.38914585e-01,  6.94209456e-01, -4.41901796e-02/)
         
fortran_layer_1_weights_15 =         (/-6.25549376e-01, -8.21106553e-01, -2.88444710e+00,&
          -1.35602772e+00, -6.17289126e-01, -8.09465230e-01,&
           1.04339194e+00, -9.46680069e-01, -2.92517853e+00,&
          -1.48377225e-01,  6.35707915e-01,  1.77979786e-02,&
          -2.82236814e+00,  1.05786443e+00, -2.82072878e+00,&
           5.28717339e-01, -9.37703013e-01,  1.75542012e-01/)
         
fortran_layer_1_weights_16 =         (/-1.21570051e-01, -1.00806868e+00, -1.83453870e+00,&
          -4.95603591e-01,  2.99847215e-01, -6.84235036e-01,&
           1.19487846e+00, -3.23593706e-01,  1.00924551e-01,&
          -2.81116888e-02,  1.37532651e+00,  3.73760518e-03,&
           2.54469980e-02,  1.77175379e+00, -4.30875830e-02,&
          -5.03891259e-02,  1.36339396e-01, -1.32950544e-01/)
         
fortran_layer_1_weights_17 =         (/ 1.54590678e+00, -1.98984432e+00,  5.97017944e-01,&
          -2.39708006e-01,  1.16201699e+00, -5.62999392e+00,&
          -5.85864723e-01, -2.63847113e-01, -1.43209979e-01,&
          -6.56966865e-02, -1.08446431e+00,  9.66987848e-01,&
           1.62798408e-02,  2.88331419e-01,  1.66104153e-01,&
           1.28961816e-01,  5.42708874e-01,  4.80869979e-01/)
         
fortran_layer_1_weights_18 =         (/ 1.79221541e-01, -2.25003079e-01,  5.73252797e-01,&
           3.15885067e-01,  1.74035937e-01,  1.51017487e-01,&
          -3.60386282e-01,  6.05808496e-01, -1.47065461e-01,&
           6.02594757e+00, -4.74337339e-01,  8.37316811e-01,&
          -3.95663857e+00, -2.60224100e-02, -5.72261155e-01,&
          -1.50730804e-01,  5.64021766e-01,  5.23075294e+00/)
  

  
fortran_layer_1_consts =         (/-0.9164963 ,  0.2924885 , -1.6785355 ,  4.0989594 , -0.93711317,&
         -0.47355673,  0.34073138,  1.5443676 , -1.2162033 ,  1.4767244 ,&
          0.779868  ,  0.41047722, -1.2183926 ,  0.02930823, -0.5462315 ,&
         -0.37763873, -1.5859716 ,  0.21565591/)

 fortran_layer_2_weights_1 =        (/-2.0382209e+00,  1.2637959e-01, -1.5550919e-01,  2.2305219e-01/)
 fortran_layer_2_weights_2 =           (/-8.3203744e-03, -3.1474564e-02,  2.5547924e+00, -1.2132616e+00/)
 fortran_layer_2_weights_3 =           (/ 6.8245423e-03,  5.6050724e-04, -1.4687733e-05,  2.2283049e+00/)
 fortran_layer_2_weights_4 =           (/ 1.9336683e-03, -4.5355142e-04,  1.1667073e-03,  3.0836620e+00/)
 fortran_layer_2_weights_5 =           (/ 2.0419343e+00, -8.2099766e-02,  9.7993597e-02, -6.3515723e-01/)
 fortran_layer_2_weights_6 =           (/ 1.3373217e-02,  3.0887997e-02, -2.5921769e+00,  2.7811360e-01/)
 fortran_layer_2_weights_7 =           (/-3.2425273e-04, -6.4034224e-04, -1.2463514e-03, -3.2584918e+00/)
 fortran_layer_2_weights_8 =           (/-4.8256596e-03, -2.7713514e-04,  7.0269057e-03, -5.4241414e+00/)
 fortran_layer_2_weights_9 =           (/-5.1082755e-03,  3.3260563e-03, -1.3122180e-02,  2.7641842e+00/)
 fortran_layer_2_weights_10 =           (/ 7.6884567e-04, -1.8884348e-03,  6.4021032e-03,  2.6751335e+00/)
 fortran_layer_2_weights_11 =           (/-4.3521887e-03, -2.6692827e-03, -1.7259895e-03,  9.9515474e-01/)
 fortran_layer_2_weights_12 =           (/ 8.5380352e-01,  8.7711555e-01,  5.4775372e-02, -5.5035746e-01/)
 fortran_layer_2_weights_13 =           (/-1.6388702e-03, -1.5424542e-02,  4.7086809e-02, -6.0008645e+00/)
 fortran_layer_2_weights_14 =           (/ 5.1045367e-03,  4.7537074e-03, -9.4384374e-03,  1.3258998e+00/)
 fortran_layer_2_weights_15 =           (/ 2.3830070e-03,  4.7654356e-03, -2.0403244e-02,  2.6448138e+00/)
 fortran_layer_2_weights_16 =           (/-8.4060774e-04, -1.8176779e-03, -3.7535343e-03,  2.8039877e+00/)
 fortran_layer_2_weights_17 =           (/ 1.0890267e-02,  6.8359268e-03, -1.7945178e-02,  1.2960835e+00/)
 fortran_layer_2_weights_18 =           (/-1.3706553e-01, -1.4154488e-01,  1.7187421e-01, -2.7636540e+00/)


fortran_layer_2_consts =    (/-0.2800149 , -0.24015978, -0.10848322, -1.6201998 /) 

      x_t(1) = dt
      x_t(2) = nilyr
      x_t(3) = hi_min
      x_t(4) = rfrac
      x_t(5) = meltt
      x_t(6) = melts
      x_t(7) = frain
      x_t(8) = Tair
      x_t(9) = fsurfn
      x_t(10) = dhs
      x_t(11) = aicen
      x_t(12) = vicen
      x_t(13) = vsnon
      x_t(14) = Tsfcn
      x_t(15) = alvl
      x_t(16) = apnd
      x_t(17) = hpnd
      x_t(18) = ipnd

      !Print *, "This is before scaling and ffrac is: ", ffrac
      !Print *, "This is before scaling and x_t(18) is: ", x_t(18)
      !Print *, "This is before scaling and x_t(18) is: ", x_t
      !Print *, "this is Tsfcn ", Tsfcn, alvl, rfrac, meltt, melts, frain, Tair, fsurfn, aicen, vsnon
      !Print *, "Unscaled Data:", x_t
      !-----------------------------------------------------------------
      ! Scale Data
      !-----------------------------------------------------------------
      do i_loop = 1, 18
        tempmyreal = fortran_scaler_stds(i_loop)
        if (tempmyreal==0.0) THEN
          x_t(i_loop) = 0.0
        else 
          x_t(i_loop) = (x_t(i_loop)-fortran_scaler_means(i_loop))/(fortran_scaler_stds(i_loop))
        endif
      end do
      x_t(3)=1.301043e-16
      
      !Print *, "Scaled Data:", x_t
      !Print *, "fortran scaler means", fortran_scaler_means
      !Print *, "fortran scaler stddevs", fortran_scaler_stds
      !Print *, "This is after scaling and x_t(18) is: ", x_t
      !-----------------------------------------------------------------
      ! Tsfcn + alvl + rfrac + meltt + melts + frain + Tair + fsurfn 
      ! + aicen + vicen + vsnon + dhs + apnd_in + hpnd_in + ipnd_in 
      ! + ffrac_in 
      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      ! Initalise arrays
      !-----------------------------------------------------------------
      do i_loop = 1, 18
          layer_1(i_loop) = 0.0
      end do

      do i_loop = 1, 4 
          layer_2(i_loop) = 0.0 
      end do
      !-----------------------------------------------------------------
      ! Calculate layers of the neural network with the layer weights 
      ! and constants for each input to each respective layer 
      !-----------------------------------------------------------------
      do i_loop = 1, 18

!            Print *, "fortran_layer_1_weights_1(i_loop) is", fortran_layer_1_weights_1(i_loop)
!            Print *, "x_t(1) is", x_t(1)
!            Print *, "fortran_layer_1_weights_1(i_loop) is", fortran_layer_1_weights_2(i_loop)
!            Print *, "x_t(2) is", x_t(2)
!            Print *, "fortran_layer_1_weights_1(i_loop) is", fortran_layer_1_weights_3(i_loop)
!            Print *, "x_t(3) is", x_t(3)
!            Print *, "fortran_layer_1_weights_1(i_loop) is", fortran_layer_1_weights_4(i_loop)
!            Print *, "x_t(4) is", x_t(4)
!            Print *, "fortran_layer_1_weights_1(i_loop) is", fortran_layer_1_weights_5(i_loop)
!            Print *, "x_t(5) is", x_t(5)
!            Print *, "fortran_layer_1_weights_1(i_loop) is", fortran_layer_1_weights_6(i_loop)
!            Print *, "x_t(6) is", x_t(6)
!            Print *, "fortran_layer_1_weights_1(i_loop) is", fortran_layer_1_weights_7(i_loop)
!            Print *, "x_t(7) is", x_t(7)
!            Print *, "fortran_layer_1_weights_1(i_loop) is", fortran_layer_1_weights_8(i_loop)
!            Print *, "x_t(8) is", x_t(8)
!            Print *, "fortran_layer_1_weights_1(i_loop) is", fortran_layer_1_weights_9(i_loop)
!            Print *, "x_t(9) is", x_t(9)
!            Print *, "fortran_layer_1_weights_1(i_loop) is", fortran_layer_1_weights_10(i_loop)
!            Print *, "x_t(10) is", x_t(10)
!            Print *, "fortran_layer_1_weights_1(i_loop) is", fortran_layer_1_weights_11(i_loop)
!            Print *, "x_t(11) is", x_t(11)
!            Print *, "fortran_layer_1_weights_1(i_loop) is", fortran_layer_1_weights_12(i_loop)
!            Print *, "x_t(12) is", x_t(12)
!            Print *, "fortran_layer_1_weights_1(i_loop) is", fortran_layer_1_weights_13(i_loop)
!            Print *, "x_t(13) is", x_t(13)
!            Print *, "fortran_layer_1_weights_1(i_loop) is", fortran_layer_1_weights_14(i_loop)
!            Print *, "x_t(14) is", x_t(14)
!            Print *, "fortran_layer_1_weights_1(i_loop) is", fortran_layer_1_weights_15(i_loop)
!            Print *, "x_t(15) is", x_t(15)
!            Print *, "fortran_layer_1_weights_1(i_loop) is", fortran_layer_1_weights_16(i_loop)
!            Print *, "x_t(16) is", x_t(16)
!            Print *, "fortran_layer_1_weights_1(i_loop) is", fortran_layer_1_weights_17(i_loop)
!            Print *, "x_t(17) is", x_t(17)
!            Print *, "fortran_layer_1_weights_1(i_loop) is", fortran_layer_1_weights_18(i_loop)
!            Print *, "x_t(18) is", x_t(18)

            layer_1(i_loop) = layer_1(i_loop) + fortran_layer_1_weights_1(i_loop) * x_t(1)
            layer_1(i_loop) = layer_1(i_loop) + fortran_layer_1_weights_2(i_loop) * x_t(2)
            layer_1(i_loop) = layer_1(i_loop) + fortran_layer_1_weights_3(i_loop) * x_t(3)
            layer_1(i_loop) = layer_1(i_loop) + fortran_layer_1_weights_4(i_loop) * x_t(4)
            layer_1(i_loop) = layer_1(i_loop) + fortran_layer_1_weights_5(i_loop) * x_t(5)
            layer_1(i_loop) = layer_1(i_loop) + fortran_layer_1_weights_6(i_loop) * x_t(6)
            layer_1(i_loop) = layer_1(i_loop) + fortran_layer_1_weights_7(i_loop) * x_t(7)
            layer_1(i_loop) = layer_1(i_loop) + fortran_layer_1_weights_8(i_loop) * x_t(8)
            layer_1(i_loop) = layer_1(i_loop) + fortran_layer_1_weights_9(i_loop) * x_t(9)
            layer_1(i_loop) = layer_1(i_loop) + fortran_layer_1_weights_10(i_loop) * x_t(10)
            layer_1(i_loop) = layer_1(i_loop) + fortran_layer_1_weights_11(i_loop) * x_t(11)
            layer_1(i_loop) = layer_1(i_loop) + fortran_layer_1_weights_12(i_loop) * x_t(12)
            layer_1(i_loop) = layer_1(i_loop) + fortran_layer_1_weights_13(i_loop) * x_t(13)
            layer_1(i_loop) = layer_1(i_loop) + fortran_layer_1_weights_14(i_loop) * x_t(14)
            layer_1(i_loop) = layer_1(i_loop) + fortran_layer_1_weights_15(i_loop) * x_t(15)
            layer_1(i_loop) = layer_1(i_loop) + fortran_layer_1_weights_16(i_loop) * x_t(16)
            layer_1(i_loop) = layer_1(i_loop) + fortran_layer_1_weights_17(i_loop) * x_t(17)
            layer_1(i_loop) = layer_1(i_loop) + fortran_layer_1_weights_18(i_loop) * x_t(18)
      end do
      
      !Print *, "fortran for python layer_1 before constants are added. Alison layer_1(1): ", layer_1(1)
      !Print *, "fortran layer 1 weights", fortran_layer_1_weights
 
      do i_loop = 1, 18
         layer_1(i_loop) = fortran_layer_1_consts(i_loop) + layer_1(i_loop)
      end do

      !Print *, "fortran for python layer_1 after constants are added. Alison : ", layer_1

      do i_loop = 1, 18
         selu_x_dummy = layer_1(i_loop) 
         !Print *, "Selu for Python"
         !Print *, "Selu X In: ", selu_x_dummy
         CALL selu(selu_x_dummy)
         !Print *, "Selu X Out: ", selu_x_dummy
         layer_1(i_loop) = selu_x_dummy       
      end do

      !Print *, "This is layer1 as with python: ", layer_1

      do i_loop = 1, 4
            layer_2(i_loop) = layer_2(i_loop) + fortran_layer_2_weights_1(i_loop)*layer_1(1) 
            layer_2(i_loop) = layer_2(i_loop) + fortran_layer_2_weights_2(i_loop)*layer_1(2) 
            layer_2(i_loop) = layer_2(i_loop) + fortran_layer_2_weights_3(i_loop)*layer_1(3) 
            layer_2(i_loop) = layer_2(i_loop) + fortran_layer_2_weights_4(i_loop)*layer_1(4) 
            layer_2(i_loop) = layer_2(i_loop) + fortran_layer_2_weights_5(i_loop)*layer_1(5) 
            layer_2(i_loop) = layer_2(i_loop) + fortran_layer_2_weights_6(i_loop)*layer_1(6) 
            layer_2(i_loop) = layer_2(i_loop) + fortran_layer_2_weights_7(i_loop)*layer_1(7) 
            layer_2(i_loop) = layer_2(i_loop) + fortran_layer_2_weights_8(i_loop)*layer_1(8) 
            layer_2(i_loop) = layer_2(i_loop) + fortran_layer_2_weights_9(i_loop)*layer_1(9) 
            layer_2(i_loop) = layer_2(i_loop) + fortran_layer_2_weights_10(i_loop)*layer_1(10) 
            layer_2(i_loop) = layer_2(i_loop) + fortran_layer_2_weights_11(i_loop)*layer_1(11) 
            layer_2(i_loop) = layer_2(i_loop) + fortran_layer_2_weights_12(i_loop)*layer_1(12) 
            layer_2(i_loop) = layer_2(i_loop) + fortran_layer_2_weights_13(i_loop)*layer_1(13) 
            layer_2(i_loop) = layer_2(i_loop) + fortran_layer_2_weights_14(i_loop)*layer_1(14) 
            layer_2(i_loop) = layer_2(i_loop) + fortran_layer_2_weights_15(i_loop)*layer_1(15) 
            layer_2(i_loop) = layer_2(i_loop) + fortran_layer_2_weights_16(i_loop)*layer_1(16) 
            layer_2(i_loop) = layer_2(i_loop) + fortran_layer_2_weights_17(i_loop)*layer_1(17) 
            layer_2(i_loop) = layer_2(i_loop) + fortran_layer_2_weights_18(i_loop)*layer_1(18) 
      end do

      do i_loop = 1, 4
         layer_2(i_loop) = layer_2(i_loop) + fortran_layer_2_consts(i_loop)
      end do

      !Print *, "fortran_layer_2_weights: ", fortran_layer_2_weights
      !Print *, "fortran_layer_2_consts: ", fortran_layer_2_consts

      do i_loop = 1, 4
         selu_x_dummy = layer_2(i_loop)
         !Print *, "Selu X In: ", selu_x_dummy
         CALL selu(selu_x_dummy)
         !Print *, "Selu X Out: ", selu_x_dummy
         layer_2(i_loop) = selu_x_dummy
      end do
      
!      Print *, "This is layer2 as with python: ", layer_2

      !-----------------------------------------------------------------
      ! Update pond variables and fraction of fsurfn.
      !-----------------------------------------------------------------
      apnd = layer_2(1)
      hpnd = layer_2(2)
      ipnd = layer_2(3)
      ffrac = layer_2(4)
!      Print *, "END OF SUBR. apnd: ", apnd 
!      Print *, "hpnd: ", hpnd
!      Print *, "ipnd: ", ipnd
!      Print *, "ffrac: ", ffrac
!      Print *, "fortran_scaler_stds: ", fortran_scaler_stds
!      Print *, "fortran_scaler_means: ", fortran_scaler_means
      !Print *, " layer_2: ", layer_2
      !Print *, "layer_1: ", layer_1
      Print *, "fortran_scaler_means(15)", fortran_scaler_means(15)
      Print *, "fortran_scaler_means(16)", fortran_scaler_means(16)
      Print *, "fortran_scaler_means(17)", fortran_scaler_means(17)
      Print *, "fortran_scaler_means(18)", fortran_scaler_means(18)
      Print *, "fortran_scaler_stds(15)", fortran_scaler_stds(15)
      Print *, "fortran_scaler_stds(16)", fortran_scaler_stds(16)
      Print *, "fortran_scaler_stds(17)", fortran_scaler_stds(17)
      Print *, "fortran_scaler_stds(18)", fortran_scaler_stds(18)
      

      Print *, "apnd before scaling and apnd scalers: ", apnd
      Print *, fortran_scaler_stds(15)
      Print *, fortran_scaler_means(15)
      Print *, "layer2", layer_2
      apnd=apnd*y_fortran_scaler_stds(1)+y_fortran_scaler_means(1)
      hpnd=hpnd*y_fortran_scaler_stds(2)+y_fortran_scaler_means(2)
      ipnd=ipnd*y_fortran_scaler_stds(3)+y_fortran_scaler_means(3)
      ffrac=ffrac*y_fortran_scaler_stds(4)+y_fortran_scaler_means(4)
!      apnd=apnd*(fortran_scaler_stds(15))+fortran_scaler_means(15) 
!      hpnd=hpnd*(fortran_scaler_stds(16))+fortran_scaler_means(16) 
!      ipnd=ipnd*(fortran_scaler_stds(17))+fortran_scaler_means(17) 
!      ffrac=ffrac*(fortran_scaler_stds(18))+fortran_scaler_means(18) 

      !Print *, "postscaling output: ", apnd, hpnd, ipnd, ffrac
      !if (apnd < 0.0) apnd = 0.0 
      !if (hpnd < 0.0) hpnd = 0.0
      !if (ipnd < 0.0) ipnd = 0.0 
      !if (ffrac < 0.0) ffrac = 0.0 

      !if (apnd > 1.0) apnd = 1.0 

      !if (ffrac > 1.0) ffrac = 1.0 

      Print *, "++++++++++++++++++++"
      Print *, "Inside compute_ponds_lvl SubRoutine (3)."
      Print *, "++++++++++++++++++++"
      Print *, "apnd (out):                   ", apnd
      Print *, "hpnd (out):                   ", hpnd
      Print *, "ipnd (out):                   ", ipnd
      Print *, "ffrac (out):                  ", ffrac
      Print *, "++++++++++++++++++++"
      Print *, "Inside compute_ponds_lvl SubRoutine (4)."
      Print *, "++++++++++++++++++++"

      end subroutine compute_ponds_lvl

!=======================================================================

!===============================================================
!
! selu activation function
! 
! This copies the keras selu activation function
! For more details see: https://keras.io/api/layers/activations/
!
!===============================================================

      subroutine selu(x)

      real (kind=dbl_kind), &
         intent(inout) :: &
         x

      ! local variables

      real (kind=dbl_kind) :: &
         selu_alpha, &    ! alpha 
         selu_scale, &    ! scale
         selu_x ! local selu_x function
      !Print *, "We are inside Selu, and x is: ", x
      selu_scale = 1.05070098
      selu_alpha = 1.67326324      
      !-----------------------------------------------------------------
      ! Compute selu activation function
      !-----------------------------------------------------------------

      if (x > 0.0) selu_x = selu_scale * x 
      if (x < 0.0) selu_x = selu_scale * selu_alpha * (exp(x) - 1.0)
      !Print *, "Inside selu and selu_x is: ", selu_x
      x = selu_x
      !Print *, "Inside selu and selu_x is: ", x

      end subroutine selu
  
!=======================================================================

      end module icepack_meltpond_lvl

!=======================================================================
