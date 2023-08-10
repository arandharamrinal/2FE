module potvars 
use molparams,only : Ncarts,Natoms,N_int,Nang,Nbonds
implicit none 
integer::Mmax,Nmax,n_eb,n_ea,n_ed, n_ob,n_oa,n_od,n_fix_tors
integer,parameter,dimension(0:1)::tors_idx =[15,20]
real(8),parameter :: sqrt_2 = dsqrt(2.d0),sqrt_2_inv = 1.d0/dsqrt(2.d0)
real(8) :: dx = 0.002,dx_inv = 500.d0
real(8),allocatable::sym_mat(:,:),sym_mat_transpose(:,:)
character(len=200)::potParamDir
!Numerical derivative array
real(8) :: five_pt_1st_deriv(0:4)
integer,allocatable::atom_ind(:,:)
real(8),parameter:: degree_to_radian = 4.d0*datan(1.d0)/180.d0,radian_to_degree = 180.d0 /4.d0/datan(1.d0)
integer,parameter,dimension(0:11)  :: sym_idx     = [3,4,5,6,10,11,12,13,16,17,18,19]
integer,parameter,dimension(0:5)  :: odd_fn_idx  = [4,6,11,13,16,18] 
integer,parameter,dimension(0:12) :: even_fn_idx = [0,1,2,3,5,7,8,9,10,12,14,17,19] 
integer,allocatable::max_id(:)
!stride_arr
integer,allocatable::stride_arr_quad_even(:,:),stride_arr_quad_odd(:,:),stride_arr_qubic_even_gt5(:,:),stride_arr_qubic_even_lt5(:,:),stride_arr_qubic_even_lt2(:,:),stride_arr_qubic_odd_gt5(:,:),stride_arr_qubic_odd_lt5(:,:),stride_arr_qubic_odd_lt2(:,:),stride_arr_quartic_even_gt5(:,:),stride_arr_quartic_even_lt5(:,:),stride_arr_quartic_even_lt2(:,:),stride_arr_quartic_odd_gt5(:,:),stride_arr_quartic_odd_lt5(:,:),stride_arr_quartic_odd_lt2(:,:) 
integer,allocatable::stride_arr_qo_even(:,:),stride_arr_qo_odd(:,:),stride_arr_qo_ang_even(:,:),stride_arr_qo_ang_odd(:,:)
!fn_order
integer,allocatable::stride_arr_pot(:,:),stride_arr_b_even(:,:),stride_arr_b_odd(:,:),stride_arr_a_even(:,:),stride_arr_a_odd(:,:),stride_arr_d_even(:,:),stride_arr_d_odd(:,:)
integer :: fn_order_pot,fn_order_bonds_even,fn_order_angs_even,fn_order_dihs_even,fn_order_bonds_odd,fn_order_angs_odd,fn_order_dihs_odd
integer :: Ncoeff_pot(0:1),Ncoeff_bonds_even(0:1),Ncoeff_bonds_odd(0:1),Ncoeff_angs_even(0:1),Ncoeff_angs_odd(0:1),Ncoeff_dihs_even(0:1),Ncoeff_dihs_odd(0:1) 
integer :: fn_order_fij_even,fn_order_fij_odd,fn_order_fijk_even_gt5,fn_order_fijk_even_lt5,fn_order_fijk_even_lt2,fn_order_fijk_odd_gt5,fn_order_fijk_odd_lt5,fn_order_fijk_odd_lt2,fn_order_fijkl_even_gt5,fn_order_fijkl_even_lt5,fn_order_fijkl_even_lt2,fn_order_fijkl_odd_gt5,fn_order_fijkl_odd_lt5,fn_order_fijkl_odd_lt2,fn_order_fijkl_even_ltpt1,fn_order_fijkl_odd_ltpt1
integer :: fn_order_qo_odd,fn_order_qo_even
integer :: fn_order_qo_ang_odd,fn_order_qo_ang_even
!n_params
integer::n_fcs_quad_even,n_fcs_quad_odd,n_fcs_qubic_even_gt5,n_fcs_qubic_even_lt5,n_fcs_qubic_even_lt2,n_fcs_qubic_odd_gt5,n_fcs_qubic_odd_lt5,n_fcs_qubic_odd_lt2,n_fcs_quartic_even_gt5,n_fcs_quartic_even_lt5,n_fcs_quartic_even_lt2,n_fcs_quartic_odd_gt5,n_fcs_quartic_odd_lt5,n_fcs_quartic_odd_lt2
integer::n_fcs_qo_even,n_fcs_qo_odd,n_fcs_quartic_even_ltpt1,n_fcs_quartic_odd_ltpt1
integer::n_fcs_qo_ang_even,n_fcs_qo_ang_odd

!Ncoeff
integer::Ncoeff_fij_even(0:1),Ncoeff_fij_odd(0:1),Ncoeff_fijk_even_gt5(0:1),Ncoeff_fijk_even_lt5(0:1),Ncoeff_fijk_even_lt2(0:1),Ncoeff_fijk_odd_gt5(0:1),Ncoeff_fijk_odd_lt5(0:1),Ncoeff_fijk_odd_lt2(0:1),Ncoeff_fijkl_even_gt5(0:1),Ncoeff_fijkl_even_lt5(0:1),Ncoeff_fijkl_even_lt2(0:1),Ncoeff_fijkl_odd_gt5(0:1),Ncoeff_fijkl_odd_lt5(0:1),Ncoeff_fijkl_odd_lt2(0:1) 
integer::Ncoeff_qo_even(0:1),Ncoeff_qo_odd(0:1)
integer::Ncoeff_qo_ang_even(0:1),Ncoeff_qo_ang_odd(0:1)
integer,allocatable::fc_idx_quad_even(:,:),fc_idx_quad_odd(:,:),fc_idx_qubic_even_gt5(:,:),fc_idx_qubic_even_lt5(:,:),fc_idx_qubic_even_lt2(:,:),fc_idx_qubic_odd_gt5(:,:),fc_idx_qubic_odd_lt5(:,:),fc_idx_qubic_odd_lt2(:,:),fc_idx_quartic_even_gt5(:,:),fc_idx_quartic_even_lt5(:,:),fc_idx_quartic_even_lt2(:,:),fc_idx_quartic_odd_gt5(:,:),fc_idx_quartic_odd_lt5(:,:),fc_idx_quartic_odd_lt2(:,:) ,fc_idx_quartic_odd_ltpt1(:,:) ,fc_idx_quartic_even_ltpt1(:,:) 
integer,allocatable::fc_idx_qo_even(:,:),fc_idx_qo_odd(:,:), fc_idx_qo_ang_even(:,:),fc_idx_qo_ang_odd(:,:)
!pot coeff
real(8),allocatable::Vt_coeff(:)
!fitted internal coordinates parameter
real(8),allocatable :: fitted_bonds_coeff_even(:,:),fitted_bonds_coeff_odd(:,:),fitted_angs_coeff_even(:,:),fitted_angs_coeff_odd(:,:),fitted_dihs_coeff_even(:,:),fitted_dihs_coeff_odd(:,:) 
!force const 
real(8),allocatable::fitted_fc_coeff_ij_even(:,:),fitted_fc_coeff_ij_odd(:,:),fitted_fc_coeff_ijk_even_gt5(:,:),fitted_fc_coeff_ijk_even_lt5(:,:),fitted_fc_coeff_ijk_even_lt2(:,:),fitted_fc_coeff_ijk_odd_gt5(:,:),fitted_fc_coeff_ijk_odd_lt5(:,:),fitted_fc_coeff_ijk_odd_lt2(:,:),fitted_fc_coeff_ijkl_even_gt5(:,:),fitted_fc_coeff_ijkl_even_lt5(:,:),fitted_fc_coeff_ijkl_even_lt2(:,:),fitted_fc_coeff_ijkl_odd_gt5(:,:),fitted_fc_coeff_ijkl_odd_lt5(:,:),fitted_fc_coeff_ijkl_odd_lt2(:,:) ,fitted_fc_coeff_ijkl_even_ltpt1(:) ,fitted_fc_coeff_ijkl_odd_ltpt1(:) 
real(8),allocatable::fitted_fc_coeff_qo_even(:,:),fitted_fc_coeff_qo_odd(:,:), fitted_fc_coeff_qo_ang_even(:,:),fitted_fc_coeff_qo_ang_odd(:,:)
real(8),allocatable:: dim_scal_factor(:),dim_scal_factor_inv(:)
real(8),allocatable::dsdx(:,:),dVds(:),dVdx(:),dVdx_num(:)

contains

subroutine initializePotvars()
    implicit none 
    integer::i,ctr,id1,id2,j
    allocate(dsdx(0:Ncarts-1,0:N_int-1),dVdx(0:Ncarts-1),dVdx_num(0:Ncarts-1))
    allocate(dVds(0:N_int-1))
    allocate(dim_scal_factor(0:N_int-1),dim_scal_factor_inv(0:N_int-1))
    allocate(sym_mat(0:N_int-1,0:N_int-1),sym_mat_transpose(0:N_int-1,0:N_int-1))
    sym_mat_transpose = 0.d0
    five_pt_1st_deriv = [1.d0,-8.d0,0.d0,8.d0,-1.d0]
    five_pt_1st_deriv = five_pt_1st_deriv/12.d0
    sym_mat      = 0.d0
    ctr = 0
    i = 0
    open(unit =13,file ="pot.param",status="old")
    read(13,*)potParamDir
    !Read_scale_factor
    open(unit =11,file ="scale_factor_avtz_ref.dat",status="old") 
    read(11,*)(dim_scal_factor(j),j=0,N_int-1)
    do i = 0,N_int-1
    dim_scal_factor_inv(i) = 1.d0 / dim_scal_factor(i)
    enddo
    !****************************************************************************************#
    !Parameters for Minimum energy path along OCCF and HOCC diherdrals
    Mmax = 8;Nmax = 6
    call get_fn_order_even(Mmax,Nmax,fn_order_pot,Ncoeff_pot)
    allocate(stride_arr_pot(0:fn_order_pot-1,0:1))
    call get_stride_arr_even(Mmax,Nmax,fn_order_pot,stride_arr_pot)
    !********************************For bonds********************************#
    !Even expansion
    Mmax = 8;Nmax = 6
    call get_fn_order_even(Mmax,Nmax,fn_order_bonds_even,Ncoeff_bonds_even)
    allocate(stride_arr_b_even(0:fn_order_bonds_even-1,0:1))
    call get_stride_arr_even(Mmax,Nmax,fn_order_bonds_even,stride_arr_b_even)
    !Odd expansion 
    Mmax = 8;Nmax = 6
    call get_fn_order_odd(Mmax,Nmax,fn_order_bonds_odd,Ncoeff_bonds_odd)
    allocate(stride_arr_b_odd(0:fn_order_bonds_odd-1,0:1))
    call get_stride_arr_odd(Mmax,Nmax,fn_order_bonds_odd,stride_arr_b_odd)
    !********************************For Angles*******************************#
    !Even expansion
    Mmax = 8;Nmax = 6
    call get_fn_order_even(Mmax,Nmax,fn_order_angs_even,Ncoeff_angs_even)
    allocate(stride_arr_a_even(0:fn_order_angs_even-1,0:1))
    call get_stride_arr_even(Mmax,Nmax,fn_order_angs_even,stride_arr_a_even)
    !Odd expansion 
    Mmax = 8;Nmax = 6
    call get_fn_order_odd(Mmax,Nmax,fn_order_angs_odd,Ncoeff_angs_odd )
    allocate(stride_arr_a_odd(0:fn_order_angs_odd-1,0:1))
    call get_stride_arr_odd(Mmax,Nmax,fn_order_angs_odd,stride_arr_a_odd)
    !********************************For Dihedrals*******************************#
    !Even expansion
    Mmax = 8;Nmax = 6
    call get_fn_order_even(Mmax,Nmax,fn_order_dihs_even,Ncoeff_dihs_even)
    allocate(stride_arr_d_even(0:fn_order_dihs_even-1,0:1))
    call get_stride_arr_even(Mmax,Nmax,fn_order_dihs_even,stride_arr_d_even)
    !Odd expansion 
    Mmax = 8;Nmax = 6
    call get_fn_order_odd(Mmax,Nmax,fn_order_dihs_odd,Ncoeff_dihs_odd)
    allocate(stride_arr_d_odd(0:fn_order_dihs_odd-1,0:1))
    call get_stride_arr_odd(Mmax,Nmax,fn_order_dihs_odd,stride_arr_d_odd)
    !**************************************************************************#
    !Quadratic force constants
    !Even 
    Mmax = 8;Nmax = 6
    call get_fn_order_even(Mmax,Nmax,fn_order_fij_even,Ncoeff_fij_even)
    allocate(stride_arr_quad_even(0:fn_order_fij_even-1,0:1))
    call get_stride_arr_even(Mmax,Nmax,fn_order_fij_even,stride_arr_quad_even)
    !Odd
    Mmax = 8;Nmax = 6
    call get_fn_order_odd(Mmax,Nmax,fn_order_fij_odd,Ncoeff_fij_odd)
    allocate(stride_arr_quad_odd(0:fn_order_fij_odd-1,0:1))
    call get_stride_arr_odd(Mmax,Nmax,fn_order_fij_odd,stride_arr_quad_odd)
    !****************************************************************************************************************************************#
    !Cubic force constants
    !Even
    !One Body 
    Mmax = 8;Nmax = 6
    call get_fn_order_even(Mmax,Nmax,fn_order_fijk_even_gt5,Ncoeff_fijk_even_gt5)
    allocate(stride_arr_qubic_even_gt5(0:fn_order_fijk_even_gt5-1,0:1))
    call get_stride_arr_even(Mmax,Nmax,fn_order_fijk_even_gt5,stride_arr_qubic_even_gt5)
    !Two Body 
    Mmax = 6;Nmax = 6
    call get_fn_order_even(Mmax,Nmax,fn_order_fijk_even_lt5,Ncoeff_fijk_even_lt5)
    allocate(stride_arr_qubic_even_lt5(0:fn_order_fijk_even_lt5-1,0:1))
    call get_stride_arr_even(Mmax,Nmax,fn_order_fijk_even_lt5,stride_arr_qubic_even_lt5)
    !Three Body
    Mmax = 5;Nmax = 5
    call get_fn_order_even(Mmax,Nmax,fn_order_fijk_even_lt2,Ncoeff_fijk_even_lt2)
    allocate(stride_arr_qubic_even_lt2(0:fn_order_fijk_even_lt2-1,0:1))
    call get_stride_arr_even(Mmax,Nmax,fn_order_fijk_even_lt2,stride_arr_qubic_even_lt2)
    !Odd
    !One Body 
    Mmax = 8;Nmax = 6
    call get_fn_order_odd(Mmax,Nmax,fn_order_fijk_odd_gt5,Ncoeff_fijk_odd_gt5)
    allocate(stride_arr_qubic_odd_gt5(0:fn_order_fijk_odd_gt5-1,0:1))
    call get_stride_arr_odd(Mmax,Nmax,fn_order_fijk_odd_gt5,stride_arr_qubic_odd_gt5)
    !Two Body 
    Mmax = 6;Nmax = 6
    call get_fn_order_odd(Mmax,Nmax,fn_order_fijk_odd_lt5,Ncoeff_fijk_odd_lt5)
    allocate(stride_arr_qubic_odd_lt5(0:fn_order_fijk_odd_lt5-1,0:1))
    call get_stride_arr_odd(Mmax,Nmax,fn_order_fijk_odd_lt5,stride_arr_qubic_odd_lt5)
    !Three Body
    Mmax = 5;Nmax = 5
    call get_fn_order_odd(Mmax,Nmax,fn_order_fijk_odd_lt2,Ncoeff_fijk_odd_lt2)
    allocate(stride_arr_qubic_odd_lt2(0:fn_order_fijk_odd_lt2-1,0:1))
    call get_stride_arr_odd(Mmax,Nmax,fn_order_fijk_odd_lt2,stride_arr_qubic_odd_lt2)
    !****************************************************************************************************************************************#
    !Quartic force constants
    !Even
    Mmax = 8;Nmax = 6
    call get_fn_order_even(Mmax,Nmax,fn_order_fijkl_even_gt5,Ncoeff_fijkl_even_gt5)
    allocate(stride_arr_quartic_even_gt5(0:fn_order_fijkl_even_gt5-1,0:1))
    call get_stride_arr_even(Mmax,Nmax,fn_order_fijkl_even_gt5,stride_arr_quartic_even_gt5)
    Mmax = 6;Nmax = 6
    call get_fn_order_even(Mmax,Nmax,fn_order_fijkl_even_lt5,Ncoeff_fijkl_even_lt5)
    allocate(stride_arr_quartic_even_lt5(0:fn_order_fijkl_even_lt5-1,0:1))
    call get_stride_arr_even(Mmax,Nmax,fn_order_fijkl_even_lt5,stride_arr_quartic_even_lt5)
    Mmax = 5;Nmax = 5
    call get_fn_order_even(Mmax,Nmax,fn_order_fijkl_even_lt2,Ncoeff_fijkl_even_lt2)
    allocate(stride_arr_quartic_even_lt2(0:fn_order_fijkl_even_lt2-1,0:1))
    call get_stride_arr_even(Mmax,Nmax,fn_order_fijkl_even_lt2,stride_arr_quartic_even_lt2)
    !Odd
    Mmax = 8;Nmax = 6
    call get_fn_order_odd(Mmax,Nmax,fn_order_fijkl_odd_gt5,Ncoeff_fijkl_odd_gt5)
    allocate(stride_arr_quartic_odd_gt5(0:fn_order_fijkl_odd_gt5-1,0:1))
    call get_stride_arr_odd(Mmax,Nmax,fn_order_fijkl_odd_gt5,stride_arr_quartic_odd_gt5)
    Mmax = 6;Nmax = 6
    call get_fn_order_odd(Mmax,Nmax,fn_order_fijkl_odd_lt5,Ncoeff_fijkl_odd_lt5)
    allocate(stride_arr_quartic_odd_lt5(0:fn_order_fijkl_odd_lt5-1,0:1))
    call get_stride_arr_odd(Mmax,Nmax,fn_order_fijkl_odd_lt5,stride_arr_quartic_odd_lt5)
    Mmax = 5;Nmax = 5
    call get_fn_order_odd(Mmax,Nmax,fn_order_fijkl_odd_lt2,Ncoeff_fijkl_odd_lt2)
    allocate(stride_arr_quartic_odd_lt2(0:fn_order_fijkl_odd_lt2-1,0:1))
    call get_stride_arr_odd(Mmax,Nmax,fn_order_fijkl_odd_lt2,stride_arr_quartic_odd_lt2)
    !****************************************************************************************************************************************#
    Mmax = 3;Nmax = 3
    call get_fn_order_odd(Mmax,Nmax,fn_order_qo_odd,Ncoeff_qo_odd)
    allocate(stride_arr_qo_odd(0:fn_order_qo_odd-1,0:1))
    call get_stride_arr_odd(Mmax,Nmax,fn_order_qo_odd,stride_arr_qo_odd)
    !**************************************************************************************************************************************#
    Mmax = 3;Nmax = 3
    call get_fn_order_even(Mmax,Nmax,fn_order_qo_even,Ncoeff_qo_even)
    allocate(stride_arr_qo_even(0:fn_order_qo_even-1,0:1))
    call get_stride_arr_even(Mmax,Nmax,fn_order_qo_even,stride_arr_qo_even)
    !****************************************************************************************************************************************#
    Mmax = 3;Nmax = 3
    call get_fn_order_even(Mmax,Nmax,fn_order_qo_ang_even,Ncoeff_qo_ang_even)
    allocate(stride_arr_qo_ang_even(0:fn_order_qo_ang_even-1,0:1))
    call get_stride_arr_even(Mmax,Nmax,fn_order_qo_ang_even,stride_arr_qo_ang_even)
    !****************************************************************************************************************************************#
    open(unit=101,file=trim(adjustl(potParamDir))//'/PotFitParams2FE/fitted_param_quad_8_6_even_fortran.dat',status='old')
    read(101,*)n_fcs_quad_even
    allocate(fitted_fc_coeff_ij_even(0:fn_order_fij_even-1,0:n_fcs_quad_even-1),fc_idx_quad_even(0:n_fcs_quad_even-1,0:3))
    open(unit=102,file=trim(adjustl(potParamDir))//'/PotFitParams2FE/fitted_param_qubic_gt5_8_6_even_fortran.dat',status='old')
    read(102,*)n_fcs_qubic_even_gt5
    allocate(fitted_fc_coeff_ijk_even_gt5(0:fn_order_fijk_even_gt5-1,0:n_fcs_qubic_even_gt5-1),fc_idx_qubic_even_gt5(0:n_fcs_qubic_even_gt5-1,0:3))
    open(unit=103,file=trim(adjustl(potParamDir))//'/PotFitParams2FE/fitted_param_qubic_lt5_6_6_even_fortran.dat',status='old')
    read(103,*)n_fcs_qubic_even_lt5
    allocate(fitted_fc_coeff_ijk_even_lt5(0:fn_order_fijk_even_lt5-1,0:n_fcs_qubic_even_lt5-1),fc_idx_qubic_even_lt5(0:n_fcs_qubic_even_lt5-1,0:3))
    open(unit=104,file=trim(adjustl(potParamDir))//'/PotFitParams2FE/fitted_param_qubic_lt2_5_5_even_fortran.dat',status='old')
    read(104,*)n_fcs_qubic_even_lt2
    allocate(fitted_fc_coeff_ijk_even_lt2(0:fn_order_fijk_even_lt2-1,0:n_fcs_qubic_even_lt2-1),fc_idx_qubic_even_lt2(0:n_fcs_qubic_even_lt2-1,0:3))
    open(unit=105,file=trim(adjustl(potParamDir))//'/PotFitParams2FE/fitted_param_quartic_gt5_8_6_even_fortran.dat',status='old')
    read(105,*)n_fcs_quartic_even_gt5
    allocate(fitted_fc_coeff_ijkl_even_gt5(0:fn_order_fijkl_even_gt5-1,0:n_fcs_quartic_even_gt5-1),fc_idx_quartic_even_gt5(0:n_fcs_quartic_even_gt5-1,0:3))
    open(unit=106,file=trim(adjustl(potParamDir))//'/PotFitParams2FE/fitted_param_quartic_lt5_6_6_even_fortran.dat',status='old')
    read(106,*)n_fcs_quartic_even_lt5
    allocate(fitted_fc_coeff_ijkl_even_lt5(0:fn_order_fijkl_even_lt5-1,0:n_fcs_quartic_even_lt5-1),fc_idx_quartic_even_lt5(0:n_fcs_quartic_even_lt5-1,0:3))
    open(unit=107,file=trim(adjustl(potParamDir))//'/PotFitParams2FE/fitted_param_quartic_lt2_5_5_even_fortran.dat',status='old')
    read(107,*)n_fcs_quartic_even_lt2
    allocate(fitted_fc_coeff_ijkl_even_lt2(0:fn_order_fijkl_even_lt2-1,0:n_fcs_quartic_even_lt2-1),fc_idx_quartic_even_lt2(0:n_fcs_quartic_even_lt2-1,0:3))
    open(unit=108,file=trim(adjustl(potParamDir))//'/PotFitParams2FE/fitted_param_quartic_ltpt1_even_fortran.dat',status='old')
    read(108,*)n_fcs_quartic_even_ltpt1
    allocate(fitted_fc_coeff_ijkl_even_ltpt1(0:n_fcs_quartic_even_ltpt1-1),fc_idx_quartic_even_ltpt1(0:n_fcs_quartic_even_ltpt1-1,0:3))
    call Read_fc_coeff(101,N_int,fn_order_fij_even,4,n_fcs_quad_even,fc_idx_quad_even,fitted_fc_coeff_ij_even)
    call Read_fc_coeff(102,N_int,fn_order_fijk_even_gt5,4,n_fcs_qubic_even_gt5,fc_idx_qubic_even_gt5,fitted_fc_coeff_ijk_even_gt5)
    call Read_fc_coeff(103,N_int,fn_order_fijk_even_lt5,4,n_fcs_qubic_even_lt5,fc_idx_qubic_even_lt5,fitted_fc_coeff_ijk_even_lt5)
    call Read_fc_coeff(104,N_int,fn_order_fijk_even_lt2,4,n_fcs_qubic_even_lt2,fc_idx_qubic_even_lt2,fitted_fc_coeff_ijk_even_lt2)
    call Read_fc_coeff(105,N_int,fn_order_fijkl_even_gt5,4,n_fcs_quartic_even_gt5,fc_idx_quartic_even_gt5,fitted_fc_coeff_ijkl_even_gt5)
    !print*,'potvars= ',fitted_fc_coeff_ijkl_even_gt5
    call Read_fc_coeff(106,N_int,fn_order_fijkl_even_lt5,4,n_fcs_quartic_even_lt5,fc_idx_quartic_even_lt5,fitted_fc_coeff_ijkl_even_lt5)
    call Read_fc_coeff(107,N_int,fn_order_fijkl_even_lt2,4,n_fcs_quartic_even_lt2,fc_idx_quartic_even_lt2,fitted_fc_coeff_ijkl_even_lt2)
    call Read_fc_coeff_ltpt1(108,N_int,4,n_fcs_quartic_even_ltpt1,fc_idx_quartic_even_ltpt1,fitted_fc_coeff_ijkl_even_ltpt1)

    open(unit=201,file=trim(adjustl(potParamDir))//'/PotFitParams2FE/fitted_param_quad_8_6_odd_fortran.dat',status='old')
    read(201,*)n_fcs_quad_odd
    allocate(fitted_fc_coeff_ij_odd(0:fn_order_fij_odd-1,0:n_fcs_quad_odd-1),fc_idx_quad_odd(0:n_fcs_quad_odd-1,0:3))

    open(unit=202,file=trim(adjustl(potParamDir))//'/PotFitParams2FE/fitted_param_qubic_gt5_8_6_odd_fortran.dat',status='old')
    read(202,*)n_fcs_qubic_odd_gt5
    allocate(fitted_fc_coeff_ijk_odd_gt5(0:fn_order_fijk_odd_gt5-1,0:n_fcs_qubic_odd_gt5-1),fc_idx_qubic_odd_gt5(0:n_fcs_qubic_odd_gt5-1,0:3))

    open(unit=203,file=trim(adjustl(potParamDir))//'/PotFitParams2FE/fitted_param_qubic_lt5_6_6_odd_fortran.dat',status='old')
    read(203,*)n_fcs_qubic_odd_lt5
    allocate(fitted_fc_coeff_ijk_odd_lt5(0:fn_order_fijk_odd_lt5-1,0:n_fcs_qubic_odd_lt5-1),fc_idx_qubic_odd_lt5(0:n_fcs_qubic_odd_lt5-1,0:3))

    open(unit=204,file=trim(adjustl(potParamDir))//'/PotFitParams2FE/fitted_param_qubic_lt2_5_5_odd_fortran.dat',status='old')
    read(204,*)n_fcs_qubic_odd_lt2
    allocate(fitted_fc_coeff_ijk_odd_lt2(0:fn_order_fijk_odd_lt2-1,0:n_fcs_qubic_odd_lt2-1),fc_idx_qubic_odd_lt2(0:n_fcs_qubic_odd_lt2-1,0:3))
    open(unit=205,file=trim(adjustl(potParamDir))//'/PotFitParams2FE/fitted_param_quartic_gt5_8_6_odd_fortran.dat',status='old')
    read(205,*)n_fcs_quartic_odd_gt5
    allocate(fitted_fc_coeff_ijkl_odd_gt5(0:fn_order_fijkl_odd_gt5-1,0:n_fcs_quartic_odd_gt5-1),fc_idx_quartic_odd_gt5(0:n_fcs_quartic_odd_gt5-1,0:3))
    open(unit=206,file=trim(adjustl(potParamDir))//'/PotFitParams2FE/fitted_param_quartic_lt5_6_6_odd_fortran.dat',status='old')
    read(206,*)n_fcs_quartic_odd_lt5
    allocate(fitted_fc_coeff_ijkl_odd_lt5(0:fn_order_fijkl_odd_lt5-1,0:n_fcs_quartic_odd_lt5-1),fc_idx_quartic_odd_lt5(0:n_fcs_quartic_odd_lt5-1,0:3))
    open(unit=207,file=trim(adjustl(potParamDir))//'/PotFitParams2FE/fitted_param_quartic_lt2_5_5_odd_fortran.dat',status='old')
    read(207,*)n_fcs_quartic_odd_lt2
    allocate(fitted_fc_coeff_ijkl_odd_lt2(0:fn_order_fijkl_odd_lt2-1,0:n_fcs_quartic_odd_lt2-1),fc_idx_quartic_odd_lt2(0:n_fcs_quartic_odd_lt2-1,0:3))
    open(unit=208,file=trim(adjustl(potParamDir))//'/PotFitParams2FE/fitted_param_quartic_ltpt1_odd_fortran.dat',status='old')
    read(208,*)n_fcs_quartic_odd_ltpt1
    allocate(fitted_fc_coeff_ijkl_odd_ltpt1(0:n_fcs_quartic_odd_ltpt1-1),fc_idx_quartic_odd_ltpt1(0:n_fcs_quartic_odd_ltpt1-1,0:3))

    call Read_fc_coeff(201,N_int,fn_order_fij_odd,4,n_fcs_quad_odd,fc_idx_quad_odd,fitted_fc_coeff_ij_odd)
    call Read_fc_coeff(202,N_int,fn_order_fijk_odd_gt5,4,n_fcs_qubic_odd_gt5,fc_idx_qubic_odd_gt5,fitted_fc_coeff_ijk_odd_gt5)
    call Read_fc_coeff(203,N_int,fn_order_fijk_odd_lt5,4,n_fcs_qubic_odd_lt5,fc_idx_qubic_odd_lt5,fitted_fc_coeff_ijk_odd_lt5)
    call Read_fc_coeff(204,N_int,fn_order_fijk_odd_lt2,4,n_fcs_qubic_odd_lt2,fc_idx_qubic_odd_lt2,fitted_fc_coeff_ijk_odd_lt2)
    call Read_fc_coeff(205,N_int,fn_order_fijkl_odd_gt5,4,n_fcs_quartic_odd_gt5,fc_idx_quartic_odd_gt5,fitted_fc_coeff_ijkl_odd_gt5)
    call Read_fc_coeff(206,N_int,fn_order_fijkl_odd_lt5,4,n_fcs_quartic_odd_lt5,fc_idx_quartic_odd_lt5,fitted_fc_coeff_ijkl_odd_lt5)
    call Read_fc_coeff(207,N_int,fn_order_fijkl_odd_lt2,4,n_fcs_quartic_odd_lt2,fc_idx_quartic_odd_lt2,fitted_fc_coeff_ijkl_odd_lt2)
    call Read_fc_coeff_ltpt1(208,N_int,4,n_fcs_quartic_odd_ltpt1,fc_idx_quartic_odd_ltpt1,fitted_fc_coeff_ijkl_odd_ltpt1)
    !****************************************************************************************#

    open(unit=109,file=trim(adjustl(potParamDir))//'/PotFitParams2FE/fitted_param_qo_3_3_even_fortran.dat',status='old')
    read(109,*)n_fcs_qo_even
    allocate(fitted_fc_coeff_qo_even(0:fn_order_qo_even-1,0:n_fcs_qo_even-1),fc_idx_qo_even(0:n_fcs_qo_even-1,0:7))


    open(unit=110,file=trim(adjustl(potParamDir))//'/PotFitParams2FE/fitted_param_qo_3_3_ang_even_fortran.dat',status='old')
    read(110,*)n_fcs_qo_ang_even
    allocate(fitted_fc_coeff_qo_ang_even(0:fn_order_qo_ang_even-1,0:n_fcs_qo_ang_even-1),fc_idx_qo_ang_even(0:n_fcs_qo_ang_even-1,0:7))

    open(unit=209,file=trim(adjustl(potParamDir))//'/PotFitParams2FE/fitted_param_qo_3_3_odd_fortran.dat',status='old')
    read(209,*)n_fcs_qo_odd
    allocate(fitted_fc_coeff_qo_odd(0:fn_order_qo_odd-1,0:n_fcs_qo_odd-1),fc_idx_qo_odd(0:n_fcs_qo_odd-1,0:7))

    call Read_fc_coeff(109,N_int,fn_order_qo_even,8,n_fcs_qo_even,fc_idx_qo_even,fitted_fc_coeff_qo_even)
    !print*,'inside =',n_fcs_qo_even,fitted_fc_coeff_qo_even
    call Read_fc_coeff(110,N_int,fn_order_qo_ang_even,8,n_fcs_qo_ang_even,fc_idx_qo_ang_even,fitted_fc_coeff_qo_ang_even)
    call Read_fc_coeff(209,N_int,fn_order_qo_odd,8,n_fcs_qo_odd,fc_idx_qo_odd,fitted_fc_coeff_qo_odd)



    !****************************************************************************************#
    !MEP Path fitted parameters
    allocate(Vt_coeff(0:fn_order_pot-1))
    call Read_Vt_coeff(trim(adjustl(potParamDir))//"/PotFitParams2FE/pot_par_8_6.dat",fn_order_pot,Vt_coeff )
    !Vt_coeff = In Hartee
    !Symmetrized internal coordinate fit
    open(unit=302,file=trim(adjustl(potParamDir))//'/PotFitParams2FE/eq_par_8_6_even_fortran.dat',status='old')
    read(302,*)n_eb,n_ea,n_ed
    allocate(fitted_bonds_coeff_even(0:fn_order_bonds_even-1,0:n_eb-1),fitted_angs_coeff_even(0:fn_order_angs_even-1,0:n_ea-1),fitted_dihs_coeff_even(0:fn_order_dihs_even-1,0:n_ed-1))
    call   Read_eq_int_coeff(302,n_eb,n_ea,n_ed,Nbonds,Nang,N_int,fn_order_bonds_even,fn_order_angs_even,fn_order_dihs_even,even_fn_idx,fitted_bonds_coeff_even,fitted_angs_coeff_even,fitted_dihs_coeff_even)
    open(unit=303,file=trim(adjustl(potParamDir))//'/PotFitParams2FE/eq_par_8_6_odd_fortran.dat',status='old')
    read(303,*)n_ob,n_oa,n_od
    allocate(fitted_bonds_coeff_odd(0:fn_order_bonds_odd-1,0:n_eb-1),fitted_angs_coeff_odd(0:fn_order_angs_odd-1,0:n_ea-1),fitted_dihs_coeff_odd(0:fn_order_dihs_odd-1,0:n_ed-1))
    call  Read_eq_int_coeff(303,n_ob,n_oa,n_od,Nbonds,Nang,N_int,fn_order_bonds_odd,fn_order_angs_odd,fn_order_dihs_odd,odd_fn_idx,fitted_bonds_coeff_odd,fitted_angs_coeff_odd,fitted_dihs_coeff_odd)
    !Fitted_b_coeff in Hartee/angstrom2(or degree2 or angstrom*degree)  
    allocate(atom_ind(0:Natoms-1,0:3))
    call Read_atom_ind("atom_ind.dat",Natoms,atom_ind)
    n_fix_tors = 2 
    allocate(max_id(0:n_fix_tors-1))	
    max_id(0) = 8;max_id(1) = 6
    sym_mat = 0.d0
    ctr = 0
    i = 0
    do while (i<N_int)
        if  (any(sym_idx==i)) then
            id1 = sym_idx(ctr);id2 = sym_idx(ctr+1)
            sym_mat(id1,id1) =  sqrt_2_inv
            sym_mat(id1,id2) =  sqrt_2_inv
            sym_mat(id2,id1) =  sqrt_2_inv
            sym_mat(id2,id2) = -sqrt_2_inv
            i = i + 2
            ctr = ctr + 2
        else
            sym_mat(i,i) = 1.d0
            i = i + 1
    endif
    enddo
    sym_mat_transpose = transpose(sym_mat)
end subroutine 
!****************************************************************************************#
! Reads the fitted force constant coefficients into  pot_coeff array.
! pot_coeff is an 5 dimensional array, first four corresponds to internal coordinates,
! and the fifth corresponds to order of the term.
subroutine Read_fc_coeff(f_n,N_int,fn_order,order,n_params,fc_idx,fitted_coeff)
    implicit none
    integer::i,j,k,ctr,idum
    integer,intent(in)::f_n,n_params,N_int,fn_order,order
    integer,intent(out),allocatable::fc_idx(:,:)
    real(8),intent(out),allocatable::fitted_coeff(:,:)
    !Input file has 5 columns. First four columns contains the index of the internal coordinates. 
    !and the last column has the corresponding fitted parameters. fc_idx stores the unique indexes 
    !and the corresponding row of fitted_coeff contains fitted coefficients for that index values
    allocate(fc_idx(0:n_params-1,0:order-1),fitted_coeff(0:fn_order-1,0:n_params-1))
    !n_params: no of coefficients in potential expansion
    fc_idx = 0;fitted_coeff = 0.d0
    ctr = 0
    do i = 1,fn_order * n_params,fn_order
        do k = 0,fn_order-1
            read(f_n,*)(fc_idx(ctr,j),j=0,order-1),idum,fitted_coeff(k,ctr)
        enddo
        ctr = ctr + 1
    enddo
    fc_idx  = fc_idx-1
    close(f_n)
endsubroutine 
!****************************************************************************************#
! Reads the fitted force constant coefficients into  pot_coeff array.
! pot_coeff is an 5 dimensional array, first four corresponds to internal coordinates,
! and the fifth corresponds to order of the term.
subroutine Read_fc_coeff_ltpt1(f_n,N_int,order,n_params,fc_idx,fitted_coeff)
    implicit none
    integer::i,j,k,ctr,idum
    integer,intent(in)::f_n,n_params,N_int,order
    integer,intent(out),allocatable::fc_idx(:,:)
    real(8),intent(out),allocatable::fitted_coeff(:)
    !Input file has 5 columns. First four columns contains the index of the internal coordinates. 
    !and the last column has the corresponding fitted parameters. fc_idx stores the unique indexes 
    !and the corresponding row of fitted_coeff contains fitted coefficients for that index values
    allocate(fc_idx(0:n_params-1,0:order-1),fitted_coeff(0:n_params-1))
    !n_params: no of coefficients in potential expansion
    fc_idx = 0;fitted_coeff = 0.d0
    do i = 0,n_params-1
        read(f_n,*)(fc_idx(i,j),j=0,order-1),idum,fitted_coeff(i)
    enddo
    fc_idx  = fc_idx-1
    close(f_n)
endsubroutine 
!****************************************************************************************#
subroutine  Read_Vt_coeff(f_name,fn_order,Vt_coeff)
    implicit none
    integer::i,idum1,idum2,idx
    real(8)::val
    character(len=*)::f_name
    integer,intent(in)::fn_order
    real(8),intent(out),dimension(0:fn_order-1)::Vt_coeff
    open(unit=2222,file=f_name,status="old")
    do i = 0,fn_order-1
        read(2222,*)idx,idum1,idum2,val
        if (i/=idx) then
            print*,"Error reading file ",f_name
            exit
        else
            Vt_coeff(idx)    = val
        endif
    enddo
    close(2222)
endsubroutine Read_Vt_coeff
!****************************************************************************************#
subroutine Read_eq_int_coeff(f_n,nb,na,nd,Nbonds,Nang,N_int,fn_order_bonds,fn_order_angs,fn_order_dihs,coordidx,fitted_bond_coeff,fitted_ang_coeff,fitted_dih_coeff)
    implicit none
    integer::i,j,idum,a0,a1,idx1,idx2,odd_b,odd_a,odd_d,b_ctr,a_ctr,d_ctr,n_params
    real(8)::val
    integer,intent(in) :: f_n,nb,na,nd,Nbonds,Nang,N_int,fn_order_bonds,fn_order_angs,fn_order_dihs,coordidx(*)
    !Dimension coordidx(*)
    real(8),intent(out),allocatable::fitted_bond_coeff(:,:),fitted_ang_coeff(:,:),fitted_dih_coeff(:,:)
    allocate(fitted_bond_coeff(0:fn_order_bonds-1,0:nb-1))
    allocate( fitted_ang_coeff(0:fn_order_angs-1,0:na-1))
    allocate( fitted_dih_coeff(0:fn_order_dihs-1,0:nd-1))
    fitted_bond_coeff = 0.d0;fitted_ang_coeff=0.d0;fitted_dih_coeff=0.d0
    n_params = nb + na  + nd 
    a_ctr = 0;b_ctr =0 ;d_ctr=0
    do i = 0,n_params-1
        read(f_n,*)a0,idum,val
        if (a0<Nbonds) then
        fitted_bond_coeff(0,b_ctr) = val
            do j = 0,fn_order_bonds-2
                read(f_n,*)a1,idum,val
                if (a1==a0) then
                    fitted_bond_coeff(j+1,b_ctr) = val
                else
                    exit
                endif
            enddo
            b_ctr = b_ctr + 1
        else  if ((a0>=Nbonds).and.(a0<Nbonds+Nang)) then
            fitted_ang_coeff(0,a_ctr) = val
            do j = 0,fn_order_angs-2
                read(f_n,*)a1,idum,val
                if (a1==a0) then
                    fitted_ang_coeff(j+1,a_ctr) = val
                else
                    exit
                endif
            enddo
            a_ctr = a_ctr + 1
        else
            fitted_dih_coeff(0,d_ctr) = val
            do j = 0,fn_order_dihs-2
                read(f_n,*)a1,idum,val
                if (a1==a0) then
                    fitted_dih_coeff(j+1,d_ctr) = val
                else
                    exit
                endif
            enddo
            d_ctr = d_ctr + 1
        endif
    enddo
    close(f_n)
end subroutine Read_eq_int_coeff
!****************************************************************************************#
subroutine Read_atom_ind(f_name,Natoms,atom_ind)
    implicit none
    integer::i,j
    character(len=*)::f_name
    integer,intent(in)::Natoms
    integer,intent(out),dimension(0:Natoms-1,0:3)::atom_ind
    open(unit=12,file=f_name,status="old")
    atom_ind = 0
    do i = 0,Natoms-1
        read(12,*)(atom_ind(i,j),j=0,3)
    enddo
endsubroutine Read_atom_ind
!****************************************************************************************#
subroutine get_stride_arr_even(Mmax,Nmax,fn_order,stride_arr)
    implicit none
    integer,intent(in):: Mmax,Nmax,fn_order
    integer,intent(Out)::stride_arr(0:fn_order-1,0:1)
    integer::ctr,i,j
    stride_arr = 0
    ctr=0
    do i = 0,Mmax
        do j = 0,Nmax
            stride_arr(ctr,:) = [i,j]
            ctr = ctr + 1
        enddo
    enddo
    do i = 1,Mmax
        do j = 1,Nmax
            stride_arr(ctr,:) = [i,j]
            ctr = ctr + 1
        enddo
    enddo
endsubroutine
!*****************************************************************************************#
subroutine get_stride_arr_odd(Mmax,Nmax,fn_order,stride_arr)
    implicit none
    integer,intent(in):: Mmax,Nmax,fn_order
    integer,intent(out)::stride_arr(0:fn_order-1,0:1)
    integer::ctr,i,j
    stride_arr = 0
    ctr=0
    do i = 0,Mmax
        do j = 1,Nmax
            stride_arr(ctr,:) = [i,j]
            ctr = ctr + 1
        enddo
    enddo
    do i = 1,Mmax
        do j = 0,Nmax
            stride_arr(ctr,:) = [i,j]
            ctr = ctr + 1
        enddo
    enddo
endsubroutine
!*****************************************************************************************#
subroutine get_fn_order_even(Mmax,Nmax,fn_order,Ncoeff)
    implicit none
    integer,intent(in)::Mmax,Nmax
    integer,intent(out)::fn_order,Ncoeff(0:1)
    integer::i,j,ctr1,ctr2
    fn_order = 0
    ctr1 =0;ctr2=0;
    do i = 0,Mmax
        do j = 0,Nmax
            fn_order = fn_order + 1
            ctr1 = ctr1 + 1
            if ((i>0).and.(j>0)) then
                fn_order = fn_order + 1
                ctr2 = ctr2 + 1
            endif
        enddo
    enddo
    Ncoeff(0) = ctr1
    Ncoeff(1) = ctr2
endsubroutine
!****************************************************************************************#
subroutine get_fn_order_odd(Mmax,Nmax,fn_order,Ncoeff)
    implicit none
    integer,intent(in)::Mmax,Nmax
    integer,intent(out)::fn_order,Ncoeff(0:1)
    integer::i,j,ctr1,ctr2
    fn_order = 0
    ctr1 =0;ctr2=0;
    do i = 0,Mmax
        do j = 0,Nmax
            if (j>0) then
                fn_order = fn_order + 1
                ctr1 = ctr1 + 1
            endif
            if (i>0) then
                fn_order = fn_order + 1
                ctr2 = ctr2 + 1
            endif
        enddo
    enddo
    Ncoeff(0) = ctr1
    Ncoeff(1) = ctr2
endsubroutine
!****************************************************************************************#
end module
