!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief routines printing various data related to POD method
!> \author Z. Futera (03.2019)
! **************************************************************************************************
MODULE et_coupling_pod_print

   USE cp_blacs_types,                  ONLY: cp_blacs_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE et_coupling_pod_mo,              ONLY: pod_mo_sum_c2
   USE et_coupling_pod_types,           ONLY: et_pod
   USE input_section_types,             ONLY: section_get_lval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind_set,&
                                              qs_kind_type
   USE qs_mo_types,                     ONLY: mo_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'et_coupling_pod_print'

   PUBLIC :: pod_print_couplings, &
             pod_print_state_fracs, &
             pod_print_states

CONTAINS

! **************************************************************************************************
!> \brief Print out specific MO coefficients
!> \param output_unit unit number of the open output stream
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param blk atomic-block ID
!> \param n_spins number of spin components
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE pod_print_state_fracs(output_unit, qs_env, ec, blk, n_spins)

      ! Routine arguments
      INTEGER                                            :: output_unit
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_pod), POINTER                              :: ec
      INTEGER                                            :: blk, n_spins

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_print_state_fracs', &
         routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: j, k, l, m, n, n_ao, n_mo, handle
      INTEGER, DIMENSION(:), POINTER                     :: list_at, list_mo
      REAL(KIND=dp)                                      :: c1, c2
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mat_w
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: block_sec, print_sec

      ! Time mark
      CALL timeset(routineN, handle)

      ! Initialization
      NULLIFY (block_sec)
      NULLIFY (print_sec)
      NULLIFY (qs_kind_set)

      ! Atomic block data
      block_sec => section_vals_get_subs_vals(qs_env%input, &
                                              'PROPERTIES%ET_COUPLING%PROJECTION%BLOCK')

      print_sec => section_vals_get_subs_vals(block_sec, 'PRINT', i_rep_section=blk)

      ! List of atoms
      CALL section_vals_val_get(print_sec, keyword_name='MO_COEFF_ATOM', n_rep_val=n)

      IF (n > 0) THEN

         IF (output_unit > 0) &
            WRITE (output_unit, '(/,T3,A/)') 'Block state fractions:'

         ! Number of AO functions
         CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
         CALL get_qs_kind_set(qs_kind_set, nsgf=n_ao)

         ! MOs in orthonormal basis set
         ALLOCATE (mat_w(n_spins))
         CPASSERT(ASSOCIATED(mat_w))
         DO j = 1, n_spins
            n_mo = ec%block(blk)%n_ao
            CALL cp_fm_create(matrix=mat_w(j), &
                              matrix_struct=ec%block(blk)%mo(j)%mo_coeff%matrix_struct, &
                              name='BLOCK MOs IN ORTHONORMAL BASIS SET')
            CALL parallel_gemm("N", "N", n_ao, n_mo, n_ao, 1.0_dp, ec%m_transf_inv, &
                               ec%block(blk)%mo(j)%mo_coeff, 0.0_dp, mat_w(j))
         END DO

         DO j = 1, n
            NULLIFY (list_at)
            CALL section_vals_val_get(print_sec, keyword_name='MO_COEFF_ATOM', &
                                      i_rep_val=j, i_vals=list_at)
            IF (ASSOCIATED(list_at)) THEN

               ! List of states
               CALL section_vals_val_get(print_sec, keyword_name='MO_COEFF_ATOM_STATE', n_rep_val=m)

               IF (m > 0) THEN

                  DO k = 1, m
                     NULLIFY (list_mo)
                     CALL section_vals_val_get(print_sec, keyword_name='MO_COEFF_ATOM_STATE', &
                                               i_rep_val=k, i_vals=list_mo)
                     IF (ASSOCIATED(list_mo)) THEN

                        IF (j > 1) THEN
                           IF (output_unit > 0) &
                              WRITE (output_unit, *)
                        END IF

                        DO l = 1, SIZE(list_mo)

                           IF (n_spins > 1) THEN
                              c1 = pod_mo_sum_c2(ec%block(blk)%atom, mat_w(1), list_mo(l), list_at)
                              c2 = pod_mo_sum_c2(ec%block(blk)%atom, mat_w(2), list_mo(l), list_at)
                              IF (output_unit > 0) &
                                 WRITE (output_unit, '(I5,A,I5,2F20.10)') j, ' /', list_mo(l), c1, c2
                           ELSE
                              c1 = pod_mo_sum_c2(ec%block(blk)%atom, mat_w(1), list_mo(l), list_at)
                              IF (output_unit > 0) &
                                 WRITE (output_unit, '(I5,A,I5,F20.10)') j, ' /', list_mo(l), c1
                           END IF

                        END DO

                     END IF
                  END DO

               END IF

            END IF
         END DO

         ! Clean memory
         CALL cp_fm_release(mat_w)

      END IF

      ! Running time
      CALL timestop(handle)

   END SUBROUTINE pod_print_state_fracs

! **************************************************************************************************
!> \brief Print out electronic states (MOs)
!> \param output_unit unit number of the open output stream
!> \param mo array of MO sets
!> \param n_spins number of spin components
!> \param label output label
!> \param mx_mo_a maximum number of alpha states to print out
!> \param mx_mo_b maximum number of beta states to print out
!> \param fermi print out Fermi level and number of electrons
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE pod_print_states(output_unit, mo, n_spins, label, mx_mo_a, mx_mo_b, fermi)

      ! Routine arguments
      INTEGER                                            :: output_unit
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mo
      INTEGER                                            :: n_spins
      CHARACTER(LEN=*)                                   :: label
      INTEGER, OPTIONAL                                  :: mx_mo_a, mx_mo_b
      LOGICAL, OPTIONAL                                  :: fermi

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_print_states', &
         routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: i, mx_a, mx_b, n, handle
      LOGICAL                                            :: prnt_fm

      ! Time mark
      CALL timeset(routineN, handle)

      prnt_fm = .FALSE.
      IF (PRESENT(fermi)) &
         prnt_fm = fermi

      IF (output_unit > 0) THEN

         WRITE (output_unit, '(/,T3,A/)') 'State energies ('//TRIM(ADJUSTL(label))//'):'

         ! Spin-polarized calculation
         IF (n_spins > 1) THEN

            mx_a = mo(1)%nmo
            IF (PRESENT(mx_mo_a)) &
               mx_a = MIN(mo(1)%nmo, mx_mo_a)
            mx_b = mo(2)%nmo
            IF (PRESENT(mx_mo_b)) &
               mx_b = MIN(mo(2)%nmo, mx_mo_b)
            n = MAX(mx_a, mx_b)

            DO i = 1, n
               WRITE (output_unit, '(T3,I10)', ADVANCE='no') i
               IF (i <= mx_a) THEN
                  WRITE (output_unit, '(2F12.4)', ADVANCE='no') &
                     mo(1)%occupation_numbers(i), mo(1)%eigenvalues(i)
               ELSE
                  WRITE (output_unit, '(A)', ADVANCE='no') '                        '
               END IF
               WRITE (output_unit, '(A)', ADVANCE='no') '     '
               IF (i <= mx_b) THEN
                  WRITE (output_unit, '(2F12.4)') &
                     mo(2)%occupation_numbers(i), mo(2)%eigenvalues(i)
               ELSE
                  WRITE (output_unit, *)
               END IF
            END DO

            IF (prnt_fm) THEN
               WRITE (output_unit, '(/,T3,I10,F24.4,I10,F19.4)') &
                  mo(1)%nelectron, mo(1)%mu, &
                  mo(2)%nelectron, mo(2)%mu
            END IF

            ! Spin-restricted calculation
         ELSE

            mx_a = mo(1)%nmo
            IF (PRESENT(mx_mo_a)) &
               mx_a = MIN(mo(1)%nmo, mx_mo_a)

            DO i = 1, mx_a
               WRITE (output_unit, '(T3,I10,2F12.4)') &
                  i, mo(1)%occupation_numbers(i), mo(1)%eigenvalues(i)
            END DO

            IF (prnt_fm) THEN
               WRITE (output_unit, '(/,T3,I10,F24.4)') &
                  mo(1)%nelectron, mo(1)%mu
            END IF

         END IF

      END IF

      ! Running time
      CALL timestop(handle)

   END SUBROUTINE pod_print_states

! **************************************************************************************************
!> \brief Print out donor-acceptor state couplings
!> \param output_unit unit number of the open output stream
!> \param logger output logger
!> \param ec_sec electronic-coupling input-file section
!> \param ec electronic coupling data structure
!> \param mo MO coefficients
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE pod_print_couplings(output_unit, logger, ec_sec, ec, mo)

      IMPLICIT NONE

      ! Routine arguments
      INTEGER                                            :: output_unit
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: ec_sec
      TYPE(et_pod), POINTER                              :: ec
      TYPE(mo_set_type), DIMENSION(:), INTENT(IN)        :: mo

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_print_couplings', &
                                     routineP = moduleN//':'//routineN

      ! Local variables
      CHARACTER(LEN=default_path_length)                 :: filename, my_pos, title
      LOGICAL                                            :: do_print, has_cpl
      INTEGER                                            :: i, j, k, l, id, handle, n_spins, unit_nr
      INTEGER                                            :: n_states(2)
      TYPE(section_vals_type), POINTER                   :: print_key
#if defined(__SCALAPACK)
     
      LOGICAL                                            :: blank_line
      INTEGER                                            :: tag_m1, tag_m2
      INTEGER                                            :: n, ib, ic, ir, ix, info
      INTEGER                                            :: ir_l_m1, ic_l_m1, ir_l_m2, ic_l_m2
      INTEGER                                            :: pr_m1, pc_m1, pr_m2, pc_m2
      INTEGER                                            :: nr_m1, nr_m2, nc_m1, nc_m2
      INTEGER                                            :: nr_b_m1, nr_b_m2
      INTEGER, DIMENSION(9)                              :: desc_m1, desc_m2
      REAL(KIND=dp), DIMENSION(:), POINTER               :: v1, v2
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: d1, d2
      TYPE(cp_fm_type), POINTER                          :: m1, m2
      TYPE(cp_blacs_type)                                :: cntxt_m1, cntxt_m2
      TYPE(mp_para_env_type), POINTER                    :: para_m1, para_m2

      INTEGER, EXTERNAL                                  :: numroc
#endif

      ! Time mark
      CALL timeset(routineN, handle)

#if defined(__SCALAPACK)
      NULLIFY (d1)
      NULLIFY (d2)
      NULLIFY (v1)
      NULLIFY (v2)

      tag_m1 = 0
      tag_m2 = 0
#endif

      n_states = 0
      DO i = 1, SIZE(mo)
         n_states(i) = mo(i)%nmo
      END DO
      n_spins = 1
      IF (n_states(2) > 0) n_spins = 2

      print_key => section_vals_get_subs_vals(section_vals=ec_sec, &
                                              subsection_name="PRINT%COUPLING_ELEMENTS")
      do_print = BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)

      ! Check whether we have any coupling values to print
      has_cpl = .FALSE.
      DO i = 1, ec%n_blocks
        IF (ASSOCIATED(ec%block(i)%hab)) THEN
           DO j = 1, ec%block(i)%n_hab_blocks
             IF (ASSOCIATED(ec%block(i)%hab(1,j)%local_data)) THEN
               has_cpl = .TRUE.
               EXIT
             END IF
           END DO
        END IF
      END DO

      IF (do_print .AND. has_cpl) THEN

         my_pos = "REWIND"
         IF (section_get_lval(print_key, "APPEND")) &
            my_pos = "APPEND"

         ! Header
         IF (output_unit > 0) &
            WRITE (output_unit, '(/,T3,A/)') 'Printing coupling elements to output files'

#if defined(__SCALAPACK)

         DO i = 1, ec%n_blocks
            id = 1
            DO j = i + 1, ec%n_blocks

               ! Alpha states
               m1 => ec%block(i)%hab(1, id)
               CALL cp_fm_get_info(m1, nrow_global=nr_m1, ncol_global=nc_m1, &
                                   nrow_block=nr_b_m1, para_env=para_m1)
               CALL cntxt_m1%gridinit(para_m1,'R',1, para_m1%num_pe)
               CALL descinit(desc_m1, nr_m1, nc_m1, nr_b_m1, nc_m1, 0, 0, &
                             cntxt_m1%get_handle(), nr_m1, info)
               CPASSERT(info == 0)
               n = numroc(nr_m1, nr_b_m1, cntxt_m1%mepos(1), 0, cntxt_m1%num_pe(1))
               ALLOCATE (d1(MAX(1, n), nc_m1))
               CALL pdgemr2d(nr_m1, nc_m1, m1%local_data(1, 1), 1, 1, &
                             m1%matrix_struct%descriptor, d1, 1, 1, desc_m1, cntxt_m1%get_handle())
               ALLOCATE (v1(nr_b_m1*nc_m1))
               v1 = HUGE(1.0_dp)

               ! Beta states
               IF (n_spins > 1) THEN
                  m2 => ec%block(i)%hab(2, id)
                  CALL cp_fm_get_info(m2, nrow_global=nr_m2, ncol_global=nc_m2, &
                                      nrow_block=nr_b_m2, para_env=para_m2)
                  CALL cntxt_m2%gridinit(para_m2,'R',1, para_m2%num_pe)
                  CALL descinit(desc_m2, nr_m2, nc_m2, nr_b_m2, nc_m2, 0, 0, &
                                cntxt_m2%get_handle(), nr_m2, info)
                  CPASSERT(info == 0)
                  n = numroc(nr_m2, nr_b_m2, cntxt_m2%mepos(1), 0, cntxt_m2%num_pe(1))
                  ALLOCATE (d2(MAX(1, n), nc_m2))
                  CALL pdgemr2d(nr_m2, nc_m2, m2%local_data(1, 1), 1, 1, &
                                m2%matrix_struct%descriptor, d2, 1, 1, desc_m2, cntxt_m2%get_handle())
                  ALLOCATE (v2(nr_b_m2*nc_m2))
                  v2 = HUGE(1.0_dp)
               END IF

               ! Open the output file
               IF (output_unit > 0) THEN

                  WRITE (filename, '(a5,I1.1,a1,I1.1)') "ET_BL_", i, "-", j
                  unit_nr = cp_print_key_unit_nr(logger, ec_sec, "PRINT%COUPLING_ELEMENTS", &
                                                 extension=".elcoup", middle_name=TRIM(filename), &
                                                 file_position=my_pos, log_filename=.FALSE.)
   
                  WRITE (title, *) 'Coupling elements [meV] between blocks:', i, j
   
                  WRITE (unit_nr, *) TRIM(title)
                  IF (n_spins > 1) THEN
                     WRITE (unit_nr, '(T3,A8,T13,A8,T28,A,A)') 'State A', 'State B', &
                           'Coupling spin 1', '  Coupling spin 2'
                  ELSE
                     WRITE (unit_nr, '(T3,A8,T13,A8,T28,A)') 'State A', 'State B', 'Coupling'
                  END IF

               END IF

               ! Row / column index
               k = 1
               l = 1

               ! Loop over blocks of rows
               IF (n_spins > 1) THEN
                  CPASSERT(nr_m1 == nr_m2)
                  CPASSERT(nc_m1 == nc_m2)
                  CPASSERT(nr_b_m1 == nr_b_m2)
               END IF

               DO ir = 1, nr_m1, MAX(nr_b_m1, 1)
   
                  ! Number of rows in the block
                  ib = MIN(nr_b_m1, nr_m1 - ir + 1)
   
                  ! Transfer alpha-spin data
                  CALL infog2l(ir, 1, desc_m1, cntxt_m1%num_pe(1), cntxt_m1%num_pe(2), &
                               cntxt_m1%mepos(1), cntxt_m1%mepos(2), ir_l_m1, ic_l_m1, pr_m1, pc_m1)
                  IF (pr_m1 == cntxt_m1%mepos(1)) THEN
                     DO ix = 1, ib
                        v1((ix - 1)*nc_m1 + 1:nc_m1*ix) = d1(ir_l_m1 + ix - 1, :)
                     END DO
                  END IF
                  IF (pr_m1 /= 0) THEN
                     IF (pr_m1 == cntxt_m1%mepos(1)) &
                        CALL para_m1%send(v1(:), 0, tag_m1)
                     IF (cntxt_m1%mepos(1) == 0) &
                        CALL para_m1%recv(v1(:), pr_m1, tag_m1)
                  END IF
 
                  ! Transfer beta-spin data
                  IF (n_spins > 1) THEN
                     CALL infog2l(ir, 1, desc_m2, cntxt_m2%num_pe(1), cntxt_m2%num_pe(2), &
                                  cntxt_m2%mepos(1), cntxt_m2%mepos(2), ir_l_m2, ic_l_m2, pr_m2, pc_m2)
                     IF (pr_m2 == cntxt_m2%mepos(1)) THEN
                        DO ix = 1, ib
                           v2((ix - 1)*nc_m2 + 1:nc_m2*ix) = d2(ir_l_m2 + ix - 1, :)
                        END DO
                     END IF
                     IF (pr_m2 /= 0) THEN
                        IF (pr_m2 == cntxt_m2%mepos(1)) &
                           CALL para_m2%send(v2(:), 0, tag_m2)
                        IF (cntxt_m2%mepos(1) == 0) &
                           CALL para_m2%recv(v2(:), pr_m2, tag_m2)
                     END IF
                  END IF

                  ! Output
                  IF (output_unit > 0) THEN

                     DO ix = 1, ib
                        blank_line = .FALSE.
                        DO ic = (ix - 1)*nc_m1 + 1, nc_m1*ix
   
                           IF (k <= MIN(ec%block(i)%n_ao, n_states(1)) .AND. &
                               l <= MIN(ec%block(j)%n_ao, n_states(1))) THEN
                              ! alpha spin
                              WRITE (unit_nr, '(T3,I5,T13,I5,T22,E20.6)', ADVANCE='no') &
                                 k, l, v1(ic)*evolt*1000.0_dp
                              ! beta spin
                              IF (n_spins > 1 .AND. k <= n_states(2) .AND. l <= n_states(2)) &
                                 WRITE (unit_nr, '(E20.6)', ADVANCE='no') v2(ic)*evolt*1000.0_dp
   
                              WRITE (unit_nr, *)
                              blank_line = .TRUE.
                           END IF
   
                           ! Row / column IDs
                           l = l + 1
                           IF (l > nc_m1) THEN
                              l = 1
                              k = k + 1
                           END IF
   
                        END DO

                        IF (blank_line) &
                           WRITE (unit_nr, *)
                     END DO

                  END IF

               END DO

               ! close the output file
               IF (output_unit > 0) &
                 CALL cp_print_key_finished_output(unit_nr, logger, ec_sec, "PRINT%COUPLING_ELEMENTS")

               ! clean memory
               call cntxt_m1%gridexit()
               IF (n_spins > 1) &
                  call cntxt_m2%gridexit()
   
               IF (ASSOCIATED(d1)) THEN
                  DEALLOCATE (d1)
                  NULLIFY (d1)
               END IF
               IF (ASSOCIATED(d2)) THEN
                  DEALLOCATE (d2)
                  NULLIFY (d2)
               END IF
               IF (ASSOCIATED(v1)) THEN
                  DEALLOCATE (v1)
                  NULLIFY (v1)
               END IF
               IF (ASSOCIATED(v2)) THEN
                  DEALLOCATE (v2)
                  NULLIFY (v2)
               END IF

               ! block index
               id = id + 1

            END DO
         END DO

#else

         IF (output_unit > 0) THEN

            DO i = 1, ec%n_blocks
               id = 1
               DO j = i + 1, ec%n_blocks

                  WRITE (filename, '(A5,I1.1,A1,I1.1)') "ET_BL_", i, "-", j
                  unit_nr = cp_print_key_unit_nr(logger, ec_sec, "PRINT%COUPLING_ELEMENTS", &
                                        extension=".elcoup", middle_name=TRIM(filename), &
                                        file_position=my_pos, log_filename=.FALSE.)
                  WRITE (title, *) 'Coupling elements [meV] between blocks:', i, j
                  WRITE (unit_nr, *) TRIM(title)
                  IF (n_spins > 1) THEN
                     WRITE (unit_nr, '(T3,A8,T13,A8,T28,A,A)') 'State A', 'State B', &
                            'Coupling spin 1', '  Coupling spin 2'
                  ELSE
                     WRITE (unit_nr, '(T3,A8,T13,A8,T28,A)') 'State A', 'State B', 'Coupling'
                  END IF

                  DO k = 1, MIN(ec%block(i)%n_ao, n_states(1))
                     DO l = 1, MIN(ec%block(j)%n_ao, n_states(1))
   
                        IF (n_spins > 1) THEN

                           WRITE (unit_nr, '(T3,I5,T13,I5,T22,E20.6)', ADVANCE='no') k, l, &
                              ec%block(i)%hab(1, id)%local_data(k, l)*evolt*1000.0_dp
                           IF ((k <= n_states(2)) .AND. (l <= n_states(2))) THEN
                              WRITE (unit_nr, '(E20.6)') &
                                     ec%block(i)%hab(2, id)%local_data(k, l)*evolt*1000.0_dp
                           ELSE
                              WRITE (unit_nr, *)
                           END IF
   
                        ELSE
   
                           WRITE (unit_nr, '(T3,I5,T13,I5,T22,E20.6)') k, l, &
                                  ec%block(i)%hab(1, id)%local_data(k, l)*evolt*1000.0_dp
   
                        END IF
   
                     END DO
                     
                     WRITE (unit_nr, *)
                  END DO

                  CALL cp_print_key_finished_output(unit_nr, logger, &
                         ec_sec, "PRINT%COUPLING_ELEMENTS")

                  id = id + 1

               END DO
            END DO

         END IF

#endif

      END IF

      ! Running time
      CALL timestop(handle)

   END SUBROUTINE pod_print_couplings

END MODULE et_coupling_pod_print
