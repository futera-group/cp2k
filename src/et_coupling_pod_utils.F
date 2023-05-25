!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief different utilities used by POD method routines
!> \author Z. Futera (03.2019)
! **************************************************************************************************
MODULE et_coupling_pod_utils

   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_diag,                      ONLY: cp_fm_power
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_generate_filename
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE et_coupling_pod_save,            ONLY: pod_save_transf_mat
   USE et_coupling_pod_types,           ONLY: et_pod
   USE input_section_types,             ONLY: section_get_lval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'et_coupling_pod_utils'

   PUBLIC :: pod_check_restart, &
             pod_get_transf_mat

CONTAINS

! **************************************************************************************************
!> \brief calculated transformation matrix for basis-set orthogonalization (S^{-1/2})
!> \param output_unit unit number of the open output stream
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param mat_w working matrix of the same dimension
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE pod_get_transf_mat(output_unit, qs_env, ec, mat_w)

      ! Routine arguments
      INTEGER                                            :: output_unit
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_pod), POINTER                              :: ec
      TYPE(cp_fm_type)                                   :: mat_w

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_get_transf_mat', &
                                     routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: handle, n_deps
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_s
      TYPE(scf_control_type), POINTER                    :: scf_cntrl

      print *,'DEBUG: pod_get_transf_mat'

      ! Time mark
      CALL timeset(routineN, handle)

      ! Initialization
      NULLIFY (mat_s)
      NULLIFY (scf_cntrl)

      print *,'DEBUG: --> call get_qs_env'
      CALL get_qs_env(qs_env, matrix_s=mat_s, scf_control=scf_cntrl)

      ! Transformation S -> S^{-1/2}
      IF (output_unit > 0) &
         WRITE (output_unit, '(/,T3,A)') &
         'Constructing forward basis-set transformation matrix'

      CPASSERT(.NOT. ASSOCIATED(ec%m_transf))
      print *,'DEBUG: --> allocate m_transf'
      ALLOCATE (ec%m_transf)
      CPASSERT(ASSOCIATED(ec%m_transf))
      print *,'DEBUG: --> call cp_fm_create'
      CALL cp_fm_create(matrix=ec%m_transf, matrix_struct=mat_w%matrix_struct, &
                        name='S^(-1/2) TRANSFORMATION MATRIX')
      print *,'DEBUG: --> call copy_dbcsr_to_fm'
      CALL copy_dbcsr_to_fm(mat_s(1)%matrix, ec%m_transf)
      print *,'DEBUG: --> call cp_fm_power'
      CALL cp_fm_power(ec%m_transf, mat_w, -0.5_dp, scf_cntrl%eps_eigval, n_deps)

      ! Sanity check
      IF (n_deps /= 0) THEN
         CALL cp_warn(__LOCATION__, &
                      "Overlap matrix exhibits linear dependencies. At least some "// &
                      "eigenvalues have been quenched.")
      END IF

      ! Save the transformation matrix to file
      CALL pod_save_transf_mat(qs_env, ec, 'F')

      ! Transformation S -> S^{+1/2}
      IF (output_unit > 0) &
         WRITE (output_unit, '(T3,A)') &
         'Constructing reverse basis-set transformation matrix'

      CPASSERT(.NOT. ASSOCIATED(ec%m_transf_inv))
      print *,'DEBUG: --> allocate m_transf_inv'
      ALLOCATE (ec%m_transf_inv)
      CPASSERT(ASSOCIATED(ec%m_transf_inv))
      print *,'DEBUG: --> call cp_fm_create'
      CALL cp_fm_create(matrix=ec%m_transf_inv, matrix_struct=mat_w%matrix_struct, &
                       name='S^(+1/2) TRANSFORMATION MATRIX')
      print *,'DEBUG: --> call copy_dbcsr_to_fm'
      CALL copy_dbcsr_to_fm(mat_s(1)%matrix, ec%m_transf_inv)
      print *,'DEBUG: --> call cp_fm_power'
      CALL cp_fm_power(ec%m_transf_inv, mat_w, +0.5_dp, scf_cntrl%eps_eigval, n_deps)

      ! Sanity check
      IF (n_deps /= 0) THEN
         CALL cp_warn(__LOCATION__, &
                      "Overlap matrix exhibits linear dependencies. At least some "// &
                      "eigenvalues have been quenched.")
      END IF

      ! Save the transformation matrix to file
      CALL pod_save_transf_mat(qs_env, ec, 'R')

      ! Running time
      CALL timestop(handle)

   END SUBROUTINE pod_get_transf_mat

! **************************************************************************************************
!> \brief check whether the POD data are read from the restart files
!> \param logger output stream data
!> \param ec_sec POD input-file section
!> \param file_tmf name of the POD restart file with forward transformation matrix
!> \param file_tmr name of the POD restart file with reverse transformation matrix
!> \param file_hab name of the POD restart file with Hamiltonian
!> \param file_wfn name of the POD restart file with MOs
!> \return ...
!> \author Z. Futera (03.2019)
! **************************************************************************************************
   FUNCTION pod_check_restart(logger, ec_sec, file_tmf, file_tmr, file_hab, file_wfn, blocks) &
      RESULT(restart)

      ! Routine arguments
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: ec_sec
      CHARACTER(LEN=default_path_length), OPTIONAL       :: file_tmf, file_tmr, file_hab, file_wfn
      LOGICAL, OPTIONAL, INTENT(OUT)                     :: blocks
      ! Returning value
      LOGICAL                                            :: restart

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_check_restart', &
         routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: handle
      LOGICAL                                            :: ok,blocks_only
      TYPE(section_vals_type), POINTER                   :: print_key

      ! Time mark
      CALL timeset(routineN, handle)

      ! Initialization 
      blocks_only = .FALSE.

      ! Check value of the RESTART option in the input file
      restart = section_get_lval(ec_sec, 'RESTART')

      ! Forward transformation matrix
      IF (PRESENT(file_tmf)) THEN
         print_key => section_vals_get_subs_vals(ec_sec, &
                                                 'PRINT%TRANSF_MATRIX_F')
         file_tmf = cp_print_key_generate_filename(logger, print_key, &
                                                   middle_name='RESTART', extension='.tmf', my_local=.FALSE.)
         IF (restart) THEN
            INQUIRE (file=file_tmf, exist=ok)
            IF (.NOT. ok) THEN
               CPWARN('restart file "'//TRIM(ADJUSTL(file_tmf))//'" not found')
               file_tmf = ''
            END IF
         END IF
      END IF

      ! Reverse transformation matrix
      IF (PRESENT(file_tmr)) THEN
         print_key => section_vals_get_subs_vals(ec_sec, &
                                                 'PRINT%TRANSF_MATRIX_R')
         file_tmr = cp_print_key_generate_filename(logger, print_key, &
                                                   middle_name='RESTART', extension='.tmr', my_local=.FALSE.)
         IF (restart) THEN
            INQUIRE (file=file_tmr, exist=ok)
            IF (.NOT. ok) THEN
               CPWARN('restart file "'//TRIM(ADJUSTL(file_tmr))//'" not found')
               file_tmr = ''
            END IF
         END IF
      END IF

      ! Hamiltonian restart file
      IF (PRESENT(file_hab)) THEN
         print_key => section_vals_get_subs_vals(ec_sec, &
                                                 'PRINT%HAMILTONIAN')
         file_hab = cp_print_key_generate_filename(logger, print_key, &
                                                   middle_name='RESTART', extension='.hab', my_local=.FALSE.)
         IF (restart) THEN
            INQUIRE (file=file_hab, exist=ok)
            IF (.NOT. ok) THEN
               CPWARN('restart file "'//TRIM(ADJUSTL(file_hab))//'" not found')
               file_hab = ''
            END IF
         END IF

         ! Blocks rather than full Hamiltonian
         blocks_only = section_get_lval(print_key, 'BLOCKS_ONLY')

      END IF

      ! MO restart file
      IF (PRESENT(file_wfn)) THEN
         print_key => section_vals_get_subs_vals(ec_sec, &
                                                 'PRINT%WAVEFUNCTION')
         file_wfn = cp_print_key_generate_filename(logger, print_key, &
                                                   middle_name='RESTART', extension='.wfn', my_local=.FALSE.)
         IF (restart) THEN
            INQUIRE (file=file_wfn, exist=ok)
            IF (.NOT. ok) THEN
               CPWARN('restart file "'//TRIM(ADJUSTL(file_wfn))//'" not found')
               file_wfn = ''
            END IF
         END IF
      END IF

      ! Check whether there is enough data to restart from
      IF (restart .and. .not.blocks_only .and. &
          (file_tmf == '' .or. file_tmr == '') .and. &
          (file_hab == '' .or. file_wfn == '')) THEN
         restart = .FALSE.
      END IF

      ! Hamiltonian blocks
      IF (PRESENT(blocks)) &
        blocks = blocks_only

      ! Running time
      CALL timestop(handle)

   END FUNCTION pod_check_restart

END MODULE et_coupling_pod_utils
