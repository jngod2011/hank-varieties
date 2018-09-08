        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep  7 20:35:46 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RTFLSP__genmod
          INTERFACE 
            SUBROUTINE RTFLSP(FUNC,LX1,LX2,XACC,FTOL,IFLAG,MAXIT)
              REAL(KIND=8) :: FUNC
              EXTERNAL FUNC
              REAL(KIND=8), INTENT(IN) :: LX1
              REAL(KIND=8), INTENT(IN) :: LX2
              REAL(KIND=8), INTENT(IN) :: XACC
              REAL(KIND=8), INTENT(IN) :: FTOL
              INTEGER(KIND=4), INTENT(OUT) :: IFLAG
              INTEGER(KIND=4), INTENT(IN) :: MAXIT
            END SUBROUTINE RTFLSP
          END INTERFACE 
        END MODULE RTFLSP__genmod
