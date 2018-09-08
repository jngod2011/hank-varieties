        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep  7 20:35:44 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RTSEC__genmod
          INTERFACE 
            SUBROUTINE RTSEC(FUNC,LX1,LX2,LFACC,LXSOL,IFLAG)
              REAL(KIND=8) :: FUNC
              EXTERNAL FUNC
              REAL(KIND=8), INTENT(IN) :: LX1
              REAL(KIND=8), INTENT(IN) :: LX2
              REAL(KIND=8), INTENT(IN) :: LFACC
              REAL(KIND=8), INTENT(OUT) :: LXSOL
              INTEGER(KIND=4), INTENT(OUT) :: IFLAG
            END SUBROUTINE RTSEC
          END INTERFACE 
        END MODULE RTSEC__genmod
