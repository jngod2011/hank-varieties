        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep  7 20:35:46 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RTBIS__genmod
          INTERFACE 
            SUBROUTINE RTBIS(FUNC,X1,X2,XACC,FTOL,XOUT)
              REAL(KIND=8) :: FUNC
              EXTERNAL FUNC
              REAL(KIND=8), INTENT(IN) :: X1
              REAL(KIND=8), INTENT(IN) :: X2
              REAL(KIND=8), INTENT(IN) :: XACC
              REAL(KIND=8), INTENT(IN) :: FTOL
              REAL(KIND=8), INTENT(OUT) :: XOUT
            END SUBROUTINE RTBIS
          END INTERFACE 
        END MODULE RTBIS__genmod
