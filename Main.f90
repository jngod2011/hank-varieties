PROGRAM Main

USE Parameters
USE Globals

CALL SetParameters
! IF (SobolSequenceExploration==1) CALL Exploration
! IF (CalibrateCostFunction==1) CALL Calibration
CALL InitialSteadyState
CALL SaveSteadyStateOutput
CALL FinalSteadyState
IF(DoImpulseResponses==1) CALL ImpulseResponses
! IF(DoFeedInPrices==1) CALL FeedInPrices

END PROGRAM Main