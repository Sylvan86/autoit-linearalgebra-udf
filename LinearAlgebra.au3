;~ #AutoIt3Wrapper_Run_AU3Check=Y
;~ #AutoIt3Wrapper_Au3Check_Parameters=-d -w 1 -w 2 -w 3 -w 4 -w 5 -w 6 -w 7
;~ ;~ #AutoIt3Wrapper_AU3Check_Stop_OnWarning=Y
;~ Opt("MustDeclareVars", 1)

#include-once
#include <Math.au3>
#include "LAPACK.au3"

; #INDEX# =======================================================================================================================
; Title .........: LinearAlgebra
; AutoIt Version : 3.16.1
; Description ...: Linear Algebra functionality based on BLAS/LAPACK
; Author(s) .....: AspirinJunkie
; Dll ...........: a BLAS/LAPACK dll
; Last changed ..: 2024-09-19
; Version .......: 0.1
; ===============================================================================================================================

; #CURRENT# =====================================================================================================================
; ---- vector/matrix creation ----
; _la_fromArray                 - converts a AutoIt array or array define string into a matrix map
; _la_fromStruct                - creates a matrix/vector map from a DllStruct as used here in the UDF
; _la_createVector              - creates new empty vector
; _la_createMatrix              - creates new empty matrix
; _la_createIdentity            - create identity matrix/vector
; _la_duplicate                 - creates an independent copy of a matrix/vector map
; _la_fromFile                  - reads a matrix or a vector from a file created by _la_toFile()
;
; ---- extraction/transforming ----
; _la_join                      - combines 2 matrices
; _la_transpose                 - transposes a matrix in-place or out-place and [optional] scaling
; _la_ReDim                     - changes the shape of a matrix by by changing the number of columns (also matrix <-> vector conversion)
; _la_getRow                    - extracts a row of a matrix as a vector
; _la_getColumn                 - extracts a column of a matrix as a vector
; _la_getDiag                   - extracts the diagonal of a matrix as a vector
; _la_getTriangle               - extract upper or lower triangle part of a matrix
; _la_VectorToDiag              - creates a diagonal matrix from a vector
;
; ---- data output ----
; _la_display                   - displays a matrix/vector map, similar to _ArrayDisplay
; _la_toArray                   - converts a matrix/vector map into an AutoIt array
; _la_toFile                    - write a matrix/vector into a file
;
; ---- scalar operations ----
; _la_rotate                    - applies a plane rotation to coordinate-pairs
;
; ---- matrix attributes ----
; _la_isPositiveDefinite        - checks whether a matrix is positive definite
; _la_isSymmetric               - checks whether a matrix is symmetrical
; _la_rank                      - determines the rank of a matrix
; _la_determinant               - calculate the determinant of a matrix
; _la_conditionNumber           - determine the condition number of a matrix
;
; ---- unary operations ----
; _la_inverse                   - calculates the inverse of a matrix
; _la_pseudoInverse             - calculate the Moore-Penrose pseudo inverse of a matrix
; _la_sum                       - calculates the sum of the elements of a matrix, vector or parts thereof
; _la_asum                      - calculate the sum of the absolute(!) values of a matrix/vector
; _la_amin                      - finds the first element having the minimum absolute(!) value
; _la_amax                      - finds the first element having the maximum absolute(!) value
; _la_norm                      - calculate the euclidian norm of a vector
; _la_mean                      - calculate the mean of a vector or parts of a matrix
;
; ---- element wise operations ----
; _la_sqrtElements              - calculates the square root of each element of a matrix/vector
; _la_squareElements            - calculates the square of each element of a matrix/vector
; _la_invElements               - forms the reciprocal (1/x) for each element of the matrix/vector
;
; ---- addition subtraction ----
; _la_sub                       - subtracts a matrix/vector B from matrix/vector A
; _la_add                       - calculate the sum of a matrix/vector/scalar mA and a matrix/vector/scalar mB
;
; ---- multiplication ----
; _la_mul                       - calculates a multiplication between a matrix/vector/scalar A and a matrix/vector/scalar B
; _la_outerproduct              - calculates the outer product ("tensor product") of two vectors
; _la_dot                       - calculate the "dot product"/"scalar product"/"inner product" of two vectors
; _la_scale                     - multiplies the elements of a matrix/vector by a scalar value
; _la_mulElementWise            - calculates the element-wise ("Hadarmard") product between two matrices/vectors
; _la_cross                     - calculates the cross product between two 3-element vectors
;
; ---- factorization / decomposition ----
; _la_LU                        - calculates the LU decomposition of a matrix
; _la_QR                        - calculates the QR decomposition of a matrix
; _la_SVD                       - calculates the singular value decomposition (SVD) of a matrix
; _la_cholesky                  - calculate the cholesky decomposition of a symmetric, positive definite matrix ( A --> L * Lᵀ or A --> U * Uᵀ )
;
; ---- eigenvalues / eigenvectors ----
; _la_eigen                     - computes for an N-by-N real matrix A, the eigenvalues and the left and/or right eigenvectors.
;
; ---- solve linear equation systems ----
; _la_solve                     - computes the solution to a system of linear equations A * X = B
;
; ---- least squares solving ----
; _la_lstsq                     - solves overdetermined or underdetermined [weighted] linear system
;
; ---- regression ----
; _la_regression                - calculates an n-dimensional linear or non-linear regression
;
; ---- adjustment ----
; _la_adjustement               - performs a least-squares adjustment calculation for a system of different [weighted] non-linear equations
; _la_adjustment_l1             - performs a adjustment calculation to L1 norm for a system of different [weighted] non-linear equations
; _la_adj_addObservation        - adds an observation to the adjustment system
;
; ---- additional helper functions ----
; _la_adj_showResult            - formats the results of _la_adj more clearly and display them in a window
; ===============================================================================================================================

; #INTERNAL_USE_ONLY# ===========================================================================================================
; ---- addition subtraction ----
; __la_addScalar                     - adds a constant value to the elements of a matrix/vector
;
; ---- least squares solving ----
; __la_lstsq_qr                      - solves overdetermined or underdetermined [weighted] linear system using QR decomposition
; __la_lstsq_svd                     - solves overdetermined or underdetermined [weighted] linear system using singular value decomposition
; __la_lstsq_cholesky                - solves overdetermined or underdetermined [weighted] linear system using cholesky decomposition
;
; ---- regression ----
; __la_regression_GaussNewton        - calculates an n-dimensional linear or non-linear regression using Gauss Newton method
; __la_regression_LevenbergMarquardt - calculates an n-dimensional linear or non-linear regression using levenberg-marquardt model
; __la_Regression_getJacobiParams    - derive the Jacobian matrix for a regression
; __la_Regression_getYfromModel      - determines the vector of model results based on the approximated parameters
; __la_regression_prepareApproxVals  - prepares the list of parameters with their initial values
; __la_regression_prepareVars        - prepares the list of variables
; __la_derivate1D                    - calculates the 1st derivative of a function
;
; ---- adjustment ----
; __la_adj_LevenbergMarquardt        - performs a adjustment calculation for a system of different [weighted] non-linear equations by using the Levenberg-Marquardt algorithm
; __la_adj_GaussNewton               - performs a adjustment calculation for a system of different [weighted] non-linear equations by using the Gauss-Newton algorithm
; __la_adj_getYfromModel             - determines the vector of model results based on the approximated parameters
; __la_adj_getJacobiParams           - calculate the jacobian matrix out of the observation formulas and the approximated parameters
; __la_adj_setApproxValue            - sets/changes the initial value of a parameter in the parameter list
; __la_adj_getWeights                - extracts the vector (or diagonal matrix) of the observation weights from the observations
; __la_adj_getParamList              - extract the parameters to be estimated from the observation equations into a map
; __la_adj_getVarComponents          - extracts a list of existing variance component groups in the system and calculates the weight of an observation
;
; ---- additional helper functions ----
; __la_getSignificantDecimals        - determines the position of the 1st significant decimal place of a number
; ===============================================================================================================================


; #CONSTANTS# ===================================================================================================================
Global Const $f_LA_FLT_EPS = _lp_lamch("e", "FLOAT"),  $f_LA_FLT_MIN = _lp_lamch("S", "FLOAT"),  $f_LA_FLT_MAX = 3.402823E+38,             $f_LA_FLT_PREC = _lp_lamch("p", "FLOAT"),  $f_LA_FLT_STEP = $f_LA_FLT_EPS^(1/3), _
             $f_LA_DBL_EPS = _lp_lamch("e", "DOUBLE"), $f_LA_DBL_MIN = _lp_lamch("S", "DOUBLE"), $f_LA_DBL_MAX = 1.79769313486231570E+308, $f_LA_DBL_PREC = _lp_lamch("p", "DOUBLE"), $f_LA_DBL_STEP = $f_LA_DBL_EPS^(1/3)

; constants/flags to control what elements should be returned by the least square solution functions
Global Enum Step *2 $__LA_LSTSQ_R = 1, $__LA_LSTSQ_R2Sum, $__LA_LSTSQ_S0, $__LA_LSTSQ_QX, $__LA_LSTSQ_SDX, $__LA_LSTSQ_QY, $__LA_LSTSQ_QYD, $__LA_LSTSQ_SDY, $__LA_LSTSQ_QR, $__LA_LSTSQ_SDR, $__LA_LSTSQ_REDUNDANCY, $__LA_LSTSQ_COND, $__LA_LSTSQ_RANK
; $__LA_LSTSQ_R		      - Residuals r:  r = y - y_d
; $__LA_LSTSQ_R2Sum       - square sum of [weighted] residuals: r2sum = rᵀ * W * v
; $__LA_LSTSQ_S0          - a posteriori standard deviation factor
; $__LA_LSTSQ_QX          - cofactor matrix of parameters
; $__LA_LSTSQ_SDX         - standard deviations of parameters: sdₓ = sqrt(s₀² * diag(Qₓ))
; $__LA_LSTSQ_QY          - a-priori cofactor matrix for the observations: Q_y = P⁻¹
; $__LA_LSTSQ_QYD         - a posteriori cofactormatrix for the adjusted observations: Q_yd = A * Qₓ * Aᵀ
; $__LA_LSTSQ_SDY         - a posteriori standard deviations for the adjusted observations
; $__LA_LSTSQ_QR          - cofactor matrix for the residuals r: Q_r = Q_y - Q_yd
; $__LA_LSTSQ_SDR         - standard deviations for the residuals: sd_r = sqrt(s₀² * diag(Q_r))
; $__LA_LSTSQ_REDUNDANCY  - R = Q_y * W
; $__LA_LSTSQ_COND        - condition number (only if "SVD" is used)
; $__LA_LSTSQ_RANK        - rank for the jacobian matrix (only if "SVD" is used)
; ===============================================================================================================================


#Region vector/matrix creation

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_fromArray()
; Description ...: converts a AutoIt array or array define string into a matrix map
; Syntax ........: _la_fromArray($aArray, [$nMode = 0, [$sType = "DOUBLE", [$iKL = 0, [$iKU = 0]]]])
; Parameters ....: aArray - [Array] AutoIt array or array define string (see examples) which should be converted into a matrix map
;                  nMode  - (Default: 0)
;                         ↳ BLAS memory layout for the matrix (see $__g_BLAS_STYPEₓXX flags)
;                  sType  - [String] (Default: "DOUBLE")
;                         ↳ data type of the elements ("DOUBLE" or "FLOAT")
;                  iKL    - [UInt] (Default: 0)
;                         ↳ If nMode contains $__g_BLAS_STYPE_BAND: number of lower band diagonals
;                  iKU    - [UInt] (Default: 0)
;                         ↳ If nMode contains $__g_BLAS_STYPE_BAND: number of upper band diagonals
; Return value ..: Success: [Map] Matrix Map
;                  Failure: Null and set @error to:
;                           | 1: invalid value for sType
;                           | 2: invalid value for aArray
;                           | 1X: error X during _blas_fromArray() (@extended: @extended from _blas_fromArray())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......: wrapper for _blas_fromArray()
; Related .......: _blas_fromArray()
; Link ..........:
; Example .......: Yes
;                  Global $aArray[][] = [[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97],[93,89,13,-53,43]]
;                  Global $mA = _la_fromArray($aArray)
;                  Global $mS = _la_fromArray('[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97],[93,89,13,-53,43]]')
;                  _la_display($mA, "from AutoIt Array")
;                  _la_display($mS, "from string definition")
; ===============================================================================================================================
Func _la_fromArray($aArray, Const $nMode = 0, Const $sType = "DOUBLE", Const $iKL = 0, Const $iKU = 0)

	; validation of the input parameters
	If $sType <> "DOUBLE" And $sType <> "FLOAT"  Then Return SetError(1, 0, Null)
	If Not (IsArray($aArray) Or IsString($aArray)) Then Return SetError(2, 0, Null)

	Local $mRet = _blas_fromArray($aArray, $nMode, $sType, $iKL, $iKU)
	Return SetError(@error + 10, @extended, $mRet)
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_fromStruct()
; Description ...: creates a matrix/vector map from a DllStruct as used here in the UDF
; Syntax ........: _la_fromStruct($tStruct, $iRows, [$iCols = 0, [$sDatatype = "DOUBLE", [$nMode = $__g_BLAS_STYPE_MATRIX]]])
; Parameters ....: tStruct   - [DllStruct] the DllStruct variable with the payload data of the matrix/vector
;                  iRows     - [UInt] Number of rows in the matrix/vector
;                  iCols     - [UInt] (Default: 0)
;                            ↳ Number of columns in the matrix/vector
;                  sDatatype - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
;                  nMode     - (Default: $__g_BLAS_STYPE_MATRIX)
;                            ↳ BLAS memory layout for the matrix (see $__g_BLAS_STYPEₓXX flags)
; Return value ..: Success: [Map] matrix map
;                  Failure: Null and set @error to:
;                           | 1: invalid value for sDatatype
;                           | 2: invalid value for iRows (@extended: iRows)
;                           | 3: invalid value for iCols (@extended: iCols)
;                           | 4: invalid value for tStruct
;                           | 1X: error X during _blas_fromStruct() (@extended: @extended from _blas_fromStruct())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......: wrapper for _blas_fromStruct()
; Related .......: _blas_fromStruct()
; Link ..........:
; Example .......: Yes
;                  Global $tStruct = DllStructCreate("DOUBLE[10]")
;                  For $i = 1 To 10
;                     DllStructSetData($tStruct, 1, $i * 2, $i)
;                  Next
;                  Global $mVector = _la_fromStruct($tStruct, 10)
;                  _la_display($mVector)
; ===============================================================================================================================
Func _la_fromStruct(ByRef $tStruct, $iRows, $iCols = 0, $sDatatype = "DOUBLE", $nMode = $__g_BLAS_STYPE_MATRIX)

	; validation of the input parameters
	If $sDatatype <> "DOUBLE" And $sDatatype <> "FLOAT"  Then Return SetError(1, 0, Null)
	If $iRows < 0 Then Return SetError(2, $iRows, Null)
	If $iCols < 0 Then Return SetError(3, $iCols, Null)
	If Not IsDllStruct($tStruct) Then Return SetError(4, 0, Null)

	Local $mRet = _blas_fromStruct($tStruct, $iRows, $iCols, $sDatatype, $nMode)
	Return SetError(@error + 10, @extended, $mRet)
EndFunc   ;==>_la_fromStruct

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_createVector()
; Description ...: creates new empty vector
; Syntax ........: _la_createVector($iN, [$sType = "DOUBLE", [$pTarget = Default]])
; Parameters ....: iN      - [UInt] number of elements in the vector
;                  sType   - [String] (Default: "DOUBLE")
;                          ↳ data type of the elements ("DOUBLE" or "FLOAT")
;                  pTarget - [Pointer] (Default: Default)
;                          ↳ pointer to a memory area in which the vector is to be created
;                            if Default a new memory area is reserved
; Return value ..: Success: [Map] vector map
;                  Failure: Null and set @error to:
;                           | 1: invalid value for sType
;                           | 2: invalid value for iN (@extended: $iN)
;                           | 1X: error X during _blas_createVector()  (@extended: @extenden from _blas_createVector())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......: wrapper for _blas_createVector()
; Related .......: _blas_createVector()
; Link ..........:
; Example .......: Yes
;                  Global $mVec = _la_createVector(10)
;                  _la_display($mVec, "new empty vector")
; ===============================================================================================================================
Func _la_createVector(Const $iN, Const $sType = "DOUBLE", Const $pTarget = Default)

	; validation of the input parameters
	If $sType <> "DOUBLE" And $sType <> "FLOAT" Then Return SetError(1, 0, Null)
	If $iN < 1 Then Return SetError(2, $iN, Null)

	Local $mRet = _blas_createVector($iN, $sType, $pTarget)
	Return SetError(@error + 10, @extended, $mRet)
EndFunc   ;==>_la_createVector

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_createMatrix()
; Description ...: creates new empty matrix
; Syntax ........: _la_createMatrix($iR, [$iC = $iR, [$sType = "DOUBLE"]])
; Parameters ....: iR    - [UInt] number of matrix rows
;                  iC    - [UInt] (Default: $iR)
;                        ↳ number of matrix columns
;                  sType - [String] (Default: "DOUBLE")
;                        ↳ data type of the elements ("DOUBLE" or "FLOAT")
; Return value ..: Success: [Map] matrix map
;                  Failure: $mRet and set @error to:
;                           | 1: invalid value for sType
;                           | 2: invalid value for iR (@extended: $iR)
;                           | 3: invalid value for iC (@extended: $iC)
;                           | 1X: error X during _blas_createMatrix()  (@extended: @extenden from _blas_createMatrix())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......: wrapper for _blas_createMatrix()
; Related .......: _blas_createMatrix()
; Link ..........:
; Example .......: Yes
;                  Global $mMat = _la_createMatrix(5, 3)
;                  _la_display($mMat, "new empty matrix")
; ===============================================================================================================================
Func _la_createMatrix(Const $iR, Const $iC = $iR, Const $sType = "DOUBLE")

	; validation of the input parameters
	If $sType <> "DOUBLE" And $sType <> "FLOAT" Then Return SetError(1, 0, Null)
	If $iR < 1 Then Return SetError(2, $iR, Null)
	If $iC < 1 Then Return SetError(3, $iC, Null)

	Local $mRet = _blas_createMatrix($iR, $iC, $sType)
	Return SetError(@error + 10, @extended, $mRet)
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_createIdentity()
; Description ...: create identity matrix/vector
; Syntax ........: _la_createIdentity($iRows, [$iCols = 0, [$sDatatype = "DOUBLE"]])
; Parameters ....: iRows     - [UInt] number of matrix rows or vector elements
;                  iCols     - [UInt] (Default: 0)
;                            ↳ 0: result is a identity vector
;                              else: result is a identity matrix
;                  sDatatype - [String] (Default: "DOUBLE")
;                            ↳ data type of the elements ("DOUBLE" or "FLOAT")
; Return value ..: Success: SetError(3, $iCols, Null)
;                  Failure: Null and set @error to:
;                           | 1X: error X during _blas_createVector/_blas_createMatrix  (@extended: @extended from call)
;                           | 1: invalid value for sDatatype
;                           | 2: invalid value for iRows (@extended: $iRows)
;                           | 3: invalid value for iCols (@extended: $iCols)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _blas_createVector, _blas_createMatrix, __blas_fillWithScalar
; Link ..........:
; Example .......: Yes
;                  Global $mIV = _la_createIdentity(10)
;                  Global $mIM = _la_createIdentity(10, 10)
;                  _la_display($mIV, "identity vector")
;                  _la_display($mIM, "identity matrix")
; ===============================================================================================================================
Func _la_createIdentity($iRows, $iCols = 0, $sDatatype = "DOUBLE")

	; validation of the input parameters
	If $sDatatype <> "DOUBLE" And $sDatatype <> "FLOAT" Then Return SetError(1, 0, Null)
	If $iRows < 1 Then Return SetError(2, $iRows, Null)

	Local $mRet
	Select
		Case $iCols = 0 ; Vector
			$mRet = _blas_createVector($iRows, $sDatatype)
			__blas_fillWithScalar($mRet, 1, 0, 1, $iRows)
			If @error Then Return SetError(@error + 10, @extended, Null)

		Case $iCols > 0 ; Matrix
			$mRet = _blas_createMatrix($iRows, $iCols, $sDatatype)
			__blas_fillWithScalar($mRet, 1, 0, $iRows + 1)
			If @error Then Return SetError(@error + 10, @extended, Null)

		Case $iCols < 0 ; error
			Return SetError(3, $iCols, Null)

	EndSelect

	Return $mRet
EndFunc   ;==>_la_createIdentity


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_duplicate()
; Description ...: creates an independent copy of a matrix/vector map
; Syntax ........: _la_duplicate($mMatrix)
; Parameters ....: mMatrix - [Map] matrix/vector as a map/array/definition string
; Return value ..: Success: [Map] the copied matrix map
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mMatrix
;                           | 1X: error X during _blas_duplicate() (@extended: @extended from _blas_duplicate())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......: wrapper for _blas_duplicate()
; Related .......: _blas_duplicate()
; Link ..........:
; Example .......: Yes
;                  Global $mA = _la_fromArray('[[1,2,3],[4,5,6],[7,8,9]]'), _
;                  $mDouble = _la_duplicate($mA)
;                  _la_display($mDouble, "copy of Matrix")
; ===============================================================================================================================
Func _la_duplicate($mMatrix)
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; validation of the input parameters
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, $mMatrix) ; valid AutoIt-BLAS/LAPACK-Map?

	Local $mRet = _blas_duplicate($mMatrix)
	Return @error ? SetError(@error + 10, @extended, Null) : $mRet
EndFunc   ;==>_la_duplicate


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_fromFile()
; Description ...: reads a matrix or a vector from a file created by _la_toFile()
; Syntax ........: _la_fromFile($sFile)
; Parameters ....: sFile - [String] file path
; Return value ..: Success: [Map] matrix/vector map
;                  Failure: Null and set @error to:
;                           | 1: error during FileOpen (@extended: @error from FileOpen())
;                           | 2: error during _blas_fromStruct() (@extended: @error from _blas_fromStruct())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func _la_fromFile($sFile)
	Local $hFile = FileOpen($sFile, 16)
	If $hFile = -1 Then Return SetError(1, @error, Null)

	Local $bWhole = FileRead($hFile)
	FileClose($hFile)
	Local $tWhole = DllStructCreate(StringFormat("BYTE[%d]", BinaryLen($bWhole)))
	Local $pWhole = DllStructGetPtr($tWhole)
	DllStructSetData($tWhole, 1, $bWhole)

	Local $tMeta = DllStructCreate('BOOLEAN type; UINT rows; UINT cols; USHORT flags; USHORT kl; USHORT ku; DOUBLE Tmp', $pWhole)

	Local $iSize = (DllStructGetSize($tWhole) - DllStructGetSize($tMeta) + 8 ) / ($tMeta.type = 0 ? $iBLAS_SIZE_DOUBLE : $iBLAS_SIZE_FLOAT)  ; +8 because of FLOAT Tmp - the align dummy

	Local $tData = DllStructCreate(StringFormat("%s[%d]", $tMeta.type = 0 ? "DOUBLE" : "FLOAT", $iSize), DllStructGetPtr($tMeta, 7))

	Local $mData =  _blas_fromStruct($tData, $tMeta.rows, $tMeta.cols, $tMeta.type = 0 ? "DOUBLE" : "FLOAT", $tMeta.flags, $tMeta.kl, $tMeta.ku)
	If @error Then Return SetError(2, @error, Null)

	Return _la_duplicate($mData)
EndFunc

#EndRegion

#Region extraction/transforming

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_join()
; Description ...: combines 2 matrices
; Syntax ........: _la_join($mA, [$mB = Default, [$sPosition = "vertical"]])
; Parameters ....: mA        - [Map] matrix A as a map/array/defintion string
;                  mB        - [Map] matrix B as a map
;                  sPosition - [String] (Default: "vertical")
;                            ↳ "horizontal": result matrix = [A, B]
;                              "vertical":   result matrix = [[A], [B]]
;                              "diagonal":   result matrix = [[A, 0], [0, B]]
; Return value ..: Success: [Map] result matrix
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mA/mB  (@extended: 1 = mA, 2 = mB )
;                           | 2: invalid dimensions between mA and mB (@extended: 1 = vertical, 2 = horizontal )
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: Yes
;                  Global $mA = _la_fromArray('[[1,2,3],[4,5,6],[7,8,9],[10,11,12],[13,14,15],[16,17,18],[19,20,21]]')
;                  Global $mB = _la_fromArray('[[1,0,5],[3,1,0],[0,2,1]]')
;                  Global $mC = _la_join($mA, $mB, "vertical")
;                  _la_display($mC, "vertical")
;                  $mC = _la_join($mA, $mB, "diagonal")
;                  _la_display($mC, "diagonal")
;                  Global $mA = _la_fromArray('[[1,2,3],[4,5,6],[7,8,9]]')
;                  Global $mB = _la_fromArray('[[1,0],[0,1],[1,0]]')
;                  Global $mC = _la_join($mA, $mB, "horizontal")
;                  _la_display($mC)
; ===============================================================================================================================
Func _la_join($mA, $mB, $sPosition = "vertical")
	; direct AutoIt-type input
	If IsArray($mA) Or IsString($mA) Then $mA = _blas_fromArray($mA)

	; validation of the input parameters
	If Not MapExists($mA, "ptr") Then Return SetError(1, 1, Null)
	If Not MapExists($mB, "ptr") Then Return SetError(1, 2, Null)

	Local $iMA = $mA.rows, $iNA = $mA.cols = 0 ? 1 : $mA.cols, $iSizeA = $iMA * $iNA, _
	      $iMB = $mB.rows, $iNB = $mB.cols = 0 ? 1 : $mB.cols, _
		  $sDataType = $mA.datatype, _
		  $dSize = ($sDataType = "DOUBLE") ? $iBLAS_SIZE_DOUBLE : $iBLAS_SIZE_FLOAT, _
		  $mC, $pC, $iMC, $iNC

	Switch $sPosition
		Case "vertical" ; B under A
			If $iNA <> $iNB Then Return SetError(2, 1, Null)

			; create new buffer
			$iMC = $iMA + $iMB
			$iNC = $iNA
			$mC = $iNA < 2 _
				? _blas_createVector($iMC * $iNC, $sDataType) _
				: _blas_createMatrix($iMC, $iNC, $sDataType)
			$pC = $mC.ptr

			; copy A
			_lp_lacpy($mA.ptr, $pC, "X", $iMA, $iNA, $iMA, $iMC, $sDataType)
			; copy B
			_lp_lacpy($mB.ptr, $pC + $dSize * $iMA, "X", $iMB, $iNB, $iMB, $iMC, $sDataType)

		Case "horizontal" ; B right side of A
			If $iMA <> $iMB Then Return SetError(2, 2, Null)

			; create new buffer
			$iMC = $iMA
			$iNC = $iNA + $iNB
			$mC = $iNA < 2 _
				? _blas_createVector($iMC * $iNC, $sDataType) _
				: _blas_createMatrix($iMC, $iNC, $sDataType)
			$pC = $mC.ptr

			; copy A
			_blas_copy($mA.ptr, 0, 1, 0, 1, $iSizeA, $pC, False, $sDataType)
			; copy B
			_blas_copy($mB.ptr, 0, 1, $iSizeA, 1, $iMB * $iNB, $pC, False, $sDataType)

		Case Else ; = "diagonal"
			; create new buffer
			$iMC = $iMA + $iMB
			$iNC = $iNA + $iNB
			$mC = $iNA < 2 _
				? _blas_createVector($iMC * $iNC, $sDataType) _
				: _blas_createMatrix($iMC, $iNC, $sDataType)
			$pC = $mC.ptr

			; copy A
			_lp_lacpy($mA.ptr, $pC, "X", $iMA, $iNA, $iMA, $iMC, $sDataType)
			; copy B
			_lp_lacpy($mB.ptr, $pC + $dSize * ($iMC * $iNA + $iMA), "X", $iMB, $iNB, $iMB, $iMC, $sDataType)

	EndSwitch

	Return $mC

EndFunc



; #FUNCTION# ====================================================================================================================
; Name ..........: _la_transpose()
; Description ...: transposes a matrix in-place or out-place and [optional] scaling
; Syntax ........: _la_transpose($mMatrix, [$fAlpha = 1, [$bInPlace = False]])
; Parameters ....: mMatrix  - [Map] matrix/vector as a map/array/definition string
;                  fAlpha   - [Float] (Default: 1)
;                           ↳ optional scaling factor
;                  bInPlace - [Bool] (Default: False)
;                           ↳ True: mMatrix is overwritten with the transposed matrix
;                             False: the transposed matrix will be returned (mMatrix remains untouched)
; Return value ..: Success: True or [Map] result matrix
;                  Failure: Null and set @error to:
;                           | 1X: error X during imatcopy/omatcopy (@extended: @extended from imatcopy/omatcopy)
;                           | 1: invalid value for mMatrix
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _blas_imatcopy(), _blas_omatcopy()
; Link ..........:
; Example .......: Yes
;                  Global $mMatrix = _la_fromArray('[[1,2,3,4,5],[6,7,8,9,10],[11,12,13,14,15],[16,17,18,19,20]]')
;                  _la_display($mMatrix, "original matrix")
;                  Global $mTransposed = _la_transpose($mMatrix)
;                  _la_display($mTransposed, "transposed")
;                  ConsoleWrite(@error & @TAB & @extended & @CRLF)
; ===============================================================================================================================
Func _la_transpose($mMatrix, $fAlpha = 1, $bInPlace = False)
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; validation of the input parameters
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, False) ; valid AutoIt-BLAS/LAPACK-Map?

	Local $iTmp
	If $bInPlace Then
		_blas_imatcopy($mMatrix, $fAlpha)
		If @error Then Return SetError(@error + 10, @extended, False)

		; switch rows <> cols
		$iTmp = $mMatrix.rows
		$mMatrix.rows = $mMatrix.cols
		$mMatrix.cols = $iTmp

		Return True

	Else
		Local $mRet = _blas_createMatrix($mMatrix.cols, $mMatrix.rows, $mMatrix.datatype)
		_blas_omatcopy($mMatrix, $fAlpha, "T", $mRet)
		Return @error ? SetError(@error + 10, @extended, Null) : $mRet
	EndIf
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_ReDim()
; Description ...: changes the shape of a matrix by by changing the number of columns (also matrix <-> vector conversion)
; Syntax ........: _la_ReDim($mMatrix, [$iC = Default])
; Parameters ....: mMatrix - [Map] matrix/vector as a map/array/definition string
;                  iC      - [UInt] (Default: Default)
;                          ↳ new number of columns (If Default: convert matrix to vector)
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: invalid value for mMatrix
;                           | 2: invalid value for iC (@extended: iC)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......: TODO: handle packed and triangular matrices. Currently for unpacked full matrices and vectors only
; Related .......:
; Link ..........:
; Example .......: Yes
;                  Global $mMatrix = _la_fromArray('[[1,2,3],[4,5,6],[7,8,9]]')
;                  _la_ReDim($mMatrix) ; to Vector
;                  _la_display($mMatrix, "as vector (column-order!)")
;                  _la_ReDim($mMatrix, 3) ; back to matrix with 3 columns
;                  _la_display($mMatrix, "matrix again")
; ===============================================================================================================================
Func _la_ReDim(ByRef $mMatrix, $iC = Default)
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; validation of the input parameters
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, False) ; valid AutoIt-BLAS/LAPACK-Map?

	If IsKeyword($iC) = 1 Then ; to Vector
		$mMatrix.rows        = $mMatrix.elements
		$mMatrix.cols        = 1
		$mMatrix.storageType = 0
	ElseIf $iC > 0 Then ; to Matrix
		$mMatrix.rows = Ceiling($mMatrix.elements / $iC)
		$mMatrix.cols = $iC
		$mMatrix.storageType = BitOR($mMatrix.storageType, $__g_BLAS_STYPE_MATRIX)
	Else
		Return SetError(2, $iC, False)
	EndIf

	Return True
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_getRow()
; Description ...: extracts a row of a matrix as a vector
; Syntax ........: _la_getRow($mMatrix, [$iRow = 0, [$pTarget = Default]])
; Parameters ....: mMatrix - [Map] matrix as a map/array/definition string
;                  iRow    - [UInt] (Default: 0)
;                          ↳ index of the row to be extracted (0-based)
;                  pTarget - [Pointer] (Default: Default)
;                          ↳ Default: new vector is beeing created and returned
;                            Else: pointer to the memory address to which the data is to be written
; Return value ..: Success: pTarget = Default ? [Map] result vector : True
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mMatrix
;                           | 2: invalid value for iRow (@extended: $iRow)
;                           | 1X: error X during _blas_copy() (@extended: @extended from _blas_copy())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _blas_copy()
; Link ..........:
; Example .......: Yes
;                  Global $mRow = _la_getRow('[[1,2,3],[4,5,6],[7,8,9]]', 1)
;                  _la_display($mRow, "2nd row")
; ===============================================================================================================================
Func _la_getRow(ByRef $mMatrix, Const $iRow = 0, $pTarget = Default)
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; validation of the input parameters
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, Null) ; valid AutoIt-BLAS/LAPACK-Map?
	If $iRow < 0 Or $iRow >= $mMatrix.rows Then Return SetError(2, $iRow, Null)

	Local $mRet = _blas_copy($mMatrix, $iRow, $mMatrix.cols, 0, 1, $mMatrix.cols, $pTarget)
	Return @error ? SetError(@error + 10, @extended, Null) : $mRet
EndFunc   ;==>_la_getRow


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_getColumn()
; Description ...: extracts a column of a matrix as a vector
; Syntax ........: _la_getColumn($mMatrix, [$iColumn = 0, [$pTarget = Default]])
; Parameters ....: mMatrix - [Map] matrix as a map/array/definition string
;                  iColumn - [UInt] (Default: 0)
;                          ↳ index of the column to be extracted (0-based)
;                  pTarget - [Pointer] (Default: Default)
;                          ↳ Default: new vector is beeing created and returned
;                            Else: pointer to the memory address to which the data is to be written
; Return value ..: Success: pTarget = Default ? [Map] result vector : True
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mMatrix
;                           | 2: invalid value for iRow (@extended: $iRow)
;                           | 1X: error X during _blas_copy() (@extended: @extended from _blas_copy())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _blas_copy()
; Link ..........:
; Example .......: Yes
;                  Global $mColumn = _la_getColumn('[[1,2,3],[4,5,6],[7,8,9]]', 1)
;                  _la_display($mColumn, "2nd column")
; ===============================================================================================================================
Func _la_getColumn(ByRef $mMatrix, Const $iColumn = 0, $pTarget = Default)
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; validation of the input parameters
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, Null) ; valid AutoIt-BLAS/LAPACK-Map?
	If $iColumn < 0 Or $iColumn >= $mMatrix.cols Then Return SetError(2, $iColumn, Null)

	Local $mRet = _blas_copy($mMatrix, $iColumn * $mMatrix.rows, 1, 0, 1, $mMatrix.rows, $pTarget)
	Return @error ? SetError(@error + 10, @extended, Null) : $mRet
EndFunc   ;==>_la_getColumn

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_getDiag()
; Description ...: extracts the diagonal of a matrix as a vector
; Syntax ........: _la_getDiag($mMatrix, [$pTarget = Default])
; Parameters ....: mMatrix - [Map] matrix as a map/array/definition string
;                  pTarget - [Pointer] (Default: Default)
;                          ↳ Default: new vector is beeing created and returned
;                            Else: pointer to the memory address to which the data is to be written
; Return value ..: Success: pTarget = Default ? [Map] result vector : True
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mMatrix
;                           | 1X: error X during _blas_copy() (@extended: @extended from _blas_copy())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......:_blas_copy
; Link ..........:
; Example .......: Yes
;                  Global $mDiag = _la_getDiag('[[1,2,3],[4,5,6],[7,8,9]]')
;                  _la_display($mDiag, "diag(A)")
; ===============================================================================================================================
Func _la_getDiag($mMatrix, $pTarget = Default)
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; validation of the input parameters
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, Null) ; valid AutoIt-BLAS/LAPACK-Map?

	Local $mRet = _blas_copy($mMatrix, 0, $mMatrix.rows + 1, 0, 1, $mMatrix.rows, $pTarget)
	Return @error ? SetError(@error + 10, @extended, Null) : $mRet
EndFunc   ;==>_la_getDiag


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_getTriangle()
; Description ...: extract upper or lower triangle part of a matrix
; Syntax ........: _la_getTriangle($mMatrix, [$cUpperLower = "U", [$bDiag = True, [$bShrink = True]]])
; Parameters ....: mMatrix     - [Map] matrix as a map/array/definition string
;                  cUpperLower - [Char] (Default: "U")
;                              ↳ "U": upper triangular part of the matrix is extracted
;                                "L": lower triangular part of the matrix is extracted
;                  bDiag       - [Bool] (Default: True)
;                              ↳ True: Diagonal is also extracted
;                                False: Diagonal is not extracted
;                  bShrink     - [Bool] (Default: False)
;                              ↳ If rows > columns then remove empty rows (True) or not (False)
; Return value ..: Success: [Map] result triangle matrix (unpacked)
;                  Failure: Null and set @error to:
;                           | 1X: error X during _lp_lacpy() (@extended: @extended from _lp_lacpy())
;                           | 2X: error X during _blas_createMatrix() (@extended: @extended from _blas_createMatrix())
;                           | 1: invalid value for mMatrix
;                           | 2: invalid value for cUpperLower
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _lp_lacpy()
; Link ..........:
; Example .......: Yes
;                  Global $mTriangle = _la_getTriangle('[[1,2,3,4,5],[6,7,8,9,10],[11,12,13,14,15],[16,17,18,19,20],[21,22,23,24,25],[26,27,28,29,30]]', "U")
;                  _la_display($mTriangle, "upper right triangle without diagonal")
; ===============================================================================================================================
Func _la_getTriangle(ByRef $mMatrix, $cUpperLower = "U", $bDiag = True, $bShrink = False)
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; validation of the input parameters
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, Null) ; valid AutoIt-BLAS/LAPACK-Map?
	If $cUpperLower <> "U" And $cUpperLower <> "L" Then Return SetError(2, 0, Null)

	Local $iRows = $bShrink ? ($mMatrix.rows > $mMatrix.cols ? $mMatrix.cols : $mMatrix.rows) : $mMatrix.rows

	; target matrix
	Local $mB = _blas_createMatrix($iRows, $mMatrix.cols, $mMatrix.datatype)
	If @error Then Return SetError(@error + 20, @extended, Null)

	; copy the triangle part to target matrix
	_lp_lacpy($mMatrix, $mB, $cUpperLower, $mMatrix.rows, $mMatrix.cols, $mMatrix.rows, $mB.rows, $mMatrix.datatype)
	If @error Then Return SetError(@error + 10, @extended, Null)

	; set diagonal to zero if needed
	If Not $bDiag Then __blas_fillWithScalar($mB, 0, 0, $mB.rows + 1)

	$mB.storageType = BitOR($mB.storageType, $__g_BLAS_STYPE_TRIANGLE + ($cUpperLower = "U" ? $__g_BLAS_STYPE_UPPER : $__g_BLAS_STYPE_LOWER))

	Return $mB
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_VectorToDiag()
; Description ...: creates a diagonal matrix from a vector
; Syntax ........: _la_VectorToDiag($mVector, [$bInPlace = False])
; Parameters ....: mVector  - [Map] vector as a map/array/definition string
;                  bInPlace - [Bool] (Default: False)
;                           ↳ overwrite variable mVector with the matrix (True) or return a new matrix (False)
; Return value ..: Success: bInPlace ? True : [Map] result diagonal matrix
;                  Failure: bInPlace ? Null : False and set @error to:
;                           | 1: invalid value for mVector
;                           | 1X: error X during _blas_copy() (@extended: @extended from _blas_copy())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _blas_copy()
; Link ..........:
; Example .......: Yes
;                  Global $mDiagMatrix = _la_VectorToDiag("[10,20,30,40,50]")
;                  _la_display($mDiagMatrix, "Diagonal matrix from vector")
;                  Global $mVector = _la_fromArray("[10,20,30,40,50]")
;                  _la_VectorToDiag($mVector, True)
;                  _la_display($mVector, "in-place converted vector")
; ===============================================================================================================================
Func _la_VectorToDiag(ByRef $mVector, $bInPlace = False)
	; direct AutoIt-type input
	If IsArray($mVector) Or IsString($mVector) Then $mVector = _blas_fromArray($mVector)

	; validation of the input parameters
	If Not MapExists($mVector, "ptr") Then Return SetError(1, 0, $bInPlace ? Null : False) ; valid AutoIt-BLAS/LAPACK-Map?

	Local $iN = $mVector.rows
	Local $mDiagMatrix = _blas_createMatrix($iN, $iN, $mVector.datatype)

	_blas_copy($mVector, 0, 1, 0, $iN + 1, $iN, $mDiagMatrix.ptr)
	If @error Then Return SetError(@error + 10, @extended, $bInPlace ? Null : False)

	If $bInPlace Then $mVector = $mDiagMatrix

	Return $bInPlace ? True : $mDiagMatrix
EndFunc   ;==>_la_VectorToDiag

#EndRegion

#Region data output

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_display()
; Description ...: displays a matrix/vector map, similar to _ArrayDisplay
; Syntax ........: _la_display($mData, [$sTitle = "", [iDecimalPlaces = 5, [$iFlags = 64]]])
; Parameters ....: mData          - [Map] matrix/vector as a map/array/definition string
;                  sTitle         - [String] (Default: "")
;                                 ↳ the window title to be displayed
;                  iDecimalPlaces - [UInt] (Default: 5)
;                                   number of decimal places to which the figures are to be rounded
;                  iFlags         - [UInt] (Default: 64)
;                                 ↳ display options - see $iFlags for _ArrayDisplay()
; Return value ..: Success: SetError(@error + 10, @extended, $mRet)
;                  Failure: $mRet and set @error to:
;                           | 1X: error X during _blas_display (@extended: @extended from _blas_display())
;                           | 1: invalid value for mData
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _blas_display()
; Link ..........:
; Example .......: Yes
;                  Global $mA = _la_fromArray('[[1,2,3],[4,5,6],[7,8,9]]')
;                  _la_display($mA, "from AutoIt Array")
; ===============================================================================================================================
Func _la_display($mData, $sTitle = "", $iDecimalPlaces = 5, $iFlags = 64)
	; direct AutoIt-type input
	If IsArray($mData) Or IsString($mData) Then $mData = _blas_fromArray($mData)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not MapExists($mData, "ptr") Then Return SetError(1, 0, False)

	Local $mRet = _blas_display($mData, $sTitle, $iDecimalPlaces, $iFlags)
	Return SetError(@error + 10, @extended, $mRet)
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_toArray()
; Description ...: converts a matrix/vector map into an AutoIt array
; Syntax ........: _la_toArray($mMatrix)
; Parameters ....: mMatrix - [Map] matrix/vector as a map/array/definition string
; Return value ..: Success: [Array] matrix/vector as a AutoIt 2D/1D array
;                  Failure: Null and set @error to:
;                           | 1: invalid value fro mMatrix
;                           | 1X: error X during _blas_toArray() (@extended: @extended from _blas_toArray())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _blas_toArray()
; Link ..........:
; Example .......: Yes
;                  Global $mA = _la_fromArray('[[1,2,3],[4,5,6],[7,8,9]]')
;                  Global $aArray = _la_toArray($mA)
;                  _ArrayDisplay($aArray, "reconstructed Array")
; ===============================================================================================================================
Func _la_toArray(Const ByRef $mMatrix)
	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, Null)

	Local $aRet = _blas_toArray($mMatrix)
	Return @error ? SetError(@error + 10, @extended, Null) : $aRet
EndFunc   ;==>_la_ToArray


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_toFile()
; Description ...: write a matrix/vector into a file
; Syntax ........: _la_toFile($mMatrix, $sFile, [$nMode = 0, [$sType = "DOUBLE", [$iKL = 0, [$iKU = 0]]]])
; Parameters ....: mMatrix - [Map] matrix/vector as a map/array/definition string
;                  sFile   - [String/Handle] file path or file handle
;                  nMode   - (Default: 0)
;                          ↳ matrix shape flags ($__g_BLAS_STYPEₓXX) in case an array or string was passed
;                  sType   - [String] (Default: "DOUBLE")
;                          ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
;                  iKL     - [UInt] (Default: 0)
;                          ↳ if nMode contains $__g_BLAS_STYPE_BAND: number of lower band diagonals
;                  iKU     - [UInt] (Default: 0)
;                          ↳ if nMode contains $__g_BLAS_STYPE_BAND: number of upper band diagonals
; Return value ..: Success: True (@extended = Number of written elements)
;                  Failure: False and set @error to:
;                           | 1: invalid value for mMatrix
;                           | 2: error during _blas_fromArray() (@extended: @error from _blas_fromArray())
;                           | 3: error during FileOpen() (@extended: @error from FileOpen())
;                           | 4: error during _blas_copy() (@extended: @error from _blas_copy())
;                           | 5: error during FileWrite() (@extended: @error from FileWrite())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _blas_fromArray(), _blas_copy()
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func _la_toFile($mMatrix, $sFile, $nMode = 0, Const $sType = "DOUBLE", $iKL = 0, $iKU = 0)
	Local $tStruct, $pStruct, $pData

	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix, $nMode, $sType, $iKL, $iKU)
	If @error Then Return SetError(2, @error, False)

	; validation of the input parameters
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, False) ; valid AutoIt-BLAS/LAPACK-Map?

	$tStruct = DllStructCreate(StringFormat("BOOLEAN type; UINT rows; UINT cols; USHORT flags; USHORT kl; USHORT ku; %s[%d]", $mMatrix.datatype, $mMatrix.size))
	$pStruct = DllStructGetPtr($tStruct)
	$pData = DllStructGetPtr($tStruct, 7)

	$tStruct.type = $mMatrix.datatype = "FLOAT" ? 1 : 0   ; 0 = "DOUBLE" / 1: "FLOAT"
	$tStruct.rows = $mMatrix.rows
	$tStruct.cols = $mMatrix.cols
	$tStruct.flags = $mMatrix.storageType
	$tStruct.kl = $mMatrix.kl
	$tStruct.ku = $mMatrix.ku

	; copy data to new structure:
	_blas_copy($mMatrix.ptr, 0, 1, 0, 1, $mMatrix.size, $pData, False, $mMatrix.datatype)
	If @error Then Return SetError(4, @error, False)

	; create output binary
	Local $tOut = DllStructCreate(StringFormat("BYTE[%d]", DllStructGetSize($tStruct)), $pStruct)
	Local $bOut = DllStructGetData($tOut, 1)

	; write to file
	Local $hFile = IsString($sFile) ? FileOpen($sFile, 2 + 16) : $sFile ; if already opened
	If @error Or $hFile = -1 Then Return SetError(3, @error, False)
	FileWrite($hFile, $bOut)
	If @error Then Return SetError(5, @error, False)
	FileClose($hFile)

	Return SetExtended($mMatrix.size, True)
EndFunc

#EndRegion

#Region scalar operations


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_rotate()
; Description ...: applies a plane rotation to coordinate-pairs
; Syntax ........: _la_rotate($mX, $mY, $fAlpha)
; Parameters ....: mX     - [Map] vector with x-coordinates as a map/array/definition string (gets overwritten)
;                  mY     - [Map] vector with y-coordinates as a map/array/definition string (gets overwritten)
;                  fAlpha - [Float] rotation angle [radian]
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: invalid value for mX (@extended = 1) or my (@extended = 2)
;                           | 1X: error X during _blas_rot() (@extended: @extended from _blas_rot())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _blas_rot()
; Link ..........:
; Example .......: Global $mX = _la_fromArray("[1,2,3,4,5]")
;                  Global $mY = _la_fromArray("[-1,-2,-3,-4,-5]")
;                  _la_rotate($mX, $mY, ACos(-1) / 4)
;                  _la_display($mX, "rotated x-coordinates")
;                  _la_display($mY, "rotated y-coordinates")
; ===============================================================================================================================
Func _la_rotate(ByRef $mX, ByRef $mY, $fAlpha)
	; direct AutoIt-type input
	If IsArray($mX) Or IsString($mX) Then $mX = _blas_fromArray($mX)
	If IsArray($mY) Or IsString($mY) Then $mY = _blas_fromArray($mY)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not (IsMap($mX) And MapExists($mX, "ptr")) Then Return SetError(1, 1, False)
	If Not (IsMap($mY) And MapExists($mY, "ptr")) Then Return SetError(1, 2, False)

	_blas_rot($mX, $mY, $fAlpha)
	Return @error ? SetError(@error + 10, @extended, False) : True
EndFunc   ;==>_la_rotate

#EndRegion


#Region matrix attributes

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_isPositiveDefinite()
; Description ...: checks whether a matrix is positive definite
; Syntax ........: _la_isPositiveDefinite($mA)
; Parameters ....: mA - [Map] matrix/vector as a map/array/definition string
; Return value ..: Success: True (positive definite) or False (not positive definite)
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mA
;                           | 1/2: error during _lp_potrf()
;                           | 1X: error X during _la_isSymmetric() (@extended: @extended from _la_isSymmetric())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _la_isSymmetric(), _lp_potrf()
; Link ..........:
; Example .......: Yes
;                  ConsoleWrite(_la_isPositiveDefinite("[[4,-1,1],[-1,4,-2],[1,-2,5]]") & @TAB & @extended & @CRLF)
; ===============================================================================================================================
Func _la_isPositiveDefinite($mA)
	; direct AutoIt-type input
	If IsArray($mA) Or IsString($mA) Then $mA = _blas_fromArray($mA)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not (IsMap($mA) And MapExists($mA, "ptr")) Then Return SetError(1, 0, Null)

	; check for quadratic matrix
	If $mA.rows <> $mA.cols Then Return SetExtended(1, False)

	; check for symmetry
	Local $bSymmetric = _la_isSymmetric($mA)
	If @error Then Return SetError(@error + 10, @extended, Null)
	If Not $bSymmetric Then Return SetExtended(2, False)

	; try cholesky factorization - if success then A must be positive definite
	_lp_potrf(_la_duplicate($mA))

	Switch @error
		Case 0
			Return True

		Case 2 ; error in INFO from potrf
			Return @extended > 0 ? SetExtended(3, False) : SetError(2, @extended, NULL)

		Case Else ; error in _lp_potrf
			Return SetError(1, @error, NULL)

	EndSwitch
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_isSymmetric()
; Description ...: checks whether a matrix is symmetrical
; Syntax ........: _la_isSymmetric($mA, [$fTolerance = Default])
; Parameters ....: mA         - [Map] matrix/vector as a map/array/definition string
;                  fTolerance - [Float] (Default: Default)
;                             ↳ threshold value above which a difference is regarded as 0 (values are regarded as equal)
;                               Default: 10× machine precision
; Return value ..: Success: True: symmetric, False: non symmetric
;                  Failure: Null and set @error to:
;                           | 1: invalid value for $mA
;                           |1X: error X during _blas_axpy() (@extended: @extended from _blas_axpy())
;                           |2X: error X during _blas_amax() (@extended: @extended from _blas_amax())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......: It would be more efficient to consider only one triangular matrix, as the differences calculated by axpy are also symmetrical.
;                  Also amax would then only have to iterate over the triangular matrix.
;                  Unfortunately, I am not aware of an efficient pure BLAS/LAPACK solution.
; Related .......: _blas_axpy(), _blas_amax()
; Link ..........:
; Example .......: Yes
;                  Global $bSym = _la_isSymmetric('[[611,197,-192,407,-8,-52,-49,29],[197,899,113,-192,-71,-43,-8,-44],[-192,113,899,196,61,49,8,52],[407,-192,196,611,8,44,59,-23],[-8,-71,61,8,411,-599,208,208],[-52,-43,49,44,-599,411,208,208],[-49,-8,8,59,208,208,99,-911],[29,-44,52,-23,208,208,-911,99]]')
;                  ConsoleWrite($bSym & @CRLF)
; ===============================================================================================================================
Func _la_isSymmetric($mA, $fTolerance = Default)
	; direct AutoIt-type input
	If IsArray($mA) Or IsString($mA) Then $mA = _blas_fromArray($mA)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not (IsMap($mA) And MapExists($mA, "ptr")) Then Return SetError(1, 0, Null)

	; calculate tolerance d of elements if not defined
	If IsKeyword($fTolerance) = 1 Then $fTolerance = 10 * ($mA.datatype = "FLOAT" ? $f_LA_FLT_PREC : $f_LA_DBL_PREC)

	; transpose A and write to new Matrix AT
	Local $mAT = _blas_omatcopy($mA, 1, "T")

	; calculate  AT <-- A - AT   (= A - Aᵀ)
	_blas_axpy($mA, $mAT, -1, 0, 0, 1, 1, $mA.elements)
	If @error Then Return SetError(10 + @error, @extended, Null)

	; calculate the euclidian norm of the result (0-matrix should return ~ 0 )
	Local $fAmax = Abs(_blas_amax($mAT))
	If @error Then Return SetError(20 + @error, @extended, Null)

	; compare to the threshold
	Return $fAmax < $fTolerance
EndFunc



; #FUNCTION# ====================================================================================================================
; Name ..........: _la_rank()
; Description ...: determines the rank of a matrix
; Syntax ........: _la_rank($mMatrix, [$fTolerance = Default, [$bInPlace = False]])
; Parameters ....: mMatrix    - [Map] matrix/vector as a map/array/definition string
;                  fTolerance - [Float] (Default: Default)
;                             ↳ Maximum absolute value up to which a number is considered to be 0 (threshold value)
;                  bInPlace   - [Bool] (Default: False)
;                             ↳ True: mMatrix gets overwritten
;                               False: mMatrix remains untouched
; Return value ..: Success: [UInt] rank of the matrix
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mMatrix
;                           |1X: error X during _lp_gesvd() (@extended: @extended from _lp_gesvd())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _lp_gesvd()
; Link ..........:
; Example .......: Yes
;                  Global $dRank = _la_rank("[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97],[93,89,13,-53,43]]")
;                  MsgBox(0,"rank", $dRank)
; ===============================================================================================================================
Func _la_rank($mMatrix, $fTolerance = Default, $bInPlace = False)
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, Null)

	; prevent overwrite if choosed
	If Not $bInPlace Then $mMatrix = _la_duplicate($mMatrix)

	; do a SVD to obtain the singular values
	Local $mSVD = _lp_gesvd($mMatrix)
	If @error Then Return SetError(10 + @error, @extended, Null)

	Local $tS = $mSVD.S.struct, $iRank = 0

	; calculate tolerance d of elements if not defined
	If IsKeyword($fTolerance) = 1 Then $fTolerance = DllStructGetData($tS, 1, 1) * _Max($mMatrix.rows, $mMatrix.cols) * ($mMatrix.datatype = "FLOAT" ? $f_LA_FLT_PREC : $f_LA_DBL_PREC) ; DllStructGetData($tS, 1, 1) = Max(S) because S is sorted

	; count the singular values above the tolerance threshold
	For $i = 1 To $mSVD.S.elements
		If DllStructGetData($tS, 1, $i) > $fTolerance Then $iRank += 1
	Next

	Return $iRank
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_determinant()
; Description ...: calculate the determinant of a matrix
; Syntax ........: _la_determinant($mMatrix, [$bInPlace = False])
; Parameters ....: mMatrix  - [Map] matrix as a map/array/definition string
;                  bInPlace - [Bool] (Default: False)
;                           ↳ True: mMatrix gets overwritten
;                             False: mMatrix remains untouched
; Return value ..: Success: [Float] determinant of mMatrix
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mMatrix
;                           |1X: error X during _lp_getrf() (@extended: @extended from _lp_getrf())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _lp_getrf()
; Link ..........:
; Example .......: Yes
;                  Global $fDet = _la_determinant('[[1,4,5,-1],[-2,3,-1,0],[2,1,1,0],[3,-1,2,1]]') ; --> -40
;                  MsgBox(0, "Determinant", $fDet)
; ===============================================================================================================================
Func _la_determinant(ByRef $mMatrix, $bInPlace = False)
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, Null)

	; prevent overwrite if choosed
	If Not $bInPlace Then $mMatrix = _la_duplicate($mMatrix)

	Local $iN = $mMatrix.rows

	; run LU decomposition
	Local $tIPIV = _lp_getrf($mMatrix)
	If @error Then Return SetError(10 + @error, @error, Null)

	; calculate the absolute value of the determinant (product of diagonal elements) and it`s sign (det(P))
	Local $fDet = 1, $iSign = 1
	Local $tMatrix = $mMatrix.struct, $iMatrix = 1
	For $i = 1 To $iN
		$fDet *= DllStructGetData($tMatrix, 1, $iMatrix)
		$iMatrix += $iN + 1
		$iSign *= DllStructGetData($tIPIV, 1, $i) <> $i ? -1 : 1
	Next

	Return $fDet * $iSign
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_conditionNumber()
; Description ...: determine the condition number of a matrix
; Syntax ........: _la_conditionNumber($mMatrix, [$cNorm = "1", [$bInPlace = False]])
; Parameters ....: mMatrix  - [Map] matrix as a map/array/definition string
;                  cNorm    - [Char] (Default: "1")
;                           ↳ "1": use 1-norm (max column sum)
;                             "I": use infinity norm (max row sum)
;                  bInPlace - [Bool] (Default: False)
;                           ↳ True: mMatrix gets overwritten
;                             False: mMatrix remains untouched
; Return value ..: Success: [Float] condition number
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mMatrix
;                           |1X: error X during _lp_lange() (@extended: @extended from _lp_lange())
;                           |2X: error X during _lp_getrf() (@extended: @extended from _lp_getrf())
;                           |3X: error X during _lp_gecon() (@extended: @extended from _lp_gecon())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _lp_lange(), _lp_getrf(), _lp_gecon()
; Link ..........:
; Example .......: Yes
;                  Global $fCond = _la_conditionNumber('[[1,4,5,-1],[-2,3,-1,0],[2,1,1,0],[3,-1,2,1]]')
;                  MsgBox(0, "Condition Number", $fCond)
; ===============================================================================================================================
Func _la_conditionNumber($mMatrix, $cNorm = "1", $bInPlace = False)
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, Null)

	; prevent overwrite if choosed
	If Not $bInPlace Then $mMatrix = _la_duplicate($mMatrix)

	Local $iM = $mMatrix.rows, $iN = $mMatrix.cols, $iMinMN = $iM < $iN ? $iM : $iN, _
	      $sDataType = $mMatrix.datatype

	; calculate the norm of A
	Local $fANORM = _lp_lange($mMatrix.ptr, $cNorm = "I" ? "I" : "1", $iM, $iN, $iM, $sDataType)
	If @error Then Return SetError(@error + 10, @extended, Null)

	; LU factorization of A
	_lp_getrf($mMatrix.ptr, $iM, $iN, $iM, $sDataType)
	If @error Then Return SetError(@error + 20, @extended, Null)

	Local $fCondInv = _lp_gecon($mMatrix.ptr, $fANORM, $cNorm = "I" ? "I" : "1", $iMinMN, $iMinMN, $sDataType)
	If @error Then Return SetError(@error + 30, @extended, Null)

	Return 1 / $fCondInv
EndFunc

#EndRegion

#Region unary operations

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_inverse()
; Description ...: calculates the inverse of a matrix
; Syntax ........: _la_inverse($mMatrix, [$bInPlace = False])
; Parameters ....: mMatrix  - [Map] matrix as a map/array/definition string
;                  bInPlace - [Bool] (Default: False)
;                           ↳ True: mMatrix gets overwritten
;                             False: mMatrix remains untouched
; Return value ..: Success: bInPlace ? True : [Map] inverse of mMatrix as a map
;                  Failure: bInPlace ? False : Null
;                           | 1: invalid value for mMatrix
;                           |1X: error X during _lp_getrf() (@extended: @extended from _lp_getrf())
;                           |2X: error X during _lp_getri() (@extended: @extended from _lp_getri())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......: Current via LU decomposition. But still offers a lot of room for optimization for special cases of certain matrix geometries.
; Related .......: _lp_getrf(), _lp_getri
; Link ..........:
; Example .......: Yes
;                  Global $mInverse = _la_inverse('[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97],[93,89,13,-53,43]]')
;                  _la_display($mInverse)
; ===============================================================================================================================
Func _la_inverse($mMatrix, $bInPlace = False)

	; ToDo: Man könnte auch einfach das Gleichungssystem  A * A⁻¹ = I nach A⁻¹ lösen.
	; Ich vermute, dass dgetri nur Sinn macht, wenn man eh eine LU-Faktorisierung (dgetrf) durchführen muss.
	; Hierbei noch schauen ob auch m x n Matrizen damit gehen.
	; Laufzeiten:
	; - dgetrf ( O(2/3 n³) ) + dgetri ( O(4/3 n³) ) = O(2 n³)				= 8/3 n³ FLOPS
	; - dgeqrf ( O(4/3 n³) ) + dgeqrs ( O(4/3 n³) ) = O(2 2/3 n³)			= 10/3 n³ FLOPS
	; - dgesvd = O(8 n³)    (4n³ SVD + 4n³ Inverse aus Ergebnis berechnen)	= 5 1/3 n³ bis 10 1/3 n³
	; - dgesv = O(2 2/3 n³)   (2/3n³ LU-Zerlegung + 2n³ System lösen)		= 8/3 n³ FLOPS
	; - dposv = O(n³)    (1/3n³ Cholesky + 2/3n³ System lösen)				= 4/3 n³ FLOPS
	; - dsytrf ( O(1/3n³) ) + dsytri ( O(2/3n³) ) = O(n³)			(gleich dsysv - jedoch mit expliziter LDLᵀ-Zerlegung)
	; - dsysv = O(n³)														= 7/3 n³ FLOPS
	; - dgbsv = 8/3 n*k²  (2/3 nk² LU-Zerlegung + 2nk² System lösen)		= 2n(2m + 1)²
	; - dtrtrs = 1/3n³														= n²

	; Die Methode per SVD:
	; - Σ⁻¹: Σ invertieren: Einfach Reziproke der einzelnen Werte, da Σ = Diagonalmatrix bzw. noch Vektor
	; - A⁻¹ = V Σ⁻¹ Uᵀ				(VΣ⁻¹ = Tmp per dtrmm, dann Tmp * Uᵀ per dgemm)

	; Pseudoinverse per SVD:
	; - ist im Grunde genau gleich, nur statt der Inverse von Σ⁻¹ muss die Pseudoinverse Σ⁺ berechnet werden.
	;   das wird genauso gemacht (Reziproke der Werte bilden) aber die Nullwerte bleiben hier unangetastet - bleiben also 0 (Grenzwerte beachten)
	; - Dann wieder genauso ausmultiplizieren:
	;   A⁺ = V Σ⁺ Uᵀ

	; Es gäbe noch die Variante dgeqrf + dgeqrs (noch zu implementieren). Diese wäre numerisch stabiler aber wohl langsamer als dgetrf + dgetri

	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, $bInPlace ? False : Null)

	; prevent overwrite if choosed
	If Not $bInPlace Then $mMatrix = _la_duplicate($mMatrix)

	; first: calculate LU decomposition
	Local $tIPIV = _lp_getrf($mMatrix)
	If @error Then Return SetError(10 + @error, @extended, $bInPlace ? False : Null)

	; calculate inverse:
	_lp_getri($mMatrix, $tIPIV)
	Return @error ? SetError(20 + @error, @extended, $bInPlace ? False : Null) : ($bInPlace ? True : $mMatrix)
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_pseudoInverse()
; Description ...: calculate the Moore-Penrose pseudo inverse of a matrix
; Syntax ........: _la_pseudoInverse($mMatrix, [$fTolerance = Default, [$bOverwrite = False]])
; Parameters ....: mMatrix    - [Map] matrix as a map/array/definition string
;                  fTolerance - [Float] (Default: Default)
;                             ↳ threshold value up to which an absolute value is regarded as 0 (singular value or not)
;                  bOverwrite - [Bool] (Default: False)
;                             ↳ True: mMatrix gets overwritten
;                               False: mMatrix remains untouched
; Return value ..: Success: [Map] pseudo inverse of mMatrix as a map
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mMatrix
;                           |1X: error X during _lp_gesvd() (@extended: @extended from _lp_gesvd() )
;                           |2X: error X during _blas_copy() (@extended: @extended from _blas_copy() )
;                           |3X: error X during 1st _blas_gemm() (@extended: @extended from _blas_gemm() )
;                           |4X: error X during 2nd _blas_gemm() (@extended: @extended from _blas_gemm() )
; Author ........: AspirinJunkie
; Modified.......: 2024-09-26
; Remarks .......: algorithm:
;                  - calculate SVD decomposition:  A = U * Σ * Vᵀ
;                  - calculate Σ⁺ by invert every element Σᵢ of Σ if |Σᵢ| > 0
;                  - calculate A⁺ = V Σ⁺ Uᵀ
; Related .......: _lp_gesvd()
; Link ..........:
; Example .......: Yes
                  Global $mInverse = _la_pseudoInverse('[[1,1,1,1],[5,7,7,9]]') ; --> [[2, -0.25], [0.25, 0], [0.25, 0], [-1.5, 0.25]]
                  _la_display($mInverse)
; ===============================================================================================================================
Func _la_pseudoInverse($mMatrix, $fTolerance = Default, $bOverwrite = False)
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, $bOverwrite ? False : Null)

	; determine datatype and shape of the matrices
	Local $sDataType = $mMatrix.datatype, _
	      $iM = $mMatrix.rows, $iN = $mMatrix.cols, _
	      $iMS = $iM, $iNS = $iN, _ ; dimension of Sigma
	      $iMU = $iM, $iNU = $iM, _ ; dimension of U
	      $iNV = $iN    ; dimension of V

	; calculate tolerance
	If IsKeyword($fTolerance) = 1 Then $fTolerance = 10 * ($sDataType = "FLOAT" ? $f_LA_FLT_PREC : $f_LA_DBL_PREC)

	; prevent overwrite if choosed
	If Not $bOverwrite Then $mMatrix = _la_duplicate($mMatrix)

	Local $mSVD = _lp_gesvd($mMatrix)
	If @error Then Return SetError(@error + 10, @extended, Null)

	; calculate Σ⁺    ( |Σᵢ| > 0 ? 1 / Σᵢ : 0  )
	Local $mS = $mSVD.S
	Local $tSigma = $mS.struct, $fValue
	; ToDo: if anyone has a idea how to do this (especially with the threshold) directly with BLAS/LAPACK - do it!
	For $i = 1 To $iM
		$fValue = DllStructGetData($tSigma, 1, $i)
		DllStructSetData($tSigma, 1, Abs($fValue) > $fTolerance ? 1.0 / $fValue : $fValue, $i)
	Next

	; Σ⁺ vector to diagonal matrix
	Local $mSigmaPlus = _blas_createMatrix($iMS, $iNS, $sDataType)
	_blas_copy($tSigma, 0, 1, 0, $iMS + 1, $iMS, $mSigmaPlus.ptr)
	If @error Then Return SetError(20 + @error, @extended, Null)

	; V * Σ⁺ --> C
	Local $mC = _blas_createMatrix($iNV, $iNS, $sDataType)
	_blas_gemm($mSVD.VT.ptr, $mSigmaPlus.ptr, $mC.ptr, 1.0, 0.0, "T", "N", $iNV, $iNS, $iMS > $iNS ? $iNS : $iMS, $iNV, $iMS, $iNV, $sDataType)
	If @error Then Return SetError(30 + @error, @extended, Null)

	; C * Uᵀ
	Local $mInverse = _blas_createMatrix($iNV, $iMU, $sDataType)
	_blas_gemm($mC.ptr, $mSVD.U.ptr, $mInverse.ptr, 1.0, 0.0, "N", "T", $iNV, $iMU, $iNS, $iNV, $iNU, $iNV, $sDataType)
	If @error Then Return SetError(40 + @error, @extended, Null)

	Return $mInverse
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_sum()
; Description ...: calculates the sum of the elements of a matrix, vector or parts thereof
; Syntax ........: _la_sum($mMatrix, [$iStart = 0, [$iInc = 1, [$iN = Default, [$sDatatype = "DOUBLE"]]]])
; Parameters ....: mMatrix   - [Map] matrix/vector as a map/array/definition string
;                  iStart    - [Int] (Default: 0)
;                            ↳ start element index (0-based)
;                  iInc      - [Int] (Default: 1)
;                            ↳ storage spacing between elements of mMatrix (can be used to handle parts of a matrix as a vector)
;                  iN        - [Int] (Default: Default)
;                            ↳ number of elements in input vector
;                  sDatatype - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: [Float] the sum {@extended: number of elements used}
;                  Failure: False and set @error to:
;                           | 1: invalid value for sDatatype
;                           | 2: invalid value for mMatrix
;                           | 3: error during DllCall of dot (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-09
; Remarks .......: Actually just a workaround. But I don't know of a more effective method using BLAS/LAPACK alone
; Related .......:
; Link ..........:
; Example .......: Yes
;                  Global $fSum = _la_sum("[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]")
;                  ConsoleWrite("Sum: " & $fSum & @CRLF & "Mean: " & $fSum / @extended & @CRLF)
; ===============================================================================================================================
Func _la_sum($mMatrix, $iStart = 0, $iInc = 1, $iN = Default, $sDatatype = "DOUBLE")
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix, 0, $sDatatype)
	; validation of the input parameters
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, Null)
	If $sDatatype <> "DOUBLE" And $sDatatype <> "FLOAT" Then Return SetError(2, 0, Null)

	If IsKeyword($iN) = 1 Then $iN = $mMatrix.elements

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d", _
	            $dSize   = ($sDataType = "DOUBLE") ? $iBLAS_SIZE_DOUBLE : $iBLAS_SIZE_FLOAT

	Local $aDLL = DllCall($__g_hBLAS_DLL, $sDataType & ":cdecl", $cPrefix & "dot", _
		"INT*",           $iN, _                             ; number of elements in vectors x and y
		"PTR",            $mMatrix.ptr + $dSize * $iStart, _ ; start ptr of x
		"INT*",           $iInc, _                           ; increment in x
		$sDataType & "*", 1.0, _                             ; start ptr of y
		"INT*",           0 _                                ; increment in y
	)
	Return @error ? SetError(3, @error, Null) : SetExtended(Ceiling($iN / $iInc) - $iStart, $aDLL[0])
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_asum()
; Description ...: calculate the sum of the absolute(!) values of a matrix/vector
; Syntax ........: _la_asum($mMatrix, [$iStart = 0, [$iStep = 1, [$iN = Default]]])
; Parameters ....: mMatrix - [Map] matrix/vector as a map/array/definition string
;                  iStart  - [Int] (Default: 0)
;                          ↳ start element index (0-based)
;                  iStep   - [UInt] (Default: 1)
;                          ↳ storage spacing between elements of mMatrix (can be used to handle parts of a matrix as a vector)
;                  iN      - [Int] (Default: Default)
;                          ↳ number of elements in input vector
; Return value ..: Success: [Float] absolute sum of the values
;                  Failure: Null and set @error to:
;                           |1X: error X during _blas_asum() (@extended: @extended from _blas_asum())
;                           | 1: invalid value for mMatrix
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _blas_asum()
; Link ..........:
; Example .......: Yes
;                  Global $fSum = _la_asum("[1,2,3,4,5]")
;                  MsgBox(0,"absolute sum", $fSum)
; ===============================================================================================================================
Func _la_asum(ByRef $mMatrix, $iStart = 0, $iStep = 1, $iN = Default)
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, Null)

	; calculate number of elements if not defined
	If IsKeyword($iN) = 1 Then $iN = Floor(($mMatrix.elements - $iStart) / $iStep)

	Local $fSum = _blas_asum($mMatrix, $iStart, $iStep, $iN)
	Return @error ? SetError(@error + 10, @extended, Null) : $fSum
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_amin()
; Description ...: finds the first element having the minimum absolute(!) value
; Syntax ........: _la_amin($mMatrix, [$iStart = 0, [$iStep = 1, [$iN = Default]]])
; Parameters ....: mMatrix - [Map] matrix/vector as a map/array/definition string
;                  iStart  - [Int] (Default: 0)
;                          ↳ start element index (0-based)
;                  iStep   - [UInt] (Default: 1)
;                          ↳ storage spacing between elements of mMatrix (can be used to handle parts of a matrix as a vector)
;                  iN      - [Int] (Default: Default)
;                          ↳ number of elements in input vector
; Return value ..: Success: [Float] minimum absolute value (@extended: index of the elements)
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mMatrix
;                           |1X: error X during _blas_amin() (@extended: @extended from _blas_amin())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _blas_amin()
; Link ..........:
; Example .......: Yes
;                  Global $fMin = _la_amin("[3,2,-3,4,5]")
;                  MsgBox(0,"absolute min", $fMin)
; ===============================================================================================================================
Func _la_amin(ByRef $mMatrix, $iStart = 0, $iStep = 1, $iN = Default)
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, Null)

	; calculate number of elements if not defined
	If IsKeyword($iN) = 1 Then $iN = Floor(($mMatrix.elements - $iStart) / $iStep)

	Local $fMin = _blas_amin($mMatrix, $iStart, $iStep, $iN)
	Return @error ? SetError(@error + 10, @extended, Null) : SetExtended(@extended, $fMin)
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_amax()
; Description ...: finds the first element having the maximum absolute(!) value
; Syntax ........: _la_amax($mMatrix, [$iStart = 0, [$iStep = 1, [$iN = Default]]])
; Parameters ....: mMatrix - [Map] matrix/vector as a map/array/definition string
;                  iStart  - [Int] (Default: 0)
;                          ↳ start element index (0-based)
;                  iStep   - [UInt] (Default: 1)
;                          ↳ storage spacing between elements of mMatrix (can be used to handle parts of a matrix as a vector)
;                  iN      - [Int] (Default: Default)
;                          ↳ number of elements in input vector
; Return value ..: Success: [Float] maximum absolute value (@extended: index of the elements)
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mMatrix
;                           |1X: error X during _blas_amax() (@extended: @extended from _blas_amax())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _blas_amax()
; Link ..........:
; Example .......: Yes
;                  Global $fMax = _la_amax(_la_fromArray("[-3,-5,-7,-4,3,2,-3]"))
;                  MsgBox(0,"absolute max", $fMax)
; ===============================================================================================================================
Func _la_amax(ByRef $mMatrix, $iStart = 0, $iStep = 1, $iN = Default)
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, Null)

	; calculate number of elements if not defined
	If IsKeyword($iN) = 1 Then $iN = Floor(($mMatrix.elements - $iStart) / $iStep)

	Local $fMax = _blas_amax($mMatrix, $iStart, $iStep, $iN)
	Return @error ? SetError(@error + 10, @extended, Null) : SetExtended(@extended, $fMax)
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_norm()
; Description ...: calculate the euclidian norm of a vector
; Syntax ........: _la_norm($mMatrix, [$iStart = 0, [$iStep = 1, [$iN = Default]]])
; Parameters ....: mMatrix - [Map] vector as a map/array/definition string
;                  iStart  - [Int] (Default: 0)
;                          ↳ start element index (0-based)
;                  iStep   - [UInt] (Default: 1)
;                          ↳ storage spacing between elements of mMatrix (can be used to handle parts of a matrix as a vector)
;                  iN      - [Int] (Default: Default)
;                          ↳ number of elements in input vector
; Return value ..: Success: [Float] euclidian norm of mMatrix
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mMatrix
;                           |1X: error X during _blas_nrm2() (@extended: @extended from _blas_nrm2())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _blas_nrm2()
; Link ..........:
; Example .......: Yes
;                  Global $fNorm = _la_norm(_la_fromArray("[3,2,7,3,4,5]"))
;                  MsgBox(0,"Norm", $fNorm)
; ===============================================================================================================================
Func _la_norm(ByRef $mMatrix, $iStart = 0, $iStep = 1, $iN = Default)
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, Null)

	; calculate number of elements if not defined
	If IsKeyword($iN) = 1 Then $iN = Floor(($mMatrix.elements - $iStart) / $iStep)

	Local $fNorm = _blas_NRM2($mMatrix, $iStart, $iStep, $iN)
	Return @error ? SetError(@error + 10, @extended, Null) : $fNorm
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_mean()
; Description ...: calculate the mean of a vector or parts of a matrix
; Syntax ........: _la_mean($mMatrix, [$iStart = 0, [$iStep = 1, [$iN = Default]]])
; Parameters ....: mMatrix - [Map] matrix/vector as a map/array/definition string
;                  iiStart  - [Int] (Default: 0)
;                          ↳ start element index (0-based)
;                  iStep   - [UInt] (Default: 1)
;                          ↳ storage spacing between elements of mMatrix (can be used to handle parts of a matrix as a vector)
;                  iN      - [Int] (Default: Default)
;                          ↳ number of elements in input vector
; Return value ..: Success: [Float] mean value
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mMatrix
;                           |1X: error X during _la_sum() (@extended: @extended from _la_sum())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _la_sum()
; Link ..........:
; Example .......: Yes
;                  Global $fMean = _la_mean("[3,2,7,3,4,5]")
;                  MsgBox(0,"Mean", $fMean)
; ===============================================================================================================================
Func _la_mean($mMatrix, $iStart = 0, $iStep = 1, $iN = Default)
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, Null)

	; calculate number of elements if not defined
	If IsKeyword($iN) = 1 Then $iN = Floor(($mMatrix.elements - $iStart) / $iStep)

	Local $fSum = _la_sum($mMatrix, $iStart, $iStep, $iN, $mMatrix.datatype)
	Return @error ? SetError(@error + 10, @extended, Null) : $fSum / @extended
EndFunc


#EndRegion

#Region element wise operations


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_sqrtElements()
; Description ...: calculates the square root of each element of a matrix/vector
; Syntax ........: _la_sqrtElements($mMatrix, [$bInPlace = False])
; Parameters ....: mMatrix  - [Map] matrix/vector as a map/array/definition string
;                  bInPlace - [Bool] (Default: False)
;                           ↳ True: mMatrix gets overwritten
;                             False: mMatrix remains untouched
; Return value ..: Success: bInplace ? True : matrix/vector with square rooted values
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mMatrix
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......: If value < 0 then sqrt(Abs(value)) is used instead (without abs() it is not defined)
; Related .......: _lp_pbtrf()
; Link ..........:
; Example .......: Yes
;                  Global $mA = _la_fromArray("[[9,0,0,0,0],[0,16,0,0,0],[0,0,25,0,0],[0,0,0,36,0],[0,0,0,0,49]]")
;                  _la_sqrtElements($mA, True)
;                  _la_display($mA)
; ===============================================================================================================================
Func _la_sqrtElements($mMatrix, $bInPlace = False)
; idea: cholesky factorization on a diagonal matrix results in the sqrt-values.
; so thread the input matrix/vector elements as a diagonal matrix
; for efficiency use the compressed band matrix form for cholesky with 0 sub/super diagonals

	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, Null)

	; duplicate input matrix to prevent overwrite if option choosed
	If Not $bInPlace Then $mMatrix = _la_duplicate($mMatrix)

	Local $iN = $mMatrix.size

	_lp_pbtrf($mMatrix.ptr, "L", 0, $iN, 1, $mMatrix.datatype)
	If @error Then ; switch to manual method
		Local $bStart = False, $tM = $mMatrix.struct
		For $i = 1 To $iN
			If Not $bStart Then
				If DllStructGetData($tM, 1, $i) <= 0 Then
					$bStart = True
					DllStructSetData($tM, 1, Sqrt(Abs(DllStructGetData($tM, 1, $i))), $i)
				EndIf
			Else
				DllStructSetData($tM, 1, Sqrt(Abs(DllStructGetData($tM, 1, $i))), $i)
			EndIf
		Next
		Return SetExtended(@error, $bInPlace ? True : $mMatrix)
	EndIf

	Return SetExtended($iN, $bInPlace ? True : $mMatrix)
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_squareElements()
; Description ...: calculates the square of each element of a matrix/vector
; Syntax ........: _la_squareElements($mMatrix, [$bInPlace = False])
; Parameters ....: mMatrix  - [Map] matrix/vector as a map/array/definition string
;                  bInPlace - [Bool] (Default: False)
;                           ↳ True: mMatrix gets overwritten
;                             False: mMatrix remains untouched
; Return value ..: Success: bInplace ? True : matrix/vector with square rooted values
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mMatrix
;                           |1X: error X during _blas_sbmv() (@extended: @extended from _blas_sbmv())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _blas_sbmv()
; Link ..........:
; Example .......: Yes
;                  Global $mA = _la_fromArray("[[1,-2,3],[4,5,6],[7,8,9]]")
;                  Global $mS = _la_squareElements($mA)
;                  _la_display($mS)
; ===============================================================================================================================
Func _la_squareElements(ByRef $mMatrix, $bInPlace = False)
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, Null)

	; duplicate input matrix to prevent overwrite if option choosed
	If Not $bInPlace And IsMap($mMatrix) Then $mMatrix = _la_duplicate($mMatrix)

	Local $iN = $mMatrix.size

	; empty (0-filled) vector Y
	Local $tTmp = DllStructCreate(StringFormat("%s[%d]", $mMatrix.datatype, $iN)), $pTmp = DllStructGetPtr($tTmp)

	; sbmv to treat the vector simply as diagonal matrix to perform a efficient element wise multiplication
	_blas_sbmv($mMatrix.ptr, $mMatrix.ptr, $pTmp, 1, 1, 0, "L", $iN, 1, 1, 1)
	If @error Then Return SetError(@error + 10, @extended, $bInPlace ? False : Null)

	If $bInPlace Then
		$mMatrix.struct = $tTmp
		$mMatrix.ptr = $pTmp
		Return True
	Else
		Local $mRet = $mMatrix.storageType = 0 ? _blas_createVector($iN, $mMatrix.datatype) : _blas_createMatrix($mMatrix.rows, $mMatrix.cols, $mMatrix.datatype)
		$mRet.struct = $tTmp
		$mRet.ptr = $pTmp
		Return $mRet
	EndIf

EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_invElements()
; Description ...: forms the reciprocal (1/x) for each element of the matrix/vector
; Syntax ........: _la_invElements($mMatrix, [$bInPlace = False])
; Parameters ....: mMatrix  - [Map] matrix/vector as a map/array/definition string
;                  bInPlace - [Bool] (Default: False)
;                           ↳ True: mMatrix gets overwritten
;                             False: mMatrix remains untouched
; Return value ..: Success: bInplace ? True : matrix/vector with inverted values
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mMatrix
;                           |1X: error X during _blas_tbsv() (@extended: @extended from _blas_tbsv())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......: idea: solve A*x = b with b = vector(n) filled with ones and A = diagonal matrix with elements from $mMatrix as diagonal
;                  use tbsv because a diag-matrix is saved as simple diag vector in band matrix form with k=0
; Related .......: _blas_tbsv()
; Link ..........:
; Example .......: Yes
;                  Global $mRet = _la_invElements("[[1,1e-308,3],[4,0,6]]")
;                  _la_display($mRet)
; ===============================================================================================================================
Func _la_invElements(ByRef $mMatrix, $bInPlace = False)
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, Null)

	; duplicate input matrix to prevent overwrite if option choosed
	If Not $bInPlace And IsMap($mMatrix) Then $mMatrix = _la_duplicate($mMatrix)

	Local $iR = $mMatrix.rows, $iC = $mMatrix.cols, $iN = $mMatrix.size, $bIsMatrix = False

	; convert matrix into vector for better handling
	If $mMatrix.storageType <> 0 Then
		$bIsMatrix           = True
		$mMatrix.rows        = $iN
		$mMatrix.cols        = 1
		$mMatrix.storageType = 0
	EndIf

	; create vector filled with ones with size as $mMatrix
	Local $mOnes = _la_createIdentity($iN)

	; use tbsv to efficiently solve a system with A = diagonal matrix
	_blas_tbsv($mMatrix.ptr, $mOnes.ptr, 0, "U", "N", "N", $iN, 1, 1)
	If @error Then Return SetError(@error + 10, @extended, $bInPlace ? False : Null)

	; reshape vector into matrix
	If $bIsMatrix Then
		$mOnes.rows        = $iR
		$mOnes.cols        = $iC
		$mOnes.storageType = $__g_BLAS_STYPE_MATRIX
	EndIf

	If $bInPlace Then
		$mMatrix = $mOnes
		Return True
	Else
		Return $mOnes
	EndIf
EndFunc

#EndRegion



#Region addition subtraction

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_sub()
; Description ...: subtracts a matrix/vector B from matrix/vector A
; Syntax ........: _la_sub($mA, $mB, [$bInPlace = False])
; Parameters ....: mA       - [Map] minuend - matrix/vector/scalar as a map/array/definition string/number
;                  mB       - [Map] subtrahend - matrix/vector/scalar as a map/array/definition string/number (may be overwritten)
;                  bInPlace - [Bool] (Default: False)
;                           ↳ True:  mB gets overwritten
;                             False: mB remains untouched
; Return value ..: Success: bInplace ? True : result as matrix/vector/scalar as a map/array/definition string/number
;                  Failure: $mRet and set @error to:
;                           | 1: invalid value for mA/mB (@extended = 1: mA, @extended = 2: mB)
;                           | 1X: error X during _la_add() (@extended: @extended from _la_add()
;                           | 2X: error X during _blas_scal() (@extended: @extended from _blas_scal()
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......: only changes the sign of mB and then executes _la_add()
; Related .......: _la_add(), _blas_scal()
; Link ..........:
; Example .......: Yes
;                  ; Matrix - Vector
;                  Global $mResult = _la_sub('[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97],[93,89,13,-53,43]]', '[11,12,13,14,15]')
;                  _la_display($mResult, "matrix - vector")
;                  ; Matrix - Scalar
;                  Global $mResult = _la_sub('[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97],[93,89,13,-53,43]]', 5)
;                  _la_display($mResult, "matrix - scalar")
; ===============================================================================================================================
Func _la_sub(ByRef $mA, ByRef $mB, $bInPlace = False)
	; direct AutoIt-type input
	If IsArray($mA) Or IsString($mA) Then $mA = _blas_fromArray($mA)
	If IsArray($mB) Or IsString($mB) Then $mB = _blas_fromArray($mB)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not (MapExists($mA, "ptr") Or IsNumber($mA)) Then Return SetError(1, 1, $bInPlace ? False : Null)
	If Not (MapExists($mB, "ptr") Or IsNumber($mB)) Then Return SetError(1, 2, $bInPlace ? False : Null)

	If Not $bInPlace And IsMap($mB) Then $mB = _la_duplicate($mB)

	; invert elements in B
	If IsMap($mB) Then
		_blas_scal($mB, -1)
		If @error Then Return SetError(@error + 20, @extended, $bInPlace ? False : Null)
	Else
		$mB = 0 - $mB
	EndIf

	Local $mRet = _la_add($mA, $mB, $bInPlace)
	Return SetError(@error + 10, @extended, $mRet)
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_add()
; Description ...: calculate the sum of a matrix/vector/scalar mA and a matrix/vector/scalar mB
; Syntax ........: _la_add($mA, $mB, [$bInPlace = False])
; Parameters ....: mA       - [Map] summand mA - matrix/vector/scalar as a map/array/definition string/number
;                  mB       - [Map] summand mB - matrix/vector/scalar as a map/array/definition string/number
;                  bInPlace - [Bool] (Default: False)
;                           ↳ True:  mB gets overwritten
;                             False: mB remains untouched
; Return value ..: Success: bInplace ? True : result as matrix/vector/scalar as a map/array/definition string/number
;                  Failure:  and set @error to:
;                           | 1: invalid value for mA/mB (@extended = 1: mA, @extended = 2: mB)
;                           | 2: invalid types for mA/mB (@extended: reference to specific location in the code)
;                           | 3: invalid dimensions between mA and mB (@extended: reference to specific location in the code)
;                           |1X: error X during _blas_ger() (@extended: @extended from _blas_ger())
;                           |2X: error X during __la_addScalar() (@extended: @extended from __la_addScalar())
;                           |3X: error X during _blas_ger() (@extended: @extended from _blas_ger())
;                           |4X: error X during _blas_axpy() (@extended: @extended from _blas_axpy())
;                           |5X: error X during __la_addScalar() (@extended: @extended from __la_addScalar())
;                           |6X: error X during __la_addScalar() (@extended: @extended from __la_addScalar())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _blas_ger(), _blas_axpy(), __la_addScalar()
; Link ..........:
; Example .......: Yes
;                  ; Vector + Vector
;                  Global $mResult = _la_add('[11,12,13,14,15]', '[11,12,13,14,15]')
;                  _la_display($mResult, "vector + vector")
;                  ; Matrix + Vector
;                  Global $mResult = _la_add('[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97],[93,89,13,-53,43]]', '[11,12,13,14,15]')
;                  _la_display($mResult, "matrix + vector")
;                  ; Vector + Matrix
;                  Global $mResult = _la_add('[11,12,13,14,15]', '[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97],[93,89,13,-53,43]]')
;                  _la_display($mResult, "vector + matrix")
;                  ; Matrix + Matrix
;                  Global $mA = '[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97],[93,89,13,-53,43]]'
;                  Global $mResult = _la_add($mA, $mA)
;                  _la_display($mResult, "Matrix + Matrix")
;                  ; Matrix + Scalar
;                  Global $mA = '[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97],[93,89,13,-53,43]]'
;                  Global $mResult = _la_add($mA, 5.0)
;                  _la_display($mResult, "Matrix + scalar")
; ===============================================================================================================================
Func _la_add(ByRef $mA, ByRef $mB, $bInPlace = False)
	; direct AutoIt-type input
	If IsArray($mA) Or IsString($mA) Then $mA = _blas_fromArray($mA)
	If IsArray($mB) Or IsString($mB) Then $mB = _blas_fromArray($mB)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not ((IsMap($mA) And MapExists($mA, "ptr")) Or IsNumber($mA)) Then Return SetError(1, 1, $bInPlace ? False : Null)
	If Not ((IsMap($mB) And MapExists($mB, "ptr")) Or IsNumber($mB)) Then Return SetError(1, 2, $bInPlace ? False : Null)

	Local $sTypeA = IsMap($mA) ? ($mA.storageType = 0 ? "V" : "M") : (IsNumber($mA) ? "S" : "")
	Local $sTypeB = IsMap($mB) ? ($mB.storageType = 0 ? "V" : "M") : (IsNumber($mB) ? "S" : "")

	; variable declarations
	Local $mY

	Switch $sTypeA
		Case "M"
			Switch $sTypeB
				Case "M"
					If $mA.rows <> $mB.rows Then Return SetError(3, 1, Null)
					If $mA.cols <> $mB.cols Then Return SetError(3, 2, Null)

					If Not $bInPlace Then $mB = _la_duplicate($mB)

					; misuse axpy because it`s an element wise addition - independent of the shape
					_blas_axpy($mA, $mB)
					Return SetError(@error + (@error ? 10 : 0), @extended, @error ? False : ($bInPlace ? True : $mB))

				Case "V"
					If $mA.rows <> $mB.elements Then Return SetError(3, 3, Null)

					If Not $bInPlace Then $mA = _la_duplicate($mA)

					; identity vector for missusing _dger_ for matrix-vector addition (incY = 0 is not possible in _blas_ger)
					$mY = _la_createIdentity($mA.cols)

					_blas_ger($mA, $mB, $mY, 1.0, $mA.rows, $mA.cols, 1, 1, $mA.rows, $mA.datatype)
					Return @error ? SetError(@error + 10, @extended, $bInPlace ? False : Null) : ($bInPlace ? True : $mA)

				Case "S"
					If Not $bInPlace Then $mA = _la_duplicate($mA)
					__la_addScalar($mA, $mB)
					Return SetError(@error + (@error ? 20 : 0), @extended, @error ? False : ($bInPlace ? True : $mA))

				Case Else
					Return SetError(2, 2, Null)

			EndSwitch

		Case "V"
			Switch $sTypeB
				Case "M"
					If $mB.rows <> $mA.elements Then Return SetError(3, 4, Null)

					; identity vector for missusing _dger_ for matrix-vector addition
					$mY = _la_createIdentity($mB.cols)

					_blas_ger($mB, $mY, $mA)
					Return @error ? SetError(@error + 30, @extended, $bInPlace ? False : Null) : ($bInPlace ? True : $mB)

				Case "V"
					If $mA.elements <> $mB.elements Then Return SetError(3, 5, Null)

					If Not $bInPlace Then $mB = _la_duplicate($mB)
					_blas_axpy($mA, $mB)
					If @error Then Return SetError(@error + 40, @extended, $bInPlace ? False : Null)
					Return $bInPlace ? True : $mB

				Case "S"
					If Not $bInPlace Then $mA = _la_duplicate($mA)
					__la_addScalar($mA, $mB)
					Return SetError(@error + (@error ? 50 : 0), @extended, @error ? False : ($bInPlace ? True : $mA))

				Case Else
					Return SetError(2, 3, Null)

			EndSwitch

		Case "S"
			Switch $sTypeB
				Case "M", "V"
					If Not $bInPlace Then $mB = _la_duplicate($mB)
					__la_addScalar($mB, $mA)
					Return SetError(@error + (@error ? 60 : 0), @extended, @error ? False : ($bInPlace ? True : $mB))

				Case "S"
					Return $mA + $mB

				Case Else
					Return SetError(2, 4, Null)

			EndSwitch

		Case Else
			Return SetError(2, 1, Null)
	EndSwitch

EndFunc   ;==>_la_add

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __la_addScalar()
; Description ...: adds a constant value to the elements of a matrix/vector
; Syntax ........: __la_addScalar($mMatrix, $fScalar, [$iStart = 0, [$iInc = 1, [$iN = Default]]])
; Parameters ....: mMatrix - [Map] matrix/vector as a map/array/definition string
;                  fScalar - [Float] the scalar value
;                  iStart  - [UInt] (Default: 0)
;                          ↳ start element index in mMatrix (0-based)
;                  iInc    - [UInt] (Default: 1)
;                          ↳ storage spacing between elements of mMatrix (can be used to handle parts of a matrix as a vector)
;                  iN      - [UInt] (Default: Default)
;                          ↳ number of elements in mMatrix
; Return value ..: Success: @error ? SetError(@error + 10, @extended, False) : True
;                  Failure: False and set @error to:
;                           | 1: invalid value for mMatrix
;                           | 2: invalid value for fScalar
;                           |1X: error X during _blas_axpy() (@extended: @extended from _blas_axpy())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: Yes
;                  Global $mMatrix = _la_fromArray('[[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]]')
;                  _la_display($mMatrix, "before")
;                  __la_addScalar($mMatrix, 3.0)
;                  _la_display($mMatrix, "after")
; ===============================================================================================================================
Func __la_addScalar(ByRef $mMatrix, $fScalar, $iStart = 0, $iInc = 1, $iN = Default)
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not (IsMap($mMatrix) And MapExists($mMatrix, "ptr")) Then Return SetError(1, 0, False)
	If Not IsNumber($fScalar) Then Return SetError(2, 0, False)

	Local Const $cPrefix = ($mMatrix.datatype = "FLOAT") ? "s" : "d"
	Local Const $dSize = ($mMatrix.datatype = "DOUBLE") ? $iBLAS_SIZE_DOUBLE : $iBLAS_SIZE_FLOAT

	If IsKeyword($iN) = 1 Then $iN = $mMatrix.size

	; struct for the scalar value
	Local $tTmp = DllStructCreate($mMatrix.datatype)
	DllStructSetData($tTmp, 1, $fScalar)

	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "axpy", _
			"INT*", $iN, _ ; number of elements in vectors x and y
			$mMatrix.datatype & "*", 1.0, _ ; the scalar a
			"PTR", DllStructGetPtr($tTmp), _ ; start ptr of x
			"INT*", 0, _ ; increment in x
			"PTR", $mMatrix.ptr + $dSize * $iStart, _ ; start ptr of y
			"INT*", $iInc _ ; increment in y
			)
	Return @error ? SetError(@error + 10, @extended, False) : True
EndFunc   ;==>__la_addScalar

#EndRegion

#Region multiplication
; Vector-Vector:
; - "scalar product"/"inner product"   a · b : _la_dot()                        --> a(n) · b(n) = scalar
; - "element wise"/"Hadarmard"         a ∘ b : _la_mulElementWise() / _la_mul() --> a(n) ∘ b(n) = c(n)
; - "cross product"                    a × b : _la_cross()                      --> a(3) × b(3) = c(3)
; - "outer product"/"tensor product"   a ⊗ b: _la_outerProduct()               --> a(m) ⊗ b(n) = C(m,n)
; Matrix-Vector:
; - "standard multiplication"          A · x : _la_mul()                        --> A(m,n) · x(n) = y(m)
; - "inverse standard multiplication"  x · A : _la_mul()                        --> x(m) · A(m,n) = y(n)
; Matrix-Matrix:
; - "standard multiplication"          A · B : _la_mul()                        --> A(m,n) · B(n,o) = C(m,o)
; - "element wise"/"Hadamard"          A ∘ B : _la_mulElementWise()             --> A(m,n) ∘ B(m,n) = C(m,n)


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_mul()
; Description ...: calculates a multiplication between a matrix/vector/scalar A and a matrix/vector/scalar B
; Syntax ........: _la_mul($mA, $mB, [$bInPlace = False, [$mC = Default]])
; Parameters ....: mA       - [Map] factor A - matrix/vector/scalar as a map/array/definition string/number
;                  mB       - [Map] factor B - matrix/vector/scalar as a map/array/definition string/number
;                  bInPlace - [Bool] (Default: False)
;                           ↳ True:  mA gets overwritten
;                             False: mA remains untouched
;                  mC       - [Map] (Default: Default)
;                           ↳ Existing matrix/vector/scalar into which the result is to be written (else: result returned by function)
; Return value ..: Success: bInplace ? True : result as matrix/vector/scalar as a map/array/definition string/number
;                  Failure: $mRet and set @error to:
;                           | 1: invalid value for mA/mB (@extended = 1: mA, @extended = 2: mB)
;                           | 2: invalid types for mA/mB (@extended: reference to specific location in the code)
;                           | 3: invalid dimensions between mA and mB (@extended: reference to specific location in the code)
;                           |1X: error X during _blas_scal() (@extended: @extended from _blas_scal())
;                           |2X: error X during _blas_gemv() (@extended: @extended from _blas_gemv())
;                           |3X: error X during _blas_gemm() (@extended: @extended from _blas_gemm())
;                           |4X: error X during _blas_gemv() (@extended: @extended from _blas_gemv())
;                           |5X: error X during _blas_scal() (@extended: @extended from _blas_scal())
;                           |6X: error X during _blas_scal() (@extended: @extended from _blas_scal())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......: ToDo: Use optimized function depending on the nature of the input matrices
;                  (e.g. _blas_symm/_trmm instead of _blas_gemm or _blas_trmv/symv/gbmv/sbmv instead of _blas_gemv)
; Related .......:
; Link ..........:
; Example .......: Yes
;                 ; Matrix * Scalar:
;                 Global $mResult = _la_mul("[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97],[93,89,13,-53,43]]", 3.0)
;                 _la_display($mResult, "matrix * scalar")
;                 ; Matrix · vector:
;                 Global $mResult = _la_mul('[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97]]', '[-67,77,-70,13,-58]')
;                 _la_display($mResult, 'A(m,n) · x(n) = y(m)')
;                 ; Vector · Matrix
;                 Global $mResult = _la_mul('[-67,77,-70,13,-58]', '[[19,-80,-55,-58],[29,-92,-67,-94],[-12,-29,77,-7],[96,96,-70,12],[93,89,13,-53]]')
;                 _la_display($mResult, 'x(m) · A(m,n) = y(n)')
;                 ; Matrix · Matrix:
;                 Global $mA = _la_fromArray("[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97]]")
;                 Global $mB = _la_fromArray("[[-63,4,25,-11,17,89],[-25,69,65,34,45,12],[-1,-63,98,46,6,71],[-31,-87,51,9,12,64],[-88,-34,11,50,1,24]]")
;                 Global $mResult = _la_mul($mA, $mB)
;                 _la_display($mResult, 'A(m,n) · B(n,o) = C(m,o)')
; ===============================================================================================================================
Func _la_mul($mA, $mB, $bInPlace = False, $mC = Default)
	; direct AutoIt-type input
	If IsArray($mA) Or IsString($mA) Then $mA = _blas_fromArray($mA)
	If IsArray($mB) Or IsString($mB) Then $mB = _blas_fromArray($mB)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not (MapExists($mA, "ptr") Or IsNumber($mA)) Then Return SetError(1, 1, $bInPlace ? False : Null)
	If Not (MapExists($mB, "ptr") Or IsNumber($mB)) Then Return SetError(1, 2, $bInPlace ? False : Null)

	; process data structure types
	Local $sTypeA = IsMap($mA) ? ($mA.storageType = 0 ? "V" : "M") : (IsNumber($mA) ? "S" : "")
	Local $sTypeB = IsMap($mB) ? ($mB.storageType = 0 ? "V" : "M") : (IsNumber($mB) ? "S" : "")

	; local variables
	Local $mRet

	Switch $sTypeA
		Case "M"
			Switch $sTypeB
				Case "S" ; scale Matrix
					; duplicate input matrix to prevent overwrite if option choosed
					If Not $bInPlace And IsMap($mA) Then $mA = _la_duplicate($mA)

					; TODO: use _lp_lascl() instead: easy multiplication for different matrix shapes

					; run BLAS function scal
					_blas_scal($mA.ptr, $mB, 0, 1, $mA.size, $mA.datatype)
					Return @error ? SetError(@error + 10, @extended, $bInPlace ? False : Null) : ($bInPlace ? True : $mA)

				Case "V" ; standard matrix-vector multiplication / dot product: A(m,n) · x(n) = y(m)
					; check for correct input dimensions
					If $mA.cols <> $mB.size Then Return SetError(3, 1, $bInPlace ? False : Null)

					$mRet = IsKeyword($mC) = 1 ? _blas_createVector($mA.rows, $mB.datatype) : $mC

					_blas_gemv($mA.ptr, $mB.ptr, $mRet.ptr, 1, 0, "N", 1, 1, $mA.rows, $mA.cols, $mA.rows, $mA.datatype)
					Return @error ? SetError(@error + 20, @extended, $bInPlace ? False : Null) : ($bInPlace ? True : $mRet)

				Case "M"
					; TODO: If StringInStr($mA.type, "S"): use _blas_symm() instead
					; TODO: If StringInStr($mA.type, "T"): use _blas_trmm() instead

					; check for correct input dimensions
					If $mA.cols <> $mB.rows Then Return SetError(3, 2, $bInPlace ? False : Null)

					$mRet = IsKeyword($mC) = 1 ? _blas_createMatrix($mA.rows, $mB.cols, $mA.datatype) : $mC

					_blas_gemm($mA.ptr, $mB.ptr, $mRet.ptr, 1, 0, "N", "N", $mA.rows, $mB.cols, $mA.cols, $mA.rows, $mB.rows, $mRet.rows, $mA.datatype)
					Return @error ? SetError(@error + 30, @extended, $bInPlace ? False : Null) : ($bInPlace ? True : $mRet)

				Case Else
					SetError(2, 2, $bInPlace ? False : Null)

			EndSwitch

		Case "V"
			Switch $sTypeB
				Case "M" ; standard matrix-vector multiplication / dot product: x(m) · A(m,n) = y(n)
					; TODO: If StringInStr($mA.type, "S"): use _blas_symm() instead
					; TODO: If StringInStr($mA.type, "T"): use _blas_trmm() instead

					; check for correct input dimensions
					If $mA.size <> $mB.rows Then Return SetError(3, 3, $bInPlace ? False : Null)

					; target vector y
					$mRet = IsKeyword($mC) = 1 ? _blas_createVector($mB.cols, $mA.datatype) : $mC

					; use transposed version of gemv to calculate x*A instead of A*x
					_blas_gemv($mB.ptr, $mA.ptr, $mRet.ptr, 1, 0, "T", 1, 1, $mB.rows, $mB.cols, $mB.rows, $mB.datatype)
					Return @error ? SetError(@error + 40, @extended, $bInPlace ? False : Null) : ($bInPlace ? True : $mRet)

				Case "V"
					$mRet = _la_mulElementWise($mA, $mB, $bInPlace)
					Return SetError(@error, @extended, $mRet)

				Case "S"
					; duplicate input matrix to prevent overwrite if option choosed
					If Not $bInPlace And IsMap($mA) Then $mA = _la_duplicate($mA)

					; run BLAS function scal
					_blas_scal($mA.ptr, $mB, 0, 1, $mA.size, $mA.datatype)
					Return @error ? SetError(@error + 50, @extended, $bInPlace ? False : Null) : ($bInPlace ? True : $mA)

				Case Else
					SetError(2, 2, $bInPlace ? False : Null)

			EndSwitch

		Case "S"
			Switch $sTypeB
				Case "M", "V"
					; duplicate input matrix to prevent overwrite if option choosed
					If Not $bInPlace And IsMap($mB) Then $mB = _la_duplicate($mB)

					; run BLAS function scal
					_blas_scal($mB.ptr, $mA, 0, 1, $mB.size, $mB.datatype)
					Return @error ? SetError(@error + 60, @extended, $bInPlace ? False : Null) : ($bInPlace ? True : $mB)

				Case "S"
					Return $mA * $mB

				Case Else
					SetError(2, 2, $bInPlace ? False : Null)

			EndSwitch

		Case Else
			Return SetError(2, 1, $bInPlace ? False : Null)

	EndSwitch
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_outerproduct()
; Description ...: calculates the outer product ("tensor product") of two vectors
; Syntax ........: _la_outerproduct($mVecA, $mVecB)
; Parameters ....: mVecA - [Map] vector A as a map/array/definition string
;                  mVecB - [Map] vector B as a map/array/definition string
; Return value ..: Success: [Map] result matrix
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mVecA/mVecB (@extended = 1: mVecA, @extended = 2: mVecB)
;                           |1X: error X during _blas_ger() (@extended: @extended from _blas_ger())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _blas_ger()
; Link ..........:
; Example .......: Yes
;                  Global $mOuter = _la_outerproduct("[1,2,3]","[4,5,6]")
;                  _la_display($mOuter, "outer product")
; ===============================================================================================================================
Func _la_outerproduct($mVecA, $mVecB)
	; direct AutoIt-type input
	If IsArray($mVecA) Or IsString($mVecA) Then $mVecA = _blas_fromArray($mVecA)
	If IsArray($mVecB) Or IsString($mVecB) Then $mVecB = _blas_fromArray($mVecB)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not MapExists($mVecA, "ptr") Then Return SetError(1, 1, Null)
	If Not MapExists($mVecB, "ptr") Then Return SetError(1, 2, Null)

	Local $mResult = _blas_createMatrix($mVecA.elements, $mVecB.elements, $mVecA.datatype)

	_blas_ger($mResult, $mVecA, $mVecB, 1.0)
	Return @error ?  SetError(@error + 10, @extended, Null) : $mResult
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_dot()
; Description ...: calculate the "dot product"/"scalar product"/"inner product" of two vectors
; Syntax ........: _la_dot($mVecA, $mVecB)
; Parameters ....: mVecA - [Map] vector A as a map/array/definition string
;                  mVecB - [Map] vector B as a map/array/definition string
; Return value ..: Success: [Float] result as a scalar value
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mVecA/mVecB (@extended = 1: mVecA, @extended = 2: mVecB)
;                           | 2: invalid dimensions between mVecA and mVecB
;                           |1X: error X during _blas_dot() (@extended: @extended from _blas_dot())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: Yes
;                  Global $fRet = _la_dot("[1,2,3]","[4,5,6]")
;                  MsgBox(0,"dot product", $fRet)
; ===============================================================================================================================
Func _la_dot($mVecA, $mVecB)
	; direct AutoIt-type input
	If IsArray($mVecA) Or IsString($mVecA) Then $mVecA = _blas_fromArray($mVecA)
	If IsArray($mVecB) Or IsString($mVecB) Then $mVecB = _blas_fromArray($mVecB)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not MapExists($mVecA, "ptr") Then Return SetError(1, 1, Null)
	If Not MapExists($mVecB, "ptr") Then Return SetError(1, 2, Null)

	If $mVecA.elements <> $mVecB.elements Then Return SetError(2, 0, Null)

	Local $fRet = _blas_dot($mVecA, $mVecB)
	Return @error ?  SetError(@error + 10, @extended, Null) : $fRet
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_scale()
; Description ...: multiplies the elements of a matrix/vector by a scalar value
; Syntax ........: _la_scale($mMatrix, $fScalar, [$bInPlace = False, [$iStart = 0, [$iInc = 1]]])
; Parameters ....: mMatrix  - [Map] matrix/vector A as a map/array/definition string (may be overwritten)
;                  fScalar  - [Float] the scalar value
;                  bInPlace - [Bool] (Default: False)
;                           ↳ True:  mMatrix gets overwritten
;                             False: mMatrix remains untouched
;                  iStart   - [UInt] (Default: 0)
;                           ↳ start index in mMatrix (0-based)
;                  iInc     - [UInt] (Default: 1)
;                           ↳ storage spacing between elements of mMatrix (can be used to handle parts of a matrix as a vector)
; Return value ..: Success: bInplace ? True : [Map] result matrix/vector
;                  Failure: bInPlace ? False : Null and set @error to:
;                           | 1: invalid value for mMatrix/fScalar (@extended = 1: mMatrix, @extended = 2: fScalar)
;                           |1X: error X during _blas_scal() (@extended: @extended from _blas_scal())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _blas_scal()
; Link ..........:
; Example .......: Yes
;                  ; vector * scalar
;                   Global $mResult = _la_scale("[-67,77,-70,13]", ACos(-1))
;                   _la_display($mResult, "after")
;                  ; matrix * scalar
;                   Global $mResult = _la_scale("[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97],[93,89,13,-53,43]]", 3.0)
;                   _la_display($mResult, "after")
; ===============================================================================================================================
Func _la_scale($mMatrix, $fScalar, $bInPlace = False, $iStart = 0, $iInc = 1)
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not MapExists($mMatrix, "ptr") Then Return SetError(1, 1, $bInPlace ? False : Null)
	If Not (IsNumber($fScalar) Or StringIsFloat($fScalar) Or StringIsInt($fScalar)) Then Return SetError(1, 2, $bInPlace ? False : Null)

	; duplicate input matrix to prevent overwrite if option choosed
	If Not $bInPlace And IsMap($mMatrix) Then $mMatrix = _la_duplicate($mMatrix)

	; run BLAS function scal
	;~ _blas_scal($mMatrix, $fScalar, $iStart, $iInc, $iN)
	_blas_scal($mMatrix.ptr, $fScalar, $iStart, $iInc, $mMatrix.size, $mMatrix.datatype)
	Return @error ? SetError(@error + 10, @extended, $bInPlace ? False : Null) : ($bInPlace ? True : $mMatrix)
EndFunc   ;==>_la_scale

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_mulElementWise()
; Description ...: calculates the element-wise ("Hadarmard") product between two matrices/vectors
; Syntax ........: _la_mulElementWise($mVecA, $mVecB, [$bInPlace = False])
; Parameters ....: mVecA - [Map] matrix/vector A as a map/array/definition string
;                  mVecB - [Map] matrix/vector B as a map/array/definition string (may be overwritten)
;                  bInPlace - [Bool] (Default: False)
;                           ↳ True:  mVecB gets overwritten
;                             False: mVecB remains untouched
; Return value ..: Success: bInplace ? True : [Map] result matrix/vector
;                  Failure: bInPlace ? False : Null and set @error to:
;                           | 1: invalid value for mVecA/mVecB (@extended = 1: mVecA, @extended = 2: mVecB)
;                           | 2: invalid dimensions between mVecA and mVecB
;                           |1X: error X during _blas_sbmv() (@extended: @extended from _blas_sbmv())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _blas_sbmv()
; Link ..........:
; Example .......: Yes
;                  ; Vector * Vector element wise
;                  Global $mResult = _la_mulElementWise('[-67,77,-70,13,-58,-94,-7]', '[-67,77,-70,13,-58,-94,-7]')
;                  _la_display($mResult)
;                  ; Matrix * Matrix element wise product
;                  Global $mResult = _la_mulElementWise("[[1,2],[3,4]]", "[[5,6],[7,8]]")
;                  _la_display($mResult, "after")
; ===============================================================================================================================
Func _la_mulElementWise(ByRef $mVecA, ByRef $mVecB, $bInPlace = False)
	; direct AutoIt-type input
	If IsArray($mVecA) Or IsString($mVecA) Then $mVecA = _blas_fromArray($mVecA)
	If IsArray($mVecB) Or IsString($mVecB) Then $mVecB = _blas_fromArray($mVecB)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not (MapExists($mVecA, "ptr") Or IsNumber($mVecA)) Then Return SetError(1, 1, $bInPlace ? False : Null)
	If Not (MapExists($mVecB, "ptr") Or IsNumber($mVecB)) Then Return SetError(1, 2, $bInPlace ? False : Null)

	; dimension check
	If $mVecB.elements < $mVecA.elements Then Return SetError(2, 1, $bInPlace ? False : Null)

	; empty (0-filled) vector Y
	Local $tTmp = DllStructCreate(StringFormat("%s[%d]", $mVecA.datatype, $mVecA.elements)), $pTmp = DllStructGetPtr($tTmp)

	; sbmv to treat the vector simply as diagonal matrix to perform a efficient element wise multiplication
	_blas_sbmv($mVecA, $mVecB, $pTmp, 1, 1, 0, "L", $mVecA.elements, 1, 1, 1)
	If @error Then Return SetError(@error + 10, @extended, $bInPlace ? False : Null)

	If $bInPlace Then
		$mVecB.struct = $tTmp
		$mVecB.ptr = $pTmp
		Return True
	Else
		Local $mRet = BitAND($mVecA.storageType, $__g_BLAS_STYPE_MATRIX) ? _blas_createMatrix($mVecA.rows, $mVecA.cols, $mVecA.datatype) : _blas_createVector($mVecA.elements)
		$mRet.struct = $tTmp
		$mRet.ptr = $pTmp
		Return $mRet
	EndIf
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_cross()
; Description ...: calculates the cross product between two 3-element vectors
; Syntax ........: _la_cross($mVecA, $mVecB)
; Parameters ....: mVecA - [Map] vector A as a map/array/definition string
;                  mVecB - [Map] vector B as a map/array/definition string
; Return value ..: Success: [Map] result 3-element vector
;                  Failure: False and set @error to:
;                           | 1: invalid value for mVecA/mVecB (@extended = 1: mVecA, @extended = 2: mVecB)
;                           | 2: invalid storage type for mVecA/mVecB (@extended = 1: mVecA, @extended = 2: mVecB)
;                           | 3: invalid dimension for mVecA ($mVecA.elements must be 3)
;                           | 4: invalid dimension for mVecB ($mVecB.elements must be 3)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: Yes
;                  Global $mCross = _la_cross("[1, 2, 0]", "[4, 5, 6]")
;                  _la_display($mCross, "cross product")
; ===============================================================================================================================
Func _la_cross($mVecA, $mVecB)
	; direct AutoIt-type input
	If IsArray($mVecA) Or IsString($mVecA) Then $mVecA = _blas_fromArray($mVecA)
	If IsArray($mVecB) Or IsString($mVecB) Then $mVecB = _blas_fromArray($mVecB)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not MapExists($mVecA, "ptr") Then Return SetError(1, 1, False)
	If Not MapExists($mVecB, "ptr") Then Return SetError(1, 2, False)

	If $mVecA.storageType <> 0 Then Return SetError(2, 1, Null)
	If $mVecB.storageType <> 0 Then Return SetError(2, 2, Null)

	If $mVecA.elements <> 3 Then Return SetError(3, $mVecA.elements, Null)
	If $mVecB.elements <> 3 Then Return SetError(4, $mVecB.elements, Null)

	Local $mRet = _la_createVector(3, $mVecA.datatype)
	Local $tA = $mVecA.struct, $tB = $mVecB.struct, $tR = $mRet.struct

	DllStructSetData($tR, 1, DllStructGetData($tA, 1, 2) * DllStructGetData($tb, 1, 3) - DllStructGetData($tA, 1, 3) * DllStructGetData($tb, 1, 2), 1)
	DllStructSetData($tR, 1, DllStructGetData($tA, 1, 3) * DllStructGetData($tb, 1, 1) - DllStructGetData($tA, 1, 1) * DllStructGetData($tb, 1, 3), 2)
	DllStructSetData($tR, 1, DllStructGetData($tA, 1, 1) * DllStructGetData($tb, 1, 2) - DllStructGetData($tA, 1, 2) * DllStructGetData($tb, 1, 1), 3)

	Return $mRet
EndFunc

#EndRegion

#Region factorization / decomposition

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_LU()
; Description ...: calculates the LU decomposition of a matrix
; Syntax ........: _la_LU($mA, [$bInPlace = False])
; Parameters ....: mA       - [Map] matrix as a map/array/definition string (may be overwritten)
;                  bInPlace - [Bool] (Default: False)
;                           ↳ True:  mA gets overwritten
;                             False: mA remains untouched
; Return value ..: Success: bInPlace ? True : [Map] map with 3 elements: {"L": lower triangular matrix L, "U": upper triangular matrix U, "P": permutation matrix P}
;                  Failure: bInPlace ? False : Null and set @error to:
;                           | 1: invalid value for mA
;                           |1X: error X during _lp_getrf() (@extended: _lp_getrf())
;                           |2X: error X during _lp_laswp() (@extended: _lp_laswp())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _lp_getrf(), _lp_laswp()
; Link ..........:
; Example .......: Yes
;                  Global $mA = _blas_fromArray('[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97]]')
;                  Global $mLU = _la_LU($mA)
;                  _la_display($mLU.L, "lower triangular matrix L")
;                  _la_display($mLU.U, "upper triangular matrix U")
;                  _la_display($mLU.P, "permutation matrix P")
; ===============================================================================================================================
Func _la_LU($mA, $bInPlace = False)
	; direct AutoIt-type input
	If IsArray($mA) Or IsString($mA) Then $mA = _blas_fromArray($mA)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not MapExists($mA, "ptr") Then Return SetError(1, 0, $bInPlace ? False : Null)

	; data structure for return matrices list
	Local $mRet[]

	; prevent overwriting if user not choosed in-place processing
	If Not $bInPlace Then $mA = _la_duplicate($mA)

	; run the lapack function
	Local $tIPIV = _lp_getrf($mA)
	If @error Then Return SetError(@error + 10, @extended, $bInPlace ? False : Null)

	; extract the U and L triangular matrices
	$mRet.U = _la_getTriangle($mA, "U", True, True)
	$mRet.L = _la_getTriangle($mA, "L", False, True)

	; reconstruct the permutation matrix P
	Local $mP = _la_createIdentity($mA.rows, $mA.cols)
	_lp_laswp($mP, $tIPIV, 1, $mA.rows, -1) ; remark: "-1" to obtain A = P*L*U, "+1" if P*A = L*U
	If @error Then Return SetError(@error + 20, @extended, $bInPlace ? False : Null)
	$mRet.P = $mP

	Return $mRet
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_QR()
; Description ...: calculates the QR decomposition of a matrix
; Syntax ........: _la_QR($mMatrix, [$bTileAlgorithm = False, [$bInPlace = False]])
; Parameters ....: mMatrix        - [Map] matrix as a map/array/definition string (may be overwritten)
;                  bTileAlgorithm - [Bool] (Default: False)
;                                 ↳ True: tiled QR algorithm is used
;                                   False: conventional / "Householder" QR algorithm is used
;                  bInPlace       - [Bool] (Default: False)
;                                 ↳ True:  mMatrix gets overwritten
;                                   False: mMatrix remains untouched
; Return value ..: Success: bInPlace ? True : [Map] map with 2 elements: {"Q": orthogonal matrix Q, "R": upper triangular matrix R}
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mMatrix
;                           |1X: error X during _lp_geqr() (@extended: _lp_geqr())
;                           |2X: error X during _lp_gemqr() (@extended: _lp_gemqr())
;                           |3X: error X during _lp_geqrf() (@extended: _lp_geqrf())
;                           |4X: error X during _lp_orgqr() (@extended: _lp_orgqr())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _lp_geqr(), _lp_gemqr(), _lp_geqrf(), _lp_orgqr()
; Link ..........:
; Example .......: Yes
;                  Global $mQR = _la_QR('[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97],[93,89,13,-53,43]]')
;                  _la_display($mQR.Q, "Q")
;                  _la_display($mQR.R, "R")
; ===============================================================================================================================
Func _la_QR($mMatrix, Const $bTileAlgorithm = False, Const $bInPlace = False)
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not MapExists($mMatrix, "ptr") Then Return SetError(1, 0, $bInPlace ? False : Null)

	; data structure for return matrices list
	Local $mRet[]

	; prevent overwriting if user not choosed in-place processing
	If Not $bInPlace Then $mMatrix = _la_duplicate($mMatrix)

	; choose QR algorithm
	If $bTileAlgorithm Then ; tiled QR algorithm

		; calculate QR decomposition
		Local $tT = _lp_geqr($mMatrix)
		If @error Then Return SetError(@error + 10, @extended, Null)
		Local $iTSIZE = @extended

		; extract R (upper triangle of A)
		$mRet.R = _la_getTriangle($mMatrix, "U", True, True)

		; initialize Q (= C) as identity matrix to calculate whole Q
		Local $mQ = _blas_createMatrix($mMatrix.rows, $mMatrix.cols, $mMatrix.datatype)
		__blas_fillWithScalar($mQ, 1, 0, $mMatrix.rows + 1)

		; calculate Q out of T
		_lp_gemqr($mMatrix, $tT, $mQ, $iTSIZE)
		If @error Then Return SetError(@error + 20, @extended, Null)
		$mRet.Q = $mQ

	Else ; conventional / "Householder" QR algorithm

		; calculate QR decomposition
		Local $tTau = _lp_geqrf($mMatrix)
		If @error Then Return SetError(@error + 30, @extended, Null)

		; extract R (upper triangle of A)
		$mRet.R = _la_getTriangle($mMatrix, "U", True)

		; calculate Q out of the householder reflectors
		_lp_orgqr($mMatrix, $tTau)
		If @error Then Return SetError(@error + 40, @extended, Null)
		$mRet.Q = $mMatrix

	EndIf

	Return $mRet
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_SVD()
; Description ...: calculates the singular value decomposition (SVD) of a matrix
; Syntax ........: _la_SVD($mMatrix, [$bInPlace = False])
; Parameters ....: mMatrix  - [Map] matrix as a map/array/definition string (may be overwritten)
;                  bInPlace - [Bool] (Default: False)
;                           ↳ True:  mMatrix gets overwritten
;                             False: mMatrix remains untouched
; Return value ..: Success: [Map] map with 3 elements: {"S": singular value vector diag(Σ), "VT": matrix Vᵀ, "U": matrix U}
;                  Failure: bInPlace ? False : Null and set @error to:
;                           | 1: invalid value for mMatrix
;                           |1X: error X during _lp_gesvd() (@extended: @error from _lp_gesvd())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _lp_gesvd()
; Link ..........:
; Example .......: Yes
;                  Global $mSVD = _la_SVD('[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97],[93,89,13,-53,43]]')
;                  _la_display($mSVD.S, "S")
;                  _la_display($mSVD.VT, "Vᵀ")
;                  _la_display($mSVD.U, "U")
; ===============================================================================================================================
Func _la_SVD($mMatrix, Const $bInPlace = False)
	; direct AutoIt-type input
	If IsArray($mMatrix) Or IsString($mMatrix) Then $mMatrix = _blas_fromArray($mMatrix)

	; prevent overwriting if user not choosed in-place processing
	If Not $bInPlace Then $mMatrix = _la_duplicate($mMatrix)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not MapExists($mMatrix, "ptr") Then Return SetError(1, 0, $bInPlace ? False : Null)

	Local $mRet = _lp_gesvd($mMatrix)
	If @error Then Return SetError(@error + 10, @extended, $bInPlace ? False : Null)
	Return $mRet
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_cholesky()
; Description ...: calculate the cholesky decomposition of a symmetric, positive definite matrix ( A --> L * Lᵀ or A --> U * Uᵀ )
; Syntax ........: _la_cholesky($mA, [$cUPLO = "L", [$bInPlace = False]])
; Parameters ....: mA       - [Map] symmetric, positive definite matrix as a map/array/definition string (may be overwritten)
;                  cUPLO    - [Char] (Default: "L") which part holds the value in the symmetric matrix
;                           ↳ "U": upper triangle of A is stored
;                             "L": lower triangle of A is stored
;                  bInPlace - [Bool] (Default: False)
;                           ↳ True:  mA gets overwritten
;                             False: mA remains untouched
; Return value ..: Success: bInPlace ? True : [Map] L-Matrix
;                  Failure: bInPlace ? False : Null and set @error to:
;                           | 1: invalid value for mMatrix
;                           |1X: error X during _lp_potrf() (@extended: @error from _lp_potrf())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _lp_potrf()
; Link ..........:
; Example .......: Yes
;                  Global $mL = _la_cholesky("[[6,2,1],[2,5,2],[1,2,4]]")
;                  _la_display($mL, "Cholesky: L-Matrix")
; ===============================================================================================================================
Func _la_cholesky(ByRef $mA, $cUPLO = "L", $bInPlace = False)
	; direct AutoIt-type input
	If IsArray($mA) Or IsString($mA) Then $mA = _blas_fromArray($mA)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not MapExists($mA, "ptr") Then Return SetError(1, 0, $bInPlace ? False : Null)

	; prevent overwrite if choosed
	If Not $bInPlace Then $mA = _la_duplicate($mA)

	; calculate the cholesky factorization
	_lp_potrf($mA, $cUPLO)
	If @error Then Return SetError(@error + 10, @extended, $bInPlace ? False : Null)

	; extract the lower triangle matrix to get L
	$mA = _la_getTriangle($mA, $cUPLO, True, True)

	Return $bInPlace ? True : $mA
EndFunc

#EndRegion

#Region eigenvalues / eigenvectors


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_eigen()
; Description ...: computes for an N-by-N real matrix A, the eigenvalues and the left and/or right eigenvectors.
; Syntax ........: _la_eigen($mA, [$bSymmetric, [$bComputeEigenvectors, [$bInPlace, [$cUPLO]]]])
; Parameters ....: mA                   - [Map] matrix as a map/array/definition string (may be overwritten)
;                  bSymmetric           - [Bool] (Default: _la_isSymmetric($mA)
;                                       ↳ True: mA is symmetric
;                                         False: mA is not symmetric
;                  bComputeEigenvectors - [Bool] (Default: False)
;                                       ↳ True: also calculate the eigenvectors
;                                         False: don`t calculate the eigenvectors
;                  bInPlace             - [Bool] (Default: False)
;                                       ↳ True:  mA gets overwritten
;                                         False: mA remains untouched
;                  cUPLO                - [Char] (Default: "U") if bSymmetric: which part holds the value in the symmetric matrix
;                                       ↳ "U": upper triangle of A is stored
;                                         "L": lower triangle of A is stored
; Return value ..: Success: [Map] map with 2 elements: {"R": eigenvalues, "VL": array of eigenvectors}
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mA
;                           |1X: error X during _lp_geev() (@extended: @error from _lp_geev())
;                           |2X: error X during _lp_syev() (@extended: @error from _lp_syev())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _lp_geev(), _lp_syev()
; Link ..........:
; Example .......: Yes
;                  Global $mEigen = _la_eigen('[[611,196,-192,407,-8,-52,-49,29],[196,899,113,-192,-71,-43,-8,-44],[-192,113,899,196,61,49,8,52],[407,-192,196,611,8,44,59,-23],[-8,-71,61,8,411,-599,208,208],[-52,-43,49,44,-599,411,208,208],[-49,-8,8,59,208,208,99,-911],[29,-44,52,-23,208,208,-911,99]]', False, True)
;                  _la_display($mEigen.R, "eigenvalues")
;                  ; display all eigen vectors
;                  For $i = 0 To UBound($mEigen.VL) - 1
;                     _la_display($mEigen.VL[$i], "vector " & $i + 1)
;                  Next
; ===============================================================================================================================
Func _la_eigen($mA, $bSymmetric = _la_isSymmetric($mA), $bComputeEigenvectors = False, $bInPlace = False, Const $cUPLO = "U")
	; direct AutoIt-type input
	If IsArray($mA) Or IsString($mA) Then $mA = _blas_fromArray($mA)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not MapExists($mA, "ptr") Then Return SetError(1, 0, $bInPlace ? False : Null)

	; duplicate input matrix to prevent overwrite if option choosed
	If Not $bInPlace Then $mA = _la_duplicate($mA)

	; calculate the eigen values
	Local $mEigen[], $iStart = 0

	If $bSymmetric Then
		$mEigen.R = _lp_syev($mA, $bComputeEigenvectors ? "V" : "N", $cUPLO)
		If @error Then Return SetError(@error + 20, @extended, Null)

		; extract the eigenvectors
		If $bComputeEigenvectors Then

			; put vectors in array of vectors:
			Local $aEigenVectors[$mA.cols], $aEigenVectorsRight[$mA.cols]

			For $i = 0 To $mA.cols - 1
				$aEigenVectors[$i] = _blas_copy($mA, $iStart, 1, 0, 1, $mA.rows)
				$iStart += $mA.rows
			Next

			$mEigen.V = $aEigenVectors

		EndIf
	Else
		$mEigen = _lp_geev($mA, $bComputeEigenvectors ? "V" : "N", $bComputeEigenvectors ? "V" : "N")
		If @error Then Return SetError(@error + 10, @extended, Null)

		; extract the eigenvectors
		If $bComputeEigenvectors Then

			; put vectors in array of vectors:
			Local $mL = $mEigen.VL, $mR = $mEigen.VR
			Local $aEigenVectorsLeft[$mA.cols], $aEigenVectorsRight[$mA.cols]

			For $i = 0 To $mA.cols - 1
				$aEigenVectorsLeft[$i] = _blas_copy($mL, $iStart, 1, 0, 1, $mA.rows)
				$aEigenVectorsRight[$i] = _blas_copy($mR, $iStart, 1, 0, 1, $mA.rows)
				$iStart += $mA.rows
			Next

			$mEigen.VL = $aEigenVectorsLeft
			$mEigen.VR = $aEigenVectorsRight
		EndIf
	EndIf

	Return SetExtended($bSymmetric, $mEigen)
EndFunc

#EndRegion

#Region solve linear equation systems


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_solve()
; Description ...: computes the solution to a system of linear equations A * X = B
; Syntax ........: _la_solve($mA, $mB, [$bIsPositive = _la_isPositiveDefinite($mA])
; Parameters ....: mA          - [Map] matrix A as a map, DllStruct or pointer (may be overwritten)
;                  mB          - [Map] vector/matrix B as a map, DllStruct or pointer (may be overwritten)
;                  bIsPositive - [Bool] (Default: _la_isPositiveDefinite($mA)
;                              ↳ True: mA is positive definite (more efficient algorithm is used)
;                              ↳ False: mA is not positive definite
;                  bInPlace    - [Bool] (Default: False)
;                              ↳ True:  mA/mB gets overwritten
;                                False: mA/mB remains untouched
; Return value ..: Success: [Map] solution vector/matrix X
;                  Failure: bInplace ? False : Null and set @error to:
;                           | 1: invalid value for mA
;                           | 2: mA is not quadratic
;                           |1X: error X during _lp_posv()/_lp_gesv()  (@extended: @extended from _lp_posv()/_lp_gesv())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......: ToDo: Use more efficient functions depending on the shape of the matrix A
;                  (e.g. _blas_tbsv/trsv/trsm or _lp_sysv/gbsv/trtrs/posv)
; Related .......: _lp_gesv(), _lp_posv()
; Link ..........:
; Example .......: Yes
;                  Global $mA = _la_fromArray("[[6,2,1],[2,5,2],[1,2,4]]")
;                  Global $mB = _la_fromArray("[[9,3],[5,2],[5,1]]") ; calculate for 2 different B`s: [9,5,5] and [3,2,1]
;                  Global $mX = _la_solve($mA, $mB)
;                  _la_display($mX, "solution vector X")
; ===============================================================================================================================
Func _la_solve($mA, $mB, $bIsPositive = _la_isPositiveDefinite($mA), $bInPlace = False)
	; direct AutoIt-type input
	If IsArray($mA) Or IsString($mA) Then $mA = _blas_fromArray($mA)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not MapExists($mA, "ptr") Then Return SetError(1, 0, $bInPlace ? False : Null)

	; dimension check
	If $mA.rows <> $mA.cols Then Return SetError(2, 0, $bInPlace ? False : Null)

	; duplicate input matrix to prevent overwrite if option choosed
	If Not $bInPlace Then
		$mA = _la_duplicate($mA)
		$mB = _la_duplicate($mB)
	EndIf

	; check if only one solution vector or multiple used
	Local $iNRHS = $mB.storageType = 0 ? 1 : $mB.cols

	If $bIsPositive Then
		_lp_posv($mA, $mB, $iNRHS)
	Else
		_lp_gesv($mA, $mB, $iNRHS)
	EndIf
	If @error Then Return SetError(@error + 10, @extended, Null)

	Return SetExtended($iNRHS, $mB)
EndFunc


#EndRegion

#Region least squares solving

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_lstsq()
; Description ...: solves overdetermined or underdetermined [weighted] linear system
; Syntax ........: _la_lstsq($mA, $mB, [$mP = Default, [$sAlgorithm = "QR", [$iFlagsResults = 0, [$fTolerance = Default, [$fTallSkinnyThreshhold = 10]]]]])
; Parameters ....: mA                    - [Map] jacobian matrix A (M × N) as a map, DllStruct or pointer
;                  mB                    - [Map] observation vector/matrix B as a map, DllStruct or pointer
;                  mP                    - [Map] (Default: Default)
;                                        ↳ weight matrix or weight vector P as a map, DllStruct or pointer (Default = unweighted)
;                  sAlgorithm            - [String] (Default: "QR") algorithm for solving the system
;                                        ↳ "QR": QR decomposition                    - good stability and performance
;                                          "SVD": singular value decomposition (SVD) - highest stability but lowest performance
;                                          "Cholesky": cholesky decomposition        - better stability and performance than QR but only for positive-definite matrices
;                  iFlagsResults         - [UInt] (Default: 0) components of the output as a BitOr combination of the following flags:
;                                        ↳ $__LA_LSTSQ_R		  - $mRet.r     = Residuals r:  r = y - y_d
;                                          $__LA_LSTSQ_R2Sum      - $mRet.r2sum = square sum of [weighted] residuals: r2sum = rᵀ * W * v
;                                          $__LA_LSTSQ_S0         - $mRet.s0    = a posteriori standard deviation factor
;                                          $__LA_LSTSQ_QX         - $mRet.Qx    = cofactor matrix of parameters
;                                          $__LA_LSTSQ_SDX        - $mRet.sdx   = standard deviations of parameters: sdₓ = sqrt(s₀² * diag(Qₓ))
;                                          $__LA_LSTSQ_QY         - $mRet.Qy0   = a-priori cofactor matrix for the observations: Q_y = P⁻¹
;                                          $__LA_LSTSQ_QYD        - $mRet.Qy    = a-posteriori cofactormatrix for the adjusted observations: Q_yd = A * Qₓ * Aᵀ
;                                          $__LA_LSTSQ_SDY        - $mRet.sdY   = a posteriori standard deviations for the adjusted observations
;                                          $__LA_LSTSQ_QR         - $mRet.Qr    = cofactor matrix for the residuals r: Q_r = Q_y - Q_yd
;                                          $__LA_LSTSQ_SDR        - $mRet.sdR   = standard deviations for the residuals: sd_r = sqrt(s₀² * diag(Q_r))
;                                          $__LA_LSTSQ_REDUNDANCY - $mRet.R     = redundandy matrix R = Q_y * W
;                                          $__LA_LSTSQ_COND       - $mRet.cond  = condition number (only if "SVD" is used)
;                                          $__LA_LSTSQ_RANK       - $mRet.rank  = rank for the jacobian matrix (only if "SVD" is used)
;                  fTolerance            - [Float] (Default: Default)
;                                        ↳ threshold value up to which singular values are regarded as zero (-1 = machine precision)
;                  fTallSkinnyThreshhold - [Float] (Default: 10)
;                                        ↳ threshold value of the ratio rows / columns from which the QR decomposition switches to the tall-skinny algorithm
; Return value ..: Success: [Map] results of the solution depending on iFlagsResults:
;                           { "x":     solution vector x,
;                             "f":     degrees of freedom,
;                             "r":     Residuals r:  r = y - y_d,
;                             "r2sum": square sum of [weighted] residuals: r2sum = rᵀ * W * v,
;                             "s0":    a posteriori standard deviation factor,
;                             "Qx":    cofactor matrix of parameters,
;                             "Qy0":   a-priori cofactor matrix for the observations: Qᵧ = P⁻¹,
;                             "Qy":    a-posteriori cofactormatrix for the adjusted observations: Qᵧd = A * Qₓ * Aᵀ,
;                             "Qr":    cofactor matrix for the residuals r: Qᵣ = Qᵧ - Qᵧd,
;                             "R":     redundandy matrix R = Qᵧ * W,
;                             "sdX":   standard deviations of parameters: sdₓ = sqrt(s₀² * diag(Qₓ)),
;                             "sdY":   a posteriori standard deviations for the adjusted observations,
;                             "sdR":   standard deviations for the residuals: sdᵣ = sqrt(s₀² * diag(Qᵣ)),
;                             "cond":  condition number (only if "SVD" is used),
;                             "rank":  rank for the jacobian matrix (only if "SVD" is used)
;                           }
;                  Failure: Null and set @error to:
;                           |1X: error X during __la_lstsq_svd/__la_lstsq_cholesky/__la_lstsq_qr (@extended: @extended from __la_lstsq_svd/__la_lstsq_cholesky/__la_lstsq_qr)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: __la_lstsq_svd(), __la_lstsq_cholesky(), __la_lstsq_qr
; Link ..........:
; Example .......: Yes
;                  Global $iFlags = 2 * $__LA_LSTSQ_RANK - 1 ; = all
;                  Global $mA = _la_fromArray("[[1,0,0,0,0],[0,1,0,0,0],[-1,0,1,0,0],[-1,0,0,1,0],[-1,0,0,0,1],[0,-1,1,0,0],[0,-1,0,1,0],[0,-1,0,0,1],[0,0,1,-1,0],[0,0,0,1,-1]]")
;                  Global $mb = _la_fromArray("[0,0,0.055,0.001,0.057,0.014,-0.002,0.051,0.015,0.052]")
;                  Global $mP = _la_fromArray("[10,10,0.7,1.0,0.7,0.7,1,0.7,1,1]")
;                  _la_VectorToDiag($mP, True)
;                  Global $mLstSq = _la_lstsq($mA, $mB, $mP, "Cholesky", $iFlags)
;                  ConsoleWrite("s0: " & $mLstSq.s0 & @CRLF)
;                  ConsoleWrite("rᵀPr: " & $mLstSq.r2sum & @CRLF)
;                  _la_display($mLstSq.x, "x")
;                  _la_display($mLstSq.r, "residuals")
;                  _la_display($mLstSq.Qx, "Qx")
;                  _la_display($mLstSq.Qy0, "Ql0")
;                  _la_display($mLstSq.Qy, "Ql")
;                  _la_display($mLstSq.Qr, "Qr")
;                  _la_display($mLstSq.sdX, "std devs(x)")
;                  _la_display($mLstSq.sdR, "std devs(r)")
;                  _la_display($mLstSq.sdY, "std devs(y)")
; ===============================================================================================================================
Func _la_lstsq($mA, $mB, $mP = Default, $sAlgorithm = "QR", $iFlagsResults = 0, $fTolerance = Default, $fTallSkinnyThreshhold = 10)
	Local $mRet, _
	      $iM = $mA.rows, $iN = $mA.cols, $iF = $iM - $iN, $sDataType = $mA.datatype

	; add dependencies for user defined flags
	If BitAND($iFlagsResults, $__LA_LSTSQ_SDY)        Then $iFlagsResults = BitOR($iFlagsResults, $__LA_LSTSQ_QR)
	If BitAND($iFlagsResults, $__LA_LSTSQ_SDY)        Then $iFlagsResults = BitOR($iFlagsResults, $__LA_LSTSQ_QYD)
	If BitAND($iFlagsResults, $__LA_LSTSQ_SDX)        Then $iFlagsResults = BitOR($iFlagsResults, $__LA_LSTSQ_QX, $__LA_LSTSQ_S0)
	If BitAND($iFlagsResults, $__LA_LSTSQ_S0)         Then $iFlagsResults = BitOR($iFlagsResults, $__LA_LSTSQ_R2Sum)
	If BitAND($iFlagsResults, $__LA_LSTSQ_REDUNDANCY) Then $iFlagsResults = BitOR($iFlagsResults, $__LA_LSTSQ_QR)
	If BitAND($iFlagsResults, $__LA_LSTSQ_QR)         Then $iFlagsResults = BitOR($iFlagsResults, $__LA_LSTSQ_QYD, $__LA_LSTSQ_QY)
	If BitAND($iFlagsResults, $__LA_LSTSQ_QYD)        Then $iFlagsResults = BitOR($iFlagsResults, $__LA_LSTSQ_QX)

	If BitAND($iFlagsResults, $__LA_LSTSQ_R2Sum + $__LA_LSTSQ_QY) Then
		Local $mPOrig = _la_duplicate($mP)
	EndIf

	Switch $sAlgorithm
		Case "SVD"
			If BitAND($iFlagsResults, $__LA_LSTSQ_R2Sum)  Then $iFlagsResults = BitOR($iFlagsResults, $__LA_LSTSQ_R)
			$mRet = __la_lstsq_svd($mA, $mB, $mP, $iFlagsResults, $fTolerance)
		Case "Cholesky"
			$mRet = __la_lstsq_cholesky($mA, $mB, $mP, $iFlagsResults)
		Case Else ; = "QR"
			$mRet = __la_lstsq_qr($mA, $mB, $mP, $iFlagsResults, $fTallSkinnyThreshhold)
	EndSwitch
	If @error Then Return SetError(10 + @error, @extended, Null)

	; calculate the residuals v
	Local $mV
	If BitAND($iFlagsResults, $__LA_LSTSQ_R) Then
		$mV = _la_duplicate($mB)
		; calculate v = A*X - b
		_blas_gemv($mA.ptr, $mRet.x.ptr, $mV.ptr, 1, -1, "N", 1, 1, $iM, $iN, $iM, $sDataType)
		$mRet.r = $mV
	EndIf

	; calculate vᵀ*P*v and s0
	If BitAND($iFlagsResults, $__LA_LSTSQ_R2Sum) And Not MapExists($mRet, "r2sum") Then
		If IsKeyword($mPOrig) <> 1 Then
			Local $tTmp = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iM))
			Local $pTmp = DllStructGetPtr($tTmp)

			; calculate vᵀ*P*v
			If BitAND($mPOrig.storageType, $__g_BLAS_STYPE_MATRIX) Then ; P-Matrix
				; use symv because P should be a symmetric matrix
				_blas_symv($mPOrig.ptr, $mV.ptr, $pTmp, 1, 0, "U", $iM, $iM, 1, 1, $sDataType)

			Else ; P-Vector
				; handle P-vector as packed diagonal band matrix to use sbmv(most efficient)
				_blas_sbmv($mPOrig.ptr, $mV.ptr, $pTmp, 1, 0, 0, "U", $iM, 1, 1, 1, $sDataType)

			EndIf

			$mRet.r2sum = _blas_dot($mV.ptr, $pTmp, 0, 0, 1, 1, $iM, $sDataType)

		Else ; non weighted problem
			$mRet.r2sum = _lp_lassq($mV.ptr)

		EndIf

	EndIf

	; calculate σ₀
	If BitAND($iFlagsResults, $__LA_LSTSQ_S0) Then $mRet.s0 = Sqrt($mRet.r2sum / $iF)

	; calculate the standard deviations for every parameter in x
	If BitAND($iFlagsResults, $__LA_LSTSQ_SDX) Then
		; extract the standard deviations for every parameter
		Local $mStdDevs = _la_getDiag($mRet.Qx)
		; multiply with σ₀² --> variance for every element
		_blas_scal($mStdDevs.ptr, $mRet.s0^2 , 0, 1, $iN, $sDataType)
		; sqrt() for every element to determine the standard deviation
		_la_sqrtElements($mStdDevs, True)

		$mRet.sdX = $mStdDevs
	EndIf

	; retrieve a-priori cofactormatrix of the observations Q_l0 (= user input through P)
	If BitAND($iFlagsResults, $__LA_LSTSQ_QY) Then
		If IsKeyword($mPOrig) <> 1 Then
			If BitAND($mPOrig.storageType, $__g_BLAS_STYPE_MATRIX) Then
			; full P-Matrix

				Local $mY = _la_createIdentity($iM, $iM, $sDatatype)
				_lp_sysv($mPOrig, $mY, $iM)
				$mRet.Qy0 = $mY

			Else
			; P-Vector = Diagonal Matrix
				Local $tP = $mPOrig.struct
				For $i = 1 To $iM
					DllStructSetData($tP, 1, 1 / DllStructGetData($tP, 1, $i), $i)
				Next
				$mRet.Qy0 = _la_VectorToDiag($mPOrig)
			EndIf
		Else
		; unweighted case
			$mRet.Qy0 = _la_createIdentity($iM, $iM, $sDataType)

		EndIf
	EndIf

	; calculate a-posteriori cofactor matrix of the adjusted observations Q_ld
	If BitAND($iFlagsResults, $__LA_LSTSQ_QYD) Then
		Local $mC = _blas_createMatrix($iM, $iN, $sDataType) ; temporary matrix
		Local $mQl = _blas_createMatrix($iM, $iM, $sDataType)

		; calc A * Qₓ --> C (use symmetry of Qₓ)
		_blas_symm($mRet.Qx, $mA, $mC, 1, 0, "U", "R", $iM, $iN, $iN, $iM, $iM, $sDataType)

		; calc C * Aᵀ
		_blas_gemm($mC.ptr, $mA.ptr, $mQl.ptr, 1, 0, "N", "T", $iM, $iM, $iN, $iM, $iM, $iM, $sDataType)

		$mRet.Qy = $mQl
	EndIf

	; calculate cofactor matrix of the residuals Q_v
	If BitAND($iFlagsResults, $__LA_LSTSQ_QR) Then
		Local $mQl0 = $mRet.Qy0
		$mRet.Qr = _la_sub($mQl0, $mQl)
	EndIf

	; calculate the standard deviations for the residuals
	If BitAND($iFlagsResults, $__LA_LSTSQ_SDY) Then
		; extract the standard deviations for every parameter
		Local $mStdDevsR = _la_getDiag($mRet.Qr)
		; multiply with σ₀² --> variance for every element
		_blas_scal($mStdDevsR.ptr, $mRet.s0^2 , 0, 1, $iM, $sDataType)
		; sqrt() for every element to determine the standard deviation
		_la_sqrtElements($mStdDevsR, True)

		$mRet.sdR = $mStdDevsR
	EndIf

	; standard deviations for the adjusted observations
	If BitAND($iFlagsResults, $__LA_LSTSQ_SDY) Then
		; extract the standard deviations for every parameter
		Local $mStdDevsL = _la_getDiag($mRet.Qy)
		; multiply with σ₀² --> variance for every element
		_blas_scal($mStdDevsL.ptr, $mRet.s0^2 , 0, 1, $iM, $sDataType)
		; sqrt() for every element to determine the standard deviation
		_la_sqrtElements($mStdDevsL, True)

		$mRet.sdY = $mStdDevsL
	EndIf

	Return SetError(@error, @extended, $mRet)

EndFunc


; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __la_lstsq_qr()
; Description ...: solves overdetermined or underdetermined [weighted] linear system using QR decomposition
; Syntax ........: __la_lstsq_qr($mAOrig, $mBOrig, [$mP = Default, [$iFlagsResults = BitOR($__LA_LSTSQ_QX, $__LA_LSTSQ_QR, $__LA_LSTSQ_R, $__LA_LSTSQ_S0, $__LA_LSTSQ_COND, $__LA_LSTSQ_R2Sum, $__LA_LSTSQ_SDX]])
; Parameters ....: mAOrig                - [Map] jacobian matrix A (M × N) as a map, DllStruct or pointer
;                  mBOrig                - [Map] observation vector/matrix as a map, DllStruct or pointer
;                  mP                    - [Map] (Default: Default)
;                                        ↳ weight matrix or weight vector P as a map, DllStruct or pointer (Default = unweighted)
;                  iFlagsResults         - [UInt] (Default: BitOR($__LA_LSTSQ_QX)
;                                        ↳ output components flags (see description in _la_lstsq() for details)
;                  fTallSkinnyThreshhold - [Float] (Default: 10)
;                                        ↳ threshold value of the ratio rows/columns from which the QR decomposition switches to the tall-skinny algorithm
; Return value ..: Success: [Map] results of the solution depending on iFlagsResults (see description in _la_lstsq() for details - not all flags used here)
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mAOrig/mBOrig (@extended = 1: mAOrig, @extended = 2: mBOrig)
;                           | 2: invalid dimension combination between mAOrig and mBOrig
;                           |1X: error X during _lp_gels/_lp_getsls() (@extended: @extended from _lp_gels/_lp_getsls())
;                           |2X: error X during _lp_lacpy()           (@extended: @extended from _lp_lacpy())
;                           |3X: error X during _lp_potrf()           (@extended: @extended from _lp_potrf())
;                           |4X: error X during _blas_trmm()          (@extended: @extended from _blas_trmm())
;                           |5X: error X during _blas_trmv/trmm()     (@extended: @extended from _blas_trmv/trmm())
;                           |6X: error X during _la_sqrtElements()    (@extended: @extended from _la_sqrtElements())
;                           |7X: error X during _blas_sbmv()          (@extended: @extended from _blas_sbmv())
;                           |8X: error X during _blas_trmm()          (@extended: @extended from _blas_trmm())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _lp_gels(), _lp_getsls()
; Link ..........:
; Example .......: Yes
;                  Global $mA = _la_fromArray("[[1,0,0,0,0],[0,1,0,0,0],[-1,0,1,0,0],[-1,0,0,1,0],[-1,0,0,0,1],[0,-1,1,0,0],[0,-1,0,1,0],[0,-1,0,0,1],[0,0,1,-1,0],[0,0,0,1,-1]]")
;                  Global $mb = _la_fromArray("[0,0,0.055,0.001,0.057,0.014,-0.002,0.051,0.015,0.052]")
;                  Global $mLstSq = __la_lstsq_qr($mA, $mB)
;                  ConsoleWrite("vᵀPv: " & $mLstSq.r2sum & @CRLF)
;                  _la_display($mLstSq.x, "x")
; ===============================================================================================================================
Func __la_lstsq_qr($mAOrig, $mBOrig, $mP = Default, $iFlagsResults = BitOR($__LA_LSTSQ_QX,$__LA_LSTSQ_QR,$__LA_LSTSQ_R,$__LA_LSTSQ_S0,$__LA_LSTSQ_COND,$__LA_LSTSQ_R2Sum, $__LA_LSTSQ_SDX), $fTallSkinnyThreshhold = 10)
	; direct AutoIt-type input
	If IsArray($mAOrig) Or IsString($mAOrig) Then $mAOrig = _blas_fromArray($mAOrig)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not MapExists($mAOrig, "ptr") Then Return SetError(1, 1, Null)
	If Not MapExists($mBOrig, "ptr") Then Return SetError(1, 2, Null)

	; dimension check
	If $mAOrig.rows <> $mBOrig.rows Then Return SetError(2, 0, Null)

	; dimension parameters
	Local $bVecB = $mBOrig.cols = 0 ? True : False, _
	      $iM    = $mAOrig.rows, _
	      $iN    = $mAOrig.cols, _
	      $iK    = $mBOrig.cols < 1 ? 1 : $mBOrig.cols

	; working variables
	Local $mResults[], _
          $sDatatype = $mAOrig.datatype, _
		  $mTmpM, $mX, $mY, $mR, _
		  $tTmp, $pTmp _

	; degrees of freedom
	Local $iF   = $iM - $iN
	$mResults.f = $iF

	; duplicate input matrix to prevent overwrite if option choosed
	Local $mA     = _la_duplicate($mAOrig), _   ; m x n
	      $mB     = _la_duplicate($mBOrig)      ; m x k

	; handle different cases for P
 	; and adjust A and B with P to continue the adjustment
	If IsKeyword($mP) <> 1 Then
		If BitAND($mP.storageType, $__g_BLAS_STYPE_MATRIX) Then
		; full P-Matrix
			; P -> L * Lᵀ   (Cholesky-factorization)
			_lp_potrf($mP.ptr, "L", $mP.rows, $mP.rows, $sDataType)
			If @error Then Return SetError(@error + 30, @extended, Null)

			; L * A = Ad
			_blas_trmm($mP.ptr, $mA.ptr, 1, "L", "L", "N", "N", $iM, $iN, $iM, $iM, $sDataType)
			If @error Then Return SetError(@error + 40, @extended, Null)

			; L * B = Bd
			If $iK = 1 Then
				_blas_trmv($mP.ptr, $mB.ptr, "L", "N", "N", $iM, 1, $iM, $sDataType)
			Else
				_blas_trmm($mP.ptr, $mB.ptr, 1, "L", "L", "N", "N", $iM, $iK, $iM, $iM, $sDataType)
			EndIf
			If @error Then Return SetError(@error + 50, @extended, Null)

		Else
		; P-Vector = Diagonal Matrix
			; calculate sqrt(i) for every element in P:
			$mTmpM = _la_sqrtElements($mP, False)
			If @error Then Return SetError(@error + 60, @extended, Null)

			; P * b --> b
			$tTmp = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iM))
			$pTmp = DllStructGetPtr($tTmp)
			_blas_sbmv($mTmpM.ptr, $mB.ptr, $pTmp, 1, 1, 0, "L", $iM, 1, 1, 1, $sDataType)
			If @error Then Return SetError(@error + 70, @extended, Null)
			$mB.struct = $tTmp
			$mB.ptr = $pTmp

			; P-vector to P- diagonal matrix
			_la_VectorToDiag($mTmpM, True)

			; P * A --> A
			_blas_trmm($mTmpM.ptr, $mA.ptr, 1, "L", "U", "N", "N", $iM, $iN, $iM, $iM, $sDataType)
			If @error Then Return SetError(@error + 80, @extended, Null)
		EndIf
	EndIf

	; call the LAPACK functions for solving least squares problems by QR factorization
	;~ If $bRRQR Then ; solve by using RRQR-algorithm
	; TODO: still not working correct: R-part in A is not the same as in dgels, so Qx-determination fails
	;~ 	Local $tJPVT = _lp_gelsy($mA.ptr, $mB.ptr, $iK, $iM, $iN, $iM, $iM, -1, $sDatatype)
	;~ 	$mResults.rank = @extended

	;~ 	; rearrange the R-Matrix columns
	;~ 	_lp_lapmt($mA.ptr, $tJPVT, False, $iM, $iN, $iM, $sDataType)

	;~ Else
	If $iM / $iN > $fTallSkinnyThreshhold Then ; tall-skinny problem - use optimized algorithm for this case
		_lp_getsls($mA.ptr, $mB.ptr, "N", $iK, $iM, $iN, $iM, $iM, $sDatatype)
	Else ; normal QR factorization
		_lp_gels($mA.ptr, $mB.ptr, "N", $iK, $iM, $iN, $iM, $iM, $sDatatype)
	EndIf
	If @error Then Return SetError(@error + 10, @extended, Null)

	; extract the first n lines of B to get the solution vector/matrix x
	$mX = $bVecB ? _blas_createVector($iN, $sDatatype) : _blas_createMatrix($iN, $iK, $sDatatype)
	_lp_lacpy($mB, $mX, "X", $iN, $iK, $iM, $iK, $sDatatype)
	If @error Then Return SetError(@error + 20, @extended, Null)
	$mResults.x = $mX

	; calculate vᵀ*P*v and s0
	If BitAND($iFlagsResults, $__LA_LSTSQ_R2Sum) Then
		; calculate the residual sum (the last elements in B)
		$mResults.r2sum = _blas_nrm2($mB, $iN, 1, $iM - $iN, $sDataType)^2
	EndIf

	; calculate cofactor matrix for the solved parameters Qₓ
	If BitAND($iFlagsResults, $__LA_LSTSQ_QX) Then
		; extract R from A:
		$mR = _la_getTriangle($mA, "U", True, True)

		; identity matrix for solving Rᵀ * R * Qₓ = I
		$mY = _la_createIdentity($iN, $iN, $sDatatype)

		; step 1: solve Rᵀ * Y = I
		_lp_trtrs($mR.ptr, $mY.ptr, "U", "N", "T", $iN, $iN, $iN, $iN, $sDatatype)
		; step 2: solve R * Qx = Y
		_lp_trtrs($mR.ptr, $mY.ptr, "U", "N", "N", $iN, $iN, $iN, $iN, $sDatatype)

		; add flags to mark the matrix as upper symmetric
		$mY.storageType = BitOR($mY.storageType, $__g_BLAS_STYPE_SYMMETRIC + $__g_BLAS_STYPE_UPPER)

		$mResults.Qx = $mY
		; Alternative:
		; - extract triangular part of R
		; - calc (Rᵀ * R) with dsyrk
		; - solve (Rᵀ*R) * Qx = I with dposv
		; BUT: this should be less efficient and less numerical stable
	EndIf

	Return $mResults

EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __la_lstsq_svd()
; Description ...: solves overdetermined or underdetermined [weighted] linear system using singular value decomposition
; Syntax ........: __la_lstsq_svd($mAOrig, $mBOrig, [$mP = Default, [$iFlagsResults = BitOR($__LA_LSTSQ_QX, $__LA_LSTSQ_QR, $__LA_LSTSQ_R, $__LA_LSTSQ_S0, $__LA_LSTSQ_COND, $__LA_LSTSQ_R2Sum, $__LA_LSTSQ_SDX, $__LA_LSTSQ_RANK]])
; Parameters ....: mAOrig           - [Map] jacobian matrix A (M × N) as a map, DllStruct or pointer
;                  mBOrig           - [Map] observation vector/matrix as a map, DllStruct or pointer
;                  mP               - [Map] (Default: Default)
;                                   ↳ weight matrix or weight vector P as a map, DllStruct or pointer (Default = unweighted)
;                  iFlagsResults    - [UInt] (Default: BitOR($__LA_LSTSQ_QX)
;                                   ↳ output components flags (see description in _la_lstsq() for details)
;                  fTolerance       - [Float] (Default: Default)
;                                   ↳ threshold value up to which singular values are regarded as zero (-1 = machine precision)
; Return value ..: Success: [Map] results of the solution depending on iFlagsResults (see description in _la_lstsq() for details - not all flags used here)
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mAOrig/mBOrig (@extended = 1: mAOrig, @extended = 2: mBOrig)
;                           | 2: invalid dimension combination between mAOrig and mBOrig
;                           |1X: error X during _lp_gelss()        (@extended: @extended from _lp_gelss())
;                           |2X: error X during _lp_lacpy()        (@extended: @extended from _lp_lacpy())
;                           |3X: error X during _lp_potrf()        (@extended: @extended from _lp_potrf())
;                           |4X: error X during _blas_trmm()       (@extended: @extended from _blas_trmm())
;                           |5X: error X during _blas_trmv/trmm()  (@extended: @extended from _blas_trmv/trmm())
;                           |6X: error X during _la_sqrtElements() (@extended: @extended from _la_sqrtElements())
;                           |7X: error X during _blas_sbmv()       (@extended: @extended from _blas_sbmv())
;                           |8X: error X during _blas_trmm()       (@extended: @extended from _blas_trmm())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _lp_gelss()
; Link ..........:
; Example .......: Yes
;                  Global $mA = _la_fromArray("[[1,0,0,0,0],[0,1,0,0,0],[-1,0,1,0,0],[-1,0,0,1,0],[-1,0,0,0,1],[0,-1,1,0,0],[0,-1,0,1,0],[0,-1,0,0,1],[0,0,1,-1,0],[0,0,0,1,-1]]")
;                  Global $mb = _la_fromArray("[0,0,0.055,0.001,0.057,0.014,-0.002,0.051,0.015,0.052]")
;                  Global $mLstSq = __la_lstsq_svd($mA, $mB)
;                  ConsoleWrite(StringFormat("\n% 20s: % 6.4f\n% 20s: %-5d\n% 20s: %-5d\n\n", "condition nr.", $mLstSq.cond, "rank", $mLstSq.rank, "degrees of freedom", $mLstSq.f))
;                  _la_display($mLstSq.x, "solution")
; ===============================================================================================================================
Func __la_lstsq_svd($mAOrig, $mBOrig, $mP = Default, $iFlagsResults = BitOR($__LA_LSTSQ_QX,$__LA_LSTSQ_QR,$__LA_LSTSQ_R,$__LA_LSTSQ_S0,$__LA_LSTSQ_COND,$__LA_LSTSQ_R2Sum, $__LA_LSTSQ_SDX, $__LA_LSTSQ_RANK), $fTolerance = Default)
	; direct AutoIt-type input
	If IsArray($mAOrig) Or IsString($mAOrig) Then $mAOrig = _blas_fromArray($mAOrig)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not MapExists($mAOrig, "ptr") Then Return SetError(1, 1, Null)
	If Not MapExists($mBOrig, "ptr") Then Return SetError(1, 2, Null)

	; dimension check
	If $mAOrig.rows <> $mBOrig.rows Then Return SetError(2, 0, Null)

	; calculate tolerance
	If IsKeyword($fTolerance) = 1 Then $fTolerance = 10 * ($mAOrig.datatype = "FLOAT" ? $f_LA_FLT_PREC : $f_LA_DBL_PREC)

	; dimension parameters
	Local $bVecB = $mBOrig.cols = 0 ? True : False, _
	      $iM    = $mAOrig.rows, _
	      $iN    = $mAOrig.cols, _
	      $iK    = $mBOrig.cols < 1 ? 1 : $mBOrig.cols

	; working variables
	Local $mResults[], _
          $sDatatype = $mAOrig.datatype, _
		  $mSVD, $tS, _
		  $mTmpM, $mX, _
		  $tTmp, $pTmp _

	; degrees of freedom
	Local $iF   = $iM - $iN
	$mResults.f = $iF

	; duplicate input matrix to prevent overwrite if option choosed
	Local $mA = _la_duplicate($mAOrig), _   ; m x n
	      $mB = _la_duplicate($mBOrig)      ; m x k

	; handle different cases for P
 	; and adjust A and B with P to continue the adjustment
	If IsKeyword($mP) <> 1 Then
		If BitAND($mP.storageType, $__g_BLAS_STYPE_MATRIX) Then
		; full P-Matrix
			; P -> L * Lᵀ   (Cholesky-factorization)
			_lp_potrf($mP.ptr, "L", $mP.rows, $mP.rows, $sDataType)
			If @error Then Return SetError(@error + 30, @extended, Null)

			; L * A = Ad
			_blas_trmm($mP.ptr, $mA.ptr, 1, "L", "L", "N", "N", $iM, $iN, $iM, $iM, $sDataType)
			If @error Then Return SetError(@error + 40, @extended, Null)

			; L * B = Bd
			If $iK = 1 Then
				_blas_trmv($mP.ptr, $mB.ptr, "L", "N", "N", $iM, 1, $iM, $sDataType)
			Else
				_blas_trmm($mP.ptr, $mB.ptr, 1, "L", "L", "N", "N", $iM, $iK, $iM, $iM, $sDataType)
			EndIf
			If @error Then Return SetError(@error + 50, @extended, Null)

		Else
		; P-Vector = Diagonal Matrix
			; calculate sqrt(i) for every element in P:
			$mTmpM = _la_sqrtElements($mP, False)
			If @error Then Return SetError(@error + 60, @extended, Null)

			; P * b --> b
			$tTmp = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iM))
			$pTmp = DllStructGetPtr($tTmp)
			_blas_sbmv($mTmpM.ptr, $mB.ptr, $pTmp, 1, 1, 0, "L", $iM, 1, 1, 1, $sDataType)
			If @error Then Return SetError(@error + 70, @extended, Null)
			$mB.struct = $tTmp
			$mB.ptr = $pTmp

			; P-vector to P- diagonal matrix
			_la_VectorToDiag($mTmpM, True)

			; P * A --> A
			_blas_trmm($mTmpM.ptr, $mA.ptr, 1, "L", "U", "N", "N", $iM, $iN, $iM, $iM, $sDataType)
			If @error Then Return SetError(@error + 80, @extended, Null)

		EndIf
	EndIf

	If BitAND($iFlagsResults, $__LA_LSTSQ_QX) Then
	; manual svd solution (because with results of gelss you cannot derive Qₓₓ)

		$mSVD = _lp_gesvd($mA, $iM, $iN, $iM, $sDataType)
		Local $mS  = $mSVD.S  ; n x n --> diag (here vector only)
		Local $mU  = $mSVD.U  ; m x m --> general
		Local $mVt = $mSVD.VT ; n x n --> general

		$mResults.cond = DllStructGetData($mS.struct, 1, 1) / DllStructGetData($mS.struct, 1, $iM < $iN ? $iM : $iN)

		; pseudo inverse σ⁺ = inv(el) of Σ
		Local $mSp = _blas_duplicate($mS)
		_la_invElements($mSp, True)
		_la_VectorToDiag($mSp, True)

		; x = V * Σ⁺ * Uᵀ * B
		Local $mTmp = _blas_createMatrix($iN, $iN, $sDataType)
		$mX   = _blas_createVector($iN, $sDataType)
		_blas_gemv($mU.ptr, $mB.ptr, $mX.ptr, 1, 0, "T", 1, 1, $iM, $iN, $iM, $sDataType)
		_blas_gemv($mSP.ptr, $mX.ptr, $mTmp.ptr, 1, 0, "N", 1, 1, $iN, $iN, $iN, $sDataType)
		_blas_gemv($mVt.ptr, $mTmp.ptr, $mX.ptr, 1, 0, "T", 1, 1, $iN, $iN, $iN, $sDataType)

		$mResults.x = $mX

		; Qₓ = V * (Σᵀ * Σ)⁻¹ * Vᵀ
		$tS = $mS.struct
		For $i = 1 To $iN
			DllStructSetData($tS, 1, 1 / (DllStructGetData($tS, 1, $i)^2), $i)
		Next
		Local $mQx = _la_VectorToDiag($mS) ; --> (Σᵀ * Σ)⁻¹
		; alternative (but maybe irregular values)
		;~ _la_squareElements($mS, True)      ; Σᵀ * Σ
		;~ Local $mQx = _la_invElements($mS)
		;~ _la_VectorToDiag($mQx, True)
		_blas_gemm($mVt.ptr, $mQx.ptr, $mTmp.ptr, 1, 0, "T", "N", $iN, $iN, $iN, $iN, $iN, $iN, $sDataType) ; V * [...]
		_blas_gemm($mTmp.ptr, $mVt.ptr, $mQx.ptr, 1, 0, "N", "N", $iN, $iN, $iN, $iN, $iN, $iN, $sDataType) ; [...] * Vᵀ

		; add flags to mark the matrix as upper symmetric
		$mQx.storageType = BitOR($mQx.storageType, $__g_BLAS_STYPE_SYMMETRIC + $__g_BLAS_STYPE_UPPER)

		$mResults.Qx = $mQx

		If BitAND($iFlagsResults, $__LA_LSTSQ_RANK) Then
		; determine rank
			Local $iRank
			For $iRank = $iM < $iN ? $iM : $iN To 1 Step -1
				If DllStructGetData($mS.struct, 1, $iRank) > $fTolerance Then
					$mResults.rank = $iRank
					ExitLoop
				EndIf
			Next
		EndIf

	Else
	; direct use of gelss

		$tS = _lp_gelss($mA.ptr, $mB.ptr, "N", $iK, $iM, $iN, $iM, $iM, -1, $sDataType)
		If @error Then Return SetError(@error + 10, @extended, Null)
		$mResults.rank = @extended
		$mResults.cond = DllStructGetData($tS, 1, 1) / DllStructGetData($tS, 1, $iM < $iN ? $iM : $iN)

		; extract the first n lines of B to get the solution vector/matrix x
		$mX = $bVecB ? _blas_createVector($iN, $sDatatype) : _blas_createMatrix($iN, $iK, $sDatatype)
		_lp_lacpy($mB, $mX, "X", $iN, $iK, $iM, $iK, $sDatatype)
		If @error Then Return SetError(@error + 20, @extended, Null)
		$mResults.x = $mX
	EndIf

	Return $mResults

EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __la_lstsq_cholesky()
; Description ...: solves overdetermined or underdetermined [weighted] linear system using cholesky decomposition
; Syntax ........: __la_lstsq_cholesky($mAOrig, $mBOrig, [$mP = Default, [$iFlagsResults = BitOR($__LA_LSTSQ_QX, $__LA_LSTSQ_QR, $__LA_LSTSQ_R, $__LA_LSTSQ_S0, $__LA_LSTSQ_COND, $__LA_LSTSQ_R2Sum, $__LA_LSTSQ_SDX]])
; Parameters ....: mAOrig                - [Map] jacobian matrix A (M × N) as a map, DllStruct or pointer
;                  mBOrig                - [Map] observation vector/matrix as a map, DllStruct or pointer
;                  mP                    - [Map] (Default: Default)
;                                        ↳ weight matrix or weight vector P as a map, DllStruct or pointer (Default = unweighted)
;                  iFlagsResults         - [UInt] (Default: BitOR($__LA_LSTSQ_QX)
;                                        ↳ output components flags (see description in _la_lstsq() for details)
; Return value ..: Success: results of the solution depending on iFlagsResults (see description in _la_lstsq() for details - not all flags used here)
;                  Failure: Null and set @error to:
;                           | 1: invalid value for mAOrig/mBOrig (@extended = 1: mAOrig, @extended = 2: mBOrig)
;                           | 2: invalid dimension combination between mAOrig and mBOrig
;                           |2X: error X during _lp_lacpy()           (@extended: @extended from _lp_lacpy())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: _lp_potrf(), _lp_potrs()
; Link ..........:
; Example .......: Yes
;                  Global $mA = _la_fromArray("[[1,0,0,0,0],[0,1,0,0,0],[-1,0,1,0,0],[-1,0,0,1,0],[-1,0,0,0,1],[0,-1,1,0,0],[0,-1,0,1,0],[0,-1,0,0,1],[0,0,1,-1,0],[0,0,0,1,-1]]")
;                  Global $mb = _la_fromArray("[0,0,0.055,0.001,0.057,0.014,-0.002,0.051,0.015,0.052]")
;                  Global $mP = _la_fromArray("[10,10,0.7,1.0,0.7,0.7,1,0.7,1,1]")
;                  Global $mLstSq = __la_lstsq_cholesky($mA, $mB, $mP, 0.001)
;                  _la_display($mLstSq.x, "solution")
; ===============================================================================================================================
Func __la_lstsq_cholesky($mAOrig, $mBOrig, $mP = Default, $iFlagsResults = BitOR($__LA_LSTSQ_QX, $__LA_LSTSQ_QR,$__LA_LSTSQ_R,$__LA_LSTSQ_S0,$__LA_LSTSQ_COND,$__LA_LSTSQ_R2Sum, $__LA_LSTSQ_SDX))
	; direct AutoIt-type input
	If IsArray($mAOrig) Or IsString($mAOrig) Then $mAOrig = _blas_fromArray($mAOrig)

	; check if Input is a valid AutoIt-BLAS/LAPACK-Map
	If Not MapExists($mAOrig, "ptr") Then Return SetError(1, 1, Null)
	If Not MapExists($mBOrig, "ptr") Then Return SetError(1, 2, Null)

	; dimension check
	If $mAOrig.rows <> $mBOrig.rows Then Return SetError(2, 0, Null)

	; dimension parameters
	Local $iM    = $mAOrig.rows, _
	      $iN    = $mAOrig.cols

	; working variables
	Local $mResults[], _
          $sDatatype = $mAOrig.datatype, _
		  $mATP, $mN, $mX, $mL, _
		  $tTmp, $iStart _

	; degrees of freedom
	Local $iF   = $iM - $iN
	$mResults.f = $iF

	; duplicate input matrix to prevent overwrite if option choosed
	Local $mB = _la_duplicate($mBOrig)      ; m x k

	If IsKeyword($mP) <> 1 Then
		If BitAND($mP.storageType, $__g_BLAS_STYPE_MATRIX) Then
		; full P-Matrix

			$mATP = _blas_createMatrix($iN, $iM, $sDataType)
			; ToDo: replace with dsymm/dtrmm (SIDE = "R")
			_blas_gemm($mAOrig.ptr, $mP.ptr, $mATP.ptr, 1, 0, "T", "N", $iN, $iM, $iM, $iM, $iM, $iN, $sDataType)

		Else
		; P-Vector = Diagonal Matrix
			; calc Aᵀ * P  (with P = diagonal matrix)
			; scale every column is much more efficient and stable then _blas_gemm for this case
			; anyway: because it`s an AutoIt-loop it still could be slower than the brute force _blas_gemm method
			$mATP = _la_transpose($mAOrig)
			$tTmp = $mP.struct
			$iStart = 0
			For $i = 1 To $iM
				_blas_scal($mATP.ptr, DllStructGetData($tTmp, 1, $i), $iStart, 1, $iN, $sDataType)
				$iStart += $iN
			Next
		EndIf

		; N = ATP * A
		$mN = _blas_createMatrix($iN, $iN, $sDataType)
		_blas_gemm($mATP.ptr, $mAOrig.ptr, $mN.ptr, 1, 0, "N", "N", $iN, $iN, $iM, $iN, $iM, $iN, $sDataType)

		; n = ATP * b
		$mL = _blas_createVector($iN, $sDataType)
		_blas_gemv($mATP.ptr, $mB.ptr, $mL.ptr, 1, 0, "N", 1, 1, $iN, $iM, $iN, $sDataType)

	Else
	; P = Identity Matrix (default case / "unweighted" / "uncorrelated" case)
		$mN = _blas_duplicate($mAOrig)
		$mL = _blas_duplicate($mBOrig)

	EndIf

	; idea of solving with cholesky:
	; normal equations: Aᵀ * P * A * x = Aᵀ * P * b
	;                 :           N * x = n
	; cholesky factorize N:       N = L * Lᵀ
	; solve   L * Lᵀ * x = n

	; cholesky factorization of N --> L * Lᵀ
	_lp_potrf($mN.ptr, "U", $iN, $iN, $sDataType) ; --> L

	; solve L * Lᵀ * x = n
	_lp_potrs($mN, $mL, "U", 1, $iN, $iN, $iN, $sDataType)  ; n x n * n = n

	; extract first N values (= solution values)
	$mX = _blas_createVector($iN, $sDatatype)
	_lp_lacpy($mL, $mX, "X", $iN, 1, $iN, 1, $sDatatype)
	If @error Then Return SetError(@error + 20, @extended, Null)
	$mResults.x = $mX

	If BitAND($iFlagsResults, $__LA_LSTSQ_QX) Then
	; calculate complete Qx
		_lp_potri($mN.ptr, "U", $iN, $iN, $sDataType) ; inverse of N
		$mN.storageType = BitOR($mN.storageType, $__g_BLAS_STYPE_SYMMETRIC + $__g_BLAS_STYPE_UPPER) ; add flags to mark the matrix as upper symmetric
		If BitAND($iFlagsResults, $__LA_LSTSQ_QX) Then $mResults.Qx = _blas_duplicate($mN)

		; to check: maybe it`s more efficient to calc x by Qx * n in this case instead of potrs
		; calc x = 	Qx * n
		;~ $mXd = _blas_createVector($iN, $sDataType)
		;~ _blas_symv($mN, $mL, $mXd, 1, 0, "L")
		;~ $mResults.x = $mXd

		; remark: it`s possible to calculate the inverse value for the diagonal elements of Qx only
		; for this you have to calculate the forward and backward substiutions (with dtrsv) for the unit vectors
		; but this would lead to AutoIt loops which are much less efficient than calculating the whole inverse by _lp_potri
	EndIf

	Return $mResults
EndFunc

#EndRegion


#Region regression

; #FUNCTION# ====================================================================================================================
; Name ..........: _la_regression()
; Description ...: calculates an n-dimensional linear or non-linear regression
; Syntax ........: _la_regression($sFunc, $mVars, $mY, [$iFlagsLstSq = 0, [$mApprox = Default, [$sSolutionMethod = "GN", [$fApproxDefault = 1.0, [$iMaxIterations = 50, [$fEpsTermination = 1e-7, [$fLambda = 0.01, [$sLstSqAlgo = "QR", [$sDeriveMethod = "Central"]]]]]]]]])
; Parameters ....: sFunc           - [String] function equation as AutoIt syntax string where (see examples):
;                                    names are parameters to be estimated
;                                    names beginning with "$" are independent function values
;                  mVars           - [Array/Map] Array: independent function values for 1D case; Map of Arrays: named independent function values especially for n-dimensional case
;                  mY              - [Array/Map] dependent function values as array or vector as a map
;                  iFlagsLstSq     - [UInt] (Default: 0)
;                                  ↳ output components flags (see description in _la_lstsq() for details)
;                  mApprox         - [Map] (Default: Default)
;                                  ↳ approximate values for the model parameters in the case of non-linear correlation of the form
;                                    {"ParamName1": fValue1, "ParamName2": fValue2, ...}
;                  sSolutionMethod - [String] (Default: "GN") iterative solution algorithm (especially for non-linear case)
;                                  ↳ "GN": Gauss-Newton method           - fast but not as robust with poor approximation values
;                                    "LM": Levenberg-Marquardt algorithm - slightly slower but more robust with less good approximation values
;                  fApproxDefault  - [Float] (Default: 1.0)
;                                  ↳ default value for the output values of the parameters if mApprox is not set
;                  iMaxIterations  - [UInt] (Default: 50)
;                                  ↳ maximum number of iterations
;                  fEpsTermination - [Float] (Default: 1e-7)
;                                  ↳ threshold value between the results of two iterations from which equality is assumed (termination criterion)
;                  fLambda         - [Float] (Default: 0.01)
;                                  ↳ initial lambda value for the Levenberg-Marquardt model
;                  sLstSqAlgo      - [String] (Default: "QR") algorithm for solving the system
;                                  ↳ "QR": QR decomposition                    - good stability and performance
;                                    "SVD": singular value decomposition (SVD) - highest stability but lowest performance
;                                    "Cholesky": cholesky decomposition        - better stability and performance than QR but only for positive-definite matrices
;                  sDeriveMethod   - [String] (Default: "Central") algorithm to determine the 1st derivative
;                                  ↳ "Central"  - central difference quotient
;                                    "Central2", "Central3", "Central4" - central difference quotient of higher error order
;                                    "Forward"  - forward difference quotient
;                                    "Backward" - backward difference quotient
;                                    "Ridder"   - Ridder`s Method
;                                    "Higham"   - Highams algorithm
; Return value ..: Success: [Map] results as returned by _la_lstsq() plus:
;                           { "R2":      coefficient of determination,
;                             "R2_corr": adjusted coefficient of determination,
;                             "pearson": Pearson correlation coefficient,
;                             "params":  array of parameter names in correct order
;                             "vars":    array of variable names in correct order
;                             "Func":    the resulting formula with determined parameters as a string
;                           }
;                  Failure: Null and set @error to:
;                           | 1: error during __la_regression_prepareApproxVals() (@extended: @error from __la_regression_prepareApproxVals())
;                           |1X: error X during __la_regression_LevenbergMarquardt()/__la_regression_GaussNewton()  (@extended: @extended from them)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......: __la_regression_LevenbergMarquardt(), __la_regression_GaussNewton(), __la_regression_prepareApproxVals()
; Link ..........:
; Example .......: Yes (+ example_linear regression.au3, example_linear surface regression.au3)
;                  Global $sFunc = "B0 + B1 * $x + B2 * $x^2"
;                  Global $aX[] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
;                  Global $aY[] = [41,58,64,85,80,87,85,79,60,54,48]
;                  Global $mRegression = _la_regression($sFunc, $aX, $aY)
;                  ConsoleWrite(StringFormat("\n% 11s: %4.2f\n% 11s: %4.2f\n% 11s: %4.2f\n% 11s: %s\n\n", "R²", $mRegression.R2, "R² adjusted", $mRegression.R2_corr, "Pearson", $mRegression.pearson, "result func", $mRegression.Func))
;                  _la_display($mRegression.x, "iterations: " & @extended)
;                  _la_display($mRegression.sdX, "std. devs for X")
; ===============================================================================================================================
Func _la_regression($sFunc, $mVars, $mY, $iFlagsLstSq = 0, $mApprox = Default, $sSolutionMethod = "GN", $fApproxDefault = 1.0, $iMaxIterations = 50, $fEpsTermination = 1e-7, $fLambda = 0.01, $sLstSqAlgo = "QR", $sDeriveMethod = "Central")
	; adjust input parameters
	$sFunc = StringUpper($sFunc) ; because map keys are case sensitive
	; add vtpv to least square solutions because it is needed for R²
	$iFlagsLstSq = BitOR($iFlagsLstSq, $__LA_LSTSQ_R2Sum)
	; convert $mY to correct form
	If IsArray($mY) Or IsString($mY) Then $mY = _blas_fromArray($mY)
	; prepare the approximation values for further use
	If Not __la_regression_prepareApproxVals($mApprox) Then Return SetError(1, @error, Null)

	Local $mLstSq
	If $sSolutionMethod = "LM" Then ; Levenberg-Marquardt algorithm
		$mLstSq = __la_regression_LevenbergMarquardt($sFunc, $mVars, $mY, $mApprox, $fApproxDefault, $fLambda, $iMaxIterations, $fEpsTermination, $sLstSqAlgo, $iFlagsLstSq, $sDeriveMethod)
	Else ; "GN" Gauß-Newton algorithm
		$mLstSq = __la_regression_GaussNewton($sFunc, $mVars, $mY, $mApprox, $fApproxDefault, $iMaxIterations, $fEpsTermination, $sLstSqAlgo, $iFlagsLstSq, $sDeriveMethod)
	EndIf
	If @error Then Return SetError(10 + @error, @extended, Null)
	Local $iIterations = @extended

	Local $aParams = $mLstSq.params, $aX = _blas_toArray($mLstSq.x)
	For $i = 0 To UBound($aParams) - 1
		$sFunc = StringRegExpReplace($sFunc, '(?i)(\$\w+(*SKIP)(*FAIL)|\b\Q' & $aParams[$i] & '\E\b(?!\())', StringFormat("%.4g", $aX[$i]))
	Next
	$mLstSq.Func = StringReplace(StringRegExpReplace($sFunc, '(?i)\$(\w+)', "\1"), '+ -', '- ', 0, 1)

	Return SetExtended($iIterations, $mLstSq)
EndFunc


; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __la_regression_GaussNewton()
; Description ...: calculates an n-dimensional linear or non-linear regression using gauss markov model
; Syntax ........: __la_regression_GaussNewton($sFunc, $mVars, $mY, [$mApprox = Default, [$fApproxDefault = 1.0, [$iMaxIterations = 50, [$fEpsTermination = 1e-7, [$sLstSqAlgo = "QR", [$iFlagsLstSq = 0, [$sDeriveMethod = "Central"]]]]]]])
; Parameters ....: sFunc           - [String] function equation as AutoIt syntax string where (see examples):
;                                    names are parameters to be estimated
;                                    names beginning with "$" are independent function values
;                  mVars           - [Array/Map] Array: independent function values for 1D case; Map of Arrays: named independent function values especially for n-dimensional case
;                  mY              - [Array/Map] dependent function values as array or vector as a map
;                  mApprox         - [Array/Map] (Default: Default)
;                                  ↳ approximate values for the model parameters in the case of non-linear correlation of the form
;                                    {"ParamName1": fValue1, "ParamName2": fValue2, ...} or as corresponding 2D-Array[n][2]
;                  fApproxDefault  - [Float] (Default: 1.0)
;                                  ↳ default value for the output values of the parameters if mApprox is not set
;                  iMaxIterations  - [UInt] (Default: 50)
;                                  ↳ maximum number of iterations
;                  fEpsTermination - [Float] (Default: 1e-7)
;                                  ↳ threshold value between the results of two iterations from which equality is assumed (termination criterion)
;                  sLstSqAlgo      - [String] (Default: "QR") algorithm for solving the system
;                                  ↳ "QR": QR decomposition                    - good stability and performance
;                                    "SVD": singular value decomposition (SVD) - highest stability but lowest performance
;                                    "Cholesky": cholesky decomposition        - better stability and performance than QR but only for positive-definite matrices
;                  iFlagsLstSq     - [UInt] (Default: 0)
;                                  ↳ output components flags (see description in _la_lstsq() for details)
;                  sDeriveMethod   - [String] (Default: "Central") algorithm to determine the 1st derivative
;                                  ↳ "Central"  - central difference quotient
;                                    "Central2", "Central3", "Central4" - central difference quotient of higher error order
;                                    "Forward"  - forward difference quotient
;                                    "Backward" - backward difference quotient
;                                    "Ridder"   - Ridder`s Method
;                                    "Higham"   - Highams algorithm
; Return value ..: Success: [Map] results as returned by _la_lstsq() plus:
;                           { "R2":      coefficient of determination,
;                             "R2_corr": adjusted coefficient of determination,
;                             "pearson": Pearson correlation coefficient,
;                             "params":  array of parameter names in correct order
;                             "vars":    array of variable names in correct order
;                             "Func":    the resulting formula with determined parameters as a string }
;                  Failure: Null and set @error to:
;                           | 1: no parameters determined in sFunc (@extended: @error from StringRegExp)
;                           | 2: no variables determined in sFunc (@extended: @error from StringRegExp)
;                           | 3: error during __la_regression_prepareVars() (@extended: @error from __la_regression_prepareVars())
;                           | 4: different number of variables between sFunc and mVars (@extended: Ubound(aVars))
;                           | 5: regression did not converge until maximum number of iterations was reached (@extended: number of iterations )
;                           |1X: error X during _la_lstsq() (@extended: @extended from _la_lstsq())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: Yes
;                  Global $sFunc = "B0 + B1 * $x + B2 * $x^2"
;                  Global $aX[] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
;                  Global $aY[] = [41,58,64,85,80,87,85,79,60,54,48]
;                  Global $mRegression = __la_regression_GaussNewton($sFunc, $aX, $aY, Default, 1, 50, 1e-7, "QR", $__LA_LSTSQ_SDX) ; SDX for .sdX
;                  ConsoleWrite(StringFormat("\n% 11s: %4.2f\n% 11s: %4.2f\n% 11s: %4.2f\n\n", "R²", $mRegression.R2, "R² adjusted", $mRegression.R2_corr, "Pearson", $mRegression.pearson))
;                  _la_display($mRegression.x, "iterations: " & @extended)
;                  _la_display($mRegression.sdX, "std. devs for X")
; ===============================================================================================================================
Func __la_regression_GaussNewton($sFunc, $mVars, $mY, $mApprox = Default, $fApproxDefault = 1.0, $iMaxIterations = 50, $fEpsTermination = 1e-7, $sLstSqAlgo = "QR", $iFlagsLstSq = 0, $sDeriveMethod = "Central")
	Local $aParams, $mParams[], _ ; model parameter
	      $mXd, $mX0, $tX0, _     ; model parameter vectors
          $aVars, _               ; model variables
          $mA, _                  ; jacobian matrix
          $tA, _                  ; corresponding data structure
          $mYd, $mR, _            ; observation vectors
          $mLstSq, _              ; the least square solver return map
		  $i, $j                  ; additional helper variables like indices etc.

  	; adjust input parameters
	$sFunc = StringUpper($sFunc) ; because map keys are case sensitive
	; add vtpv to least square solutions because it is needed for R²
	$iFlagsLstSq = BitOR($iFlagsLstSq, $__LA_LSTSQ_R2Sum)

	; extract unknown params and assign their approximation value
	$aParams = StringRegExp($sFunc, '(?ix)(\$\w+(*SKIP)(*FAIL)|\b[a-z]\w*\b(?!\())', 3)
	If UBound($aParams) < 1 Then Return SetError(1, @error, Null)
	$aParams = _ArrayUnique($aParams, 0, 0, 1, 0)
	For $sParam In $aParams
		; use given approximation value if exists or the default value instead
		$mParams[$sParam] = MapExists($mApprox, $sParam) ? $mApprox[$sParam] : $fApproxDefault
	Next

	; extract the variables
	$aVars = StringRegExp($sFunc, '\$\K\w+', 3)
	If UBound($aVars) < 1 Then Return SetError(2, @error, Null)
	$aVars = _ArrayUnique($aVars, 0, 0, 1, 0)

	; prepare $mVars for further use
	If Not __la_regression_prepareVars($mVars, $aVars[0]) Then Return SetError(3, @error, Null)

	; $mVars and detected variable number should correspond
	If UBound($mVars) <> UBound($aVars) Then Return SetError(4, UBound($aVars), Null)

	Local $iM = UBound($mVars[$aVars[0]]), _
		  $iN = UBound($aParams)

	; declare the model solution vector yd = f(vars, X0)
	$mYd = _blas_createVector($iM, "DOUBLE")

	; declare the jacobian Matrix A (m x n) for the model parameters
	$mA = _blas_createMatrix($iM, $iN, "DOUBLE")
	$tA = $mA.struct

	; create approximation parameter vector X0
	$tX0 = DllStructCreate(StringFormat("DOUBLE[%d]", $iN))
	$i = 1
	For $sParam In $aParams
		DllStructSetData($tX0, 1, $mParams[$sParam], $i)
		$i += 1
	Next
	$mX0 = _blas_fromStruct($tX0, $iN, 0, "DOUBLE")


	For $i = 1 To $iMaxIterations - 1
		; retrieve the jacobian matrix A
		__la_Regression_getJacobiParams($sFunc, $mVars, $mParams, $tA, $iM, $iN, $sDeriveMethod)

		; derive the model solution vector yd to calculate the residuals to y
		__la_Regression_getYfromModel($sFunc, $mVars, $mParams, $mYd.struct, $iM)
		; shortened solution vector
		$mR = _la_sub($mY, $mYd)

		; least square solution of A * dX = r
		$mLstSq = _la_lstsq($mA, $mR, Default, $sLstSqAlgo, $iFlagsLstSq)
		If @error Then Return SetError(@error + 10, 0, Null)
		$mXd = $mLstSq.x

		; improved parameter vector
		_la_add($mXd, $mX0, True)
		$mLstSq.x = $mX0

		; if changes to the previous iteration are below the threshold value - stop
		If _la_norm($mXd) < $fEpsTermination Then ExitLoop

		; adjust the map of approximate values for the parameters
		$j = 1
		For $sParam In $aParams
			$mParams[$sParam] = DllStructGetData($tX0, 1, $j)
			$j += 1
		Next

	Next

	; throw error if result does not converge
	If $i >= $iMaxIterations Then Return SetError(5, $i, Null)

	;------ calculate R²
	; calc RSS = square sum of each difference Y_i with mean(Y)

	; calculate mean(Y)
	Local $fSumY = _la_sum($mY)
	Local $fMeanY = $fSumY / @extended

	; calculate differences Y - mean(Y)  --> Y
	__la_addScalar($mY, - $fMeanY, 0, 1, $iM)

	; TSS = square sum of diffs
	Local $fTSS = _lp_lassq($mY)

	; RSS = vᵀ * v
	Local $fRSS = $mLstSq.r2sum

	; calculate both coefficients of determination
	$mLstSq.R2      = 1 - $fRSS / $fTSS
	$mLstSq.R2_corr = 1 - (1 - $mLstSq.R2) * (($iM - 1) / ($iM - $iN - 1))

	; also calculate the pearson correlation between x and y
	$mLstSq.pearson = Sqrt($mLstSq.R2)

	$mLstSq.params = $aParams
	$mLstSq.vars   = $aVars
	Return SetExtended($i, $mLstSq)
EndFunc


; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __la_regression_LevenbergMarquardt()
; Description ...: calculates an n-dimensional linear or non-linear regression using levenberg-marquardt model
; Syntax ........: __la_regression_LevenbergMarquardt($sFunc, $mVars, $mY, [$mApprox = Default, [$fApproxDefault = 1.0, [$fLambda = 0.01, [$iMaxIterations = 50, [$fEpsTermination = 1e-7, [$sLstSqAlgo = "QR", [$iFlagsLstSq = 0, [$sDeriveMethod = "Central"]]]]]]]])
; Parameters ....: sFunc           - [String] function equation as AutoIt syntax string where (see examples):
;                                    names are parameters to be estimated
;                                    names beginning with "$" are independent function values
;                  mVars           - [Array/Map] Array: independent function values for 1D case; Map of Arrays: named independent function values especially for n-dimensional case
;                  mY              - [Array/Map] dependent function values as array or vector as a map
;                  mApprox         - [Map] (Default: Default)
;                                  ↳ approximate values for the model parameters in the case of non-linear correlation of the form
;                                    {"ParamName1": fValue1, "ParamName2": fValue2, ...}
;                  fApproxDefault  - [Float] (Default: 1.0)
;                                  ↳ default value for the output values of the parameters if mApprox is not set
;                  fLambda         - [Float] (Default: 0.01)
;                                  ↳ initial lambda value for the Levenberg-Marquardt model
;                  iMaxIterations  - [UInt] (Default: 50)
;                                  ↳ maximum number of iterations
;                  fEpsTermination - [Float] (Default: 1e-7)
;                                  ↳ threshold value between the results of two iterations from which equality is assumed (termination criterion)
;                  sLstSqAlgo      - [String] (Default: "QR") algorithm for solving the system
;                                  ↳ "QR": QR decomposition                    - good stability and performance
;                                    "SVD": singular value decomposition (SVD) - highest stability but lowest performance
;                                    "Cholesky": cholesky decomposition        - better stability and performance than QR but only for positive-definite matrices
;                  iFlagsLstSq     - [UInt] (Default: 0)
;                                  ↳ output components flags (see description in _la_lstsq() for details)
;                  sDeriveMethod   - [String] (Default: "Central") algorithm to determine the 1st derivative
;                                  ↳ "Central"  - central difference quotient
;                                    "Central2", "Central3", "Central4" - central difference quotient of higher error order
;                                    "Forward"  - forward difference quotient
;                                    "Backward" - backward difference quotient
;                                    "Ridder"   - Ridder`s Method
;                                    "Higham"   - Highams algorithm
; Return value ..: Success: [Map] results as returned by _la_lstsq() plus:
;                           { "R2":      coefficient of determination,
;                             "R2_corr": adjusted coefficient of determination,
;                             "pearson": Pearson correlation coefficient,
;                             "params":  array of parameter names in correct order
;                             "vars":    array of variable names in correct order
;                             "Func":    the resulting formula with determined parameters as a string }
;                  Failure: Null and set @error to:
;                           | 1: error during __la_regression_prepareApproxVals() (@extended: @error from __la_regression_prepareApproxVals())
;                           | 2: no parameters determined in sFunc (@extended: @error from StringRegExp)
;                           | 3: no variables determined in sFunc (@extended: @error from StringRegExp)
;                           | 4: error during __la_regression_prepareVars() (@extended: @error from __la_regression_prepareVars())
;                           | 5: different number of variables between sFunc and mVars (@extended: Ubound(aVars))
;                           | 6: regression did not converge until maximum number of iterations was reached (@extended: number of iterations )
;                           |1X: error X during _la_lstsq() (@extended: @extended from _la_lstsq())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: Yes
;                  Global $sFunc = "B0 + B1 * $x + B2 * $x^2"
;                  Global $aX[] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
;                  Global $aY[] = [41,58,64,85,80,87,85,79,60,54,48]
;                  Global $mRegression = __la_regression_LevenbergMarquardt($sFunc, $aX, $aY, Default, 1, 0.1, 50, 1e-7, "QR", $__LA_LSTSQ_SDX) ; SDX for .sdX
;                  ConsoleWrite(StringFormat("\n% 11s: %4.2f\n% 11s: %4.2f\n% 11s: %4.2f\n\n", "R²", $mRegression.R2, "R² adjusted", $mRegression.R2_corr, "Pearson", $mRegression.pearson))
;                  _la_display($mRegression.x, "iterations: " & @extended)
; ===============================================================================================================================
Func __la_regression_LevenbergMarquardt($sFunc, $mVars, $mY, $mApprox = Default, $fApproxDefault = 1.0, $fLambda = 0.01, $iMaxIterations = 50, $fEpsTermination = 1e-7, $sLstSqAlgo = "QR", $iFlagsLstSq = 0, $sDeriveMethod = "Central")
	Local $aParams, $mParams[], $mParams1[], _ ; model parameter
	      $mXd, $mX0, $tX0, $mX1, $tX1, _      ; model parameter vectors
          $aVars, _                            ; model variables
          $mA, $mLambdaI, $mJ, _               ; jacobian matrices
          $tA, _                               ; corresponding data structure
          $mYd, $mR, $mL, _                    ; observation vectors
          $fNorm, $fNormR_old, _               ; decision values for the Levenberg-Marquardt method
		  $mLstSq, _                           ; the least square solver return map
		  $i, $j, $k                           ; additional helper variables like indices etc.

	; adjust input parameters
	$sFunc = StringUpper($sFunc) ; because map keys are case sensitive
	; add vtpv to least square solutions because it is needed for R²
	$iFlagsLstSq = BitOR($iFlagsLstSq, $__LA_LSTSQ_R2Sum)

	; convert $mY to correct form
	If IsArray($mY) Or IsString($mY) Then $mY = _blas_fromArray($mY)

	; prepare the approximation values for further use
	If Not __la_regression_prepareApproxVals($mApprox) Then Return SetError(1, @error, Null)

	; extract unknown params and assign their approximation value
	$aParams = StringRegExp($sFunc, '(?ix)(\$\w+(*SKIP)(*FAIL)|\b[a-z]\w*\b(?!\())', 3)
	If UBound($aParams) < 1 Then Return SetError(2, @error, Null)
	$aParams = _ArrayUnique($aParams, 0, 0, 1, 0)
	For $sParam In $aParams
		; use given approximation value if exists or the default value instead
		$mParams[$sParam] = MapExists($mApprox, $sParam) ? $mApprox[$sParam] : $fApproxDefault
	Next

	; extract the variables
	$aVars = StringRegExp($sFunc, '\$\K\w+', 3)
	If UBound($aVars) < 1 Then Return SetError(3, @error, Null)
	$aVars = _ArrayUnique($aVars, 0, 0, 1, 0)

	; prepare $mVars for further use
	If Not __la_regression_prepareVars($mVars, $aVars[0]) Then Return SetError(4, @error, Null)

	; $mVars and detected variable number should correspond
	If UBound($mVars) <> UBound($aVars) Then Return SetError(5, UBound($aVars), Null)

	Local $iM = UBound($mVars[$aVars[0]]), _
	      $iN = UBound($aParams)

	; declare the model solution vector yd = f(vars, X0)
	$mYd = _blas_createVector($iM, "DOUBLE")

	; declare the jacobian Matrix A (m x n) for the model parameters
	$mA = _blas_createMatrix($iM, $iN, "DOUBLE")
	$tA = $mA.struct

	; declare the levenberg-marquardt extension for the jacobian:  A --> [[A], [sqrt(Lambda) * I]]
	$mLambdaI = _blas_createMatrix($iN, $iN, "DOUBLE")

	; create approximation parameter vector X0
	$tX0 = DllStructCreate(StringFormat("DOUBLE[%d]", $iN))
	$i = 1
	For $sParam In $aParams
		DllStructSetData($tX0, 1, $mParams[$sParam], $i)
		$i += 1
	Next
	$mX0 = _blas_fromStruct($tX0, $iN, 0, "DOUBLE")


	For $i = 1 To  $iMaxIterations - 1
		; retrieve the jacobian matrix A
		__la_Regression_getJacobiParams($sFunc, $mVars, $mParams, $tA, $iM, $iN, $sDeriveMethod)
		; append the lambda-part for doing the gradient method
		__blas_fillWithScalar($mLambdaI, Sqrt($fLambda), 0, $iN + 1)
		$mJ = _la_join($mA, $mLambdaI, "vertical")

		; derive the model solution vector yd to calculate the residuals to y
		__la_Regression_getYfromModel($sFunc, $mVars, $mParams, $mYd.struct, $iM)
		; shortened solution vector
		$mR = _la_sub($mY, $mYd)
		$fNormR_old = _la_norm($mR)
		; Levenberg-Marquardt extension to the residual vector  r --> [r, 0]
		$mL = _blas_createVector($iM + $iN, "DOUBLE")
		_blas_copy($mYd.struct, 0, 1, 0, 1, $iM, $mL.ptr, False, "DOUBLE")

		; least square solution of J * dX = l
		$mLstSq = _la_lstsq($mJ, $mL, Default, $sLstSqAlgo, $iFlagsLstSq)
		If @error Then Return SetError(@error + 10, 0, Null)
		$mXd = $mLstSq.x

		; if changes to the previous iteration are below the threshold value - stop
		If _la_norm($mXd) < $fEpsTermination Then
			$mLstSq.x = _la_add($mXd, $mX0, False)
			ExitLoop
		EndIf

		; improved parameter vector
		$mX1 = _la_add($mXd, $mX0, False)
		$tX1 = $mX1.struct

		; adjust the map of approximate values for the parameters
		$j = 1
		For $sParam In $aParams
			$mParams1[$sParam] = DllStructGetData($tX1, 1, $j)
			$j += 1
		Next

		; derive the model solution vector yd to calculate the residuals to y
		__la_Regression_getYfromModel($sFunc, $mVars, $mParams1, $mYd.struct, $iM)
		; shortened solution vector
		$mR = _la_sub($mY, $mYd)
		$fNorm = _la_norm($mR)

		If $fNorm < $fNormR_old Then
			$mParams = $mParams1
			$k = 1
			For $sParam In $aParams
				DllStructSetData($tX0, 1, $mParams[$sParam], $k)
				$k += 1
			Next

			$fNormR_old = $fNorm

			$fLambda /= 10

		Else
			$fLambda *= 10

		EndIf
	Next

; throw error if result does not converge
	If $i >= $iMaxIterations Then Return SetError(6, $i, Null)

	;------ calculate R²
	; calc RSS = square sum of each difference Y_i with mean(Y)

	; calculate mean(Y)
	Local $fSumY = _la_sum($mY)
	Local $fMeanY = $fSumY / @extended

	; calculate differences Y - mean(Y)  --> Y
	__la_addScalar($mY, - $fMeanY, 0, 1, $iM)

	; TSS = square sum of diffs
	Local $fTSS = _lp_lassq($mY)

	; RSS = vᵀ * v
	Local $fRSS = $mLstSq.r2sum

	; calculate both coefficients of determination
	$mLstSq.R2      = 1 - $fRSS / $fTSS
	$mLstSq.R2_corr = 1 - (1 - $mLstSq.R2) * (($iM - 1) / ($iM - $iN - 1))

	; also calculate the pearson correlation between x and y
	$mLstSq.pearson = Sqrt($mLstSq.R2)

	$mLstSq.params = $aParams
	$mLstSq.vars   = $aVars
	Return SetExtended($i, $mLstSq)

EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __la_Regression_getJacobiParams()
; Description ...: derive the Jacobian matrix for a regression
; Syntax ........: __la_Regression_getJacobiParams($sFunc, $mVars, $mParams, $tJacobian, $iM, $iN, [$sDeriveMethod = "Central"])
; Parameters ....: sFunc         - [String] function equation as AutoIt syntax string where (see examples):
;                                    names are parameters to be estimated
;                                    names beginning with "$" are independent function values
;                  mVars         - [Map] Map of Arrays: named independent function values
;                  mParams       - [Map] approximate values for the model parameters of the form
;                                  {"ParamName1": fValue1, "ParamName2": fValue2, ...}
;                  tJacobian     - [DllStruct] target area in the memory where the Jacobian matrix is to be written
;                  iM            - [UInt] number of observations (rows of matrix)
;                  iN            - [UInt] number of unknown parameters (columns of matrix)
;                  sDeriveMethod   - [String] (Default: "Central") algorithm to determine the 1st derivative
;                                  ↳ "Central"  - central difference quotient
;                                    "Central2", "Central3", "Central4" - central difference quotient of higher error order
;                                    "Forward"  - forward difference quotient
;                                    "Backward" - backward difference quotient
;                                    "Ridder"   - Ridder`s Method
;                                    "Higham"   - Highams algorithm
; Return value ..: iN (@extended: iM)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: Yes
;                  Global $sFunc = StringUpper("B0 + B1 * $X + B2 * $x^2")
;                  Global $aX[] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
;                  Global $mVars[]
;                  $mVars["X"] = $aX
;                  Global $mModelParams[]
;                  $mModelParams["B0"] = 42
;                  $mModelParams["B1"] = 17
;                  $mModelParams["B2"] = -1.7
;                  Global $iM = UBound($aX)
;                  Global $iN = UBound($mModelParams)
;                  Global $mJacobi = _blas_createMatrix($iM, $iN, "DOUBLE")
;                  __la_Regression_getJacobiParams($sFunc, $mVars, $mModelParams, $mJacobi.struct, $iM, $iN, "Higham")
;                  _la_display($mJacobi, "Jacobi")
; ===============================================================================================================================
Func __la_Regression_getJacobiParams($sFunc, ByRef Const $mVars, ByRef Const $mParams, $tJacobian, $iM, $iN, $sDeriveMethod = "Central")
	Local $sFuncOrig = $sFunc, $sFuncTmp, $sFuncTmp2, $sParam, $iIndex = 1, _
	      $aParams = MapKeys($mParams), _
		  $aVars   = MapKeys($mVars)

	; go through every parameter-variable combination
	For $j = 0 To $iN - 1
		$sParam = $aParams[$j]
		$sFunc = $sFuncOrig

		; replace parameters with their approximate value except the current one
		For $y = 0 To $iN - 1
			If $y = $j Then ContinueLoop ; leave current parameter as is
			$sFunc = StringRegExpReplace($sFunc, "\b\Q" & $aParams[$y] & "\E\b", " " & StringFormat("%.16g", $mParams[$aParams[$y]]) & " ")
		Next

		; go through all observations
		$sFuncTmp = $sFunc
		For $i = 0 To $iM - 1
			$sFuncTmp2 = $sFuncTmp

			; set the current value for each variable
			For $sVar In $aVars
				$sFuncTmp2 = StringReplace($sFuncTmp2, "$" & $sVar, StringFormat("%.16g", ($mVars[$sVar])[$i]), 0, 1)
			Next

			; form the 1st derivative at the approximate value of the current parameter (linearize the function)
			DllStructSetData($tJacobian, 1, __la_derivate1D($sFuncTmp2, $sParam, $mParams[$sParam], $sDeriveMethod), $iIndex)
			$iIndex += 1

		Next
	Next

	Return SetExtended($iN, $iM)
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __la_Regression_getYfromModel()
; Description ...: determines the vector of model results based on the approximated parameters
; Syntax ........: __la_Regression_getYfromModel($sFunc, $mVars, $mParams, $tYModel, $iM)
; Parameters ....: sFunc   - [String] function equation as AutoIt syntax string where (see examples):
;                                    names are parameters to be estimated
;                                    names beginning with "$" are independent function values
;                  mVars   - [Map] Map of Arrays: named independent function values
;                  mParams - [Map] approximate values for the model parameters of the form
;                            {"ParamName1": fValue1, "ParamName2": fValue2, ...}
;                  tYModel - [DllStruct] target area in the memory where the value should be written
;                  iM      - [UInt] number of observations (rows of matrix)
; Return value ..: iM
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: Yes
;                  Global $sFunc = StringUpper("B0 + B1 * $X + B2 * $x * $x")
;                  Global $aX[] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
;                  Global $mVars[]
;                  $mVars["X"] = $aX
;                  Global $mY = _blas_createVector(UBound($aX), "DOUBLE")
;                  Global $mModelParams[]
;                  $mModelParams["B0"] = 42
;                  $mModelParams["B1"] = 17
;                  $mModelParams["B2"] = -1.7
;                  Global $iM = UBound($aX)
;                  Global $mYModel = _blas_createVector($iM, "DOUBLE")
;                  __la_Regression_getYfromModel($sFunc, $mVars, $mModelParams, $mYModel.struct, UBound($aX))
;                  _la_display($mYModel, "Y")
; ===============================================================================================================================
Func __la_Regression_getYfromModel($sFunc, ByRef Const $mVars, ByRef Const $mParams, Const $tYModel, $iM)
	Local $sFuncTmp, $sParam, _
		  $aVars   = MapKeys($mVars)

	; replace parameters with their approximate value except the current one
	For $sParam In MapKeys($mParams)
		$sFunc = StringRegExpReplace($sFunc, "\b\Q" & $sParam & "\E\b", " " & StringFormat("%.16g", $mParams[$sParam]) & " ")
	Next

	; go through all observations

	For $i = 0 To $iM - 1
		$sFuncTmp = $sFunc

		; set the current value for each variable
		For $sVar In $aVars
			$sFuncTmp = StringReplace($sFuncTmp, "$" & $sVar, StringFormat("%.16g", ($mVars[$sVar])[$i]), 0, 1)
		Next

		; write the calculated value for y into the buffer
		DllStructSetData($tYModel, 1, Execute($sFuncTmp), $i + 1)
	Next

	Return $iM
EndFunc


; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __la_regression_prepareApproxVals()
; Description ...: prepares the list of parameters with their initial values
; Syntax ........: __la_regression_prepareApproxVals($mApprox)
; Parameters ....: mApprox - [Array/Map] parameter list
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: invalid value for mApprox
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......: approximation values can be passed either as a 2D array or map but for further use it has to be a map
; Related .......:
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __la_regression_prepareApproxVals(ByRef $mApprox)
	Local $mTmp[]

	If (UBound($mApprox, 0)) = 2 And (UBound($mApprox, 2) = 2) Then
		For $i = 0 To UBound($mApprox) - 1
			$mTmp[$mApprox[$i][0]] = $mApprox[$i][1]
		Next
		$mApprox = $mTmp

	ElseIf IsKeyword($mApprox) = 1 Then
		$mApprox = $mTmp

	ElseIf Not IsMap($mApprox) Then
		Return SetError(1, 0, False)

	EndIf

	Return True
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __la_regression_prepareVars()
; Description ...: prepares the list of variables
; Syntax ........: __la_regression_prepareVars($mVars, [$sFirstName = "X"])
; Parameters ....: mVars      - [Array/Map] variable list
;                  sFirstName - [String] (Default: "X")
;                             ↳ default name for the first variable
; Return value ..: Success: SetError(1,0, False)
;                  Failure: False and set @error to:
;                           | 1: invalid value for mVars
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......: variable values can be passed as BLAS-object, 1D-Array or already in target form
;                  for further use we need a map with arrays as values und variable names as keys
; Related .......:
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __la_regression_prepareVars(ByRef $mVars, $sFirstName = "X")
	Local $mTmp[]

	If IsMap($mVars) Then
		If MapExists($mVars, "ptr") Then
			$mTmp[$sFirstName] = _blas_toArray($mVars)
			$mVars = $mTmp
		; no else - $mVars seems to have the correct form in this case
		EndIf

	ElseIf IsArray($mVars) And UBound($mVars, 0) = 1 Then
		Local $mTmp[]
		$mTmp[$sFirstName] = $mVars
		$mVars = $mTmp

	Else
		Return SetError(1,0, False)

	EndIf

	Return True
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __la_derivate1D()
; Description ...: calculates the 1st derivative of a function
; Syntax ........: __la_derivate1D($sFunc, $sParam, $fValue, [$sMethod = "Central", [$fH = $f_LA_DBL_STEP, [$iIterationLimit = 10, [$fInitialH_Ridder = 1e-5]]]])
; Parameters ....: sFunc            - [String] function as a string in AutoIt syntax with one parameter (see example)
;                  sParam           - [String] the name of the parameter to which is to be derived
;                  fValue           - [Float] the parameter value at which to derive
;                  sMethod          - [String] (Default: "Central") algorithm to determine the 1st derivative
;                                   ↳ "Central"  - central difference quotient
;                                     "Central2", "Central3", "Central4" - central difference quotient of higher error order
;                                     "Forward"  - forward difference quotient
;                                     "Backward" - backward difference quotient
;                                     "Ridder"   - Ridder`s Method
;                                     "Higham"   - Highams algorithm
;                  fH               - [Float] (Default: $f_LA_DBL_STEP)
;                                   ↳ step size
;                  iIterationLimit  - [UInt] (Default: 10)
;                                   ↳ maximum number of iterations for "Ridder" and "Higham"
;                  fInitialH_Ridder - [Float] (Default: 1e-5)
;                                   ↳ initial step size for "Ridder"
; Return value ..: Success: the 1st derivative of the function sFunc to sParam at the location fValue
;                  Failure: Null and set @error to:
;                           | 1: invalid value for sMethod
;                           | 2: did not converge
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: Yes
;                  Global $sFunc = 'SQRT(( 67  -  46 )^2 + ( 37  -  78.3 )^2) - R'
;                  ConsoleWrite("Forward:  " & __la_derivate1D($sFunc, "R", 1, "Forward") & @CRLF & _
;                               "Backward: " & __la_derivate1D($sFunc, "R", 1, "Backward") & @CRLF & _
;                               "Central:  " & __la_derivate1D($sFunc, "R", 1, "Central") & @CRLF & _
;                               "Central2: " & __la_derivate1D($sFunc, "R", 1, "Central2") & @CRLF & _
;                               "Central3: " & __la_derivate1D($sFunc, "R", 1, "Central3") & @CRLF & _
;                               "Central4: " & __la_derivate1D($sFunc, "R", 1, "Central4") & @CRLF & _
;                               "Ridder:   " & __la_derivate1D($sFunc, "R", 1, "Ridder") & @CRLF & _
;                               "Higham:   " & __la_derivate1D($sFunc, "R", 1, "Higham") & @CRLF)
; ===============================================================================================================================
Func __la_derivate1D(Const $sFunc, Const $sParam, Const $fValue, $sMethod = "Central", $fH = $f_LA_DBL_STEP, $iIterationLimit = 10, $fInitialH_Ridder = 1e-5)
	;~ $fH = target accuracy if "Ridder"
	Local $fA, $fB, _
		  $fStep

	Switch $sMethod
		Case "Central", "Central1"
			$fA = Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue - $fH)))
			$fB = Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue + $fH)))
			Return ($fB - $fA) / (2 * $fH)

		Case "Central2"
			Return (Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue - 2 * $fH))) / 12 _
			- 2/3 * Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue - $fH))) _
			+ 2/3 * Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue + $fH))) _
			- Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue + 2 * $fH))) / 12) _
			/ $fH

		Case "Central3"
			Return ( _
			-     Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue - 3 * $fH))) / 60 _
			+ 3 * Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue - 2 * $fH))) / 20 _
			- 3 * Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue - $fH)))     / 4 _
			+ 3 * Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue + $fH)))     / 4 _
			- 3 * Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue + 2 * $fH))) / 20 _
			+     Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue + 3 * $fH))) / 60 _
			) / ($fH)

		Case "Central4"
			Return ( _
			      Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue - 4 * $fH))) / 280 _
			- 4 * Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue - 3 * $fH))) / 105 _
			+     Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue - 2 * $fH))) / 5 _
			- 4 * Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue - $fH)))     / 5 _
			+ 4 * Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue + $fH)))     / 5 _
			-     Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue + 2 * $fH))) / 5 _
			+ 4 * Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue + 3 * $fH))) / 105 _
			-     Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue + 4 * $fH))) / 280 _
			) / ($fH)

		Case "Forward"
			$fA = Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue)))
			$fB = Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue + $fH)))
			Return ($fB - $fA) / $fH

		Case "Backward"
			$fA = Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue - $fH)))
			$fB = Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue)))
			Return ($fB - $fA) / $fH

		Case "Ridder"
			Local $fD1, $fD2, $fD3, $fErr, $fErrOld = 1e100

			$fStep = $fInitialH_Ridder

			; initial difference quotient
			$fA = Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue + $fStep)))
			$fB = Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue - $fStep)))
			$fD1 = ($fA - $fB) / (2 * $fStep)

			For $i = 0 To $iIterationLimit
				; half step size
				$fStep /= 2

				; calculate new quotient with halfed step
				$fA = Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue + $fStep)))
				$fB = Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue - $fStep)))
				$fD2 = ($fA - $fB) / (2 * $fStep)

				; extrapolate the new value:
				$fD3 = (4 * $fD2 - $fD1) / 3

				; estimate the error
				$fErr = Abs($fD3 - $fD2) / 3  ; = |D3 - D2| / (2^O - 1)); "O" is typically 2 for the central difference

				; leave if not converge
				If $fErr > $fErrOld Then Return SetError(2, $fD1)
				$fErrOld = $fErr

				$fD1 = $fD3

				; leave if target accuracy is reached
				If $fErr < $fH Then Return $fD1
			Next

		Case "Higham"
			; initial parameters
			Local $fTargetError = $f_LA_DBL_EPS^(6/7) ; 6/7 correct digits for double type
			Local $fError = 1e100
			$fStep = $fH

			; first derivation
			$fA = Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue + $fStep)))
			$fB = Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue - $fStep)))
			Local $fD = ($fA - $fB) / (2 * $fStep)

			; iterative refinement
			Local $fDprev = $fD
			For $i = 0 To $iIterationLimit
				$fStep = $fStep / 2

				$fA = Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue + $fStep)))
				$fB = Execute(StringRegExpReplace($sFunc, '\b\Q' & $sParam & '\E\b', StringFormat("%.16g", $fValue - $fStep)))
				$fD = ($fA - $fB) / (2 * $fStep)

				$fError = Abs($fD - $fDprev)
				$fDprev = $fD

				If $fError > $fTargetError Then ExitLoop
			Next

			Return $fD

		Case Else
			Return SetError(1, 0, Null)

	EndSwitch

EndFunc
#EndRegion


#Region adjustment


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_adjustment()
; Description ...: performs a least-squares adjustment calculation for a system of different [weighted] non-linear equations
; Syntax ........: _la_adjustment($mObservations, [$mParams = Default, [$sSolutionMethod = Default, [$sNumSolutionMethode = "QR", [$iFlagsLstSq = 0, [$fEpsTermination = 1e-9, [$sDeriveMethod = "Central", [$fLambda = 1, [$iMaxIterations = 50, [$bIterativeVariance = False, [$sDataType = "DOUBLE"]]]]]]]]]])
; Parameters ....: mObservations       - [Map] set of observations as builded with _la_adj_addObservation()
;                  mParams             - [Map] (Default: Default)
;                                      ↳ Map of parameters with their upper-case names as key and their initial value as value
;                                        Default: parameters are determined using the observation equations and are automatically assigned the initial value 1.0
;                  sSolutionMethod     - [String] (Default: Default)
;                                      ↳ solution algorithm to solve the non-linear least squares problem:
;                                        "GN": Gauss-Newton algorithm - rapid convergence, but not as robust against poor initial values
;                                        "LM": Levenberg–Marquardt algorithm (Mixture of Gauss-Newton and gradient descent method)
;                                              - slightly longer convergence than GN but more robust against poor initial values
;                  sNumSolutionMethode - [String] (Default: "QR") algorithm for numerical solving the system
;                                      ↳ "QR": QR decomposition                    - good stability and performance
;                                        "SVD": singular value decomposition (SVD) - highest stability but lowest performance
;                                        "Cholesky": cholesky decomposition        - better stability and performance than QR but only for positive-definite matrices
;                  iFlagsLstSq         - [UInt] (Default: 0)
;                                      ↳ output components flags (see description in _la_lstsq() for details)
;                  fEpsTermination     - [Float] (Default: 1e-9)
;                                      ↳ threshold value between the results of two iterations from which equality is assumed (termination criterion)
;                  sDeriveMethod       - [String] (Default: "Central") algorithm to determine the 1st derivative
;                                      ↳ "Central":  central difference quotient
;                                        "Central2", "Central3", "Central4": central difference quotient of higher error order
;                                        "Forward":  forward difference quotient
;                                        "Backward": backward difference quotient
;                                        "Ridder":   Ridder`s Method
;                                        "Higham":   Highams algorithm
;                  fLambda             - [Float] (Default: 1)
;                                      ↳ initial lambda value for the Levenberg-Marquardt model
;                  iMaxIterations      - [UInt] (Default: 50)
;                                      ↳ maximum number of iterations
;                  bIterativeVariance  - [Bool] (Default: False)
;                                      ↳ If True: compare the a priori variance with the a posteriori variance and iterate until both are approximately equal.
;                                        In case that variance groups are specified in _la_adj_addObservation(), this is obsolete,
;                                        as a complete variance component estimate is performed.
;                  sDataType           - [String] (Default: "DOUBLE")
;                                      ↳ data type of the elements ("DOUBLE" or "FLOAT")
; Return value ..: Success: [Map] results of the solution depending on iFlagsLstSq {extended: number of iterations}:
;                           { "x":     solution vector x,
;                             "f":     degrees of freedom,
;                             "r":     Residuals r:  r = y - y_d,
;                             "r2sum": square sum of [weighted] residuals: r2sum = rᵀ * W * v,
;                             "s0":    a posteriori standard deviation factor,
;                             "Qx":    cofactor matrix of parameters,
;                             "Qy0":   a-priori cofactor matrix for the observations: Qᵧ = P⁻¹,
;                             "Qy":    a-posteriori cofactormatrix for the adjusted observations: Qᵧd = A * Qₓ * Aᵀ,
;                             "Qr":    cofactor matrix for the residuals r: Qᵣ = Qᵧ - Qᵧd,
;                             "R":     redundandy matrix R = Qᵧ * W,
;                             "sdX":   standard deviations of parameters: sdₓ = sqrt(s₀² * diag(Qₓ)),
;                             "sdY":   a posteriori standard deviations for the adjusted observations,
;                             "sdR":   standard deviations for the residuals: sdᵣ = sqrt(s₀² * diag(Qᵣ)),
;                             "cond":  condition number (only if "SVD" is used),
;                             "rank":  rank for the jacobian matrix (only if "SVD" is used)
;                           }
;                  Failure: Null and set @error to:
;                           | 1: error during __la_adj_getVarComponents() (@extended: @error from __la_adj_getVarComponents())
;                           |1X: error X during __la_adj_GaussNewton() (@extended: @extended from __la_adj_GaussNewton())
;                           |2X: error X during __la_adj_LevenbergMarquardt() (@extended: @extended from __la_adj_LevenbergMarquardt())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......: This function basically implements the Gauss-Markov model.
;                  Results for the Gauss-Helmert model or TLS solution can be achieved with the Gauss-Markov model through the clever use of pseudo-observations.
; Related .......:
; Link ..........:
; Example .......: example_niv_net.au3, example_PointCalculation.au3, example_tachynetz.au3, example_circle.au3, example_tls.au3
; ===============================================================================================================================
Func _la_adjustment($mObservations, $mParams = Default, $sSolutionMethod = Default, $sNumSolutionMethode = "QR", $iFlagsLstSq = 0, $fEpsTermination = 1e-9, $sDeriveMethod = "Central", $fLambda = 1, $iMaxIterations = 50, $bIterativeVariance = False, $sDataType = "DOUBLE")
	Local $sVarComp, $aObs, $aTmp, $fStdDev, $mLstSq, $j = 1, $s0

	If IsKeyword($mParams) = 1 Then $mParams = __la_adj_getParamList($mObservations)
	If IsKeyword($sSolutionMethod) = 1 Then $sSolutionMethod = "GN"

	Local $mVarComps = __la_adj_getVarComponents($mObservations)
	If @error Then Return SetError(1, @error, Null)

	If UBound($mVarComps) < 2 Then
	; adjustment with a single variance component for all observations
		For $j = 1 To 10
			If $sSolutionMethod = "GN" Then
				$mLstSq = __la_adj_GaussNewton($mObservations, $mParams, $sNumSolutionMethode, BitOr($iFlagsLstSq, $__LA_LSTSQ_R, $__LA_LSTSQ_REDUNDANCY), $fEpsTermination, $sDeriveMethod, $iMaxIterations, $sDataType)
				If @error Then Return SetError(10 + @error, @extended, Null)
			Else
				$mLstSq = __la_adj_LevenbergMarquardt($mObservations, $mParams, $fLambda, $sNumSolutionMethode, BitOr($iFlagsLstSq, $__LA_LSTSQ_R, $__LA_LSTSQ_REDUNDANCY), $fEpsTermination, $sDeriveMethod, $iMaxIterations, $sDataType)
				If @error Then Return SetError(20 + @error, @extended, Null)
			EndIf

			$s0 = $mLstSq.s0

			If (Not $bIterativeVariance) Or (Abs(1 - $mLstSq.s0) < 0.4) Then ExitLoop ; ToDo: let user choose the significant border

			For $iKey in MapKeys($mObservations)
				$aTmp = $mObservations[$iKey]

				$aTmp[5] *= $s0 ; adjust the standard deviation with the empirical s0
				$aTmp[2] = 1 / ($aTmp[5])^2

				$mObservations[$iKey] = $aTmp
			Next
		Next

	Else
	; adjustment with variance component estimation

		Local $mV, $tV, $mR, $tR, $i, $bTestPassed

		For $j = 1 To 10 ; ToDo: let the user decide which value

			If $sSolutionMethod = "GN" Then
				$mLstSq = __la_adj_GaussNewton($mObservations, $mParams, $sNumSolutionMethode, BitOr($iFlagsLstSq, $__LA_LSTSQ_R, $__LA_LSTSQ_REDUNDANCY), $fEpsTermination, $sDeriveMethod, $iMaxIterations, $sDataType)
				If @error Then Return SetError(10 + @error, @extended, Null)
			Else
				$mLstSq = __la_adj_LevenbergMarquardt($mObservations, $mParams, $fLambda, $sNumSolutionMethode, BitOr($iFlagsLstSq, $__LA_LSTSQ_R, $__LA_LSTSQ_REDUNDANCY), $fEpsTermination, $sDeriveMethod, $iMaxIterations, $sDataType)
				If @error Then Return SetError(20 + @error, @extended, Null)
			EndIf

			$mV = $mLstSq.r
			$tV = $mV.struct
			$mR = _la_getDiag($mLstSq.R)
			$tR = $mR.struct

			; set var comps to zero
			For $sVarComp In MapKeys($mVarComps)
				Dim $aTmp[3] = [($mVarComps[$sVarComp])[0], 0, 0]  ; [a-priori, vtpv, r]
				$mVarComps[$sVarComp] = $aTmp
			Next

			$i = 1
			For $aObs In $mObservations
				$sVarComp = $aObs[4]
				$fStdDev = $aObs[5]

				$aTmp =	$mVarComps[$sVarComp]

				; share of improvement square sum
				$aTmp[1] += DllStructGetData($tV, 1, $i)^2 * (1/($aTmp[0]^2 * $fStdDev^2))
				$aTmp[2] += DllStructGetData($tR, 1, $i)

				$mVarComps[$sVarComp] = $aTmp

				$i += 1
			Next

			$bTestPassed = True
			For $sVarComp In MapKeys($mVarComps)
				$aTmp = $mVarComps[$sVarComp]

				$aTmp[0] = Sqrt($aTmp[1] / $aTmp[2])
				$mVarComps[$sVarComp] = $aTmp

				If Abs(1 - $aTmp[0]) > 0.4 Then $bTestPassed = False ; ToDo: let user choose the significant border
			Next

			If $bTestPassed Then ExitLoop

			; adjust weights with variance components
			For $iObs In MapKeys($mObservations)
				$aObs = $mObservations[$iObs]
				$aTmp = $mVarComps[$aObs[4]]

				; calculate new weight
				$aObs[2] /= $aTmp[0]^2

				$mObservations[$iObs] = $aObs
			Next
		Next
	EndIf

	; add the results of the variance components to the result map
	$mLstSq.VarComps = $mVarComps

	Return SetExtended($j, $mLstSq)
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_adjustment_l1()
; Description ...: performs a adjustment calculation to L1 norm for a system of different [weighted] non-linear equations
; Syntax ........: _la_adjustment_l1($mObs, [$mParams = Default, [$sSolutionMethod = Default, [$sNumSolutionMethod = "QR", [$iFlagsLstSq = 0, [$fEpsTermination = 1e-9, [$fEpsL1Iteration = 0.001, [$sDeriveMethod = "Central", [$fLambda = 1, [$iMaxIterations = 50, [$bIterativeVariance = False, [$sDataType = "DOUBLE"]]]]]]]]]]])
; Parameters ....: mObs               - [Map] set of observations as builded with _la_adj_addObservation()
;                  mParams            - [Map] (Default: Default)
;                                     ↳ Map of parameters with their upper-case names as key and their initial value as value
;                                       Default: parameters are determined using the observation equations and are automatically assigned the initial value 1.0
;                  sSolutionMethod    - [String] (Default: Default)
;                                     ↳ solution algorithm to solve the non-linear least squares problem:
;                                       "GN": Gauss-Newton algorithm - rapid convergence, but not as robust against poor initial values
;                                       "LM": Levenberg–Marquardt algorithm (Mixture of Gauss-Newton and gradient descent method)
;                                             - slightly longer convergence than GN but more robust against poor initial values
;                  sNumSolutionMethod - [String] (Default: "QR") algorithm for numerical solving the system
;                                     ↳ "QR": QR decomposition                    - good stability and performance
;                                       "SVD": singular value decomposition (SVD) - highest stability but lowest performance
;                                       "Cholesky": cholesky decomposition        - better stability and performance than QR but only for positive-definite matrices
;                  iFlagsLstSq        - [UInt] (Default: 0)
;                                     ↳ output components flags (see description in _la_lstsq() for details)
;                  fEpsTermination    - [Float] (Default: 1e-9)
;                                     ↳ threshold value between the results of two iterations from which equality is assumed (termination criterion)
;                  fEpsL1Iteration    - [Float] (Default: 0.001)
;                                     ↳ Threshold value for whether the results of two iterations are considered equal
;                                       (termination  criterion for iterative reweighting to determine the L1 solution)
;                  sDeriveMethod      - [String] (Default: "Central") algorithm to determine the 1st derivative
;                                     ↳ "Central":  central difference quotient
;                                       "Central2", "Central3", "Central4": central difference quotient of higher error order
;                                       "Forward":  forward difference quotient
;                                       "Backward": backward difference quotient
;                                       "Ridder":   Ridder`s Method
;                                       "Higham":   Highams algorithm
;                  fLambda            - [Float] (Default: 1)
;                                     ↳ initial lambda value for the Levenberg-Marquardt model
;                  iMaxIterations     - [UInt] (Default: 50)
;                                     ↳ maximum number of iterations
;                  bIterativeVariance - [Bool] (Default: False)
;                                     ↳ If True: compare the a priori variance with the a posteriori variance and iterate until both are approximately equal.
;                                       In case that variance groups are specified in _la_adj_addObservation(), this is obsolete,
;                                       as a complete variance component estimate is performed.
;                  sDataType          - [String] (Default: "DOUBLE")
;                                     ↳ data type of the elements ("DOUBLE" or "FLOAT")
; Return value ..: Success: [Map] results of the solution depending on iFlagsLstSq {extended: number of iterations} (see _la_adjustment() for details)
;                  Failure: Null and set @error to:
;                           |1XX: error X during _la_adjustment() (@extended: @extended from _la_adjustment())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......: Roughly speaking, the L1 norm corresponds to a median solution (instead of mean for L2 norm).
;                  It is therefore robust against gross errors (up to 50% share) and can therefore be used for error detection.
;                  Computationally, it corresponds to the normal least squares method in which, however,
;                  the weights of the observations are adjusted iteratively in the direction of the L1 solution.
;                  This method can easily be adapted to other estimators (M-estimator, min-max, etc.)
; Related .......:
; Link ..........:
; Example .......: example_adjustment_l1.au3
; ===============================================================================================================================
Func _la_adjustment_l1($mObs, $mParams = Default, $sSolutionMethod = Default, $sNumSolutionMethod = "QR", $iFlagsLstSq = 0, $fEpsTermination = 1e-9, $fEpsL1Iteration = 0.001, $sDeriveMethod = "Central", $fLambda = 1, $iMaxIterations = 50, $bIterativeVariance = False, $sDataType = "DOUBLE")
	Local $mObsOrig = $mObs
	Local $mLstSq, $mLstSqOld, $aObsKeys, $aObs, $aV, $i, $j, $fEps
	Local $mX, $mX_old = _blas_createVector(UBound($mParams), $sDataType)

	$iFlagsLstSq = BitOr($iFlagsLstSq, $__LA_LSTSQ_R)

	For $i = 1 To $iMaxIterations
		$mLstSq = _la_adjustment($mObs, $mParams, $sSolutionMethod, $sNumSolutionMethod, $iFlagsLstSq, $fEpsTermination, $sDeriveMethod, $fLambda, $iMaxIterations, $bIterativeVariance, $sDataType)
		If @error Then Return SetError(100 + @error, @extended, $mLstSqOld)

		$mX = $mLstSq.x
		$fEps = _blas_nrm2(_la_sub($mX, $mX_old))
		If $fEps < $fEpsL1Iteration Then ExitLoop

		$mX_old = $mX
		$mObs = $mObsOrig

		$aV = _blas_toArray($mLstSq.r)
		$aObsKeys = MapKeys($mObs)

		For $j = 0 To UBound($aV) - 1
			$aObs = $mObs[$aObsKeys[$j]]

			$aObs[2] /= Abs($aV[$j])
			$aObs[5] = 1 / Sqrt($aObs[2])

			$mObs[$aObsKeys[$j]] = $aObs
		Next
		$mLstSqOld = $mLstSq
	Next

	Return SetExtended($i, $mLstSq)
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __la_adj_LevenbergMarquardt()
; Description ...: performs a adjustment calculation for a system of different [weighted] non-linear equations by using the Levenberg-Marquardt algorithm
; Syntax ........: __la_adj_LevenbergMarquardt($mObservations, $mParams, [$fLambda = 1, [$sLstSqAlgo = "QR", [$iFlagsLstSq = 0, [$fEpsTermination = 1e-7, [$sDeriveMethod = "Central", [$iMaxIterations = 50, [$sDataType = "DOUBLE"]]]]]]])
; Parameters ....: mObservations   - [Map] set of observations as builded with _la_adj_addObservation()
;                  mParams         - [Map] (Default: Default)
;                                  ↳ Map of parameters with their upper-case names as key and their initial value as value
;                                    Default: parameters are determined using the observation equations and are automatically assigned the initial value 1.0
;                  fLambda         - [Float] (Default: 1)
;                                  ↳ initial lambda value (damping parameter)
;                  sLstSqAlgo      - [String] (Default: "QR") algorithm for numerical solving the system
;                                  ↳ "QR": QR decomposition                    - good stability and performance
;                                    "SVD": singular value decomposition (SVD) - highest stability but lowest performance
;                                    "Cholesky": cholesky decomposition        - better stability and performance than QR but only for positive-definite matrices
;                  iFlagsLstSq     - [UInt] (Default: 0)
;                                     ↳ output components flags (see description in _la_lstsq() for details)
;                  fEpsTermination - [Float] (Default: 1e-7)
;                                  ↳ threshold value between the results of two iterations from which equality is assumed (termination criterion)
;                  sDeriveMethod   - [String] (Default: "Central") algorithm to determine the 1st derivative
;                                  ↳ "Central":  central difference quotient
;                                    "Central2", "Central3", "Central4": central difference quotient of higher error order
;                                    "Forward":  forward difference quotient
;                                    "Backward": backward difference quotient
;                                    "Ridder":   Ridder`s Method
;                                    "Higham":   Highams algorithm
;                  iMaxIterations  - [UInt] (Default: 50)
;                                  ↳ maximum number of iterations
;                  sDataType       - [String] (Default: "DOUBLE")
;                                  ↳ data type of the elements ("DOUBLE" or "FLOAT")
; Return value ..: Success: [Map] see _la_adjustment() for details
;                  Failure: Null and set @error to:
;                           | 1: error during _la_sqrtElements() (@extended: @error from _la_sqrtElements())
;                           | 2: error during _blas_sbmv() (@extended: @error from _blas_sbmv())
;                           | 3: error during _blas_trmm() (@extended: @error from _blas_trmm())
;                           | 4: error during _blas_sbmv() (@extended: @error from _blas_sbmv())
;                           | 5: max number of iterations reached (@extended: number of iterations)
;                           |1X: error X during _la_lstsq() (@extended: @extended from _la_lstsq())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __la_adj_LevenbergMarquardt($mObservations, $mParams, $fLambda = 1, $sLstSqAlgo = "QR", $iFlagsLstSq = 0, $fEpsTermination = 1e-7, $sDeriveMethod = "Central", $iMaxIterations = 50, $sDataType = "DOUBLE")
	Local $iM = UBound($mObservations), _
	      $iN = UBound($mParams), _
		  $i, $j, $k, _     ; indices
		  $mLstSq, $mXd, _
		  $fValue, $sParam, _
		  $mJ, $mL, $mX1, $tX1, $fNormR, $fNormR_old, $mParams1[]    ; extensions for Levenberg-Marquardt

	Local $aParams = MapKeys($mParams)

	; declare the observation vector y
	Local $mY = _blas_createVector($iM, $sDataType), $tY = $mY.struct
	; declare vector for the model predicted observation values y0
	Local $mY0 = _blas_createVector($iM, $sDataType), $tY0 = $mY0.struct
	; declare the jacobian matrix a for the parameters
	Local $mA = _blas_createMatrix($iM, $iN, $sDataType), $tA = $mA.struct
	; declare parameter approximation vector x0
	Local $mX0 = _blas_createVector($iN, $sDataType), $tX0 = $mX0.struct
	; declare the levenberg-marquardt extension for the jacobian:  A --> [[A], [sqrt(Lambda) * I]]
	Local $mLambdaI = _blas_createMatrix($iN, $iN, $sDataType)

	; derive weights vector (or default if unweighted)
	Local $mP = __la_adj_getWeights($mObservations, $sDataType)

	; build approximation parameter vector X0
	$i = 1
	For $fValue In $mParams
		DllStructSetData($tX0, 1, $fValue, $i)
		$i += 1
	Next

	; build observation vector y
	$i = 1
	For $aObs In $mObservations
		DllStructSetData($tY, 1, $aObs[1], $i)
		$i += 1
	Next

	; prepare weight vector P
	Local $mPsqrt, $pPsqrt, $tTmp, $pTmp, $mPmat
	If IsKeyword($mP) <> 1 Then
		$mPsqrt = _la_sqrtElements($mP, False)
		If @error Then Return SetError(1, @extended, Null)
		$pPsqrt = $mPsqrt.ptr

		; P-vector to P- diagonal matrix
		$mPmat = _la_VectorToDiag($mPsqrt, False)
	EndIf


	For $i = 1 To $iMaxIterations - 1

		; retrieve the jacobian matrix A for the parameters
		__la_adj_getJacobiParams($mObservations, $mParams, $tA, $sDeriveMethod)

		; derive model predicted observation values y0
		__la_adj_getYfromModel($mObservations, $mParams, $tY0)

		; calculate residual vector r = y - y0  --> y0
		_blas_scal($mY0.ptr, -1, 0, 1, $iM)
		_blas_axpy($mY.ptr, $mY0.ptr, 1, 0, 0, 1, 1, $iM)

		; weighted case (adjust A and r with sqrt(P))
		If IsKeyword($mP) <> 1 Then
			; P * r --> r
			$tTmp = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iM))
			$pTmp = DllStructGetPtr($tTmp)
			_blas_sbmv($pPsqrt, $mY0.ptr, $pTmp, 1, 1, 0, "L", $iM, 1, 1, 1)
			If @error Then Return SetError(2, @extended, Null)

			$mY0.struct = $tTmp
			$mY0.ptr = $pTmp
			$tY0 = $tTmp

			; P * A --> A
			_blas_trmm($mPmat.ptr, $mA.ptr, 1, "L", "U", "N", "N", $iM, $iN, $iM, $iM, $sDataType)
			If @error Then Return SetError(3, @extended, Null)
		EndIf

		$fNormR_old = _la_norm($mY0)

		; append the lambda-part for doing the gradient method
		__blas_fillWithScalar($mLambdaI, Sqrt($fLambda), 0, $iN + 1)
		$mJ = _la_join($mA, $mLambdaI, "vertical")

		; Levenberg-Marquardt extension to the residual vector  r --> [r, 0]
		$mL = _blas_createVector($iM + $iN, $sDataType)
		_blas_copy($tY0, 0, 1, 0, 1, $iM, $mL.ptr, False, $sDataType)

		; least square solution of A * dX = r
		$mLstSq = _la_lstsq($mJ, $mL, Default, $sLstSqAlgo, 0)
		If @error Then Return SetError(10 + @error, 0, Null)
		$mXd = $mLstSq.x

		; if changes to the previous iteration are below the threshold value - stop
		If _la_norm($mXd) < $fEpsTermination Then
			$mLstSq.x = _la_add($mXd, $mX0, False)
			ExitLoop
		EndIf

		; improved parameter vector
		$mX1 = _la_add($mXd, $mX0, False)
		$tX1 = $mX1.struct

		; adjust the map of approximate values for the parameters
		$j = 1
		For $sParam In $aParams
			$mParams1[$sParam] = DllStructGetData($tX1, 1, $j)
			$j += 1
		Next

		; derive model predicted observation values y0
		__la_adj_getYfromModel($mObservations, $mParams1, $tY0)

		; calculate residual vector r = y - y0  --> y0
		_blas_scal($mY0.ptr, -1, 0, 1, $iM)
		_blas_axpy($mY.ptr, $mY0.ptr, 1, 0, 0, 1, 1, $iM)

		; weighted case --> adjust r with P
		If IsKeyword($mP) <> 1 Then
			; P * r --> r
			$tTmp = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iM))
			$pTmp = DllStructGetPtr($tTmp)
			_blas_sbmv($pPsqrt, $mY0.ptr, $pTmp, 1, 1, 0, "L", $iM, 1, 1, 1)
			If @error Then Return SetError(4, @extended, Null)
			$mY0.struct = $tTmp
			$mY0.ptr = $pTmp
			$tY0 = $tTmp
		EndIf

		$fNormR = _la_norm($mY0)

		If $fNormR < $fNormR_old Then
			$mParams = $mParams1
			$k = 1
			For $sParam In $aParams
				DllStructSetData($tX0, 1, $mParams[$sParam], $k)
				$k += 1
			Next

			$fNormR_old = $fNormR
			$fLambda /= 10

		Else
			$fLambda *= 10

		EndIf
	Next

	; throw error if result does not converge
	If $i >= $iMaxIterations Then Return SetError(5, $i, Null)

	; last adjustment to get correct results without the enhancements for A and l (only necessary if vᵀPv, s0, sd_l, sdₓ, sd_v is needed)
	__la_adj_getJacobiParams($mObservations, $mParams, $tA, $sDeriveMethod)
	__la_adj_getYfromModel($mObservations, $mParams, $tY0)
	_blas_scal($mY0.ptr, -1, 0, 1, $iM)	; Y - Y0 --> Y0
	_blas_axpy($mY.ptr, $mY0.ptr, 1, 0, 0, 1, 1, $iM)
	$mLstSq = _la_lstsq($mA, $mY0, $mP, $sLstSqAlgo, $iFlagsLstSq)

	; improved parameter vector
	$mXd = $mLstSq.x
	_blas_axpy($mXd.ptr, $mX0.ptr, 1, 0, 0, 1, 1, $iN)
	$mLstSq.x = $mX0

	$mLstSq.params = $aParams

	; calculate the matrix of the redundancy shares
	If BitAND($iFlagsLstSq, $__LA_LSTSQ_REDUNDANCY) Then
		If IsKeyword($mP) <> 1 Then
			If BitAND($mP.storageType, $__g_BLAS_STYPE_MATRIX) Then
			; full P-Matrix
				Local $mC = _blas_createMatrix($iM, $iM, $sDataType) ; temporary matrix

				_blas_symm($mLstSq.Qr, $mP, $mC, 1, 0, "U", "L", $iM, $iM, $iM, $iM, $iM, $sDataType)
				$mLstSq.R = $mC ;_la_mul($mLstSq.Qr, $mP)
			Else
			; P-Vector = Diagonal Matrix

				; P-vector to P- diagonal matrix
				_la_VectorToDiag($mP, True)

				; Qv * P
				Local $mR = _blas_duplicate($mLstSq.Qr)
				_blas_trmm($mP.ptr, $mR.ptr, 1, "L", "U", "N", "N", $iM, $iM, $iM, $iM, $sDataType)
				$mLstSq.R = $mR

			EndIf
		Else
		; unweighted case
			$mLstSq.R = $mLstSq.Qr
		EndIf
	EndIf

	Return SetExtended($i, $mLstSq)
EndFunc


; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __la_adj_GaussNewton()
; Description ...: performs a adjustment calculation for a system of different [weighted] non-linear equations by using the Gauss-Newton algorithm
; Syntax ........: __la_adj_GaussNewton($mObservations, $mParams, [$sLstSqAlgo = "QR", [$iFlagsLstSq = 0, [$fEpsTermination = 1e-9, [$sDeriveMethod = "Central", [$iMaxIterations = 50, [$sDataType = "DOUBLE"]]]]]])
; Parameters ....: mObservations   - [Map] set of observations as builded with _la_adj_addObservation()
;                  mParams         - [Map] (Default: Default)
;                                  ↳ Map of parameters with their upper-case names as key and their initial value as value
;                                    Default: parameters are determined using the observation equations and are automatically assigned the initial value 1.0
;                  sLstSqAlgo      - [String] (Default: "QR") algorithm for numerical solving the system
;                                  ↳ "QR": QR decomposition                    - good stability and performance
;                                    "SVD": singular value decomposition (SVD) - highest stability but lowest performance
;                                    "Cholesky": cholesky decomposition        - better stability and performance than QR but only for positive-definite matrices
;                  iFlagsLstSq     - [UInt] (Default: 0)
;                                     ↳ output components flags (see description in _la_lstsq() for details)
;                  fEpsTermination - [Float] (Default: 1e-7)
;                                  ↳ threshold value between the results of two iterations from which equality is assumed (termination criterion)
;                  sDeriveMethod   - [String] (Default: "Central") algorithm to determine the 1st derivative
;                                  ↳ "Central":  central difference quotient
;                                    "Central2", "Central3", "Central4": central difference quotient of higher error order
;                                    "Forward":  forward difference quotient
;                                    "Backward": backward difference quotient
;                                    "Ridder":   Ridder`s Method
;                                    "Higham":   Highams algorithm
;                  iMaxIterations  - [UInt] (Default: 50)
;                                  ↳ maximum number of iterations
;                  sDataType       - [String] (Default: "DOUBLE")
;                                  ↳ data type of the elements ("DOUBLE" or "FLOAT")
; Return value ..: Success: [Map] see _la_adjustment() for details
;                  Failure: Null and set @error to:
;                           | 1: max number of iterations reached (@extended: number of iterations)
;                           |1X: error X during _la_lstsq() (@extended: @extended from _la_lstsq())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __la_adj_GaussNewton($mObservations, $mParams, $sLstSqAlgo = "QR", $iFlagsLstSq = 0, $fEpsTermination = 1e-9, $sDeriveMethod = "Central", $iMaxIterations = 50, $sDataType = "DOUBLE")
	Local $iM = UBound($mObservations), _
	      $iN = UBound($mParams), _
		  $i, $j, _     ; indices
		  $mLstSq, $mXd, _
		  $fValue, $sParam

	; declare the observation vector y
	Local $mY = _blas_createVector($iM, $sDataType), $tY = $mY.struct
	; declare vector for the model predicted observation values y0
	Local $mY0 = _blas_createVector($iM, $sDataType), $tY0 = $mY0.struct
	; declare the jacobian matrix a for the parameters
	Local $mA = _blas_createMatrix($iM, $iN, $sDataType), $tA = $mA.struct
	; declare parameter approximation vector x0
	Local $mX0 = _blas_createVector($iN, $sDataType), $tX0 = $mX0.struct

	; derive weights vector (or default if unweighted)
	Local $mP = __la_adj_getWeights($mObservations, $sDataType)

	; build approximation parameter vector X0
	$i = 1
	For $fValue In $mParams
		DllStructSetData($tX0, 1, $fValue, $i)
		$i += 1
	Next

	; build observation vector y
	$i = 1
	For $aObs In $mObservations
		DllStructSetData($tY, 1, $aObs[1], $i)
		$i += 1
	Next

	For $i = 1 To $iMaxIterations - 1
		; retrieve the jacobian matrix A for the parameters
		__la_adj_getJacobiParams($mObservations, $mParams, $tA, $sDeriveMethod)
		; derive model predicted observation values y0
		__la_adj_getYfromModel($mObservations, $mParams, $tY0)

		; calculate residual vector r = y - y0  --> y0
		_blas_scal($mY0.ptr, -1, 0, 1, $iM)
		_blas_axpy($mY.ptr, $mY0.ptr, 1, 0, 0, 1, 1, $iM)

		; least square solution of A * dX = r
		$mLstSq = _la_lstsq($mA, $mY0, $mP, $sLstSqAlgo, $iFlagsLstSq)
		If @error Then Return SetError(@error + 20, 0, Null)
		$mXd = $mLstSq.x

		; improved parameter vector
		_blas_axpy($mXd.ptr, $mX0.ptr, 1, 0, 0, 1, 1, $iN)
		$mLstSq.x = $mX0

		; if changes to the previous iteration are below the threshold value - stop
		If _la_norm($mXd) < $fEpsTermination Then ExitLoop

		; adjust the map of approximate values for the parameters
		$j = 1
		For $sParam In MapKeys($mParams)
			$mParams[$sParam] = DllStructGetData($tX0, 1, $j)
			$j += 1
		Next
	Next

	; throw error if result does not converge
	If $i >= $iMaxIterations Then Return SetError(1, $i, Null)

	$mLstSq.params = MapKeys($mParams)

	; calculate the matrix of the redundancy shares
	If BitAND($iFlagsLstSq, $__LA_LSTSQ_REDUNDANCY) Then
		If IsKeyword($mP) <> 1 Then
			If BitAND($mP.storageType, $__g_BLAS_STYPE_MATRIX) Then
			; full P-Matrix
				Local $mC = _blas_createMatrix($iM, $iM, $sDataType) ; temporary matrix

				_blas_symm($mLstSq.Qr, $mP, $mC, 1, 0, "U", "L", $iM, $iM, $iM, $iM, $iM, $sDataType)
				$mLstSq.R = $mC ;_la_mul($mLstSq.Qr, $mP)
			Else
			; P-Vector = Diagonal Matrix

				; P-vector to P- diagonal matrix
				_la_VectorToDiag($mP, True)

				; Qv * P
				Local $mR = _blas_duplicate($mLstSq.Qr)
				_blas_trmm($mP.ptr, $mR.ptr, 1, "L", "U", "N", "N", $iM, $iM, $iM, $iM, $sDataType)
				$mLstSq.R = $mR

			EndIf
		Else
		; unweighted case
			$mLstSq.R = $mLstSq.Qr
		EndIf
	EndIf

	Return SetExtended($i, $mLstSq)
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _la_adj_addObservation()
; Description ...: adds an observation to the adjustment system
; Syntax ........: _la_adj_addObservation($mObservations, $sFunc, $fValue, [$fStdDev = Default, [$sS0Component = "s0"]])
; Parameters ....: mObservations - [Map] Integer key map which holds the observation arrays. (An empty map on the 1st call)
;                  sFunc         - [String] Observation equation as a string in AutoIt syntax.
;                                  Parameters to be estimated are inserted as normal words. (best to look at the examples)
;                  fValue        - [Float] Measured value of the observation (observation value)
;                  fStdDev       - [Float] (Default: Default)
;                                ↳ Estimated standard deviation of the observation as a measure of the accuracy of the observation.
;                                 Is used to weight the observations against each other:
;                                 Wᵢ = 1 / sᵢ² --> so to specify a specific weight instead, enter 1 / sqrt(Wᵢ) instead
;                  sS0Component  - [String] (Default: "s0")
;                                ↳ Name of the observation group.
;                                 If a value is entered here, observations belonging to the same group are summarized
;                                 and a variance component estimate is calculated for these groups.
;                                 This allows different observation types to be examined and optimized with regard to their actual accuracy potential.
; Return value ..: -
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: example_niv_net.au3, example_PointCalculation.au3, example_tachynetz.au3, example_circle.au3, example_tls.au3 
; ===============================================================================================================================
Func _la_adj_addObservation(ByRef $mObservations, $sFunc, $fValue, $fStdDev = Default, $sS0Component = "s0")
	If Not IsMap($mObservations) Then
		Local $mTmpObs[]
		$mObservations = $mTmpObs
	EndIf

	$sFunc = StringUpper($sFunc) ; because Map-Keys are case sensitive

	; extract unknown params and assign their approximation value
	Local $aParams = StringRegExp($sFunc, '(\$\w+(*SKIP)(?!)|\b[A-Z]\w*\b(?!\h*\())', 3)
	If UBound($aParams) < 1 Or @error Then
		Dim $aParams[0]
	Else
		$aParams = _ArrayUnique($aParams, 0, 0, 1, 0)
	EndIf

	Local $fWeight = IsKeyword($fStdDev) = 1 ? Default : Number( 1 / ($fStdDev)^2, 3)

	Local $aObservation[6] = [$sFunc, Number($fValue, 3), $fWeight, $aParams, $sS0Component, $fStdDev] ; [function, value, weight, param list, variance component name, std dev.]
	MapAppend($mObservations, $aObservation)

EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __la_adj_getYfromModel()
; Description ...: determines the vector of model results based on the approximated parameters
; Syntax ........: __la_adj_getYfromModel($mObservations, $mParams, $tStruct)
; Parameters ....: mObservations - [Map] set of observations as builded with _la_adj_addObservation()
;                  mParams       - [Map] Map of parameters with their upper-case names as key and their initial value as value
;                  tStruct       - [DllStruct] target memory area to write the values
; Return value ..: -
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __la_adj_getYfromModel($mObservations, $mParams, ByRef $tStruct)
	Local $iIndex = 1, $sFunc, $aParams, $aObs, $sP

	For $aObs In $mObservations
		$aParams = $aObs[3]
		$sFunc = $aObs[0]

		; replace parameters with their approximate value
		For $sP In $aParams
			$sFunc = StringRegExpReplace($sFunc, "\b\Q" & $sP & "\E\b", " " & StringFormat("%.16g", $mParams[$sP]) & " ")
		Next

		; form the 1st derivative at the approximate value of the current parameter (linearize the function)
		DllStructSetData($tStruct, 1, Execute($sFunc), $iIndex)

		$iIndex += 1
	Next
EndFunc


; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __la_adj_getJacobiParams()
; Description ...: calculate the jacobian matrix out of the observation formulas and the approximated parameters
; Syntax ........: __la_adj_getJacobiParams($mObservations, $mParams, $tStruct, [$sDeriveMethod = "Central"])
; Parameters ....: mObservations - [Map] set of observations as builded with _la_adj_addObservation()
;                  mParams       - [Map] Map of parameters with their upper-case names as key and their initial value as value
;                  tStruct       - [DllStruct] target memory area to write the values
;                  sDeriveMethod - [String] (Default: "Central") algorithm to determine the 1st derivative
;                                ↳ "Central":  central difference quotient
;                                  "Central2", "Central3", "Central4": central difference quotient of higher error order
;                                  "Forward":  forward difference quotient
;                                  "Backward": backward difference quotient
;                                  "Ridder":   Ridder`s Method
;                                  "Higham":   Highams algorithm
; Return value ..: -
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __la_adj_getJacobiParams(ByRef $mObservations, ByRef $mParams, ByRef $tStruct, $sDeriveMethod = "Central")
	Local $iIndex = 1 , $sParam, $aParams, $sP, $sFunc

	For $sParam In MapKeys($mParams)

		For $aObs In $mObservations
			$aParams = $aObs[3]
			$sFunc = $aObs[0]

			; skip if current parameter is not used in this observation
			If _ArraySearch($aParams, $sParam) = -1 Then
				$iIndex += 1
				ContinueLoop
			EndIf

			; replace parameters with their approximate value except the current one
			For $sP In $aParams
				If $sParam = $sP Then ContinueLoop ; leave current parameter as is

				$sFunc = StringRegExpReplace($sFunc, "\b\Q" & $sP & "\E\b", " " & StringFormat("%.16g", $mParams[$sP]) & " ")
			Next

			; form the 1st derivative at the approximate value of the current parameter (linearize the function)
			DllStructSetData($tStruct, 1, __la_derivate1D($sFunc, $sParam, $mParams[$sParam], $sDeriveMethod), $iIndex)
			$iIndex += 1
		Next
	Next
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __la_adj_setApproxValue()
; Description ...: sets/changes the initial value of a parameter in the parameter list
; Syntax ........: __la_adj_setApproxValue($mParams, $sParam, $fValue)
; Parameters ....: mParams - [Map] Map of parameters with their upper-case names as key and their initial value as value
;                  sParam  - [String] identifier (the name) of the parameter in the list
;                  fValue  - [Float] new initial value
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: parameter does not exist in the parameter list
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __la_adj_setApproxValue(ByRef $mParams, $sParam, $fValue)
	$sParam = StringUpper($sParam)
	If Not MapExists($mParams, $sParam) Then Return SetError(1, 0, False)
	$mParams[$sParam] = Number($fValue, 3)
	Return True
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __la_adj_getWeights()
; Description ...: extracts the vector (or diagonal matrix) of the observation weights from the observations
; Syntax ........: __la_adj_getWeights($mObservations, [$sDataType = "DOUBLE"])
; Parameters ....: mObservations - [Map] set of observations as builded with _la_adj_addObservation()
;                  sDataType     - [String] (Default: "DOUBLE")
;                                ↳ data type of the elements ("DOUBLE" or "FLOAT")
; Return value ..: Success: Default: observations are completely unweighted
;                           Else: [Map] vector of weights
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __la_adj_getWeights(ByRef $mObservations, $sDataType = "DOUBLE")
	Local $bUnweighted = True
	Local $aObs, $i = 1
	Local $mP = _blas_createVector(UBound($mObservations, 1), $sDataType), $tP = $mP.struct

	For $aObs In $mObservations
		If IsKeyword($aObs[2]) = 1 Then
			DllStructSetData($tP, 1, 1, $i)
		Else
			DllStructSetData($tP, 1, $aObs[2], $i)
			$bUnweighted = False
		EndIf
		$i += 1
	Next

	Return $bUnweighted ? Default : $mP
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __la_adj_getParamList()
; Description ...: extract the parameters to be estimated from the observation equations into a map
; Syntax ........: __la_adj_getParamList($mObservations, [$fApproxDefault = 1])
; Parameters ....: mObservations  - [Map] set of observations as builded with _la_adj_addObservation()
;                  fApproxDefault - [Float] (Default: 1)
;                                 ↳ initial value for the parameter
; Return value ..: Success: [Map] {"param name": initial value, ....}
;                  Failure: Null and set @error to:
;                           | 1: invalid observation type detected (@extended: )
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __la_adj_getParamList(ByRef $mObservations, $fApproxDefault = 1.0)
	Local $mParams[], $aObs, $iObs, $i, $aParams

	For $iObs In MapKeys($mObservations)
		$aObs = $mObservations[$iObs]
		If UBound($aObs, 1) < 6 Or UBound($aObs, 0) <> 1 Then Return SetError(1, $iObs, Null)

		$aParams = $aObs[3]
		For $i = 0 To UBound($aParams) - 1
			If Not MapExists($mParams, $aParams[$i]) Then $mParams[$aParams[$i]] = $fApproxDefault
		Next
	Next

	Return $mParams

EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __la_adj_getVarComponents()
; Description ...: extracts a list of existing variance component groups in the system and calculates the weight of an observation
; Syntax ........: __la_adj_getVarComponents($mObservations)
; Parameters ....: mObservations - [Map] set of observations as builded with _la_adj_addObservation()
; Return value ..: [Map] list of variance component groups
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __la_adj_getVarComponents(ByRef $mObservations)
	Local $aObs, $iObs, $aTmp[3] = [1,0,0]

	Local $mVarComps[]
	For $iObs In MapKeys($mObservations)
		$aObs = $mObservations[$iObs]

		If Not MapExists($mVarComps, $aObs[4]) Then	$mVarComps[$aObs[4]] = $aTmp

		$aObs[2] = 1 / ($aObs[5]^2)

		$mObservations[$iObs] = $aObs
	Next

	Return $mVarComps
EndFunc

#EndRegion


#Region additional helper functions
; #FUNCTION# ====================================================================================================================
; Name ..........: _la_adj_showResult()
; Description ...: formats the results of _la_adj more clearly and display them in a window
; Syntax ........: _la_adj_showResult($mLstSq, $sWhat, [$sTitle = Default])
; Parameters ....: mLstSq - [Map] result of a adjustment as returned by _la_adjustment()
;                  sWhat  - [String] element of the result map, which should displayed
;                  sTitle - [String] (Default: Default)
;                         ↳ title of the windows
; Return value ..: Success: shows a window with the results
;                  Failure: Null and set @error to:
;                           | 1: sWhat does not exist in mLstSq
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......: combines solution values with their parameter names and (if exist) also show their standard deviation
; Related .......:
; Link ..........:
; Example .......: example_circle.au3
; ===============================================================================================================================
Func _la_adj_showResult($mLstSq, $sWhat, $sTitle = Default)
	If Not MapExists($mLstSq, $sWhat) Then Return SetError(1, 0, Null)
	Local $mData = $mLstSq[$sWhat], $iM = $mData.rows
	Local $aData = _blas_toArray($mData)

	Local $aParams = MapExists($mLstSq, "params") ? $mLstSq.params : Default
	Local $sHeader = ""
	Local $iNSigDec, $aStdDevs, $aShow

	If MapExists($mLstSq, "params") Then
	Switch $sWhat
		Case "x"
			Dim $aShow[$iM][(MapExists($mLstSq, "sdX") ? 3 : 2)]
			If MapExists($mLstSq, "sdX") Then
				$aStdDevs = _blas_toArray($mLstSq["sdX"])
				For $i = 0 To $iM -1
					$iNSigDec = __la_getSignificantDecimals($aStdDevs[$i]) + 1
					$aShow[$i][0] = $aParams[$i]
					$aShow[$i][1] = StringFormat("%." & $iNSigDec & "f", $aData[$i])
					$aShow[$i][2] = StringFormat("%." & $iNSigDec & "f", $aStdDevs[$i])
				Next
			Else
				For $i = 0 To $iM -1
					$aShow[$i][0] = $aParams[$i]
					$aShow[$i][1] = $aData[$i]
				Next
			EndIf
			$sHeader = "parameter|value|σ (std. dev.)"
			$sTitle = IsKeyword($sTitle) ? "solution vector" : $sTitle

		Case "r"
			Dim $aShow[$iM][(MapExists($mLstSq, "sdX") ? 2 : 1)]
			If MapExists($mLstSq, "sdR") Then
				$aStdDevs = _blas_toArray($mLstSq["sdR"])
			EndIf

			For $i = 0 To $iM -1
				$iNSigDec = MapExists($mLstSq, "sdR") ? __la_getSignificantDecimals($aStdDevs[$i]) + 1 : 0
				$aShow[$i][0] = StringFormat("%." & $iNSigDec & "f", $aData[$i])
				If MapExists($mLstSq, "sdR") Then $aShow[$i][1] = StringFormat("%." & $iNSigDec & "f", $aStdDevs[$i])
			Next
			$sHeader = "value|σ (std. dev.)"
			$sTitle = IsKeyword($sTitle) ? "residuals" : $sTitle

		Case "Qx"
			For $sP In $aParams
				$sHeader &= $sP & "|"
			Next
			$sHeader = StringTrimRight($sHeader, 1)
			$aShow = _blas_toArray($mLstSq["Qx"])
			$sTitle = IsKeyword($sTitle) ? "cofactor matrix Qx for parameters x" : $sTitle

		Case Else
			_la_display($mData, $sTitle)
			Return
	EndSwitch

	Else
		_la_display($mData, $sTitle)
		Return
	EndIf

	Return _ArrayDisplay($aShow, $sTitle , "", 64, "|", $sHeader)
EndFunc


; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __la_getSignificantDecimals()
; Description ...: determines the position of the 1st significant decimal place of a number
; Syntax ........: __la_getSignificantDecimals($fNumber)
; Parameters ....: fNumber - [Float] t number
; Return value ..: number of zeros before the first significant decimal place accures
; Author ........: AspirinJunkie
; Modified.......: 2024-09-05
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __la_getSignificantDecimals($fNumber)
	Local $aRegEx = StringRegExp(String($fNumber), '0\.(0*)', 1)
	If @error Then Return SetExtended(@error, 0)
	Return StringLen($aRegEx[0]) + 1
EndFunc

#EndRegion

