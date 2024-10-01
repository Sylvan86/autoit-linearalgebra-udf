;~ #AutoIt3Wrapper_Run_AU3Check=Y
;~ #AutoIt3Wrapper_Au3Check_Parameters=-d -w 1 -w 2 -w 3 -w 4 -w 5 -w 6 -w 7
;~ ;~ #AutoIt3Wrapper_AU3Check_Stop_OnWarning=Y
;~ Opt("MustDeclareVars", 1)

#include-once

#include "BLAS.au3"

; #INDEX# =======================================================================================================================
; Title .........: LAPACK
; AutoIt Version : 3.3.16.1
; Description ...: Wrapper UDF for interaction with LAPACK-compatible libraries (dlls)
;                  Provides basic functionalities for linear algebra.
;                  Is based on the BLAS-UDF
; Author(s) .....: AspirinJunkie
; Dll ...........: the user defined BLAS-DLL
; ===============================================================================================================================

; #CURRENT# =====================================================================================================================
; ---- Inverse ----
; _lp_getri  - computes the inverse of a matrix using the LU factorization computed by _lp_getrf()
; _lp_potri  - computes the inverse of a real symmetric positive definite matrix A using the Cholesky factorization computed by _lp_potrf()
;
; ---- factorization ----
; _lp_geqr   - computes a QR factorization of a real matrix A using the tiled (tall skinny) QR algorithm
; _lp_gemqr  - reconstructs the matrix Q from the results of _lp_geqr()
; _lp_geqrf  - computes a QR factorization of a real matrix A using the Householder QR algorithm
; _lp_orgqr  - reconstructs the matrix Q from the results of _lp_geqrf()
; _lp_potrf  - computes the Cholesky factorization of a real symmetric positive definite matrix
; _lp_pbtrf  - computes the Cholesky factorization of a packed real symmetric positive definite band matrix
; _lp_getrf  - computes an LU factorization of a general M-by-N matrix A
; _lp_gesvd  - computes the singular value decomposition (SVD) of a real M-by-N matrix A = U * SIGMA * Vᵀ
;
; ---- eigen values ----
; _lp_geev   - computes for an N-by-N real nonsymmetric matrix A, the eigenvalues and the left and/or right eigenvectors.
; _lp_syev   - computes all eigenvalues and eigenvectors of a real symmetric matrix A.
;
; ---- solve linear systems ----
; _lp_gesv   - computes the solution to a system of linear equations A * X = B where A is a general N × N matrix by using LU decomposition
; _lp_posv   - computes the solution to a system of linear equations A * X = B where A is an symmetric positive definite N × N matrix by using cholesky decomposition
; _lp_sysv   - computes the solution to a system of linear equations A * X = B where A is an symmetric N × N matrix by using LDL decomposition
; _lp_gbsv   - computes the solution to a system of linear equations A * X = B where A is a packed(!) band matrix by using LU decomposition
; _lp_gtsv   - computes the solution to a system of linear equations A * X = B where A is a tridiagonal matrix by gaussian elimination
; _lp_trtrs  - solves a triangular system of the form A * X = B or Aᵀ * X = B where A is a triangular matrix
; _lp_potrs  - computes the solution to a system of linear equations A * X = B out of the results of _lp_potrf() (solves Uᵀ * U * X = B Or L * Lᵀ * X = B)
;
; ---- least squares ----
; _lp_gels   - solves overdetermined or underdetermined linear system A * X = B using QR/LQ factorization
; _lp_getsls - solves overdetermined or underdetermined linear system A * X = B using tall-skinny/short-wide QR/LQ factorization
; _lp_gelss  - solves overdetermined or underdetermined linear system A * X = B using SVD factorization
; _lp_gelsy  - solves overdetermined or underdetermined linear system A * X = B using QR decomposition with column pivoting
; _lp_geqrs  - solves overdetermined or underdetermined linear system A * X = B using the results of the QR decomposition from _lp_geqrf()
;
; ---- auxiliary functions ----
; _lp_lassq  - calculate the sum of squares of elements of a matrix/vector
; _lp_lasrt  - sort the values of a matrix/vector
; _lp_rscl   - divide matrix/vector elements with a scalar a ( = 1/a * X)
; _lp_lange  - calculates the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of a general rectangular matrix
; _lp_gecon  - estimates the reciprocal of the condition number of a general real matrix A
; _lp_lacpy  - copies all or part of a two-dimensional matrix A to another matrix B
; _lp_laset  - initializes the off-diagonal elements and the diagonal elements of a matrix to given values
; _lp_laswp  - performs a series of row interchanges on a general rectangular matrix
; _lp_lapmt  - rearranges the columns of the M by N matrix X as specified by a permutation
; _lp_lauum  - computes the product U * Uᵀ or Lᵀ * L, where U or L are triangular matrices (chooses if blocked or unblocked variant is used)
; _lp_lauu2  - computes the product U * Uᵀ or Lᵀ * L, where U or L are triangular matrices (unblocked variant)
; _lp_lamch  - determines float/double precision machine parameters
; ===============================================================================================================================

#Region Inverse

; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_getri()
; Description ...: computes the inverse of a matrix using the LU factorization computed by _lp_getrf()
; Syntax ........: _lp_getri($mMatrix, $tIPIV, [$iN = Default, [$iLDA = Default, [$sDataType = "DOUBLE"]]])
; Parameters ....: mMatrix   - [Map] matrix as returned by _lp_getrf()
;                  tIPIV     - [DllStruct] pivot indices as returned by _lp_getrf()
;                  iN        - [Int] (Default: Default)
;                            ↳ number of rows of matrix A
;                  iLDA      - [Int] (Default: Default)
;                            ↳ leading dimension of the matrix A (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during first DllCall of getri (@extended: @error from DllCall)
;                           | 2: error inside first call of getri (@extended: INFO-value from getri)
;                           | 3: determined size of LWORK is not valid (@extended: determined LWORK)
;                           | 4: error during second DllCall of getri (@extended: @error from DllCall)
;                           | 5: error inside second call of getri (@extended: INFO-value from getri)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......: _lp_dgetri()
; Link ..........: https://www.netlib.org/lapack/explore-html/da/d28/group__getri_ga8b6904853957cfb37eba1e72a3f372d2.html#ga8b6904853957cfb37eba1e72a3f372d2
; Example .......: Yes
;                  Global $mA = _blas_fromArray('[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97],[93,89,13,-53,43]]')
;                  Global $tIPIV = _lp_getrf($mA) ; first: calculate LU decomposition
;                  _lp_getri($mA, $tIPIV)
;                  _blas_display($mA, "inverse(A)")
; ===============================================================================================================================
Func _lp_getri($mMatrix, $tIPIV, $iN = Default, $iLDA = Default, $sDataType = "DOUBLE")
	Local $pA ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrix)
			$sDataType = $mMatrix.datatype
			If IsKeyword($iN)   = 1 Then $iN   = $mMatrix.rows
			If IsKeyword($iLDA) = 1 Then $iLDA = $iN
			$pA = $mMatrix.ptr
		Case IsPtr($mMatrix)
			$pA = $mMatrix
		Case IsDllStruct($mMatrix)
			$pA = DllStructGetPtr($mMatrix)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; first run: determine optimal size of WORK
	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "getri", _
		"INT*",           $iN, _     ; N - order of A
		"PTR",            0, _       ; A
		"INT*",           $iLDA, _   ; lda
		"PTR",            0, _       ; IPIV
		$sDataType & "*", 0, _       ; WORK buffer (here 1 element because we only determine the size)
		"INT*",           -1, _      ; LWORK
		"INT*",           0 _        ; INFO
	)
	If @error Then Return SetError(1, @error, False)
	If $aDLL[7] <> 0 Then Return SetError(2, $aDLL[7], Null)

	; declare working buffers
	Local $iLWork = $aDLL[5]
	If $iLWork < 1 Then Return SetError(3, $iLWork, Null)
	Local $tWork = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iLWork))

	; second run: determine the inverse
	$aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "getri", _
		"INT*", $iN, _                      ; N - order of A
		"PTR",  $pA, _                      ; A
		"INT*", $iLDA, _                    ; lda
		"PTR",  DllStructGetPtr($tIPIV), _  ; IPIV
		"PTR",  DllStructGetPtr($tWork), _  ; WORK buffer
		"INT*", $iLWork, _                  ; LWORK
		"INT*", 0 _ 			            ; INFO
	)
	If @error Then Return SetError(4, @error, False)
	Return $aDLL[7] = 0 ? SetExtended($iLWork, True) : SetError(5, $aDLL[7], Null)
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_potri()
; Description ...: computes the inverse of a real symmetric positive definite matrix A using the Cholesky factorization computed by _lp_potrf()
; Syntax ........: _lp_potri($mMatrix, [$cUPLO = "L", [$iN = Default, [$iLDA = Default, [$sDataType = "DOUBLE"]]]])
; Parameters ....: mMatrix   - [Map] matrix as returned by _lp_potrf()
;                  cUPLO     - [Char] (Default: "L")
;                            ↳ "U": upper triangular parts of A are used
;                              "L": lower triangular parts of A are used
;                  iN        - [Int] (Default: Default)
;                            ↳ number of rows of matrix A
;                  iLDA      - [Int] (Default: Default)
;                            ↳ leading dimension of the matrix A (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of potri (@extended: @error from DllCall)
;                           | 2: error inside call of potri (@extended: INFO-value from potri)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......: _lp_dpotri()
; Link ..........: https://www.netlib.org/lapack/explore-html/d4/d12/group__potri_gaced925926290482532c433f73b30b9f1.html#gaced925926290482532c433f73b30b9f1
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[12.4,0,-0.7,-1,-0.7],[0,12.4,-0.7,-1,-0.7],[-0.7,-0.7,2.4,-1,0],[-1,-1,-1,4,-1],[-0.7,-0.7,0,-1,2.4]]")
;                  _lp_potrf($mA, "U") ; do a cholesky factorization
;                  _lp_potri($mA, "U") ; calc the inverse out of the cholesky result
;                  $mA.storageType = BitOR($mA.storageType, $__g_BLAS_STYPE_SYMMETRIC + $__g_BLAS_STYPE_UPPER) ; A is symmetric and only upper part is stored here
;                  _blas_display($mA, "Transposed A")
; ===============================================================================================================================
Func _lp_potri($mMatrix, $cUPLO = "L", $iN = Default, $iLDA = Default, $sDataType = "DOUBLE")
	Local $pA ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrix)
			$sDataType = $mMatrix.datatype
			If IsKeyword($iN)   = 1 Then $iN   = $mMatrix.rows
			If IsKeyword($iLDA) = 1 Then $iLDA = $iN
			$pA = $mMatrix.ptr
		Case IsPtr($mMatrix)
			$pA = $mMatrix
		Case IsDllStruct($mMatrix)
			$pA = DllStructGetPtr($mMatrix)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; set char buffers and input healing
	DllStructSetData($tBLASCHAR1, 1, $cUPLO = "L" ? "L" : "U")

	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "potri", _
		"PTR", $pBLASCHAR1, _               ; UPLO
		"INT*", $iN, _                      ; N - order of A
		"PTR",  $pA, _                      ; A
		"INT*", $iLDA, _                    ; lda
		"INT*", 0 _ 			            ; INFO
	)
	If @error Then Return SetError(1, @error, False)
	Return $aDLL[5] = 0 ? True : SetError(2, $aDLL[5], False)
EndFunc

#EndRegion

#Region factorization

; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_geqr()
; Description ...: computes a QR factorization of a real matrix A using the tiled (tall skinny) QR algorithm
; Syntax ........: _lp_geqr($mA, [$iM = Default, [$iN = Default, [$iLDA = $iM, [$sDataType = "DOUBLE"]]]])
; Parameters ....: mA        - [Map] matrix A as a map, DllStruct or pointer (will be overwritten)
;                            ↳ on exit, contain the upper triangular matrix R and in the lower parts, values to retrieve Q using _lp_gemqr()
;                  iM        - [Int] (Default: Default)
;                            ↳ number of rows of matrix A
;                  iN        - [Int] (Default: Default)
;                            ↳ number of columns of matrix A
;                  iLDA      - [Int] (Default: $iM)
;                            ↳ leading dimension of the matrix A (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: [DllStruct] dllstruct containing the pivot indices for reconstruct matrix Q from processed A using _lp_gemqr()
;                  Failure: Null and set @error to:
;                           | 1: error during first DllCall of geqr (@extended: @error from DllCall)
;                           | 2: error inside first call of geqr (@extended: INFO-value from geqr)
;                           | 3: determined size of WORK is not valid (@extended: determined LWORK)
;                           | 4: determined size of T is not valid (@extended: determined TSIZE)
;                           | 5: error during second DllCall of geqr (@extended: @error from DllCall)
;                           | 6: error inside second call of geqr (@extended: INFO-value from geqr)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......: implements the tiled QR algorithm, for Householder QR algorithm see _lp_geqrf()
; Related .......: _lp_gemqr()
; Link ..........: https://www.netlib.org/lapack/explore-html/dc/d28/group__geqr_ga2fc43ead3296c6f7bc7b8b8e53d1fedb.html#ga2fc43ead3296c6f7bc7b8b8e53d1fedb
; Example .......: Yes
;                  Global $mTest = _blas_fromArray("[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97],[93,89,13,-53,43]]")
;                  Global $tT = _lp_geqr($mTest)
;                  $mTest.storageType = BitOR($mTest.storageType, $__g_BLAS_STYPE_TRIANGLE + $__g_BLAS_STYPE_UPPER) ; A is triangular
;                  _blas_display($mTest, "R (upper right part)")
; ===============================================================================================================================
Func _lp_geqr(ByRef $mA, $iM = Default, $iN = Default, $iLDA = $iM, $sDataType = "DOUBLE")
	Local $pA ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iM)   = 1 Then $iM   = $mA.rows
			If IsKeyword($iN)   = 1 Then $iN   = $mA.cols
			If IsKeyword($iLDA) = 1 Then $iLDA = $iM
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; initialize T with 5 elements (minimum size - we determine now the new size)
	Local $tT = DllStructCreate(StringFormat("%s[%d]", $sDataType, 5))

	; first call to determine LWORK and TSIZE
	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "geqr", _
		"INT*",           $iM, _                   ; M
		"INT*",           $iN, _                   ; N
		"PTR",            0, _                     ; the general matrix A
		"INT*",           $iLDA, _                 ; lda
		"PTR",            DllStructGetPtr($tT), _  ; T (here 1 element because we first determine it`s best size)
		"INT*",           -1, _                    ; TSIZE (-1 to determine best size)
		$sDataType & "*", 0, _                     ; WORK (here 1 element because we determine the best size)
		"INT*",           -1, _                    ; LWORK (-1 to determine best size)
		"INT*",           0 _ 			           ; INFO
	)
	If @error Then Return SetError(1, @error, Null)
	If $aDLL[9] <> 0 Then Return SetError(2, $aDLL[9], Null)

	; declare WORK-Buffer
	Local $iLWork = $aDLL[7]
	If $iLWork < 1 Then Return SetError(3, $iLWork, Null)
	Local $tWork = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iLWork))

	; declare T-Buffer
	Local $dTSIZE = DllStructGetData($tT, 1, 1)
	If $dTSIZE < 1 Then Return SetError(4, $dTSIZE, Null)
	$tT = DllStructCreate(StringFormat("%s[%d]", $sDataType, $dTSIZE))

	; final run with optimal dimensioned buffers
	$aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "geqr", _
		"INT*", $iM, _                     ; M
		"INT*", $iN, _                     ; N
		"PTR",  $pA, _                     ; A
		"INT*", $iLDA, _                   ; lda
		"PTR",  DllStructGetPtr($tT), _    ; T
		"INT*", $dTSIZE, _                 ; TSIZE
		"PTR",  DllStructGetPtr($tWork), _ ; WORK
		"INT*", $iLWork, _                 ; LWORK
		"INT*", 0 _ 			           ; INFO
	)
	If @error Then Return SetError(5, @error, Null)
	Return $aDLL[9] = 0 ? SetExtended($dTSIZE, $tT) : SetError(6, $aDLL[9], Null)

EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_gemqr()
; Description ...: reconstructs the matrix Q from the results of _lp_geqr()
; Syntax ........: _lp_gemqr($mA, $tT, $mC, $iTSIZE, [$cSIDE = "L", [$cTRANS = "N", [$iM = Default, [$iN = Default, [$iK = Default, [$iLDA = Default, [$iLDC = Default, [$sDataType = "DOUBLE"]]]]]]]])
; Parameters ....: mA        - [Map] matrix A as processed by _lp_geqr()
;                  tT        - [DllStruct] pivot indices as returned by _lp_geqr()
;                  mC        - [Map] target Matrix C (m × n)
;                  iTSIZE    - [Int] size of tT
;                  cSIDE     - [Char] (Default: "L")
;                            ↳ "L": TRANS="N": Q * C else: C * Q
;                              "R": TRANS="N": Qᵀ * C else: C * Qᵀ
;                  cTRANS    - [Char] (Default: "N")
;                            ↳ "N":  Q
;                              "T":  Qᵀ
;                  iM        - [Int] (Default: Default)
;                            ↳ number of rows of matrix A
;                  iN        - [Int] (Default: Default)
;                            ↳ number of columns of matrix A
;                  iK        - [Int] (Default: Default)
;                            ↳ number of elementary reflectors
;                  iLDA      - [Int] (Default: Default)
;                            ↳ leading dimension of the matrix A
;                  iLDC      - [Int] (Default: Default)
;                            ↳ leading dimension of the matrix C
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True (@extended: size of WORK)
;                  Failure: False and set @error to:
;                           | 1: error during first DllCall of gemqr (@extended: @error from DllCall)
;                           | 2: error inside first call of gemqr (@extended: INFO-value from gemqr)
;                           | 3: determined size of LWORK is not valid (@extended: determined LWORK)
;                           | 4: error during second DllCall of gemqr (@extended: @error from DllCall)
;                           | 5: error inside second call of gemqr (@extended: INFO-value from gemqr)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......: _lp_geqr()
; Link ..........: https://www.netlib.org/lapack/explore-html/de/d55/group__gemqr_ga485cfbe680f2f9b17768dbf7255558f2.html#ga485cfbe680f2f9b17768dbf7255558f2
; Example .......: Yes
;                  Global $mTest = _blas_fromArray("[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97],[93,89,13,-53,43]]")
;                  Global $tT = _lp_geqr($mTest)
;                  Global $mQ = _blas_createMatrix($mTest.rows, $mTest.cols, $mTest.datatype) ; initialize Q (= C in gemqr) as identity matrix to calculate whole Q
;                  __blas_fillWithScalar($mQ, 1, 0, $mTest.rows + 1)
;                  _lp_gemqr($mTest, $tT, $mQ, @extended)
;                  $mTest.storageType = BitOR($mTest.storageType, $__g_BLAS_STYPE_TRIANGLE + $__g_BLAS_STYPE_UPPER)
;                  _blas_display($mQ, "Q")
;                  _blas_display($mTest, "R (upper right triangle)")
; ===============================================================================================================================
Func _lp_gemqr(ByRef $mA, ByRef $tT, ByRef $mC, Const $iTSIZE, $cSIDE = "L", $cTRANS = "N", $iM = Default, $iN = Default, $iK = Default, $iLDA = Default, $iLDC = Default, $sDataType = "DOUBLE")
	Local $pA, $pC ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iM)   = 1 Then $iM   = $mA.rows
			If IsKeyword($iN)   = 1 Then $iN   = $mA.cols

			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect
	Select
		Case IsMap($mC)
			$pC = $mC.ptr
		Case IsPtr($mC)
			$pC = $mC
		Case IsDllStruct($mC)
			$pC = DllStructGetPtr($mC)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; set char buffers and input healing
	DllStructSetData($tBLASCHAR1, 1, $cSIDE = "L" ? "L" : "R")
	DllStructSetData($tBLASCHAR2, 1, $cTRANS = "T" ? "T" : "N")

	; determine dimension parameters
	If IsKeyword($iLDA) = 1 Then $iLDA = $cSIDE = "L" ? $iM : $iN
	If IsKeyword($iLDC) = 1 Then $iLDC = $iM
	If IsKeyword($iK)   = 1 Then $iK   = $cSIDE = "L" ? $iM : $iN

	; first call to determine LWORK
	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "gemqr", _
		"PTR",            $pBLASCHAR1, _           ; SIDE
		"PTR",            $pBLASCHAR2, _           ; TRANS
		"INT*",           $iM, _                   ; M
		"INT*",           $iN, _                   ; N
		"INT*",           $iK, _                   ; K
		"PTR",            0, _                     ; the general matrix A
		"INT*",           $iLDA, _                 ; lda
		"PTR",            DllStructGetPtr($tT), _  ; T
		"INT*",           $iTSIZE, _               ; TSIZE
		"PTR",            0, _                     ; C
		"INT*",           $iLDC, _                 ; LDC
		$sDataType & "*", 0, _                     ; WORK (here 1 element because we determine the best size)
		"INT*",           -1, _                    ; LWORK
		"INT*",           0 _ 			           ; INFO
	)
	If @error Then Return SetError(1, @error, False)
	If $aDLL[14] <> 0 Then Return SetError(2, $aDLL[14], False)

	; declare working buffers
	Local $iLWork = $aDLL[12]
	If $iLWork < 1 Then Return SetError(3, $iLWork, False)
	Local $tWork = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iLWork))

	; second run: calculate Q
	$aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "gemqr", _
		"PTR",  $pBLASCHAR1, _               ; SIDE
		"PTR",  $pBLASCHAR2, _               ; TRANS
		"INT*", $iM, _                       ; M
		"INT*", $iN, _                       ; N
		"INT*", $iK, _                       ; K
		"PTR",  $pA, _                       ; the general matrix A
		"INT*", $iLDA, _                     ; lda
		"PTR",  DllStructGetPtr($tT), _      ; T
		"INT*", $iTSIZE, _                   ; TSIZE
		"PTR",  $pC, _                       ; C
		"INT*", $iLDC, _                     ; LDC
		"PTR",  DllStructGetPtr($tWork), _   ; WORK (here 1 element because we determine the best size)
		"INT*", $iLWork, _                   ; LWORK
		"INT*", 0 _ 			             ; INFO
	)
	If @error Then Return SetError(4, @error, False)
	Return $aDLL[14] = 0 ? SetExtended($iLWORK, True) : SetError(5, $aDLL[14], False)
EndFunc



; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_geqrf()
; Description ...: computes a QR factorization of a real matrix A using the Householder QR algorithm
; Syntax ........: _lp_geqrf($mA, [$iM = Default, [$iN = Default, [$iLDA = $iM, [$sDataType = "DOUBLE"]]]])
; Parameters ....: mA        - [Map] matrix A as a map, DllStruct or pointer (will be overwritten)
;                            ↳ on exit, contain the upper triangular matrix R and in the lower parts values to retrieve Q using _lp_orgqr()
;                  iM        - [Int] (Default: Default)
;                            ↳ number of rows of matrix A
;                  iN        - [Int] (Default: Default)
;                            ↳ number of columns of matrix A
;                  iLDA      - [Int] (Default: $iM)
;                            ↳ leading dimension of the matrix A (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: [DllStruct] dllstruct containing the elementary reflectors for reconstruct matrix Q from processed A using _lp_orgqr()
;                  Failure: Null and set @error to:
;                           | 1: error during first DllCall of geqrf (@extended: @error from DllCall)
;                           | 2: error inside first call of geqrf (@extended: INFO-value from geqrf)
;                           | 3: determined size of WORK is not valid (@extended: determined LWORK)
;                           | 4: error during second DllCall of geqrf (@extended: @error from DllCall)
;                           | 5: error inside second call of geqrf (@extended: INFO-value from geqrf)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......: _lp_orgqr()
; Link ..........: https://www.netlib.org/lapack/explore-html/d0/da1/group__geqrf_gade26961283814bb4e62183d9133d8bf5.html#gade26961283814bb4e62183d9133d8bf5
; Example .......: Yes
;                  Global $mTest = _blas_fromArray("[[4,-1,1],[-1,4,-2],[1,-2,5]]")
;                  Global $tTau = _lp_geqrf($mTest)
;                  _blas_display($mTest)
; ===============================================================================================================================
Func _lp_geqrf($mA, $iM = Default, $iN = Default, $iLDA = $iM, $sDataType = "DOUBLE")
	Local $pA ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iM)   = 1 Then $iM   = $mA.rows
			If IsKeyword($iN)   = 1 Then $iN   = $mA.cols
			If IsKeyword($iLDA) = 1 Then $iLDA = $iM
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; first call to determine LWORK
	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "geqrf", _
		"INT*",           $iM, _   ; M
		"INT*",           $iN, _   ; N
		"PTR",            0, _     ; the general matrix A
		"INT*",           $iLDA, _ ; lda
		"PTR",            0, _     ; TAU
		$sDataType & "*", 0, _     ; WORK (here 1 element because we determine the best size)
		"INT*",           -1, _    ; LWORK
		"INT*",           0 _      ; INFO
	)
	If @error Then Return SetError(1, @error, Null)
	If $aDLL[8] <> 0 Then Return SetError(2, $aDLL[8], Null)

	; declare working buffers
	Local $iLWork = $aDLL[6]
	If $iLWork < 1 Then Return SetError(3, $iLWork, Null)
	Local $tWork = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iLWork))
	Local $iK = $iM > $iN ? $iN : $iM
	Local $tTau = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iK))

	$aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "geqrf", _
		"INT*", $iM, _                     ; M
		"INT*", $iN, _                     ; N
		"PTR",  $pA, _                     ; the general matrix A
		"INT*", $iLDA, _                   ; lda
		"PTR",  DllStructGetPtr($tTau), _  ; Tau
		"PTR",  DllStructGetPtr($tWork), _ ; WORK
		"INT*", $iLWork, _                 ; LWORK
		"INT*", 0 _ 			           ; INFO
	)
	If @error Then Return SetError(4, @error, Null)
	Return $aDLL[8] = 0 ? SetExtended($iLWork, $tTau) : SetError(5, $aDLL[8], Null)

EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_orgqr()
; Description ...: reconstructs the matrix Q from the results of _lp_geqrf()
; Syntax ........: _lp_orgqr($mA, $tTau, [$iM = Default, [$iN = Default, [$iK = Default, [$iLDA = $iM, [$sDataType = "DOUBLE"]]]]])
; Parameters ....: mA        - [Map] matrix A as processed by _lp_geqrf()
;                  tTau      - [DllStruct] elementary reflectors as returned by _lp_geqrf()
;                  iM        - [Int] (Default: Default)
;                            ↳ number of rows of matrix Q
;                  iN        - [Int] (Default: Default)
;                            ↳ number of columns of matrix Q
;                  iK        - [Int] (Default: Default)
;                            ↳ number of elementary reflectors
;                  iLDA      - [Int] (Default: $iM)
;                            ↳ leading dimension of the matrix A (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True (@extended: size of WORK)
;                  Failure: False and set @error to:
;                           | 1: error during first DllCall of orgqr (@extended: @error from DllCall)
;                           | 2: error inside first call of orgqr (@extended: INFO-value from orgqr)
;                           | 3: determined size of WORK is not valid (@extended: determined LWORK)
;                           | 4: error during second DllCall of orgqr (@extended: @error from DllCall)
;                           | 5: error inside second call of orgqr (@extended: INFO-value from orgqr)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......: _lp_geqrf()
; Link ..........: https://www.netlib.org/lapack/explore-html/d4/dfc/group__ungqr_ga9f4abcfe1543a5d6d90fcd4dd21f12f0.html#ga9f4abcfe1543a5d6d90fcd4dd21f12f0
; Example .......: Yes
;                  Global $mTest = _blas_fromArray('[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97],[93,89,13,-53,43]]')
;                  Global $tTau = _lp_geqrf($mTest)
;                  _lp_orgqr($mTest, $tTau) ; calculate Q out of the householder reflectors
;                  _blas_display($mTest, "Q-Matrix of QR")
; ===============================================================================================================================
Func _lp_orgqr(ByRef $mA, ByRef $tTau, $iM = Default, $iN = Default, $iK = Default, $iLDA = $iM, $sDataType = "DOUBLE")
	Local $pA ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iM)   = 1 Then $iM   = $mA.rows
			If IsKeyword($iN)   = 1 Then $iN   = $mA.cols
			If IsKeyword($iLDA) = 1 Then $iLDA = $iM
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	If IsKeyword($iK) = 1 Then $iK = $iM > $iN ? $iN : $iM

	; first call to determine LWORK
	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "orgqr", _
		"INT*",           $iM, _   ; M
		"INT*",           $iN, _   ; N
		"INT*",           $iK, _   ; K (don't have to be real size of tau because dgeqr stores additional data inside T)
		"PTR",            0, _     ; A
		"INT*",           $iLDA, _ ; lda
		"PTR",            0, _     ; TAU
		$sDataType & "*", 0, _     ; WORK (here 1 element because we determine the best size)
		"INT*",           -1, _    ; LWORK
		"INT*",           0 _ 	   ; INFO
	)
	If @error Then Return SetError(1, @error, False)
	If $aDLL[9] <> 0 Then Return SetError(2, $aDLL[9], False)

	; declare working buffers
	Local $iLWork = $aDLL[7]
	If $iLWork < 1 Then Return SetError(3, $iLWork, False)
	Local $tWork = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iLWork))

	; second run: calculate Q
	$aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "orgqr", _
		"INT*", $iM, _                      ; M
		"INT*", $iN, _                      ; N
		"INT*", $iK, _                      ; K (don't have to be real size of tau because dgeqr stores additional data inside T)
		"PTR",  $pA, _                      ; A
		"INT*", $iLDA, _                    ; lda
		"PTR",  DllStructGetPtr($tTau), _   ; TAU
		"PTR",  DllStructGetPtr($tWork), _  ; WORK (here 1 element because we determine the best size)
		"INT*", $iLWork, _                  ; LWORK
		"INT*", 0 _ 			            ; INFO
	)
	If @error Then Return SetError(4, @error, False)
	Return $aDLL[9] = 0 ? SetExtended($iLWork, True) : SetError(5, $aDLL[9], False)
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_potrf()
; Description ...: computes the Cholesky factorization of a real symmetric positive definite matrix
; Syntax ........: _lp_potrf($mA, [$cUPLO = "L", [$iN = Default, [$iLDA = $iN, [$sDataType = "DOUBLE"]]]])
; Parameters ....: mA        - [Map] matrix A as a map, DllStruct or pointer (will be overwritten)
;                            ↳ on exit, contain the result U or L matrix (depends on cUPLO)
;                  cUPLO     - [Char] (Default: "L")
;                            ↳ "U": upper triangle of A is stored
;                              "L": lower triangle of A is stored
;                  iN        - [Int] (Default: Default)
;                            ↳ order of matrix A (rows)
;                  iLDA      - [Int] (Default: Default)
;                            ↳ leading dimension of the matrix A (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during first DllCall of potrf (@extended: @error from DllCall)
;                           | 2: error inside first call of potrf (@extended: INFO-value from potrf)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d2/d09/group__potrf_ga84e90859b02139934b166e579dd211d4.html#ga84e90859b02139934b166e579dd211d4
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[6,2,1],[2,5,2],[1,2,4]]")
;                  _lp_potrf($mA, "U")
;                  _blas_display($mA, "Result U")
; ===============================================================================================================================
Func _lp_potrf($mA, $cUPLO = "L", $iN = Default, $iLDA = $iN, $sDataType = "DOUBLE")
	Local $pA ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iN)   = 1 Then $iN   = $mA.rows
			If IsKeyword($iLDA) = 1 Then $iLDA = $iN
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; set char buffers and input healing
	DllStructSetData($tBLASCHAR1, 1, $cUPLO = "U" ? "U" : "L")

	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "potrf", _
		"PTR",  $pBLASCHAR1, _  ; UPLO
		"INT*", $iN, _          ; N
		"PTR",  $pA, _          ; A
		"INT*", $iLDA, _        ; lda
		"INT*", 0 _ 			; INFO
	)
	If @error Then Return SetError(1, @error, False)
	Return $aDLL[5] = 0 ? True : SetError(2, $aDLL[5], False)
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_pbtrf()
; Description ...: computes the Cholesky factorization of a packed real symmetric positive definite band matrix
; Syntax ........: _lp_pbtrf($mA, [$cUPLO = "L", [$iKD = 0, [$iN = Default, [$iLDA = $iN, [$sDataType = "DOUBLE"]]]]])
; Parameters ....: mA        - [Map] packed symmetric band matrix A as a map, DllStruct or pointer (will be overwritten)
;                            ↳ on exit, contain the result U or L matrix (depends on cUPLO)
;                  cUPLO     - [Char] (Default: "L")
;                            ↳ "U": upper triangle of A is stored
;                              "L": lower triangle of A is stored
;                  iKD       - [Int] (Default: 0)
;                            ↳ number of superdiagonals (UPLO="U") or subdiagonals (UPLO="L") of the matrix A
;                  iN        - [Int] (Default: Default)
;                            ↳ order of matrix A (rows)
;                  iLDA      - [Int] (Default: $iN)
;                            ↳ leading dimension of the matrix A (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during first DllCall of pbtrf (@extended: @error from DllCall)
;                           | 2: error inside first call of pbtrf (@extended: INFO-value from pbtrf)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/dc/d58/group__pbtrf_gafd2c3ba375ecee6dba1fba62696993e0.html#gafd2c3ba375ecee6dba1fba62696993e0
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[6,2,1],[2,5,2],[1,2,4]]", $__g_BLAS_STYPE_BAND + $__g_BLAS_STYPE_PACKED + $__g_BLAS_STYPE_SYMMETRIC + $__g_BLAS_STYPE_UPPER, "DOUBLE", 2)
;                  _lp_pbtrf($mA, "U", 2)
;                  _blas_display($mA, "Result U")
; ===============================================================================================================================
Func _lp_pbtrf($mA, $cUPLO = "L", $iKD = 0, $iN = Default, $iLDA = $iN, $sDataType = "DOUBLE")
	Local $pA ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iN)   = 1 Then $iN   = $mA.rows
			If IsKeyword($iLDA) = 1 Then $iLDA = $iKD + 1
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; set char buffers and input healing
	DllStructSetData($tBLASCHAR1, 1, $cUPLO = "U" ? "U" : "L")

	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "pbtrf", _
		"PTR",  $pBLASCHAR1, _  ; UPLO
		"INT*", $iN, _          ; N
		"INT*", $iKD, _         ; KD
		"PTR",  $pA, _          ; A
		"INT*", $iLDA, _        ; lda
		"INT*", 0 _ 			; INFO
	)
	If @error Then Return SetError(1, @error, False)
	Return $aDLL[5] = 0 ? True : SetError(2, $aDLL[5], False)
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_getrf()
; Description ...: computes an LU factorization of a general M-by-N matrix A
; Syntax ........: _lp_getrf($mA, [$iM = Default, [$iN = Default, [$iLDA = $iM, [$sDataType = "DOUBLE"]]]])
; Parameters ....: mA        - [Map] matrix A as a map, DllStruct or pointer (will be overwritten)
;                            ↳ on exit, contain the upper triangular matrix U and in the lower parts (without diagonal) the lower triangle matrix L
;                  iM        - [Int] (Default: Default)
;                            ↳ number of rows of matrix A
;                  iN        - [Int] (Default: Default)
;                            ↳ number of columns of matrix A
;                  iLDA      - [Int] (Default: $iM)
;                            ↳ leading dimension of the matrix A (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: Dllstruct containing the pivot indices for reconstruct the permutation matrix P using _lp_laswp() (@extended = iN)
;                  Failure: Null and set @error to:
;                           | 1: error during first DllCall of getrf (@extended: @error from DllCall)
;                           | 2: error inside first call of getrf (@extended: INFO-value from getrf)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......: - algorithm using partial pivoting with row interchanges
;                  - calculates the pivoted LU decomposition P * A = L * U, so not A = L * U like e.g. Maple
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/db/d04/group__getrf_gaea332d65e208d833716b405ea2a1ab69.html#gaea332d65e208d833716b405ea2a1ab69
; Example .......: Yes
;                  Global $mA = _blas_fromArray('[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97],[93,89,13,-53,43]]')
;                  Global $tIPIV = _lp_getrf($mA.ptr, 5, 5)
;                  ; create identity matrix
;                  Global $mP = _blas_createMatrix($mA.rows, $mA.cols, $mA.datatype)
;                  _lp_laset($mP, "U", 0, 1)
;                  ; tIPIV --> Permutation Matrix P
;                  _lp_laswp($mP, $tIPIV, 1, 5, -1)
;                  _blas_display($mA, "U/L Matrices")
;                  _blas_display($mP, "P")
; ===============================================================================================================================
Func _lp_getrf($mA, $iM = Default, $iN = Default, $iLDA = $iM, $sDataType = "DOUBLE")
	Local $pA ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iM)   = 1 Then $iM   = $mA.rows
			If IsKeyword($iN)   = 1 Then $iN   = $mA.cols
			If IsKeyword($iLDA) = 1 Then $iLDA = $iM
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"
	Local $tIPIV = DllStructCreate("INT[" & ($iM < $iN ? $iM : $iN) & "]")

	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "getrf", _
		"INT*", $iM, _                      ; number of rows in A
		"INT*", $iN, _                      ; number of columns A
		"ptr",  $pA, _                      ; the general matrix A
		"INT*", $iLDA, _                    ; lda
		"ptr",  DllStructGetPtr($tIPIV), _  ; pivot indices
		"INT*", 0 _ 			            ; INFO
	)
	If @error Then Return SetError(1, @error, Null)
	Return $aDLL[6] = 0 ? SetExtended($iN, $tIPIV) : SetError(2, $aDLL[6], Null)
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_gesvd()
; Description ...: computes the singular value decomposition (SVD) of a real M-by-N matrix A = U * SIGMA * Vᵀ
; Syntax ........: _lp_gesvd($mA, [$cTransposed = "N", [$iM = Default, [$iN = Default, [$iLDA = Default, [$sDataType = "DOUBLE"]]]]])
; Parameters ....: mA        - [Map] matrix A as a map, DllStruct or pointer (will be overwritten)
;                            ↳ on exit, contain the the first rows of Vᵀ
;                  iM        - [Int] (Default: Default)
;                            ↳ number of rows of matrix A
;                  iN        - [Int] (Default: Default)
;                            ↳ number of columns of matrix A
;                  iLDA      - [Int] (Default: $iM)
;                            ↳ leading dimension of the matrix A (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: [Map] {"S": S matrix, VT: Vᵀ matrix, "U": U matrix} {@extended = size of WORK}
;                  Failure: Null and set @error to:
;                           | 1: error during first DllCall of gesvd (@extended: @error from DllCall)
;                           | 2: error inside first call of gesvd (@extended: INFO-value from gesvd)
;                           | 3: determined size of WORK is not valid (@extended: determined LWORK)
;                           | 4: error during second DllCall of gesvd (@extended: @error from DllCall)
;                           | 5: error inside second call of gesvd (@extended: INFO-value from gesvd)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d1/d7f/group__gesvd_gac6bd5d4e645049e49bb70691180abf07.html#gac6bd5d4e645049e49bb70691180abf07
; Example .......: Yes
;                  Global $mA = _blas_fromArray('[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97],[93,89,13,-53,43]]')
;                  Global $mSVD = _lp_gesvd($mA)
;                  _blas_display($mSVD.S, "S")
;                  _blas_display($mSVD.VT, "Vᵀ")
;                  _blas_display($mSVD.U, "U")
; ===============================================================================================================================
Func _lp_gesvd($mA, $iM = Default, $iN = Default, $iLDA = Default, $sDataType = "DOUBLE")
	Local $pA ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iM)   = 1 Then $iM = $mA.rows
			If IsKeyword($iN)   = 1 Then $iN = $mA.cols
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	Local $iMinMN = $iM > $iN ? $iN : $iM

	; leading dimensions
	If IsKeyword($iLDA) = 1 Then $iLDA = $iM
	Local $iLDU  = $iM, _
	      $iLDVT = $iN

	; output buffers
	Local $tS  = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iMinMN))
	Local $tU  = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iLDU * $iM))
	Local $tVT = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iLDVT * $iN))

	; first call to determine LWORK
	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "gesvd", _
		"STR",            'A', _     ; jobu (A, S, O, N)
		"STR",            'A', _     ; jobvt (A, S, O, N)
		"INT*",           $iM, _     ; M
		"INT*",           $iN, _     ; N
		"PTR",            0, _       ; the general matrix A
		"INT*",           $iLDA, _   ; lda
		"PTR",            0, _       ; S (singular values of A)
		"PTR",            0, _       ; U
		"INT*",           $iLDU, _   ; ldu
		"PTR",            0, _       ; Vᵀ
		"INT*",           $iLDVT, _  ; ldvt
		$sDataType & "*", 0, _       ; WORK buffer (here = 1 element because we determine the best size)
		"INT*",           -1, _      ; LWORK
		"INT*",           0 _        ; INFO
	)
	If @error Then Return SetError(1, @error, Null)
	If $aDLL[14] <> 0 Then Return SetError(2, $aDLL[14], Null)

	; declare working buffers
	Local $iLWork = $aDLL[12]
	If $iLWork < 1 Then Return SetError(3, $iLWork, Null)
	Local $tWork = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iLWork))

	$aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "gesvd", _
		"STR", 'A', _                      ; jobu (A, S, O, N)
		"STR", 'A', _                      ; jobvt (A, S, O, N)
		"INT*", $iM, _                     ; M
		"INT*", $iN, _                     ; N
		"PTR",  $pA, _                     ; the general matrix A
		"INT*", $iLDA, _                   ; lda
		"PTR",  DllStructGetPtr($tS), _    ; S (singular values of A)
		"PTR",  DllStructGetPtr($tU), _    ; U
		"INT*", $iLDU, _                   ; ldu
		"PTR",  DllStructGetPtr($tVT), _   ; Vᵀ
		"INT*", $iLDVT, _                  ; ldvt
		"PTR",  DllStructGetPtr($tWork), _ ; WORK buffer
		"INT*", $iLWork, _                 ; LWORK
		"INT*", 0 _ 			           ; INFO
	)
	If @error Then Return SetError(4, @error, Null)
	If $aDLL[14] <> 0 Then Return SetError(5, $aDLL[14], Null)

	Local $mRet[]
	$mRet.U  = _blas_fromStruct($tU,  $iLDU,   $iM, $sDataType)
	$mRet.VT = _blas_fromStruct($tVT, $iLDVT,  $iN, $sDataType)
	$mRet.S  = _blas_fromStruct($tS,  $iMinMN, 0,   $sDataType)

	Return SetExtended($iLWork, $mRet)
EndFunc



#EndRegion


#Region eigen values

; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_geev()
; Description ...: computes for an N-by-N real nonsymmetric matrix A, the eigenvalues and the left and/or right eigenvectors.
; Syntax ........: _lp_geev($mA, [$cJOBVL = "N", [$cJOBVR = "N", [$iN = Default, [$iLDA = $iN, [$sDataType = "DOUBLE"]]]]])
; Parameters ....: mA        - [Map] matrix A as a map, DllStruct or pointer (will be overwritten)
;                  cJOBVL    - [Char] (Default: "N")
;                            ↳ "N": left eigenvectors of A are not computed
;                              "V": left eigenvectors of A are computed
;                  cJOBVR    - [Char] (Default: "N")
;                            ↳ "N": right eigenvectors of A are not computed
;                              "V": right eigenvectors of A are computed
;                  iN        - [Int] (Default: Default)
;                            ↳ order of matrix A (rows)
;                  iLDA      - [Int] (Default: $iM)
;                            ↳ leading dimension of the matrix A (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: [Map] {"R": real part, "I": imaginary part, "VL": left eigenvectors, "VR": right eigenvectors} {@extended = size of WORK}
;                  Failure: Null and set @error to:
;                           | 1: error during first DllCall of geev (@extended: @error from DllCall)
;                           | 2: error inside first call of geev (@extended: INFO-value from geev)
;                           | 3: determined size of WORK is not valid (@extended: determined LWORK)
;                           | 4: error during second DllCall of geev (@extended: @error from DllCall)
;                           | 5: error inside second call of geev (@extended: INFO-value from geev)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d4/d68/group__geev_ga7d8afe93d23c5862e238626905ee145e.html#ga7d8afe93d23c5862e238626905ee145e
; Example .......: Yes
;                  Global $mA = _blas_fromArray('[[611,196,-192,407,-8,-52,-49,29],[196,899,113,-192,-71,-43,-8,-44],[-192,113,899,196,61,49,8,52],[407,-192,196,611,8,44,59,-23],[-8,-71,61,8,411,-599,208,208],[-52,-43,49,44,-599,411,208,208],[-49,-8,8,59,208,208,99,-911],[29,-44,52,-23,208,208,-911,99]]')
;                  Global $mEigen = _lp_geev($mA, "V")
;                  _blas_display($mEigen.R, "real part of eigenvalues")
;                  _blas_display($mEigen.I, "imaginary part of eigenvalues")
;                  _blas_display($mEigen.VL, "left eigenvectors")
; ===============================================================================================================================
Func _lp_geev($mA, $cJOBVL = "N", $cJOBVR = "N", $iN = Default, $iLDA = $iN, $sDataType = "DOUBLE")
	Local $pA ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iN)   = 1 Then $iN   = $mA.rows
			If IsKeyword($iLDA) = 1 Then $iLDA = $iN
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; set char buffers and input healing:
	DllStructSetData($tBLASCHAR1, 1, $cJOBVL <> "N" ? "V" : $cJOBVL)
	DllStructSetData($tBLASCHAR2, 1, $cJOBVR <> "N" ? "V" : $cJOBVR)

	; result buffers
	Local $iLDVL = $cJOBVL = "N" ? 1 : $iN ; "V" as alternative
	Local $iLDVR = $cJOBVR = "N" ? 1 : $iN ; "V" as alternative
	Local $tWR = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iN))
	Local $tWI = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iN))
	Local $tVL = $cJOBVL = "N" ? 0 : DllStructCreate(StringFormat("%s[%d]", $sDataType, $iLDVL * $iN))
	Local $tVR = $cJOBVR = "N" ? 0 : DllStructCreate(StringFormat("%s[%d]", $sDataType, $iLDVR * $iN))

	; first call to determine LWORK
	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "geev", _
		"PTR",            $pBLASCHAR1, _  ; JOBVL
		"PTR",            $pBLASCHAR2, _  ; JOBVR
		"INT*",           $iN, _          ; N
		"PTR",            0, _            ; A
		"INT*",           $iLDA, _        ; lda
		"PTR",            0, _            ; WR
		"PTR",            0, _            ; WI
		"PTR",            0, _            ; VL
		"INT*",           $iLDVL, _       ; LDVL
		"PTR",            0, _            ; VR
		"INT*",           $iLDVR, _       ; LDVR
		$sDataType & "*", 0, _            ; WORK (here 1 element because we determine the best size)
		"INT*",           -1, _           ; LWORK
		"INT*",           0 _ 			  ; INFO
	)
	If @error Then Return SetError(1, @error, Null)
	If $aDLL[14] <> 0 Then Return SetError(2, $aDLL[14], Null)

	; declare working buffers
	Local $iLWork = $aDLL[12]
	If $iLWork < 1 Then Return SetError(3, $iLWork, Null)
	Local $tWork = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iLWork))

	$aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "geev", _
		"PTR",  $pBLASCHAR1, _                                ; JOBVL
		"PTR",  $pBLASCHAR2, _                                ; JOBVR
		"INT*", $iN, _                                        ; N
		"PTR",  $pA, _                                        ; the general matrix A
		"INT*", $iN, _                                        ; lda
		"PTR",  DllStructGetPtr($tWR), _                      ; WR (real part of eigenvalues)
		"PTR",  DllStructGetPtr($tWI), _                      ; WI (imaginary part of eigenvalues)
		"PTR",  $cJOBVL = "N" ? 0 : DllStructGetPtr($tVL), _  ; VL
		"INT*", $iLDVL, _                                     ; LDVL
		"PTR",  $cJOBVR = "N" ? 0 : DllStructGetPtr($tVR), _  ; VR
		"INT*", $iLDVR, _                                     ; LDVR
		"PTR",  DllStructGetPtr($tWork), _                    ; WORK
		"INT*", $iLWork, _                                    ; LWORK
		"INT*", 0 _ 			                              ; INFO
	)
	If @error Then Return SetError(4, @error, Null)
	If $aDLL[14] <> 0 Then Return SetError(5, $aDLL[14], Null)

	Local $mRet[]
	$mRet.R = _blas_fromStruct($tWR, $iN, 0, $sDataType)
	$mRet.I = _blas_fromStruct($tWI, $iN, 0, $sDataType)

	If $cJOBVL = "V" Then $mRet.VL = _blas_fromStruct($tVL, $iLDVL, $iN, $sDataType)
	If $cJOBVR = "V" Then $mRet.VR = _blas_fromStruct($tVR, $iLDVR, $iN, $sDataType)

	Return SetExtended($iLWork, $mRet)
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_syev()
; Description ...: computes all eigenvalues and eigenvectors of a real symmetric matrix A.
; Syntax ........: _lp_syev($mA, [$cJOBZ = "N", [$cUPLO = "U", [$iN = Default, [$iLDA = $iN, [$sDataType = "DOUBLE"]]]]])
; Parameters ....: mA        - [Map] symmetric matrix A as a map, DllStruct or pointer (will be overwritten)
;                              on exit, contains the orthonormal eigenvectors (if cJOBZ = "V")
;                  cJOBZ     - [Char] (Default: "N")
;                            ↳ "N": Compute eigenvalues only
;                              "V": Compute eigenvalues and eigenvectors
;                  cUPLO     - [Char] (Default: "U")
;                            ↳ "U": upper triangle of A is stored
;                              "L": lower triangle of A is stored
;                  iN        - [Int] (Default: Default)
;                            ↳ order of matrix A (rows)
;                  iLDA      - [Int] (Default: $iM)
;                            ↳ leading dimension of the matrix A (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: vector of eigenvalues
;                  Failure: Null and set @error to:
;                           | 1: error during first DllCall of syev (@extended: @error from DllCall)
;                           | 2: error inside first call of syev (@extended: INFO-value from syev)
;                           | 3: determined size of WORK is not valid (@extended: determined LWORK)
;                           | 4: error during second DllCall of syev (@extended: @error from DllCall)
;                           | 5: error inside second call of syev (@extended: INFO-value from syev)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d8/d1c/group__heev_ga8995c47a7578fef733189df3490258ff.html#ga8995c47a7578fef733189df3490258ff
; Example .......: Yes
;                  Global $mMatrix = _blas_fromArray('[[611,196,-192,407,-8,-52,-49,29],[196,899,113,-192,-71,-43,-8,-44],[-192,113,899,196,61,49,8,52],[407,-192,196,611,8,44,59,-23],[-8,-71,61,8,411,-599,208,208],[-52,-43,49,44,-599,411,208,208],[-49,-8,8,59,208,208,99,-911],[29,-44,52,-23,208,208,-911,99]]')
;                  Global $mEigenValues = _lp_syev($mMatrix, "V")
;                  _blas_display($mEigenValues, "eigen values")
;                  _blas_display($mMatrix, "eigen vectors")
; ===============================================================================================================================
Func _lp_syev($mA, $cJOBZ = "N", $cUPLO = "U", $iN = Default, $iLDA = $iN, $sDataType = "DOUBLE")
	Local $pA ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iN)   = 1 Then $iN   = $mA.rows
			If IsKeyword($iLDA) = 1 Then $iLDA = $iN
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; set char buffers and input healing:
	DllStructSetData($tBLASCHAR1, 1, $cJOBZ <> "N" ? "V" : $cJOBZ)
	DllStructSetData($tBLASCHAR2, 1, $cUPLO <> "U" ? "L" : $cUPLO)

	; first call to determine LWORK
	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "syev", _
		"PTR",            $pBLASCHAR1, _  ; JOBZ
		"PTR",            $pBLASCHAR2, _  ; UPLO
		"INT*",           $iN, _          ; N
		"PTR",            0, _            ; the symmetric matrix A
		"INT*",           $iLDA, _        ; lda
		"PTR",            0, _            ; W
		$sDataType & "*", 0, _            ; WORK (here 1 element because we determine the best size)
		"INT*",           -1, _           ; LWORK
		"INT*",           0 _ 			  ; INFO
	)
	If @error Then Return SetError(1, @error, Null)
	If $aDLL[9] <> 0 Then Return SetError(2, $aDLL[9], Null)

	; declare working buffers
	Local $iLWork = $aDLL[7]
	If $iLWork < 1 Then Return SetError(3, $iLWork, Null)
	Local $tWork = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iLWork))

	; result buffers
	Local $tW = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iN))

	$aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "syev", _
		"PTR",  $pBLASCHAR1, _             ; JOBZ
		"PTR",  $pBLASCHAR2, _             ; UPLO
		"INT*", $iN, _                     ; N
		"PTR",  $pA, _                     ; the general matrix A
		"INT*", $iN, _                     ; lda
		"PTR",  DllStructGetPtr($tW), _    ; W
		"PTR",  DllStructGetPtr($tWork), _ ; WORK
		"INT*", $iLWork, _                 ; LWORK
		"INT*", 0 _ 			           ; INFO
	)
	If @error Then Return SetError(4, @error, Null)
	Return $aDLL[9] = 0 ? SetExtended($iLWork, _blas_fromStruct($tW, $iN, 0, $sDataType)) : SetError(5, $aDLL[9], Null)
EndFunc

#EndRegion

#Region solve linear systems


; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_gesv()
; Description ...: computes the solution to a system of linear equations A * X = B
;                  where A is a general N × N matrix by using LU decomposition
; Syntax ........: _lp_gesv($mA, $mB, [$iNRHS = 1, [$iN = Default, [$iLDA = $iN, [$iLDB = $iN, [$sDataType = "DOUBLE"]]]]])
; Parameters ....: mA        - [Map] matrix A as a map, DllStruct or pointer (will be overwritten)
;                            ↳ on exit, contain the factors from LU decomposition
;                  mB        - [Map] vector/matrix N × NRHS A as a map, DllStruct or pointer (will be overwritten)
;                            ↳ on exit, contain the solution values X
;                  iNRHS     - [Int] (Default: 1)
;                            ↳ number of right hand sides, i.e., the number of columns of the matrix B (to solve multiple systems at once)
;                  iN        - [Int] (Default: Default)
;                            ↳ order of matrix A (rows = number of linear equations)
;                  iLDA      - [Int] (Default: $iM)
;                            ↳ leading dimension of the matrix A (rows)
;                  iLDB      - [Int] (Default: $iN)
;                            ↳ leading dimension of the matrix B (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: pivot indices as DllStruct vector (@extend = $N)
;                  Failure: Null and set @error to:
;                           | 1: error during first DllCall of gesv (@extended: @error from DllCall)
;                           | 2: error inside first call of gesv (@extended: INFO-value from gesv)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d8/da6/group__gesv_ga831ce6a40e7fd16295752d18aed2d541.html#ga831ce6a40e7fd16295752d18aed2d541
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[1,3,-2],[3,5,6],[2,4,3]]")
;                  Global $mB = _blas_fromArray("[5,7,8]")
;                  _lp_gesv($mA, $mB) ; should be: [-15, 8, 2]
;                  _blas_display($mB, "solution vector/matrix x")
; ===============================================================================================================================
Func _lp_gesv($mA, $mB, $iNRHS = 1, $iN = Default, $iLDA = $iN, $iLDB = $iN, $sDataType = "DOUBLE")
	Local $pA, $pB ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iN) = 1 Then $iN = $mA.rows
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect
	Select
		Case IsMap($mB)
			If IsKeyword($iN) = 1 Then $iN = $mB.rows
			$pB = $mB.ptr
		Case IsPtr($mB)
			$pB = $mB
		Case IsDllStruct($mB)
			$pB = DllStructGetPtr($mB)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	If IsKeyword($iLDA) = 1 Then $iLDA = $iN
	If IsKeyword($iLDB) = 1 Then $iLDB = $iN

	; buffer for pivot indices
	Local $tIPIV = DllStructCreate("INT[" & $iN & "]")

	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "gesv", _
		"INT*", $iN, _                     ; N
		"INT*", $iNRHS, _                  ; NRHS - 0: no solving - only LU factorization, >0: number of systems with different b`s to be solved (every b = columns in $mB)
		"PTR",  $pA, _                     ; A
		"INT*", $iLDA, _                   ; lda
		"PTR",  DllStructGetPtr($tIPIV), _ ; IPIV (pivot indices of permutation matrix P)
		"PTR",  $pB, _                     ; B
		"INT*", $iLDB, _                   ; LDB
		"INT*", 0 _ 			           ; INFO
	)
	If @error Then Return SetError(1, @error, Null)
	Return $aDLL[8] = 0 ? SetExtended($iN, $tIPIV) : SetError(2, $aDLL[8], Null)
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_posv()
; Description ...: computes the solution to a system of linear equations A * X = B
;                  where A is an symmetric positive definite N × N matrix by using cholesky decomposition
; Syntax ........: _lp_posv($mA, $mB, [$iNRHS = 1, [$cUPLO = "L", [$iN = Default, [$iLDA = $iN, [$iLDB = $iN, [$sDataType = "DOUBLE"]]]]]])
; Parameters ....: mA        - [Map] symmetric positive definite matrix A as a map, DllStruct or pointer (will be overwritten)
;                            ↳ on exit, contain the factors U or L from the cholesky factorization
;                  mB        - [Map] vector/matrix N × NRHS A as a map, DllStruct or pointer (will be overwritten)
;                            ↳ on exit, contain the solution values X
;                  iNRHS     - [Int] (Default: 1)
;                            ↳ number of right hand sides, i.e., the number of columns of the matrix B (to solve multiple systems at once)
;                  cUPLO     - [Char] (Default: "L")
;                            ↳ "U": upper triangle of A is stored
;                              "L": lower triangle of A is stored
;                  iN        - [Int] (Default: Default)
;                            ↳ order of matrix A (rows = number of linear equations)
;                  iLDA      - [Int] (Default: $iM)
;                            ↳ leading dimension of the matrix A (rows)
;                  iLDB      - [Int] (Default: $iN)
;                            ↳ leading dimension of the matrix B (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during first DllCall of posv (@extended: @error from DllCall)
;                           | 2: error inside first call of posv (@extended: INFO-value from posv)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/de/d6c/group__posv_ga4844053bd30fe88a17a7e08d93bbae4b.html#ga4844053bd30fe88a17a7e08d93bbae4b
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[6,2,1],[2,5,2],[1,2,4]]")
;                  Global $mB = _blas_fromArray("[9,5,5]")
;                  _lp_posv($mA, $mB) ; should be: [1.313, 0.133, 0.855]
;                  _blas_display($mB, "solution vector/matrix x")
; ===============================================================================================================================
Func _lp_posv($mA, $mB, $iNRHS = 1, $cUPLO = "L", $iN = Default, $iLDA = $iN, $iLDB = $iN, $sDataType = "DOUBLE")
	Local $pA, $pB ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iN) = 1 Then $iN = $mA.rows
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect
	Select
		Case IsMap($mB)
			If IsKeyword($iN) = 1 Then $iN = $mB.rows
			$pB = $mB.ptr
		Case IsPtr($mB)
			$pB = $mB
		Case IsDllStruct($mB)
			$pB = DllStructGetPtr($mB)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	If IsKeyword($iLDA) = 1 Then $iLDA = $iN
	If IsKeyword($iLDB) = 1 Then $iLDB = $iN

	; set char buffers and input healing
	DllStructSetData($tBLASCHAR1, 1, $cUPLO = "U" ? "U" : "L")

	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "posv", _
		"PTR",  $pBLASCHAR1, _ ; UPLO
		"INT*", $iN, _         ; N
		"INT*", $iNRHS, _      ; NRHS - 0: no solving - only LU factorization, >0: number of systems with different b to be solved (every b = columns in $mB)
		"PTR",  $pA, _         ; A
		"INT*", $iLDA, _       ; lda
		"PTR",  $pB, _         ; B
		"INT*", $iLDB, _       ; LDB
		"INT*", 0 _            ; INFO
	)
	If @error Then Return SetError(1, @error, False)
	Return $aDLL[8] = 0 ? True : SetError(2, $aDLL[8], False)
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_sysv()
; Description ...: computes the solution to a system of linear equations A * X = B
;                  where A is an symmetric N × N matrix by using LDL decomposition
; Syntax ........: _lp_sysv($mA, $mB, [$iNRHS = 1, [$cUPLO = "L", [$iN = Default, [$iLDA = $iN, [$iLDB = $iN, [$sDataType = Default]]]]]])
; Parameters ....: mA        - [Map] symmetric matrix A as a map, DllStruct or pointer (will be overwritten)
;                  mB        - [Map] vector/matrix N × NRHS A as a map, DllStruct or pointer (will be overwritten)
;                            ↳ on exit, contain the solution values X
;                  iNRHS     - [Int] (Default: 1)
;                            ↳ number of right hand sides, i.e., the number of columns of the matrix B (to solve multiple systems at once)
;                  cUPLO     - [Char] (Default: "L")
;                            ↳ "U": upper triangle of A is stored
;                              "L": lower triangle of A is stored
;                  iN        - [Int] (Default: Default)
;                            ↳ order of matrix A (rows = number of linear equations)
;                  iLDA      - [Int] (Default: $iM)
;                            ↳ leading dimension of the matrix A (rows)
;                  iLDB      - [Int] (Default: $iN)
;                            ↳ leading dimension of the matrix B (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: pivot indices as DllStruct vector (@extend = $N)
;                  Failure: Null and set @error to:
;                           | 1: error during first DllCall of sysv (@extended: @error from DllCall)
;                           | 2: error inside first call of sysv (@extended: INFO-value from sysv)
;                           | 3: determined size of WORK is not valid (@extended: determined LWORK)
;                           | 4: error during second DllCall of sysv (@extended: @error from DllCall)
;                           | 5: error inside second call of sysv (@extended: INFO-value from sysv)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d8/ddb/group__hesv_ga480a0e2d5cd7e1e85e674762ea2426f8.html#ga480a0e2d5cd7e1e85e674762ea2426f8
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[6,2,1],[2,5,2],[1,2,4]]")
;                  Global $mB = _blas_fromArray("[9,5,5]")
;                  _lp_sysv($mA, $mB) ; should be: [1.313, 0.133, 0.855]
;                  _blas_display($mB, "solution vector/matrix x")
; ===============================================================================================================================
Func _lp_sysv($mA, $mB, $iNRHS = Default, $cUPLO = "L", $iN = Default, $iLDA = $iN, $iLDB = $iN, $sDataType = Default)
	Local $pA, $pB ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iN) = 1 Then $iN = $mA.rows
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect
	Select
		Case IsMap($mB)
			If IsKeyword($iN)    = 1 Then $iN    = $mB.rows
			If IsKeyword($iNRHS) = 1 Then $iNRHS = $mB.cols
			$pB = $mB.ptr
		Case IsPtr($mB)
			$pB = $mB
		Case IsDllStruct($mB)
			$pB = DllStructGetPtr($mB)
	EndSelect

	If IsKeyword($sDataType) = 1 Then $sDataType = "DOUBLE"
	If IsKeyword($iLDA)      = 1 Then $iLDA = $iN
	If IsKeyword($iLDB)      = 1 Then $iLDB = $iN

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; set char buffers and input healing
	DllStructSetData($tBLASCHAR1, 1, $cUPLO = "U" ? "U" : "L")

	; first run: determine optimal size of WORK
	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "sysv", _
		"PTR",            $pBLASCHAR1, _ ; UPLO
		"INT*",           $iN, _         ; N
		"INT*",           $iNRHS, _      ; NRHS - 0: no solving - only LU factorization, >0: number of systems with different b to be solved (every b = columns in $mB)
		"PTR",            0, _           ; A
		"INT*",           $iLDA, _       ; lda
		"PTR",            0, _           ; IPIV
		"PTR",            0, _           ; B
		"INT*",           $iLDB, _       ; LDB
		$sDataType & "*", 0, _           ; WORK buffer (here 1 element because we only determine the size)
		"INT*",           -1, _          ; LWORK
		"INT*",           0 _            ; INFO
	)
	If @error Then Return SetError(1, @error, Null)
	If $aDLL[11] <> 0 Then Return SetError(2, $aDLL[11], Null)

	; declare working buffers
	Local $iLWork = $aDLL[9]
	If $iLWork < 1 Then Return SetError(3, $iLWork, Null)
	Local $tWork = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iLWork))

	; buffer for pivot indices
	Local $tIPIV = DllStructCreate("INT[" & $iN & "]")

	$aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "sysv", _
		"PTR",            $pBLASCHAR1, _              ; UPLO
		"INT*",           $iN, _                      ; N
		"INT*",           $iNRHS, _                   ; NRHS - 0: no solving - only LU factorization, >0: number of systems with different b to be solved (every b = columns in $mB)
		"PTR",            $pA, _                      ; A
		"INT*",           $iLDA, _                    ; lda
		"PTR",            DllStructGetPtr($tIPIV), _  ; IPIV
		"PTR",            $pB, _                      ; B
		"INT*",           $iLDB, _                    ; LDB
		$sDataType & "*", DllStructGetPtr($tWork), _  ; WORK buffer (here 1 element because we only determine the size)
		"INT*",           $iLWORK, _                  ; LWORK
		"INT*",           0 _                         ; INFO
	)
	If @error Then Return SetError(4, @error, Null)
	Return $aDLL[11] = 0 ? SetExtended($iN, $tIPIV) : SetError(5, $aDLL[11], Null)
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_gbsv()
; Description ...: computes the solution to a system of linear equations A * X = B
;                  where A is a packed(!) band matrix by using LU decomposition
; Syntax ........: _lp_gbsv($mAB, $mB, [$iKL = 0, [$iKU = 0, [$iN = Default, [$iNRHS = 1, [$iLDAB = 2 * $iKL + $iKU + 1, [$iLDB = Default, [$sDataType = Default]]]]]]])
; Parameters ....: mAB       - [Map] packed(!) band matrix A as a map, DllStruct or pointer (will be overwritten)
;                  mB        - [Map] vector/matrix N × NRHS b as a map, DllStruct or pointer (will be overwritten)
;                            ↳ on exit, contain the solution values X
;                  iKL       - [Int] (Default: 0)
;                            ↳ number of subdiagonals within the band of A
;                  iKU       - [Int] (Default: 0)
;                            ↳ umber of superdiagonals within the band of A
;                  iN        - [Int] (Default: Default)
;                            ↳ order of matrix A (rows = number of linear equations)
;                  iNRHS     - [Int] (Default: 1)
;                            ↳ number of right hand sides, i.e., the number of columns of the matrix B (to solve multiple systems at once)
;                  iLDAB     - [Int] (Default: 2 * $iKL + $iKU + 1)
;                            ↳ leading dimension of the array AB
;                  iLDB      - [Int] (Default: $iN)
;                            ↳ leading dimension of the matrix B (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during first DllCall of gbsv (@extended: @error from DllCall)
;                           | 2: error inside first call of gbsv (@extended: INFO-value from gbsv)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/db/df8/group__gbsv_gaff55317eb3aed2278a85919a488fec07.html#gaff55317eb3aed2278a85919a488fec07
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[1,0,0,0,0],[0,2,0,0,0],[0,0,3,0,0],[0,0,0,4,0],[0,0,0,0,5]]", $__g_BLAS_STYPE_BAND + $__g_BLAS_STYPE_PACKED + $__g_BLAS_STYPE_MATRIX, "DOUBLE", 0, 0)
;                  Global $mB = _blas_fromArray("[9,5,5,6,2]")
;                  _lp_gbsv($mA, $mB, 0, 0)
;                  _blas_display($mB, "solution vector/matrix x")
; ===============================================================================================================================
Func _lp_gbsv($mAB, $mB, $iKL = 0, $iKU = 0, $iN = Default, $iNRHS = 1, $iLDAB = 2 * $iKL + $iKU + 1, $iLDB = Default, $sDataType = Default)
	Local $pAB, $pB ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mAB)
			$sDataType = $mAB.datatype
			If IsKeyword($iN) = 1 Then $iN = $mAB.rows
			$pAB = $mAB.ptr
		Case IsPtr($mAB)
			$pAB = $mAB
		Case IsDllStruct($mAB)
			$pAB = DllStructGetPtr($mAB)
	EndSelect
	Select
		Case IsMap($mB)
			If IsKeyword($iN)    = 1 Then $iN    = $mB.rows
			If IsKeyword($iNRHS) = 1 Then $iNRHS = $mB.cols
			$pB = $mB.ptr
		Case IsPtr($mB)
			$pB = $mB
		Case IsDllStruct($mB)
			$pB = DllStructGetPtr($mB)
	EndSelect
	
	If IsKeyword($sDataType) = 1 Then $sDataType = "DOUBLE"
	If IsKeyword($iLDB)      = 1 Then $iLDB      = $iN
				
	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; buffer for pivot indices
	Local $tIPIV = DllStructCreate("INT[" & $iN & "]")

	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "gbsv", _
		"INT*", $iN, _                     ; N
		"INT*", $iKL, _                    ; KL
		"INT*", $iKU, _                    ; KU
		"INT*", $iNRHS, _                  ; NRHS
		"PTR",  $pAB, _                    ; A
		"INT*", $iLDAB, _                  ; LDAB
		"PTR",  DllStructGetPtr($tIPIV), _ ; IPIV
		"PTR",  $pB, _                     ; B
		"INT*", $iLDB, _                   ; LDB
		"INT*", 0 _                        ; INFO
	)
	If @error Then Return SetError(1, @error, False)
	Return $aDLL[10] = 0 ? True : SetError(2, $aDLL[10], False)
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_gtsv()
; Description ...: computes the solution to a system of linear equations A * X = B
;                  where A is a tridiagonal matrix by gaussian elimination
; Syntax ........: _lp_gtsv($mD, $mB, [$mDU = Default, [$mDL = Default, [$iNRHS = Default, [$iN = Default, [$iLDB = $iN, [$sDataType = "DOUBLE"]]]]]])
; Parameters ....: mD        - [Map] Either: vector holding the elements of the main diagonal as a map, DllStruct or pointer
;                                    Or: full matrix A as a map
;                  mB        - [Map] vector/matrix N × NRHS b as a map, DllStruct or pointer (will be overwritten)
;                            ↳ on exit, contain the solution values X
;                  mDU       - [Map] (Default: Default)
;                            ↳ vector holding the elements of the super diagonal as a map, DllStruct or pointer
;                              If mD = matrix: values are extracted from mD
;                  mDL       - [Map] (Default: Default)
;                            ↳ vector holding the elements of the supper diagonal as a map, DllStruct or pointer
;                              If mD = matrix: values are extracted from mD
;                  iNRHS     - [Int] (Default: Default)
;                            ↳ number of right hand sides, i.e., the number of columns of the matrix B (to solve multiple systems at once)
;                  iN        - [Int] (Default: Default)
;                            ↳ order of matrix A (rows of A or elements in vector mD)
;                  iLDB      - [Int] (Default: $iN)
;                            ↳ leading dimension of the matrix B (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during first DllCall of gtsv (@extended: @error from DllCall)
;                           | 2: error inside first call of gtsv (@extended: INFO-value from gtsv)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-19
; Remarks .......: to solve Aᵀ * X = B instead: interchanging the order of the arguments mDU and mDL.
; Related .......:
; Link ..........:
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[4,1,0],[1,4,1],[0,1,4]]")
;                  Global $mB = _blas_fromArray("[2.4, 10.2, -9.6]")
;                  _lp_gtsv($mA, $mB)
;                  _blas_display($mB, "solution")
; ===============================================================================================================================
Func _lp_gtsv($mD, $mB, $mDU = Default, $mDL = Default, $iNRHS = Default, $iN = Default, $iLDB = $iN, $sDataType = "DOUBLE")
	Local $pD, $pDL, $pDU, $pB ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mD)
			$sDataType = $mD.datatype
			If IsKeyword($iN) = 1 Then $iN = $mD.rows

			If $mD.storageType <> 0 Then ; $mD = matrix - extract D, DU, DL from the matrix
				$mDL = _blas_createVector($iN - 1, $sDataType)
				_blas_copy($mD, 1, $iN + 1, 0, 1, $iN - 1, $mDL, False, $sDataType)

				$mDU = _blas_createVector($iN - 1, $sDataType)
				_blas_copy($mD, $iN, $iN + 1, 0, 1, $iN - 1, $mDU, False, $sDataType)

				Local $mDtmp = _blas_createVector($iN, $sDataType)
				_blas_copy($mD, 0, $iN + 1, 0, 1, $iN, $mDtmp, False, $sDataType)
				$mD = $mDtmp
			EndIf
			$pD = $mD.ptr
		Case IsPtr($mD)
			$pD = $mD
		Case IsDllStruct($mD)
			$pD = DllStructGetPtr($mD)
	EndSelect
	Select
		Case IsMap($mB)
			If IsKeyword($iN) = 1 Then $iN = $mB.rows
			$pB = $mB.ptr
			If IsKeyword($iNRHS) = 1 And $mD.storageType <> 0 Then $iNRHS = $mD.cols
		Case IsPtr($mB)
			$pB = $mB
		Case IsDllStruct($mB)
			$pB = DllStructGetPtr($mB)
	EndSelect
	Select
		Case IsMap($mDU)
			$pDU = $mDU.ptr
		Case IsPtr($mDU)
			$pDU = $mDU
		Case IsDllStruct($mDU)
			$pDU = DllStructGetPtr($mDU)
	EndSelect
	Select
		Case IsMap($mDL)
			$pDL = $mDL.ptr
		Case IsPtr($mDL)
			$pDL = $mDL
		Case IsDllStruct($mDL)
			$pDL = DllStructGetPtr($mDL)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	If IsKeyword($iLDB) = 1 Then $iLDB = $iN
	If IsKeyword($iNRHS) = 1 Then $iNRHS = 1

	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "gtsv", _
		"INT*", $iN,    _ ; N
		"INT*", $iNRHS, _ ; NRHS
		"PTR",  $pDL,   _ ; DL
		"PTR",  $pD,    _ ; D
		"PTR",  $pDU,   _ ; DU
		"PTR",  $pB,    _ ; B
		"INT*", $iLDB,  _ ; LDB
		"INT*", 0       _ ; INFO
	)
	If @error Then Return SetError(1, @error, False)
	Return $aDLL[8] = 0 ? True : SetError(2, $aDLL[8], False)
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_trtrs()
; Description ...: solves a triangular system of the form A * X = B  or  Aᵀ * X = B
;                  where A is a triangular matrix
; Syntax ........: _lp_trtrs($mA, $mB, [$cUPLO = "U", [$cDIAG = "N", [$cTRANS = "N", [$iNRHS = Default, [$iN = Default, [$iLDA = Default, [$iLDB = Default, [$sDataType = "DOUBLE"]]]]]]]])
; Parameters ....: mA        - [Map] triangular matrix A as a map, DllStruct or pointer (will be overwritten)
;                  mB        - [Map] vector/matrix N × NRHS A as a map, DllStruct or pointer (will be overwritten)
;                            ↳ on exit, contain the solution values X
;                  cUPLO     - [Char] (Default: "U")
;                            ↳ "U": upper triangle of A is stored
;                              "L": lower triangle of A is stored
;                  cDIAG     - [Char] (Default: "N")
;                            ↳ "N": A is non-unit triangular
;                              "U": A is unit triangular (diagonal values = 1)
;                  cTRANS    - [Char] (Default: "N")
;                            ↳ "N": A  * X = B
;                              "T": Aᵀ * X = B
;                              "C": Aᴴ * X = B
;                  iNRHS     - [Int] (Default: Default)
;                            ↳ number of right hand sides, i.e., the number of columns of the matrix B (to solve multiple systems at once)
;                  iN        - [Int] (Default: Default)
;                            ↳ order of matrix A (rows = number of linear equations)
;                  iLDA      - [Int] (Default: Default)
;                            ↳ leading dimension of the matrix A (rows)
;                  iLDB      - [Int] (Default: $iN)
;                            ↳ leading dimension of the matrix B (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during first DllCall of trtrs (@extended: @error from DllCall)
;                           | 2: error inside first call of trtrs (@extended: INFO-value from trtrs)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d4/dc1/group__trtrs_gab0b6a7438a7eb98fe2ab28e6c4d84b21.html#gab0b6a7438a7eb98fe2ab28e6c4d84b21
; Example .......: Yes
;                  Global $mA = _blas_FromArray("[[4,2,1],[0,3,5],[0,0,2]]", $__g_BLAS_STYPE_TRIANGLE + $__g_BLAS_STYPE_UPPER) ; "upper triangle" not necessary - just little speedup
;                  Global $mX = _blas_FromArray("[[11,16,4],[5,9,12],[19,12,14]]")
;                  _lp_trtrs($mA, $mX)
;                  _blas_display($mX)
; ===============================================================================================================================
Func _lp_trtrs($mA, $mB, $cUPLO = "U", $cDIAG = "N", $cTRANS = "N", $iNRHS = Default, $iN = Default, $iLDA = Default, $iLDB = Default, $sDataType = "DOUBLE")
	Local $pA, $pB ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iN)   = 1 Then $iN   = $cTRANS = "N" ? $mA.rows : $mA.cols
			If IsKeyword($iLDA) = 1 Then $iLDA = $iN
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect
	Select
		Case IsMap($mB)
			$sDataType = $mB.datatype
			If IsKeyword($iNRHS) = 1 Then $iNRHS = $mB.cols
			If IsKeyword($iLDB) = 1 Then $iLDB = $iN
			$pB = $mB.ptr
		Case IsPtr($mB)
			$pB = $mB
		Case IsDllStruct($mB)
			$pB = DllStructGetPtr($mB)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; set char buffers and input healing
	DllStructSetData($tBLASCHAR1, 1, $cUPLO  = "U" ? "U" : "L")
	DllStructSetData($tBLASCHAR2, 1, $cTRANS = "N" ? "N" : ($cTRANS = "T" ? "T" :"C"))
	DllStructSetData($tBLASCHAR3, 1, $cDIAG  = "N" ? "N" : "U")

	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "trtrs", _
		"PTR",    $pBLASCHAR1, _     ; UPLO  -> 'U': A = upper triangular, 'L': A = lower triangular
		"PTR",    $pBLASCHAR2, _     ; TRANS -> 'N': x = A*x, 'T': x = Aᵀ*x, 'C': x = Aᵀ*x
		"PTR",    $pBLASCHAR3, _     ; DIAG  -> 'U': A = unit triangular, 'N': A = other elements on diagonal
		"INT*",   $iN, _             ; N
		"INT*",   $iNRHS, _          ; NRHS
		"ptr",    $pA, _             ; A
		"INT*",   $iLDA, _           ; LDA
		"ptr",    $pB, _             ; B
		"INT*",   $iLDB, _           ; LDB
		"INT*",   0 _                ; INFO
	)
	If @error Then Return SetError(1, @error, False)
	Return $aDLL[10] = 0 ? True : SetError(2, $aDLL[10], False)
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_potrs()
; Description ...: computes the solution to a system of linear equations A * X = B
;                  out of the results of _lp_potrf()
;                  (solves Uᵀ * U * X = B Or L * Lᵀ * X = B)
; Syntax ........: _lp_potrs($mA, $mB, [$cUPLO = "U", [$iNRHS = Default, [$iN = Default, [$iLDA = Default, [$iLDB = Default, [$sDataType = "DOUBLE"]]]]]])
; Parameters ....: mA        - [Map] matrix A as returned by _lp_potrf()
;                  mB        - [Map] vector/matrix N × NRHS A as a map, DllStruct or pointer (will be overwritten)
;                            ↳ on exit, contain the solution values X
;                  cUPLO     - [Char] (Default: "U")
;                            ↳ "U": upper triangle of A is stored
;                              "L": lower triangle of A is stored
;                  iNRHS     - [Int] (Default: Default)
;                            ↳ number of right hand sides, i.e., the number of columns of the matrix B (to solve multiple systems at once)
;                  iN        - [Int] (Default: Default)
;                            ↳ order of matrix A (rows = number of linear equations)
;                  iLDA      - [Int] (Default: Default)
;                            ↳ leading dimension of the matrix A (rows)
;                  iLDB      - [Int] (Default: $iN)
;                            ↳ leading dimension of the matrix B (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during first DllCall of potrs (@extended: @error from DllCall)
;                           | 2: error inside first call of potrs (@extended: INFO-value from potrs)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d3/dc8/group__potrs_ga70a04d13ff2123745a26b1e236212cf7.html#ga70a04d13ff2123745a26b1e236212cf7
; Example .......: Yes
;                  Global $mA = _blas_FromArray("[[4,2,1],[0,3,5],[0,0,2]]", $__g_BLAS_STYPE_TRIANGLE + $__g_BLAS_STYPE_UPPER) ; "upper triangle" not necessary - just little speedup
;                  Global $mX = _blas_FromArray("[11,16,4]")
;                  _lp_potrs($mA, $mX)
;                  _blas_display($mX)
; ===============================================================================================================================
Func _lp_potrs($mA, $mB, $cUPLO = "U", $iNRHS = Default, $iN = Default, $iLDA = Default, $iLDB = Default, $sDataType = "DOUBLE")
	Local $pA, $pB ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iN)   = 1 Then $iN   = $mA.rows
			If IsKeyword($iLDA) = 1 Then $iLDA = $iN
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect
	Select
		Case IsMap($mB)
			$sDataType = $mB.datatype
			If IsKeyword($iNRHS) = 1 Then $iNRHS = $mB.cols
			If IsKeyword($iLDB) = 1 Then $iLDB = $iN
			$pB = $mB.ptr
		Case IsPtr($mB)
			$pB = $mB
		Case IsDllStruct($mB)
			$pB = DllStructGetPtr($mB)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; set char buffers and input healing
	DllStructSetData($tBLASCHAR1, 1, $cUPLO  = "U" ? "U" : "L")

	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "potrs", _
		"PTR",    $pBLASCHAR1, _     ; UPLO  -> 'U': A = upper triangular, 'L': A = lower triangular
		"INT*",   $iN, _             ; N
		"INT*",   $iNRHS, _          ; NRHS
		"ptr",    $pA, _             ; A
		"INT*",   $iLDA, _           ; LDA
		"ptr",    $pB, _             ; B
		"INT*",   $iLDB, _           ; LDB
		"INT*",   0 _                ; INFO
	)
	If @error Then Return SetError(1, @error, False)
	Return $aDLL[8] = 0 ? True : SetError(2, $aDLL[8], False)
EndFunc

#EndRegion


#Region least squares

; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_gels()
; Description ...: solves overdetermined or underdetermined linear system A * X = B
;                  using QR/LQ factorization
; Syntax ........: _lp_gels($mA, $mB, [$cTRANS = "N", [$iNRHS = Default, [$iM = Default, [$iN = Default, [$iLDA = $iM, [$iLDB = $iN, [$sDataType = "DOUBLE"]]]]]]])
; Parameters ....: mA        - [Map] matrix A (M × N) as a map, DllStruct or pointer (will be overwritten)
;                              on exit, contain the result of _lp_geqrf() (useful for deriving the cofactor matrix of the parameters)
;                  mB        - [Map] vector/matrix N × NRHS A as a map, DllStruct or pointer (will be overwritten)
;                            ↳ on exit, contain the solution values X in the first N elements and residual sum vᵀv in the last elements
;                  cTRANS    - [Char] (Default: "N")
;                            ↳ "N": A  * X = B
;                              "T": Aᵀ * X = B
;                  iNRHS     - [Int] (Default: Default)
;                            ↳ number of right hand sides, i.e., the number of columns of the matrix B (to solve multiple systems at once)
;                  iM        - [Int] (Default: Default)
;                            ↳ number of rows of the matrix A
;                  iN        - [Int] (Default: Default)
;                            ↳ number of columns of the matrix A
;                  iLDA      - [Int] (Default: $iM)
;                            ↳ leading dimension of the array A (rows)
;                  iLDB      - [Int] (Default: $iN)
;                            ↳ leading dimension of the array B (max(M,N))
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True (@extended = size of WORK)
;                  Failure: False and set @error to:
;                           | 1: error during first DllCall of gels (@extended: @error from DllCall)
;                           | 2: error inside first call of gels (@extended: INFO-value from gels)
;                           | 3: determined size of WORK is not valid (@extended: determined LWORK)
;                           | 4: error during second DllCall of gels (@extended: @error from DllCall)
;                           | 5: error inside second call of gels (@extended: INFO-value from gels)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d8/d83/group__gels_gaa65298f8ef218a625e40d0da3c95803c.html#gaa65298f8ef218a625e40d0da3c95803c
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[1.0,2.0,0.5],[2.0,1.0,1.0],[1.5,1.5,1.5],[1.0,1.0,2.0],[0.5,2.5,1.0],[2.5,0.5,1.5]]")
;                  Global $mb = _blas_fromArray("[5.5,7.0,8.5,8.0,7.5,8.0]")
;                  _lp_gels($mA, $mB)
;                  $mB.elements = $mA.cols
;                  $mB.size     = $mA.cols
;                  _blas_display($mB, "solution vector/matrix x")
; ===============================================================================================================================
Func _lp_gels($mA, $mB, $cTRANS = "N", $iNRHS = Default, $iM = Default, $iN = Default, $iLDA = $iM, $iLDB = $iN, $sDataType = "DOUBLE")
	Local $pA, $pB ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iM) = 1 Then $iM = $cTRANS = "N" ? $mA.rows : $mA.cols
			If IsKeyword($iN) = 1 Then $iN = $cTRANS = "N" ? $mA.cols : $mA.rows
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect
	Select
		Case IsMap($mB)
			$pB = $mB.ptr
		Case IsPtr($mB)
			$pB = $mB
		Case IsDllStruct($mB)
			$pB = DllStructGetPtr($mB)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	If IsKeyword($iLDA)  = 1 Then $iLDA = $iM
	If IsKeyword($iLDB)  = 1 Then $iLDB = ($iM > $iN ? $iM : $iN)
	If IsKeyword($iNRHS) = 1 Then $iNRHS = $mB.cols < 1 ? 1 : $mB.cols

	; set char buffers and input healing:
	DllStructSetData($tBLASCHAR1, 1, $cTRANS <> "N" ? "T" : "N")

	; first run: determine optimal size of WORK
	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "gels", _
		"PTR",            $pBLASCHAR1, _ ; TRANS
		"INT*",           $iM, _         ; M
		"INT*",           $iN, _         ; N
		"INT*",           $iNRHS, _      ; NRHS - 0: no solving - only LU factorization, >0: number of systems with different b`s to be solved (every b = columns in $mB)
		"PTR",            0, _           ; A
		"INT*",           $iLDA, _       ; lda
		"PTR",            0, _           ; B
		"INT*",           $iLDB, _       ; LDB
		$sDataType & "*", 0, _           ; WORK buffer (here 1 element because we only determine the size)
		"INT*",           -1, _          ; LWORK
		"INT*",           0 _ 			 ; INFO
	)
	If @error Then Return SetError(1, @error, False)
	If $aDLL[11] <> 0 Then Return SetError(2, $aDLL[11], False)

	; declare working buffers
	Local $iLWork = $aDLL[9]
	If $iLWork < 1 Then Return SetError(3, $iLWork, False)
	Local $tWork = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iLWork))

	$aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "gels", _
		"PTR",  $pBLASCHAR1, _              ; TRANS
		"INT*", $iM, _                      ; M
		"INT*", $iN, _                      ; N
		"INT*", $iNRHS, _                   ; NRHS - 0: no solving - only LU factorization, >0: number of systems with different b`s to be solved (every b = columns in $mB)
		"PTR",  $pA, _                      ; A
		"INT*", $iLDA, _                    ; lda
		"PTR",  $pB, _                      ; B
		"INT*", $iLDB, _                    ; LDB
		"PTR",  DllStructGetPtr($tWork) , _ ; WORK buffer (here 1 element because we only determine the size)
		"INT*", $iLWORK, _                  ; LWORK
		"INT*", 0 _ 			            ; INFO
	)
	If @error Then Return SetError(4, @error, False)
	Return $aDLL[11] = 0 ? SetExtended($iLWORK, True) : SetError(5, $aDLL[11], False)

EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_getsls()
; Description ...: solves overdetermined or underdetermined linear system A * X = B
;                  using tall-skinny/short-wide QR/LQ factorization
; Syntax ........: _lp_getsls($mA, $mB, [$cTRANS = "N", [$iNRHS = Default, [$iM = Default, [$iN = Default, [$iLDA = $iN, [$iLDB = $iN, [$sDataType = "DOUBLE"]]]]]]])
; Parameters ....: mA        - [Map] matrix A (M × N) as a map, DllStruct or pointer (will be overwritten)
;                              on exit, contain the result of _lp_geqrf() (useful for deriving the cofactor matrix of the parameters)
;                  mB        - [Map] vector/matrix N × NRHS A as a map, DllStruct or pointer (will be overwritten)
;                            ↳ on exit, contain the solution values X in the first N elements and residual sum vᵀv in the last elements
;                  cTRANS    - [Char] (Default: "N")
;                            ↳ "N": A  * X = B
;                              "T": Aᵀ * X = B
;                  iNRHS     - [Int] (Default: Default)
;                            ↳ number of right hand sides, i.e., the number of columns of the matrix B (to solve multiple systems at once)
;                  iM        - [Int] (Default: Default)
;                            ↳ number of rows of the matrix A
;                  iN        - [Int] (Default: Default)
;                            ↳ number of columns of the matrix A
;                  iLDA      - [Int] (Default: $iM)
;                            ↳ leading dimension of the array A (rows)
;                  iLDB      - [Int] (Default: $iN)
;                            ↳ leading dimension of the array B (max(M,N))
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True (@extended = size of WORK)
;                  Failure: False and set @error to:
;                           | 1: error during first DllCall of gels (@extended: @error from DllCall)
;                           | 2: error inside first call of gels (@extended: INFO-value from gels)
;                           | 3: determined size of WORK is not valid (@extended: determined LWORK)
;                           | 4: error during second DllCall of gels (@extended: @error from DllCall)
;                           | 5: error inside second call of gels (@extended: INFO-value from gels)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/dd/dc8/group__getsls_ga2d8769d20f80cde1c8e8ceeab7a1cb7d.html#ga2d8769d20f80cde1c8e8ceeab7a1cb7d
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[1.0,2.0,0.5],[2.0,1.0,1.0],[1.5,1.5,1.5],[1.0,1.0,2.0],[0.5,2.5,1.0],[2.5,0.5,1.5]]")
;                  Global $mb = _blas_fromArray("[5.5,7.0,8.5,8.0,7.5,8.0]")
;                  _lp_getsls($mA, $mB)
;                  $mB.elements = $mA.cols
;                  $mB.size     = $mA.cols
;                  _blas_display($mB, "solution vector/matrix x")
; ===============================================================================================================================
Func _lp_getsls($mA, $mB, $cTRANS = "N", $iNRHS = Default, $iM = Default, $iN = Default, $iLDA = $iN, $iLDB = $iN, $sDataType = "DOUBLE")
	Local $pA, $pB ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iM) = 1 Then $iM = $cTRANS = "N" ? $mA.rows : $mA.cols
			If IsKeyword($iN) = 1 Then $iN = $cTRANS = "N" ? $mA.cols : $mA.rows
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect
	Select
		Case IsMap($mB)
			$pB = $mB.ptr
		Case IsPtr($mB)
			$pB = $mB
		Case IsDllStruct($mB)
			$pB = DllStructGetPtr($mB)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	If IsKeyword($iLDA)  = 1 Then $iLDA = $iM
	If IsKeyword($iLDB)  = 1 Then $iLDB = ($iM > $iN ? $iM : $iN)
	If IsKeyword($iNRHS) = 1 Then $iNRHS = $mB.cols < 1 ? 1 : $mB.cols

	; set char buffers and input healing:
	DllStructSetData($tBLASCHAR1, 1, $cTRANS <> "N" ? "T" : "N")

	; first run: determine optimal size of WORK
	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "getsls", _
		"PTR",            $pBLASCHAR1, _ ; TRANS
		"INT*",           $iM, _         ; M
		"INT*",           $iN, _         ; N
		"INT*",           $iNRHS, _      ; NRHS - 0: no solving - only LU factorization, >0: number of systems with different b`s to be solved (every b = columns in $mB)
		"PTR",            0, _           ; A
		"INT*",           $iLDA, _       ; lda
		"PTR",            0, _           ; B
		"INT*",           $iLDB, _       ; LDB
		$sDataType & "*", 0, _           ; WORK buffer (here 1 element because we only determine the size)
		"INT*",           -1, _          ; LWORK
		"INT*",           0 _ 			 ; INFO
	)
	If @error Then Return SetError(1, @error, False)
	If $aDLL[11] <> 0 Then Return SetError(2, $aDLL[11], False)

	; declare working buffers
	Local $iLWork = $aDLL[9]
	If $iLWork < 1 Then Return SetError(3, $iLWork, False)
	Local $tWork = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iLWork))

	$aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "getsls", _
		"PTR",  $pBLASCHAR1, _              ; TRANS
		"INT*", $iM, _                      ; M
		"INT*", $iN, _                      ; N
		"INT*", $iNRHS, _                   ; NRHS - 0: no solving - only LU factorization, >0: number of systems with different b`s to be solved (every b = columns in $mB)
		"PTR",  $pA, _                      ; A
		"INT*", $iLDA, _                    ; lda
		"PTR",  $pB, _                      ; B
		"INT*", $iLDB, _                    ; LDB
		"PTR",  DllStructGetPtr($tWork) , _ ; WORK buffer (here 1 element because we only determine the size)
		"INT*", $iLWORK, _                  ; LWORK
		"INT*", 0 _ 			            ; INFO
	)
	If @error Then Return SetError(1, @error, False)
	Return $aDLL[11] = 0 ? SetExtended($iLWORK, True) : SetError(2, $aDLL[11], False)

EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_gelss()
; Description ...: solves overdetermined or underdetermined linear system A * X = B
;                  using SVD factorization
; Syntax ........: _lp_gelss($mA, $mB, [$cTRANS = "N", [$iNRHS = Default, [$iM = Default, [$iN = Default, [$iLDA = $iN, [$iLDB = $iN, [$fRCOND = -1, [$sDataType = "DOUBLE"]]]]]]]])
; Parameters ....: mA        - [Map] matrix A (M × N) as a map, DllStruct or pointer (will be overwritten)
;                              on exit, contain the result of _lp_geqrf() (useful for deriving the cofactor matrix of the parameters)
;                  mB        - [Map] vector/matrix N × NRHS A as a map, DllStruct or pointer (will be overwritten)
;                            ↳ on exit, contain the solution values X in the first N elements and residual sum vᵀv in the last elements
;                  cTRANS    - [Char] (Default: "N")
;                            ↳ "N": A  * X = B
;                              "T": Aᵀ * X = B
;                  iNRHS     - [Int] (Default: Default)
;                            ↳ number of right hand sides, i.e., the number of columns of the matrix B (to solve multiple systems at once)
;                  iM        - [Int] (Default: Default)
;                            ↳ number of rows of the matrix A
;                  iN        - [Int] (Default: Default)
;                            ↳ number of columns of the matrix A
;                  iLDA      - [Int] (Default: $iM)
;                            ↳ leading dimension of the array A (rows)
;                  iLDB      - [Int] (Default: $iN)
;                            ↳ leading dimension of the array B (max(M,N))
;                  fRCOND    - [Float] (Default: -1)
;                            ↳ threshold value up to which singular values are regarded as zero (-1 = machine precision)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: matrix S as Dllstruct (@extended: rank of matrix A)
;                  Failure: Null and set @error to:
;                           | 1: error during first DllCall of gelss (@extended: @error from DllCall)
;                           | 2: error inside first call of gelss (@extended: INFO-value from gelss)
;                           | 3: determined size of WORK is not valid (@extended: determined LWORK)
;                           | 4: error during second DllCall of gelss (@extended: @error from DllCall)
;                           | 5: error inside second call of gelss (@extended: INFO-value from gelss)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......: slower than _lp_gels but with higher numerical stability and more robust for poorly conditioned problems.
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/da/d55/group__gelss_gac6159de3953ae0386c2799294745ac90.html#gac6159de3953ae0386c2799294745ac90
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[1.0,2.0,0.5],[2.0,1.0,1.0],[1.5,1.5,1.5],[1.0,1.0,2.0],[0.5,2.5,1.0],[2.5,0.5,1.5]]")
;                  Global $mb = _blas_fromArray("[5.5,7.0,8.5,8.0,7.5,8.0]")
;                  Global $tS = _lp_gelss($mA, $mB)
;                  $mB.elements = $mA.cols
;                  $mB.size     = $mA.cols
;                  _blas_display($mB, "solution vector/matrix x")
; ===============================================================================================================================
Func _lp_gelss($mA, $mB, $cTRANS = "N", $iNRHS = Default, $iM = Default, $iN = Default, $iLDA = $iN, $iLDB = $iN, $fRCOND = -1, $sDataType = "DOUBLE")
	Local $pA, $pB ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iM) = 1 Then $iM = $mA.rows
			If IsKeyword($iN) = 1 Then $iN = $mA.cols
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect
	Select
		Case IsMap($mB)
			$pB = $mB.ptr
		Case IsPtr($mB)
			$pB = $mB
		Case IsDllStruct($mB)
			$pB = DllStructGetPtr($mB)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	If IsKeyword($iLDA) = 1 Then $iLDA = $iM
	If IsKeyword($iLDB) = 1 Then $iLDB = ($iM > $iN ? $iM : $iN)
	If IsKeyword($iNRHS) = 1 Then $iNRHS = $mB.cols < 1 ? 1 : $mB.cols

	; set char buffers and input healing:
	DllStructSetData($tBLASCHAR1, 1, $cTRANS <> "N" ? "T" : "N")

	; output buffers
	Local $iMinMN = ($iM < $iN ? $iM : $iN)
	Local $tS  = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iMinMN))

	; first run: determine optimal size of WORK
	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "gelss", _
		"INT*",           $iM, _                   ; M
		"INT*",           $iN, _                   ; N
		"INT*",           $iNRHS, _                ; NRHS - 0: no solving - only LU factorization, >0: number of systems with different b`s to be solved (every b = columns in $mB)
		"PTR",            0, _                     ; A
		"INT*",           $iLDA, _                 ; lda
		"PTR",            0, _                     ; B
		"INT*",           $iLDB, _                 ; LDB
		"PTR",            0, _                     ; S
		$sDataType & "*", 0, _                     ; RCOND
		"INT*",           0, _                     ; RANK
		$sDataType & "*", 0, _                     ; WORK buffer (here 1 element because we only determine the size)
		"INT*",           -1, _                    ; LWORK
		"INT*",           0 _ 			           ; INFO
	)
	If @error Then Return SetError(1, @error, Null)
	If $aDLL[13] <> 0 Then Return SetError(2, $aDLL[13], Null)

	; declare working buffers
	Local $iLWork = $aDLL[11]
	If $iLWork < 1 Then Return SetError(3, $iLWork, Null)
	Local $tWork = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iLWork))

	$aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "gelss", _
		"INT*",           $iM, _                      ; M
		"INT*",           $iN, _                      ; N
		"INT*",           $iNRHS, _                   ; NRHS - 0: no solving - only LU factorization, >0: number of systems with different b`s to be solved (every b = columns in $mB)
		"PTR",            $pA, _                      ; A
		"INT*",           $iLDA, _                    ; lda
		"PTR",            $pB, _                      ; B
		"INT*",           $iLDB, _                    ; LDB
		"PTR",            DllStructGetPtr($tS), _     ; S
		$sDataType & "*", $fRCOND, _                  ; RCOND
		"INT*",           0, _                        ; RANK
		"PTR",            DllStructGetPtr($tWork), _  ; WORK buffer (here 1 element because we only determine the size)
		"INT*",           $iLWork, _                  ; LWORK
		"INT*",           0 _ 			              ; INFO
	)
	If @error Then Return SetError(4, @error, Null)
	If $aDLL[13] <> 0 Then Return SetError(5, $aDLL[13], Null)

	Return SetExtended($aDLL[10], $tS)
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_gelsy()
; Description ...: solves overdetermined or underdetermined linear system A * X = B
;                  using QR decomposition with column pivoting
; Syntax ........: _lp_gelsy($mA, $mB, [$iNRHS = Default, [$iM = Default, [$iN = Default, [$iLDA = $iM, [$iLDB = $iN, [$fRCOND = -1, [$sDataType = "DOUBLE"]]]]]]])
; Parameters ....: mA        - [Map] matrix A (M × N) as a map, DllStruct or pointer (will be overwritten)
;                              on exit, contain the results of its complete orthogonal factorization
;                  mB        - [Map] vector/matrix N × NRHS A as a map, DllStruct or pointer (will be overwritten)
;                            ↳ on exit, contain the solution values X in the first N elements and residual sum vᵀv in the last elements
;                  iNRHS     - [Int] (Default: Default)
;                            ↳ number of right hand sides, i.e., the number of columns of the matrix B (to solve multiple systems at once)
;                  iM        - [Int] (Default: Default)
;                            ↳ number of rows of the matrix A
;                  iN        - [Int] (Default: Default)
;                            ↳ number of columns of the matrix A
;                  iLDA      - [Int] (Default: $iM)
;                            ↳ leading dimension of the array A (rows)
;                  iLDB      - [Int] (Default: $iN)
;                            ↳ leading dimension of the array B (max(M,N))
;                  fRCOND    - [Float] (Default: -1)
;                            ↳ used as threshold to determine the effective rank of A
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: permutation vector JPVT as DllStruct (@extended = rank of matrix A)
;                  Failure: Null and set @error to:
;                           | 1: error during first DllCall of gelsy (@extended: @error from DllCall)
;                           | 2: error inside first call of gelsy (@extended: INFO-value from gelsy)
;                           | 3: determined size of WORK is not valid (@extended: determined LWORK)
;                           | 4: error during second DllCall of gelsy (@extended: @error from DllCall)
;                           | 5: error inside second call of gelsy (@extended: INFO-value from gelsy)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......: numerically more stable than _lp_gels() and useful if the matrix A has numerical instabilities or rank deficiencies.
;                  for retrieve cofactor matrix for the parameters you have firstly rearange A with _lp_lapmt with bForward = False
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/dc/d8b/group__gelsy_ga6d1d46ead18df76e993cd4eda6dc1bbb.html#ga6d1d46ead18df76e993cd4eda6dc1bbb
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[1.0,2.0,0.5],[2.0,1.0,1.0],[1.5,1.5,1.5],[1.0,1.0,2.0],[0.5,2.5,1.0],[2.5,0.5,1.5]]")
;                  Global $mb = _blas_fromArray("[5.5,7.0,8.5,8.0,7.5,8.0]")
;                  _lp_gelsy($mA, $mB)
;                  $mB.elements = $mA.cols
;                  $mB.size     = $mA.cols
;                  _blas_display($mB, "rank:" & @extended)
; ===============================================================================================================================
Func _lp_gelsy($mA, $mB, $iNRHS = Default, $iM = Default, $iN = Default, $iLDA = $iM, $iLDB = $iN, $fRCOND = 1e5, $sDataType = "DOUBLE")
	Local $pA, $pB ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iM) = 1 Then $iM = $mA.rows
			If IsKeyword($iN) = 1 Then $iN = $mA.cols
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect
	Select
		Case IsMap($mB)
			$pB = $mB.ptr
		Case IsPtr($mB)
			$pB = $mB
		Case IsDllStruct($mB)
			$pB = DllStructGetPtr($mB)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	If IsKeyword($iLDA)  = 1 Then $iLDA = $iM
	If IsKeyword($iLDB)  = 1 Then $iLDB = ($iM > $iN ? $iM : $iN)
	If IsKeyword($iNRHS) = 1 Then $iNRHS = $mB.cols < 1 ? 1 : $mB.cols

	; first run: determine optimal size of WORK
	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "gelsy", _
		"INT*",           $iM, _         ; M
		"INT*",           $iN, _         ; N
		"INT*",           $iNRHS, _      ; NRHS - 0: no solving - only LU factorization, >0: number of systems with different b`s to be solved (every b = columns in $mB)
		"PTR",            0, _           ; A
		"INT*",           $iLDA, _       ; lda
		"PTR",            0, _           ; B
		"INT*",           $iLDB, _       ; LDB
		"PTR",            0, _           ; JPVT
		$sDataType & "*", 0, _           ; RCOND
		"INT*",           0, _           ; RANK
		$sDataType & "*", 0, _           ; WORK buffer (here 1 element because we only determine the size)
		"INT*",           -1, _          ; LWORK
		"INT*",           0 _ 			 ; INFO
	)
	If @error Then Return SetError(1, @error, Null)
	If $aDLL[13] <> 0 Then Return SetError(2, $aDLL[13], Null)

	; set JPVT Buffer
	Local $tJPVT = DllStructCreate("INT[" & $iN & "]")

	; declare working buffers
	Local $iLWork = $aDLL[11]
	If $iLWork < 1 Then Return SetError(3, $iLWork, Null)
	Local $tWork = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iLWork))

	$aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "gelsy", _
		"INT*",           $iM, _                      ; M
		"INT*",           $iN, _                      ; N
		"INT*",           $iNRHS, _                   ; NRHS - 0: no solving - only LU factorization, >0: number of systems with different b`s to be solved (every b = columns in $mB)
		"PTR",            $pA, _                      ; A
		"INT*",           $iLDA, _                    ; lda
		"PTR",            $pB, _                      ; B
		"INT*",           $iLDB, _                    ; LDB
		"PTR",            DllStructGetPtr($tJPVT), _  ; JPVT
		$sDataType & "*", $fRCOND, _                  ; RCOND
		"INT*",           0, _                        ; RANK
		"PTR",            DllStructGetPtr($tWork), _  ; WORK
		"INT*",           $iLWork, _                  ; LWORK
		"INT*",           0 _ 			              ; INFO
	)
	If @error Then Return SetError(4, @error, Null)
	Return $aDLL[13] = 0 ? SetExtended($aDLL[10], $tJPVT) : SetError(5, $aDLL[13], Null)

EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_geqrs()
; Description ...: solves overdetermined or underdetermined linear system A * X = B
;                  using the results of the QR decomposition from _lp_geqrf()
; Syntax ........: _lp_geqrs($mA, $mB, $tTau, [$iLWork = Default, [$iM = Default, [$iN = Default, [$iNRHS = Default, [$iLDA = Default, [$iLDB = Default, [$sDataType = "DOUBLE"]]]]]]])
; Parameters ....: mA        - [Map] Details of the QR factorization in Matrix A from _lp_geqrf()
;                  mB        - [Map] ector/matrix N × NRHS A as a map, DllStruct or pointer (will be overwritten)
;                            ↳ on exit, contain the solution values X in the first N elements and residual sum vᵀv in the last elements
;                  tTau      - [DllStruct] Details of the orthogonal matrix Q as returned by _lp_geqrf()
;                  iLWork    - [Int] (Default: Default)
;                            ↳ length of the array WORK as returned in the @extended-macro from _lp_geqrf()
;                  iM        - [Int] (Default: Default)
;                            ↳ number of rows in mA
;                  iN        - [Int] (Default: Default)
;                            ↳ number of cols in mA
;                  iNRHS     - [Int] (Default: Default)
;                            ↳ number of right hand sides, i.e., the number of columns of the matrix B (to solve multiple systems at once)
;                  iLDA      - [Int] (Default: $iM)
;                            ↳ leading dimension of the array A (rows)
;                  iLDB      - [Int] (Default: $iN)
;                            ↳ leading dimension of the array B (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True (@extended = LWORK)
;                  Failure: Null and set @error to:
;                           | 1: value for LWORK is not valid (@extended:LWORK)
;                           | 2: error during second DllCall of geqrs (@extended: @error from DllCall)
;                           | 3: error inside second call of geqrs (@extended: INFO-value from geqrs)
; Author ........: AspirinJunkie
; Modified.......: 2024-10-01
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d3/d17/dgeqrs_8f_a271f8720489177ac8963f6d8fc61923d.html#a271f8720489177ac8963f6d8fc61923d
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[1,3,-2],[3,5,6],[2,4,3]]")
;                  Global $mB = _blas_fromArray("[5,7,8]")
;                  Global $tTau = _lp_geqrf($mA)
;                  _lp_geqrs($mA, $mB, $tTau, @extended) ; should be: [-15, 8, 2]
;                  _blas_display($mB, "X")
; ===============================================================================================================================
Func _lp_geqrs($mA, $mB, $tTau, $iLWork = Default, $iM  = Default, $iN  = Default, $iNRHS = Default, $iLDA = $iM, $iLDB = $iN, $sDataType = "DOUBLE")
	Local $pA, $pB ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iM)   = 1 Then $iM   = $mA.rows
			If IsKeyword($iN)   = 1 Then $iN   = $mA.cols
			If IsKeyword($iLDA) = 1 Then $iLDA = $iM
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect
	Select
		Case IsMap($mB)
			$sDataType = $mB.datatype
			If IsKeyword($iNRHS) = 1 Then $iNRHS = $mB.cols
			If IsKeyword($iLDB) = 1 Then $iLDB = $iM
			$pB = $mB.ptr
		Case IsPtr($mB)
			$pB = $mB
		Case IsDllStruct($mB)
			$pB = DllStructGetPtr($mB)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; declare working buffers
	If IsKeyword($iLWork) = 1 Then $iLWORK = $iNRHS
	If $iLWork < 1 Then Return SetError(1, $iLWork, False)
	Local $tWork = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iLWork))

	$aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "geqrs", _
		"INT*", $iM, _                     ; M
		"INT*", $iN, _                     ; N
		"INT*", $iNRHS, _                  ; NRHS
		"PTR",  $pA, _                     ; A
		"INT*", $iLDA, _                   ; LDA
		"PTR",  DllStructGetPtr($tTau), _  ; TAU
		"PTR",  $pB, _                     ; B
		"INT*", $iLDB, _                   ; LDB
		"PTR",  DllStructGetPtr($tWork), _ ; WORK buffer
		"INT*", $iLWork, _                 ; LWORK
		"INT*", 0 _ 			           ; INFO
	)
	If @error Then Return SetError(2, @error, False)
	Return $aDLL[11] = 0 ? SetExtended($iLWork, True) : SetError(3, $aDLL[11], False)
EndFunc

#EndRegion

#Region auxiliary functions


; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_lassq()
; Description ...: calculate the sum of squares of elements of a matrix/vector
; Syntax ........: _lp_lassq($mMatrix, [$iStart = 0, [$iInc = 1, [$iN = Default, [$fScale = 1.0, [$fSumSq = 0.0, [$sDataType = "DOUBLE"]]]]]])
; Parameters ....: mMatrix   - [Map] matrix/vector as a map, DllStruct or pointer
;                  iStart    - [Int] (Default: 0)
;                            ↳ start index
;                  iInc      - [Int] (Default: 1)
;                            ↳ storage spacing between elements
;                  iN        - [Int] (Default: Default)
;                            ↳ number of elements in input vector
;                  fScale    - [Float] (Default: 1.0)
;                            ↳ Scaling factor with which large and small values can be adjusted so that there is no overflow (1 is a good start for this)
;                  fSumSq    - [Float] (Default: 0.0)
;                            ↳ The previous sum of squares to which the result is to be added
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: [Float] square sum (@extended: updated Scale value (only Integer values possible))
;                  Failure: Null and set @error to:
;                           | 1: error during DllCall of lassq (@extended: @error from DllCall )
; Author ........: AspirinJunkie
; Modified.......: 2024-09-10
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d8/d76/group__lassq_gae8f40b0a34771b4f2d9c863de3af7be5.html#gae8f40b0a34771b4f2d9c863de3af7be5
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[1,2,3,4,5]")
;                  ConsoleWrite(_lp_lassq($mA) & @CRLF)
; ===============================================================================================================================
Func _lp_lassq($mMatrix, $iStart = 0, $iInc = 1, $iN = Default, $fScale = 1.0, $fSumSq = 0.0, $sDataType = "DOUBLE")
	Local $pM ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrix)
			$sDataType = $mMatrix.datatype
			; calculate number of elements if not defined
			If IsKeyword($iN) = 1 Then $iN = $mMatrix.elements - $iStart
			$pM = $mMatrix.ptr
		Case IsPtr($mMatrix)
			$pM = $mMatrix
		Case IsDllStruct($mMatrix)
			$pM = DllStructGetPtr($mMatrix)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d", _
	            $dSize   = ($sDataType = "DOUBLE") ? $iBLAS_SIZE_DOUBLE : $iBLAS_SIZE_FLOAT

	; call the function
	Local $aDLL = DllCall($__g_hBLAS_DLL, $sDataType & ":cdecl", $cPrefix & "lassq", _
		"INT*", $iN, _                    ; number of elements to check
		"PTR",  $pM + $dSize * $iStart, _ ; D (start ptr)
		"INT*", $iInc, _                  ; the stride for vector x
		$sDataType & "*", $fScale, _      ; SCALE
		$sDataType & "*", $fSumSq _       ; SUMSQ
	)
	If @error Then Return SetError(1, @error, Null)
	Return SetExtended($aDLL[4], $aDLL[5])
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_lasrt()
; Description ...: sort the values of a matrix/vector
; Syntax ........: _lp_lasrt($mMatrix, [$cOrder = "I", [$iStart = 0, [$iN = Default, [$sDataType = "DOUBLE"]]]])
; Parameters ....: mMatrix   - [Map] matrix/vector as a map, DllStruct or pointer (will be overwritten)
;                  cOrder    - [Char] (Default: "I")
;                            ↳ "D": decreasing order
;                              "I": increasing order
;                  iStart    - [Int] (Default: 0)
;                            ↳ start index
;                  iN        - [Int] (Default: Default)
;                            ↳ number of elements in input vector
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of lasrt (@extended: @error from DllCall)
;                           | 2: error inside call of lasrt (@extended: INFO-value from lasrt)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-10
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d3/d73/group__lasrt_ga45b6a3513c2190c1056ce557826f330b.html#ga45b6a3513c2190c1056ce557826f330b
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[4,-5,2,3,1]")
;                  _lp_lasrt($mA)
;                  _blas_display($mA, "sorted vector")
; ===============================================================================================================================
Func _lp_lasrt($mMatrix, $cOrder = "I", $iStart = 0, $iN = Default, $sDataType = "DOUBLE")
	Local $pM ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrix)
			$sDataType = $mMatrix.datatype
			; calculate number of elements if not defined
			If IsKeyword($iN) = 1 Then $iN = $mMatrix.elements - $iStart
			$pM = $mMatrix.ptr
		Case IsPtr($mMatrix)
			$pM = $mMatrix
		Case IsDllStruct($mMatrix)
			$pM = DllStructGetPtr($mMatrix)
	EndSelect

	DllStructSetData($tBLASCHAR1, 1, $cOrder = "D" ? "D" : "I")

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d", _
	            $dSize   = ($sDataType = "DOUBLE") ? $iBLAS_SIZE_DOUBLE : $iBLAS_SIZE_FLOAT

	; call the function
	Local $aDLL = DllCall($__g_hBLAS_DLL, $sDataType & ":cdecl", $cPrefix & "lasrt", _
		"PTR",  $pBLASCHAR1, _            ; ID -> 'I': increasing, 'D': decreasing
		"INT*", $iN, _                    ; number of elements to check
		"PTR",  $pM + $dSize * $iStart, _ ; D (start ptr)
		"INT*", 0 _                       ; INFO
	)
	If @error Then Return SetError(1, @error, False)
	If $aDLL[4] <> 0 Then Return SetError(2, $aDLL[4], False)
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_rscl()
; Description ...: divide matrix/vector elements with a scalar a ( = 1/a * X)
; Syntax ........: _lp_rscl($mVector, $fScale, [$iStart = 0, [$iInc = 1, [$iN = Default, [$sDataType = "DOUBLE"]]]])
; Parameters ....: mVector   - [Map] matrix/vector as a map, DllStruct or pointer (will be overwritten)
;                  fScale    - [Float] the divisor scalar
;                  iStart    - [Int] (Default: 0)
;                            ↳ start index
;                  iInc      - [Int] (Default: 1)
;                            ↳ storage spacing between elements
;                  iN        - [Int] (Default: Default)
;                            ↳ number of elements in input vector
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of rscl (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-10
; Remarks .......: You could also use _blas_scal() by taking 1/a as a factor. However, this function here is numerically more stable and optimized for this case.
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/dd/dc7/group__rscl_ga51f6c8a55ac4b25e6e39e9b4996829aa.html#ga51f6c8a55ac4b25e6e39e9b4996829aa
; Example .......: Yes
;                  Global $mVector = _blas_fromArray("[1,2,3,4]")
;                  _lp_rscl($mVector, 3) ; divide elements with 3
;                  _blas_display($mVector, "reciprocal scaled vector")
; ===============================================================================================================================
Func _lp_rscl($mVector, $fScale, $iStart = 0, $iInc = 1, $iN = Default, $sDataType = "DOUBLE")
	Local $pV ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mVector)
			$sDataType = $mVector.datatype
			; calculate number of elements if not defined
			If IsKeyword($iN) = 1 Then $iN = Floor(($mVector.elements - $iStart) / $iInc)
			$pV = $mVector.ptr
		Case IsPtr($mVector)
			$pV = $mVector
		Case IsDllStruct($mVector)
			$pV = DllStructGetPtr($mVector)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d", _
	            $dSize   = ($sDataType = "DOUBLE") ? $iBLAS_SIZE_DOUBLE : $iBLAS_SIZE_FLOAT

	DllCall($__g_hBLAS_DLL, $sDataType & ":cdecl", $cPrefix & "rscl", _
		"INT*",           $iN, _                    ; number of elements to check
		$sDataType & "*", $fScale, _                ; scalar
		"PTR",            $pV + $dSize * $iStart, _ ; start ptr to read
		"INT*",           $iInc _                   ; increment
	)
	Return @error ? SetError(1, @error, False) : True
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_lange()
; Description ...: calculates the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of a general rectangular matrix
; Syntax ........: _lp_lange($mA, [$cNORM = "1", [$iM = Default, [$iN = Default, [$iLDA = $iM, [$sDataType = "DOUBLE"]]]]])
; Parameters ....: mA        - [Map] matrix A (M × N) as a map, DllStruct or pointer
;                  cNORM     - [Char] (Default: "1")
;                            ↳ "1": one norm (max column sum)
;                              "I": infinity norm (max row sum)
;                              "F": Frobenius norm (square root of sum of squares)
;                              "M": largest absolute value
;                  iM        - [Int] (Default: Default)
;                            ↳ number of rows in matrix A
;                  iN        - [Int] (Default: Default)
;                            ↳ number of columns in matrix A
;                  iLDA      - [Int] (Default: $iM)
;                            ↳ leading dimension of the array A (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: [Float] the norm value
;                  Failure: Null and set @error to:
;                           | 1: error during DllCall of lange (@extended: @error from DllCall)
;                           | 2: invalid value for cNORM
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d8/d2e/group__lange_ga8581d687290b36c6e24fe76b3be7caa3.html#ga8581d687290b36c6e24fe76b3be7caa3
; Example .......: Yes
;                  Global $mA = _blas_fromArray('[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97],[93,89,13,-53,43]]')
;                  ConsoleWrite("1-norm: " & _lp_lange($mA, "1") & @CRLF)
; ===============================================================================================================================
Func _lp_lange($mA, $cNORM = "1", $iM = Default, $iN = Default, $iLDA = $iM, $sDataType = "DOUBLE")
	Local $pA ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iM)   = 1 Then $iM   = $mA.rows
			If IsKeyword($iN)   = 1 Then $iN   = $mA.cols
			If IsKeyword($iLDA) = 1 Then $iLDA = $iM
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	Local $pWORK = 0
	Switch $cNORM
		Case '1', 'O', 'o', 'F', 'f', 'E', 'e', 'M', 'm'
		; one norm, frobenius norm, largest absolute value

		Case 'I', 'i'
		; infinity norm
			Local $tWork = DllStructCreate(StringFormat("%s[%s]", $sDataType, $iM))
			$pWORK = DllStructGetPtr($tWork)
		Case Else
			Return SetError(2, 0, Null)

	EndSwitch

	; set char buffer
	DllStructSetData($tBLASCHAR1, 1, $cNORM)

	Local $aDLL = DllCall($__g_hBLAS_DLL, $sDataType & ":cdecl", $cPrefix & "lange", _
	    "PTR",  $pBLASCHAR1, _ ; NORM
		"INT*", $iM, _         ; M
		"INT*", $iN, _         ; N
		"PTR",  $pA, _         ; A
		"INT*", $iLDA, _       ; LDA
		"PTR",  $pWORK _       ; WORK
	)
	Return @error ? SetError(1, @error, Null) : $aDLL[0]
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_gecon()
; Description ...: estimates the reciprocal of the condition number of a general real matrix A
; Syntax ........: _lp_gecon($mA, $fANORM, [$cNORM = "1", [$iN = Default, [$iLDA = $iN, [$sDataType = "DOUBLE"]]]])
; Parameters ....: mA        - [Map] matrix A (M × N) as a map, DllStruct or pointer
;                  fANORM    - [Float] norm of matrix A (calc with _lp_lange())
;                  cNORM     - [Char] (Default: "1")
;                            ↳ "1": one norm (max column sum)
;                              "I": infinity norm (max row sum)
;                  iN        - [Int] (Default: Default)
;                            ↳ order matrix A
;                  iLDA      - [Int] (Default: $iM)
;                            ↳ leading dimension of the array A (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: [Float] the condition number of matrix A
;                  Failure: Null and set @error to:
;                           | 1: error during first DllCall of gecon (@extended: @error from DllCall)
;                           | 2: error inside first call of gecon (@extended: INFO-value from gecon)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d4/daf/group__gecon_ga4f9b830e19e12c7f082ddb497a57af18.html#ga4f9b830e19e12c7f082ddb497a57af18
; Example .......: Yes
;                  Global $mA = _blas_fromArray('[[19,-80,-55,-58,21],[29,-92,-67,-94,-25],[-12,-29,77,-7,40],[96,96,-70,12,97],[93,89,13,-53,43]]')
;                  Global $fANORM = _lp_lange($mA, "1")
;                  _lp_getrf($mA) ; LU factorization of A
;                  Global $fCond = 1 / _lp_gecon($mA, $fANORM, "1")
;                  ConsoleWrite("condition number: " & $fCond & @CRLF)
; ===============================================================================================================================
Func _lp_gecon($mA, $fANORM, $cNORM = "1", $iN = Default, $iLDA = $iN, $sDataType = "DOUBLE")
	Local $pA ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iN)   = 1 Then $iN   = $mA.rows < $mA.cols ? $mA.rows : $mA.cols  ; min(m,n) - because dgetrf saves LU this way in A
			If IsKeyword($iLDA) = 1 Then $iLDA = $iN
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d", _
	            $sType = $sDataType & "*"

	; set char buffer
	DllStructSetData($tBLASCHAR1, 1, $cNORM = "I" ? "I" : "1")

	Local $tWORK = DllStructCreate(StringFormat("%s[%d]", $sDataType, 4 * $iN))
	Local $tIWORK = DllStructCreate(StringFormat("INT[%d]", $iN))

	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "gecon", _
	    "PTR",  $pBLASCHAR1, _              ; NORM
		"INT*", $iN, _                      ; N
		"PTR",  $pA, _                      ; A
		"INT*", $iLDA, _                    ; LDA
		$sType, $fANORM, _                  ; ANORM
		$sType, 0, _                        ; RCOND
		"PTR",  DllStructGetPtr($tWORK), _  ; WORK
		"PTR",  DllStructGetPtr($tIWORK), _ ; IWORK
		"INT*", 0 _                         ; INFO
	)
	If @error Then Return SetError(1, @error, Null)
	Return $aDLL[9] = 0 ? $aDLL[6] : SetError(2, $aDLL[9], Null)
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_lacpy()
; Description ...: copies all or part of a two-dimensional matrix A to another matrix B
; Syntax ........: _lp_lacpy($mA, $mB, [$cUPLO = "X", [$iM = Default, [$iN = Default, [$iLDA = $iM, [$iLDB = $iM, [$sDataType = "DOUBLE"]]]]]])
; Parameters ....: mA        - [Map] matrix A (M × N) as a map, DllStruct or pointer
;                  mB        - [Map] matrix B as a map, DllStruct or pointer (gets overwritten)
;                  cUPLO     - [Char] (Default: "X")
;                            ↳ "U": upper triangle of A to be copied to B
;                              "L": lower triangle of A to be copied to B
;                              "X": all of matrix A to be copied to B
;                  iM        - [Int] (Default: Default)
;                            ↳ number of rows of the matrix A
;                  iN        - [Int] (Default: Default)
;                            ↳ number of columns of the matrix A
;                  iLDA      - [Int] (Default: $iM)
;                            ↳ leading dimension of the array A (rows)
;                  iLDB      - [Int] (Default: $iM)
;                            ↳ leading dimension of the array B (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of lacpy (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d0/d9e/group__lacpy_gaba7ee02955a93bf8af4a432c98734e65.html#gaba7ee02955a93bf8af4a432c98734e65
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[11,12,13,14,15,16],[21,22,23,24,25,26],[31,32,33,34,35,36],[41,42,43,44,45,46],[51,52,53,54,55,56],[61,62,63,64,65,66],[71,72,73,74,75,76],[81,82,83,84,85,86]]")
;                  Global $iNewLines = 3
;                  Global $mB = _blas_createMatrix($iNewLines, $mA.cols, $mA.datatype)
;                  _lp_lacpy($mA, $mB, "X", $iNewLines, $mA.cols, $mA.rows, $mB.rows, $mA.datatype)
;                  _blas_display($mB, "extracted first " & $iNewLines & " lines of mA")
; ===============================================================================================================================
Func _lp_lacpy($mA, $mB, $cUPLO = "X", $iM = Default, $iN = Default, $iLDA = $iM, $iLDB = $iM, $sDataType = "DOUBLE")
	Local $pA, $pB ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iM)   = 1 Then $iM   = $mA.rows
			If IsKeyword($iN)   = 1 Then $iN   = $mA.cols
			If IsKeyword($iLDA) = 1 Then $iLDA = $iM
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect
	Select
		Case IsMap($mB)
			$sDataType = $mB.datatype
			If IsKeyword($iLDB) = 1 Then $iLDB = $mB.rows
			$pB = $mB.ptr
		Case IsPtr($mB)
			$pB = $mB
		Case IsDllStruct($mB)
			$pB = DllStructGetPtr($mB)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	If IsKeyword($iLDB) = 1 Then $iLDB = $iM

	; set char buffer
	DllStructSetData($tBLASCHAR1, 1, $cUPLO)

	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "lacpy", _
	    "PTR",  $pBLASCHAR1, _ ; UPLO
		"INT*", $iM, _         ; M
		"INT*", $iN, _         ; N
		"PTR",  $pA, _         ; A
		"INT*", $iLDA, _       ; LDA
		"PTR",  $pB, _         ; B
		"INT*", $iLDB _        ; LDB
	)
	Return @error ? SetError(1, @error, False) : True

EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_laset()
; Description ...: initializes the off-diagonal elements and the diagonal elements of a matrix to given values
; Syntax ........: _lp_laset($mA, [$cUPLO = "X", [$fAlpha = 0, [$fBeta = 1, [$iM = Default, [$iN = Default, [$iLDA = $iM, [$sDataType = "DOUBLE"]]]]]]])
; Parameters ....: mA        - [Map] matrix A (M × N) as a map, DllStruct or pointer
;                  cUPLO     - [Char] (Default: "X")
;                            ↳ "U": upper triangle of A to be set
;                              "L": lower triangle of A to be set
;                              "X": all of matrix A to be set
;                  fAlpha    - [Float] (Default: 0)
;                            ↳ constant to which the offdiagonal elements are to be set
;                  fBeta     - [Float] (Default: 1)
;                            ↳ constant to which the diagonal elements are to be set
;                  iM        - [Int] (Default: Default)
;                            ↳ number of rows of the matrix A
;                  iN        - [Int] (Default: Default)
;                            ↳ number of columns of the matrix A
;                  iLDA      - [Int] (Default: $iM)
;                            ↳ leading dimension of the array A (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of laswp (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d0/de5/group__laset_gad8051330f20413bd2a4ee0bccaf54ec8.html#gad8051330f20413bd2a4ee0bccaf54ec8
; Example .......: Yes
;                  Global $mMatrix = _blas_createMatrix(25,25)
;                  _lp_laset($mMatrix, "U", 2, 1)
;                  _blas_display($mMatrix)
; ===============================================================================================================================
Func _lp_laset($mA, $cUPLO = "X", $fAlpha = 0, $fBeta = 1, $iM = Default, $iN = Default, $iLDA = $iM, $sDataType = "DOUBLE")
	Local $pA ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iM)   = 1 Then $iM   = $mA.rows
			If IsKeyword($iN)   = 1 Then $iN   = $mA.cols
			If IsKeyword($iLDA) = 1 Then $iLDA = $iN
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d", _
	            $sType   = $sDataType & "*"

	; set char buffer
	DllStructSetData($tBLASCHAR1, 1, $cUPLO)

	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "laset", _
	    "PTR",  $pBLASCHAR1, _ ; UPLO
		"INT*", $iM, _         ; M
		"INT*", $iN, _         ; N
		$sType, $fAlpha, _     ; Alpha
		$sType, $fBeta, _      ; Beta
		"PTR",  $pA, _         ; the general matrix A
		"INT*", $iLDA _        ; LDA
	)
	Return @error ? SetError(1, @error, False) : True
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_laswp()
; Description ...: performs a series of row interchanges on a general rectangular matrix
; Syntax ........: _lp_laswp($mA, $tIPIV, [$iK1 = 1, [$iK2 = Default, [$iIncX = 1, [$iN = Default, [$iLDA = Default, [$sDataType = "DOUBLE"]]]]]])
; Parameters ....: mA        - [Map] matrix A (LDA × N) as a map, DllStruct or pointer
;                  tIPIV     - [DllStruct]
;                  iK1       - [Int] (Default: 1)
;                            ↳ first element of IPIV for which a row interchange will be done
;                  iK2       - [Int] (Default: Default)
;                            ↳ number of elements of IPIV for which a row interchange will be done
;                  iIncX     - [Int] (Default: 1)
;                            ↳ increment between successive values of IPIV
;                              if negative, the pivots are applied in reverse order
;                  iN        - [Int] (Default: Default)
;                            ↳ number of columns of the matrix A
;                  iLDA      - [Int] (Default: Default)
;                            ↳ leading dimension of the array A (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of laswp (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d1/d7e/group__laswp_ga5d3ea3e3cb61e32750bf062a2446aa33.html#ga5d3ea3e3cb61e32750bf062a2446aa33
; Example .......: Yes
;                  Global $mMatrix = _blas_fromArray("[[6,2,1,1],[2,5,2,2],[1,2,4,3]]")
;                  Global $tIPIV = DllStructCreate("INT[3]")
;                  DllStructSetData($tIPIV, 1, 3, 1)
;                  DllStructSetData($tIPIV, 1, 2, 2)
;                  DllStructSetData($tIPIV, 1, 3, 3)
;                  _lp_laswp($mMatrix, $tIPIV)
;                  _blas_display($mMatrix, "swapped rows")
; ===============================================================================================================================
Func _lp_laswp($mA, $tIPIV, $iK1 = 1, $iK2 = Default, $iIncX = 1, $iN = Default, $iLDA = Default, $sDataType = "DOUBLE")
	Local $pA ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iN)   = 1 Then $iN   = $mA.cols
			If IsKeyword($iLDA) = 1 Then $iLDA = $mA.rows
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect

	If IsKeyword($iK2) = 1 Then $iK2 = $iLDA

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "laswp", _
		"INT*", $iN, _                     ; N
		"PTR",  $pA, _                     ; A
		"INT*", $iLDA, _                   ; LDA
		"INT*", $iK1, _                    ; K1
		"INT*", $iK2, _                    ; K2
	    "PTR",  DllStructGetPtr($tIPIV), _ ; IPIV
		"INT*", $iIncX _                   ; LDA
	)
	Return @error ? SetError(1, @error, False) : True
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_lapmt()
; Description ...: rearranges the columns of the M by N matrix X as specified by a permutation
; Syntax ........: _lp_lapmt($mX, [$tJPVT = True, [$bForward = True, [$iM = Default, [$iN = Default, [$iLDX = $iM, [$sDataType = "DOUBLE"]]]]]])
; Parameters ....: mX        - [Map] matrix X (LDX × N) as a map, DllStruct or pointer
;                  tJPVT     - [DllStruct] the permutation vector as DllStruct (INT[N])
;                  bForward  - [Bool] (Default: True)
;                            ↳ True: forward permutation
;                              False: backward permutation
;                  iM        - [Int] (Default: Default)
;                            ↳ number of rows of the matrix X
;                  iN        - [Int] (Default: Default)
;                            ↳ number of columns of the matrix X
;                  iLDX      - [Int] (Default: $iM)
;                            ↳ leading dimension of the array X (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of lapmt (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d0/dcb/group__lapmt_ga3f8f8091894a18247775ddef5da5f817.html#ga3f8f8091894a18247775ddef5da5f817
; Example .......: Yes
;                  Global $mMatrix = _blas_fromArray("[[6,2,1,1],[2,5,2,2],[1,2,4,3]]")
;                  _blas_display($mMatrix)
;                  Global $tIPIV = DllStructCreate("INT[4]")
;                  DllStructSetData($tIPIV, 1, 3, 1)
;                  DllStructSetData($tIPIV, 1, 2, 2)
;                  DllStructSetData($tIPIV, 1, 3, 3)
;                  DllStructSetData($tIPIV, 1, 2, 4)
;                  _lp_lapmt($mMatrix, $tIPIV)
;                  _blas_display($mMatrix, "swapped cols")
; ===============================================================================================================================
Func _lp_lapmt($mX, $tJPVT, $bForward = True, $iM = Default, $iN = Default, $iLDX = $iM, $sDataType = "DOUBLE")
	Local $pX ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mX)
			$sDataType = $mX.datatype
			If IsKeyword($iM)   = 1 Then $iM   = $mX.rows
			If IsKeyword($iN)   = 1 Then $iN   = $mX.cols
			If IsKeyword($iLDX) = 1 Then $iLDX = $iM
			$pX = $mX.ptr
		Case IsPtr($mX)
			$pX = $mX
		Case IsDllStruct($mX)
			$pX = DllStructGetPtr($mX)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "lapmt", _
		"BOOLEAN*", $bForward, _              ; FORWRD
		"INT*",     $iM, _                    ; M
		"INT*",     $iN, _                    ; N
		"PTR",      $pX, _                    ; X
		"INT*",     $iLDX, _                  ; LDX
	    "PTR",      DllStructGetPtr($tJPVT) _ ; IPIV
	)
	Return @error ? SetError(1, @error, False) : True
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_lauum()
; Description ...: computes the product U * Uᵀ or Lᵀ * L, where U or L are triangular matrices (chooses if blocked or unblocked variant is used)
; Syntax ........: _lp_lauum($mA, [$cUPLO = "U", [$iN = Default, [$iLDA = $iN, [$sDataType = "DOUBLE"]]]])
; Parameters ....: mA        - [Map] matrix A (LDA × N) as a map, DllStruct or pointer (gets overwritten)
;                  cUPLO     - [Char] (Default: "U")
;                            ↳ "U": upper triangle of A is used
;                              "L": lower triangle of A is used
;                  iN        - [Int] (Default: Default)
;                            ↳ order of the triangular factor U or L (rows)
;                  iLDA      - [Int] (Default: $iN)
;                            ↳ leading dimension of the array (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during first DllCall of lauum (@extended: @error from DllCall)
;                           | 2: error inside first call of lauum (@extended: INFO-value from lauum)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d5/d18/group__lauum_ga109a49941c86e4258f842ceb82952fb1.html#ga109a49941c86e4258f842ceb82952fb1
; Example .......: Yes
;                  Global $mU = _blas_fromArray("[[1,2,3],[0,4,5],[0,0,6]]")
;                  _lp_lauum($mU)
;                  _blas_display($mU)
; ===============================================================================================================================
Func _lp_lauum($mA, $cUPLO = "U", $iN = Default, $iLDA = $iN, $sDataType = "DOUBLE")
	Local $pA ; pointer to the data in memory

	; Triangle * Triangleᵀ = symmetric

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iN) = 1 Then $iN = $mA.rows
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	If IsKeyword($iLDA) = 1 Then $iLDA = $iN

	; set char buffers and input healing
	DllStructSetData($tBLASCHAR1, 1, $cUPLO = "U" ? "U" : "L")

	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "lauum", _
		"PTR",  $pBLASCHAR1, _  ; UPLO
		"INT*", $iN, _          ; N
		"PTR",  $pA, _          ; A
		"INT*", $iLDA, _        ; lda
		"INT*", 0 _ 			; INFO
	)
	If @error Then Return SetError(1, @error, False)
	Return $aDLL[5] = 0 ? True : SetError(2, $aDLL[5], False)

EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_lauu2()
; Description ...: computes the product U * Uᵀ or Lᵀ * L, where U or L are triangular matrices (unblocked variant)
; Syntax ........: _lp_lauu2($mA, [$cUPLO = "U", [$iN = Default, [$iLDA = $iN, [$sDataType = "DOUBLE"]]]])
; Parameters ....: mA        - [Map] matrix A (LDA × N) as a map, DllStruct or pointer (gets overwritten)
;                  cUPLO     - [Char] (Default: "U")
;                            ↳ "U": upper triangle of A is used
;                              "L": lower triangle of A is used
;                  iN        - [Int] (Default: Default)
;                            ↳ order of the triangular factor U or L (rows)
;                  iLDA      - [Int] (Default: $iN)
;                            ↳ leading dimension of the array (rows)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during first DllCall of lauu2 (@extended: @error from DllCall)
;                           | 2: error inside first call of lauu2 (@extended: INFO-value from lauu2)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d3/d8f/group__lauu2_ga52dc4db093f8513ca2f74195f3aaf689.html#ga52dc4db093f8513ca2f74195f3aaf689
; Example .......: Yes
;                  Global $mU = _blas_fromArray("[[1,2,3],[0,4,5],[0,0,6]]")
;                  _lp_lauu2($mU)
;                  _blas_display($mU)
; ===============================================================================================================================
Func _lp_lauu2($mA, $cUPLO = "U", $iN = Default, $iLDA = $iN, $sDataType = "DOUBLE")
	Local $pA ; pointer to the data in memory

	; Triangle * Triangleᵀ = symmetric

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iN) = 1 Then $iN = $mA.rows
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	If IsKeyword($iLDA) = 1 Then $iLDA = $iN

	; set char buffers and input healing
	DllStructSetData($tBLASCHAR1, 1, $cUPLO = "U" ? "U" : "L")

	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "lauu2", _
		"PTR",  $pBLASCHAR1, _  ; UPLO
		"INT*", $iN, _          ; N
		"PTR",  $pA, _          ; A
		"INT*", $iLDA, _        ; lda
		"INT*", 0 _ 			; INFO
	)
	If @error Then Return SetError(1, @error, False)
	Return $aDLL[5] = 0 ? True : SetError(2, $aDLL[5], False)

EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _lp_lamch()
; Description ...: determines float/double precision machine parameters
; Syntax ........: _lp_lamch($cMach, [$sDataType = "DOUBLE"])
; Parameters ....: cMach     - [Char] specifies the value to be returned
;                              "E": eps - relative machine precision
;                              "S": sfmin - safe minimum, such that 1/sfmin does not overflow
;                              "B": base - base of the machine
;                              "P": prec = eps*base
;                              "N": t - number of (base) digits in the mantissa
;                              "R": rnd - 1.0 when rounding occurs in addition, 0.0 otherwise
;                              "M": emin - minimum exponent before (gradual) underflow
;                              "U": rmin - underflow threshold - base**(emin-1)
;                              "L": emax - largest exponent before overflow
;                              "O": rmax - overflow threshold  - (base**emax)*(1-eps)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: [Float] the value
;                  Failure: Null and set @error to:
;                           | 1: error during first DllCall of lamch (@extended: @error from DllCall)
;                           | 2: error inside first call of lamch (@extended: INFO-value from lamch)
; Author ........: AspirinJunkie
; Modified.......: 2024-09-02
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d4/d86/group__lamch_gaeab255e77cbd3b0f31aea74ed0ce099e.html#gaeab255e77cbd3b0f31aea74ed0ce099e
; Example .......: Yes
;                  ConsoleWrite(_lp_lamch("p", "FLOAT") & @CRLF)
; ===============================================================================================================================
Func _lp_lamch($cMach, $sDataType = "DOUBLE")
	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	If Not StringRegExp($cMach, '(?i)^[ESBPNRMULO]$') Then Return SetError(1, @error, Null)

	; set char buffer
	DllStructSetData($tBLASCHAR1, 1, $cMach)

	Local $aDLL = DllCall($__g_hBLAS_DLL, $sDataType & ":cdecl", $cPrefix & "lamch", _
	    "PTR",  $pBLASCHAR1 _  ; CMACH
	)
	If @error Then Return SetError(2, @error, Null)
	Return $aDLL[0]
EndFunc

#EndRegion
