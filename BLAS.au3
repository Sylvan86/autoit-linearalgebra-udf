
;~ #AutoIt3Wrapper_Run_AU3Check=Y
;~ #AutoIt3Wrapper_Au3Check_Parameters=-d -w 1 -w 2 -w 3 -w 4 -w 5 -w 6 -w 7
;~ #AutoIt3Wrapper_AU3Check_Stop_OnWarning=Y
;~ Opt("MustDeclareVars", 1)

#include-once

#include <WinAPIError.au3> ; for _WinAPI_GetLastError()
#include <Array.au3> ; for _ArrayDisplay()

; #INDEX# =======================================================================================================================
; Title .........: BLAS
; AutoIt Version : 3.3.16.1
; Description ...: Wrapper UDF for interaction with BLAS-compatible libraries (dlls)
;                  Provides basic functionalities for linear algebra.
; Author(s) .....: AspirinJunkie
; Dll ...........: ntdll.dll, the user defined BLAS-DLL
; ===============================================================================================================================

; #CURRENT# =====================================================================================================================
; ---- Library administration ----
; _blas_LoadBlasDll  - loads the DLL with the BLAS-compatible interface
;
; ---- AutoIt-BLAS interface ----
; _blas_createVector - creates new empty vector
; _blas_createMatrix - creates new empty matrix
; _blas_toArray      - converts a BLAS matrix map into an AutoIt array
; _blas_fromArray    - converts a AutoIt array or array define string into BLAS matrix map
; _blas_fromStruct   - creates a matrix/vector map from a DllStruct as used here in the UDF.
; _blas_duplicate    - creates an independent copy of a matrix/vector map
; _blas_display      - displays a matrix/vector map, similar to _ArrayDisplay
;
; ---- BLAS Level 1 (Vector Operations) ----
; _blas_axpy         - calculate a * X + Y and store it in Y (if a = 1 then X + Y)
; _blas_dot          - calculate the dot product of two vectors X and Y
; _blas_rot          - applies a plane rotation to coordinate-pairs
; _blas_rotg         - constructs a plane rotation
; _blas_swap         - interchanges two vectors (useful for matrix manipulations)
; _blas_copy         - copies a vector x to a vector y (useful for extracting parts of a matrix to other areas)
; _blas_scal         - scales a vector with a scalar (X = a * X)
; _blas_nrm2         - calculate the euclidean norm of a vector
; _blas_asum         - calculate the sum of the absolute(!) values of a matrix/vector
; _blas_amax         - finds the first element having the maximum absolute(!) value
; _blas_amin         - finds the first element having the minimum absolute(!) value
;
; ---- BLAS Level 2 (Matrix-Vector Operations) ----
; _blas_gemv         - matrix-vector multiplication: y := alpha*A*x + beta*y, or y := alpha*Aᵀ*x + beta*y where A is a general unpacked matrix
; _blas_trmv         - matrix-vector multiplication: x := A*x, or x := Aᵀ*x where A is a unpacked upper or lower triangular matrix
; _blas_symv         - matrix-vector multiplication: y := alpha*A*x + beta*y where A is an unpacked symmetric matrix
; _blas_gbmv         - matrix-vector multiplication: y := alpha*A*x + beta*y, or y := alpha*A**T*x + beta*y where A is an packed(!) banded matrix
; _blas_sbmv         - matrix-vector multiplication: alpha*A*x + beta*y --> y where A is an packed(!) symmetric band matrix
; _blas_tbsv         - solves one of the systems of equations: A * x = b, or Aᵀ * x = b Where A is an packed(!) lower or upper triangular band matrix
; _blas_trsv         - solves one of the systems of equations: A * x = b, or Aᵀ * x = b Where A is an unpacked lower or upper triangular matrix
; _blas_ger          - calculate the rank 1 operation: A := alpha * x * yᵀ + A
;
; ---- BLAS Level 3 (Matrix-Matrix Operations) ----
; _blas_gemm         - matrix multiplication: C := alpha * op(A) * op(B) + beta * C
; _blas_symm         - matrix multiplication: C := alpha * A * B + beta * C or C := alpha * B * A + beta * C where A is an unpacked symmetric matrix, B and C are general m × n matrices
; _blas_trmm         - matrix multiplication: B := alpha * op(A) * B or B := alpha * B * op(A) where op(A) = A or Aᵀ where A is an unpacked lower or upper triangular matrix, B is a general m × n matrix
; _blas_trsm         - solves one of the matrix equations: op(A) * X = alpha * B, or X * op(A) = alpha * B where op(A) = A or Aᵀ where A is an unpacked lower or upper triangular matrix, B is a general m × n matrix
; _blas_syrk         - symmetric rank k operations: C := alpha * A * Aᵀ + beta * C or alpha * Aᵀ * A + beta * C where A is an unpacked symmetric matrix, C is a general m × n matrix
;
; ---- additional functions (depends on library) ----
; _blas_tpttr        - copies a triangular matrix from the standard packed format (TP) to the standard full format (TR)
; _blas_trttp        - copies a triangular matrix from the standard full format (TR) to the standard packed format (TP)
; _blas_imatcopy     - performs scaling and in-place transposition/copying of matrices
; _blas_omatcopy     - performs scaling and out-place transposition/copying of matrices.
; ===============================================================================================================================

; #INTERNAL_USE_ONLY# ===========================================================================================================
; ---- Library administration ----
; __blas_error           - shows an error window with relevant information and terminates the script
;
; ---- helper functions ----
; __blas_ArrayFromString - creates an array from an array definition in AutoIt syntax, which is passed as a string
; __blas_fillWithScalar  - fills a matrix, a vector or parts thereof with a specific value
; __blas_GBSfromArray    - converts an AutoIt array into a banded matrix, which is stored in "General-Band Storage Mode"
; ===============================================================================================================================


; #VARIABLES# ===================================================================================================================
; DLL-Handle to the BLAS-Library (no Const because User can link to his own BLAS-DLL)
If Not IsDeclared("__g_hBLAS_DLL") Then Global $__g_hBLAS_DLL = _blas_LoadBlasDll()
; ===============================================================================================================================

; #CONSTANTS# ===================================================================================================================
; system depended sizes in Bytes for the supported data types
Global Const $iBLAS_SIZE_DOUBLE = 8
Global Const $iBLAS_SIZE_FLOAT = 4

; constants used as flags to describe the internal storage type of the matrix
Global Enum Step *2 $__g_BLAS_STYPE_MATRIX = 1, $__g_BLAS_STYPE_TRIANGLE, $__g_BLAS_STYPE_SYMMETRIC, $__g_BLAS_STYPE_BAND, _
                    $__g_BLAS_STYPE_TRIDIAGONAL, $__g_BLAS_STYPE_LOWER, $__g_BLAS_STYPE_UPPER, $__g_BLAS_STYPE_NONUNIT, _
					$__g_BLAS_STYPE_PACKED, $__g_BLAS_STYPE_POSITIVE_DEFINITE, $__g_BLAS_STYPE_DIAGONAL

;~ Single-char structures: In principle, "str" also works in DllCall, but this is inefficient because 65,536 characters are always allocated. These structures are therefore declared once in order to prevent continuous re-generation.
Global Const $tBLASCHAR1 = DllStructCreate("CHAR"), $pBLASCHAR1 = DllStructGetPtr($tBLASCHAR1)
Global Const $tBLASCHAR2 = DllStructCreate("CHAR"), $pBLASCHAR2 = DllStructGetPtr($tBLASCHAR2)
Global Const $tBLASCHAR3 = DllStructCreate("CHAR"), $pBLASCHAR3 = DllStructGetPtr($tBLASCHAR3)
Global Const $tBLASCHAR4 = DllStructCreate("CHAR"), $pBLASCHAR4 = DllStructGetPtr($tBLASCHAR4)
; ===============================================================================================================================


; #REMARKS# =====================================================================================================================
; Principles of this UDF level:
; BLAS functions are low-level functions that are primarily not used by the end user but by higher abstracted functions.
; Their primary task is the high-performance connection to the BLAS library.
; For this reason, condition checks are largely dispensed with and this task is left to the next UDF levels.
; The functions are also largely offered as they are. As a rule, this means that input data is often overwritten.
; If overwriting is to be prevented, it is therefore the task of the higher UDF levels to duplicate the data accordingly beforehand.
;
; Maps are used as the main exchange object, which hold the elements of the matrices/vectors in a DllStruct, as well as additional meta information in other elements.
; The specific structure is documented here:
;
; ---------------------- Matrix object structure -------------------------------
; .elements    : [int] number of elements in the matrix (rows * cols)
; .size        : [int] number of elements in the storage representation (for unpacked: = .elements)
; .rows        : [int] number of matrix/vector rows
; .cols        : [int] number of matrix rows (vector: = 1 or 0)
; .storageType : [int] flags describing the storage type (vector = 0); see $__g_BLAS_STYPE_xxx constants
; .datatype    : [string] the element datatype ("DOUBLE" Or "FLOAT")
; .struct      : [struct] the AutoIt DllStruct structure
; .ptr         : [ptr] pointer to the DllStruct structure
; .kl          : [int] number of sub diagonals if .storageType = band matrix
; .ku          : [int] number of super diagonals if .storageType = band matrix
; ===============================================================================================================================


#Region Library administration

; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_LoadBlasDll()
; Description ...: loads the DLL with the BLAS-compatible interface
; Syntax ........: _blas_LoadBlasDll([$sDllPath = @ScriptDir & "libopenblas.dll"])
; Parameters ....: sDllPath - [String] (Default: @ScriptDir & "libopenblas.dll")
;                           ↳ path to the dll file
; Return value ..: Success: Dll handle
;                  Failure: $sDllPath and set @error to:
;                           | 1: file in path not exist
;                           | 2: error during DllOpen (@extended: @error from DllOpen)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......: The user can change the BLAS/LAPACK DLL used by setting the following BEFORE(!) the #include:
;                  Global $__g_hBLAS_DLL = DllOpen(...)
; Related .......:
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func _blas_LoadBlasDll($sDllPath = @ScriptDir & "\libopenblas.dll")
	If Not FileExists($sDllPath) Then Return SetError(1, 0, $sDllPath)

	; necessary because some dlls need to find other module in other dll-files
	EnvSet("PATH", EnvGet("PATH") & ";" & $sDllPath)

	Local $hDLL = DllOpen($sDllPath)
	If $hDLL = -1 Then Return SetError(2, _WinAPI_GetLastError(), _WinAPI_GetLastErrorMessage())

	Return $hDLL
EndFunc   ;==>_blas_LoadBlasDll

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __blas_error()
; Description ...: shows an error window with relevant information and terminates the script
; Syntax ........: __blas_error($vRet, [$sText = "", [$vError = @error, [$iExtended = @extended, [$iLine = @ScriptLineNumber]]]])
; Parameters ....: vRet      - [Variant] return value of the function that generated the error
;                  sText     - [String] (Default: "")
;                            ↳ Text to be displayed in the error window
;                  vError    - [Variant](Default: @error)
;                            ↳ @error value of the function that caused the error
;                  iExtended - [Int] (Default: @extended)
;                            ↳ @extended value of the function that caused the error
;                  iLine     - [Int] (Default: @ScriptLineNumber)
;                            ↳ line number in which the error occurred
; Return value ..: -
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func __blas_error($vRet, $sText = "", $vError = @error, $iExtended = @extended, $iLine = @ScriptLineNumber)
	MsgBox(48, "error", "error during " & $sText & @CRLF & "@error: " & $vError & @CRLF & "@extended: " & $iExtended & @CRLF & "return value: " & $vRet & @CRLF & "line: " & $iLine & @CRLF & "_WinAPI_GetLastError(): " & _WinAPI_GetLastError() & @CRLF & "_WinAPI_GetLastErrorMessage(): " & _WinAPI_GetLastErrorMessage())
	Exit
EndFunc   ;==>__blas_error

#EndRegion

#Region AutoIt-BLAS interface

; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_createVector()
; Description ...: creates new empty vector
; Syntax ........: _blas_createVector($iN, [$sType = "DOUBLE", [$pTarget = Default]])
; Parameters ....: iN      - [Int] number of elements in the vector
;                  sType   - [String] (Default: "DOUBLE")
;                          ↳ data type of the elements ("DOUBLE" or "FLOAT")
;                  pTarget - [Int] (Default: Default)
;                          ↳ pointer to a memory area in which the vector is to be created
;                            if Default a new memory area is reserved
; Return value ..: Success: [Map] BLAS Matrix Map as it is structured in this UDF
;                  Failure: Null and set @error to:
;                           | 1: error during DllStructCreate (@extended: @error from DllStructCreate)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: Yes
;                  Global $mVec = _blas_createVector(10)
;                  _blas_display($mVec, "new empty vector")
; ===============================================================================================================================
Func _blas_createVector(Const $iN, Const $sType = "DOUBLE", Const $pTarget = Default)
	Local $tStruct = IsKeyword($pTarget) = 1 ? DllStructCreate(StringFormat("%s[%d]", $sType, $iN)) : DllStructCreate(StringFormat("%s[%d]", $sType, $iN), $pTarget)
	If @error Then Return SetError(1, @error, Null)

	Local $mRet[]
	$mRet.struct      = $tStruct
	$mRet.ptr         = DllStructGetPtr($tStruct)
	$mRet.elements    = $iN
	$mRet.size        = $iN
	$mRet.rows        = $iN
	$mRet.cols        = 0
	$mRet.storageType = 0
	$mRet.datatype    = $sType

	Return $mRet
EndFunc   ;==>_blas_createVector


; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_createMatrix()
; Description ...: creates new empty matrix
; Syntax ........: _blas_createMatrix($iR, [$iC = $iR, [$sType = "DOUBLE", [$nMode = $__g_BLAS_STYPE_MATRIX]]])
; Parameters ....: iR    - [Int] number of matrix rows
;                  iC    - [Int] (Default: $iR)
;                        ↳ number of matrix columns
;                  sType - [String] (Default: "DOUBLE")
;                        ↳ data type of the elements ("DOUBLE" or "FLOAT")
;                  nMode - (Default: $__g_BLAS_STYPE_MATRIX)
;                        ↳ BLAS memory layout for the matrix (see $__g_BLAS_STYPE_XXX flags)
; Return value ..: Success: [Map] BLAS Matrix Map as it is structured in this UDF
;                  Failure: Null and set @error to:
;                           | 1: error during DllStructCreate (@extended: @error from DllStructCreate)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: Yes
;                  Global $mMat = _blas_createMatrix(5, 3)
;                  _blas_display($mMat, "new empty matrix")
; ===============================================================================================================================
Func _blas_createMatrix(Const $iR, Const $iC = $iR, Const $sType = "DOUBLE", $nMode = $__g_BLAS_STYPE_MATRIX)
	Local $tStruct = DllStructCreate(StringFormat("%s[%d]", $sType, $iR * $iC))
	If @error Then Return SetError(1, @error, Null)

	Local $mRet[]
	$mRet.struct      = $tStruct
	$mRet.ptr         = DllStructGetPtr($tStruct)
	$mRet.elements    = $iR * $iC
	$mRet.size        = $mRet.elements
	$mRet.rows        = $iR
	$mRet.cols        = $iC
	$mRet.storageType = $nMode
	$mRet.datatype    = $sType

	Return $mRet
EndFunc   ;==>_la_createMatrix


; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_toArray()
; Description ...: converts a BLAS matrix map into an AutoIt array
; Syntax ........: _blas_toArray($mMatrix)
; Parameters ....: mMatrix - [Map] the BLAS Matrix Map as it is structured in this UDF
; Return value ..: Success: matrix/vector as AutoIt 2D/1D arraySetExtended($iN, $aV)
;                  Failure: Null and set @error to:
;                           | 1: Memory layout is currently not yet supported
;                           | 2: invalid combination in the layout flags
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.ibm.com/docs/en/essl/6.3?topic=matrices-matrix-storage-representation
; Example .......: Yes
;                  Global $mTest = _blas_fromArray("[[1.3, 1.2, 1.1, 0, 0],[0, 2.3, 2.2, 2.1, 0],[0, 0, 3.3, 3.2, 3.1],[0,0,0,4.3,4.2],[0,0,0,0,5.3]]")
;                  _blas_display($mTest, "unpacked")
;                  Global $mA = _blas_fromArray("[[1.3, 1.2, 1.1, 0, 0],[0, 2.3, 2.2, 2.1, 0],[0, 0, 3.3, 3.2, 3.1],[0,0,0,4.3,4.2],[0,0,0,0,5.3]]", $__g_BLAS_STYPE_BAND + $__g_BLAS_STYPE_PACKED + $__g_BLAS_STYPE_UPPER, "DOUBLE", 2)
;                  _blas_display($mA, "packed upper band matrix")
; ===============================================================================================================================
Func _blas_toArray(Const ByRef $mMatrix)
	Local $iRows = $mMatrix.rows, $iCols = $mMatrix.cols
	Local $tStruct = $mMatrix.struct
	Local $nMode = $mMatrix.storageType

	If $mMatrix.storageType = 0 Then ; Vector
		Local $iN = $mMatrix.elements
		Local $aV[$iN]
		For $i = 0 To $iN - 1
			$aV[$i] = DllStructGetData($tStruct, 1, $i + 1)
		Next

		Return SetExtended($iN, $aV)

	EndIf

	Local $aRet[$iRows][$iCols], _
	      $iX = 1, _
		  $iTmp, $bTmp, _
		  $iKU, $iKL, _
		  $iMin, $iMax, _
		  $iLDA


	If BitAND($nMode, $__g_BLAS_STYPE_PACKED) Then
	; packed matrices
		Select
			Case BitAND($nMode, $__g_BLAS_STYPE_BAND)
				$bTmp = BitAND($nMode, $__g_BLAS_STYPE_SYMMETRIC) <> 0

				$iKL  = $mMatrix.kl
				$iKU  = $mMatrix.ku
				$iLDA = $iKL + $iKU + 1	; size of rows in packed band matrix array ( A(m,n) --> ASB(LDA,n) )

				For $j = 0 To $iCols - 1
					$iMin = ($iRows - 1) > ($j + $iKL) ? $j + $iKL : $iRows - 1
					$iMax = 0 > ($j - $iKU) ? 0 : $j - $iKU
				    For $i = $iMax To $iMin
						$aRet[$i][$j] = DllStructGetData($tStruct, 1, ($iKU - $j + $i) + $j * $iLDA + 1)
						If $bTmp And $i <> $j And $i < $iCols And $j < $iRows Then $aRet[$j][$i] = $aRet[$i][$j] ; for symmetric matrices
				    Next
				Next

			Case BitAND($nMode, $__g_BLAS_STYPE_TRIDIAGONAL)
			; not implemented yet
				Return SetError(2,0, Null)

			Case BitAND($nMode, $__g_BLAS_STYPE_TRIANGLE + $__g_BLAS_STYPE_SYMMETRIC) ; triangle and symmetric matrices are stored the same way - only the flag differs

				$bTmp = BitAND($nMode, $__g_BLAS_STYPE_SYMMETRIC) <> 0

				If BitAND($nMode, $__g_BLAS_STYPE_LOWER) Then
				; triangular: https://www.ibm.com/docs/en/essl/6.3?topic=representation-lower-triangular-packed-storage-mode
				; symmetric:  https://www.ibm.com/docs/en/essl/6.3?topic=representation-lower-packed-storage-mode

					For $j = 0 To $iCols - 1
						For $i = $j To $iRows - 1
							$aRet[$i][$j] = DllStructGetData($tStruct, 1, $iX)
							If $bTmp And $i <> $j And $i < $iCols And $j < $iRows Then $aRet[$j][$i] = $aRet[$i][$j] ; for symmetric matrices
							$iX += 1
						Next
					Next

				Else ; = $__g_BLAS_STYPE_UPPER (or default)
				; triangular: https://www.ibm.com/docs/en/essl/6.3?topic=representation-upper-triangular-packed-storage-mode
				; symmetric:  https://www.ibm.com/docs/en/essl/6.3?topic=representation-upper-packed-storage-mode

					For $j = 0 To $iCols - 1
						For $i = 0 To ($j >= $iRows ? $iRows - 1 : $j)
							$aRet[$i][$j] = DllStructGetData($tStruct, 1, $iX)
							If $bTmp And $i <> $j And $i < $iCols And $j < $iRows Then $aRet[$j][$i] = $aRet[$i][$j] ; for symmetric matrices
							$iX += 1
						Next
					Next

				EndIf

			Case BitAND($nMode, $__g_BLAS_STYPE_DIAGONAL) ; diagonal matrix as vector
				$iN = $iRows < $iCols ? $iRows : $iCols
				For $i = 0 To $iN - 1
					$aRet[$i][$i] = DllStructGetData($tStruct, 1, $i + 1)
				Next

			Case Else
			; there is no storage representation for packed matrices which are not symmetric, triangular, tridiagonal or banded
				Return SetError(1, $nMode, Null)

		EndSelect

	Else ; unpacked matrices
	; https://www.ibm.com/docs/en/essl/6.3?topic=matrices-matrix-storage-representation

		Select
			Case BitAND($nMode, $__g_BLAS_STYPE_BAND)
				$iKL = $mMatrix.kl
				$iKU = $mMatrix.ku

				For $j = 0 To $iCols - 1
					$iMin = ($iRows - 1) > ($j + $iKL) ? $j + $iKL : $iRows - 1
					$iMax = 0 > ($j - $iKU) ? 0 : $j - $iKU
					For $i = $iMax To $iMin
						$aRet[$i][$j] = DllStructGetData($tStruct, 1, $iX + $i)
					Next
					$iX += $iRows
				Next

			Case BitAND($nMode, $__g_BLAS_STYPE_TRIDIAGONAL)
			; not implemented yet
				Return SetError(2,0,Null)

			Case BitAND($nMode, $__g_BLAS_STYPE_TRIANGLE + $__g_BLAS_STYPE_SYMMETRIC) ; triangle and symmetric matrices are stored the same way - only the flag differs
			; symmetric:  https://www.ibm.com/docs/en/essl/6.3?topic=representation-lower-storage-mode
			; triangular: https://www.ibm.com/docs/en/essl/6.3?topic=representation-lower-triangular-storage-mode

				$bTmp = BitAND($nMode, $__g_BLAS_STYPE_SYMMETRIC) <> 0

				If BitAND($nMode, $__g_BLAS_STYPE_LOWER) Then

					; with or without main diagonal
					$iTmp = BitAND($nMode, $__g_BLAS_STYPE_NONUNIT) ? 1 : 0

					; read lower left triangle matrix in column-major order (Fortran-order)
					For $j = 0 To $iCols - 1
						$iX += $j + $iTmp
						For $i = $j + $iTmp To $iRows - 1
							$aRet[$i][$j] = DllStructGetData($tStruct, 1, $iX)
							If $bTmp And $i <> $j And $i < $iCols And $j < $iRows Then $aRet[$j][$i] = $aRet[$i][$j] ; for symmetric matrices
							$iX += 1
						Next
					Next

				ElseIf BitAND($nMode, $__g_BLAS_STYPE_UPPER) Then
				; symmetric:  https://www.ibm.com/docs/en/essl/6.3?topic=representation-upper-storage-mode#am5gr_utsm
				; triangular: https://www.ibm.com/docs/en/essl/6.3?topic=representation-upper-triangular-storage-mode

					; with or without main diagonal
					$iTmp = BitAND($nMode, $__g_BLAS_STYPE_NONUNIT) ? 1 : 0

					; read upper right triangle matrix in column-major order (Fortran-order)
					For $j = 0 To $iCols - 1
						For $i = 0 To (($j - $iTmp) >= $iRows ? $iRows - 1 : $j - $iTmp)
							$aRet[$i][$j] = DllStructGetData($tStruct, 1, $iX)
							If $bTmp And $i <> $j And $i < $iCols And $j < $iRows Then $aRet[$j][$i] = $aRet[$i][$j] ; for symmetric matrices
							$iX += 1
						Next
						$iX += $iRows - $j - (1 - $iTmp)
					Next

				Else
					ContinueCase ; handle as full matrix if no upper/lower defined

				EndIf

			Case Else
			; full unpacked representation of a matrix - process all elements
			; https://www.ibm.com/docs/en/essl/6.3?topic=matrices-matrix-storage-representation

				; read matrix in column-major order (Fortran-order)
				For $j = 0 To $iCols - 1
					For $i = 0 To $iRows - 1
						$aRet[$i][$j] = DllStructGetData($tStruct, 1, $iX)
						$iX += 1
					Next
				Next

		EndSelect
	EndIf

	Return $aRet

EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_fromArray()
; Description ...: converts a AutoIt array or array define string into BLAS matrix map
; Syntax ........: _blas_fromArray($aArray, [$nMode = 0, [$sType = "DOUBLE", [$iKL = 0, [$iKU = 0]]]])
; Parameters ....: aArray - [Array/String] AutoIt array or array define string (see examples) which should be converted into BLAS map
;                  nMode  - (Default: 0)
;                         ↳ BLAS memory layout for the matrix (see $__g_BLAS_STYPE_XXX flags)
;                  sType  - [String] (Default: "DOUBLE")
;                         ↳ data type of the elements ("DOUBLE" or "FLOAT")
;                  iKL    - [Int] (Default: 0)
;                         ↳ If nMode contains $__g_BLAS_STYPE_BAND: number of lower band diagonals
;                  iKU    - [Int] (Default: 0)
;                         ↳ If nMode contains $__g_BLAS_STYPE_BAND: number of upper band diagonals
; Return value ..: Success: [Map] BLAS Matrix Map as it is structured in this UDF
;                  Failure: Null and set @error to:
;                           | 1: invalid value for $sType
;                           | 2: invalid number of dimensions in $aArray (@extended: number of dimensions of $aArray)
;                           | 3: error during __blas_ArrayFromString() (@extended: @error from __blas_ArrayFromString)
;                           | 4: error during DllStructCreate() (@extended: @error from DllStructCreate)
;                           | 5: invalid combination in the layout flags
;                           | 6: Memory layout is currently not yet supported
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.ibm.com/docs/en/essl/6.3?topic=matrices-matrix-storage-representation
; Example .......: Yes
;                  Global $mMatrix = _blas_fromArray("[[4,-1,1],[-1,4,-2],[1,-2,5]]")
;                  _blas_display($mMatrix, "general Matrix created from string definition")
;
;                  ;upper symmetric band matrix
;                  Global $mMatrix = _blas_fromArray("[[11,12,13,14,0,0],[12,22,23,24,25,0],[13,23,33,34,35,36],[14,24,34,44,45,46],[0,25,35,45,55,56],[0,0,36,46,56,66]]", $__g_BLAS_STYPE_SYMMETRIC + $__g_BLAS_STYPE_BAND + $__g_BLAS_STYPE_UPPER, "DOUBLE", 3)
;                  _blas_display($mMatrix, "symmetric upper band matrix")
; ===============================================================================================================================
;~ [[19,-80,-55,0,0],[29,-92,-67,-94,0],[0,-29,77,-7,40],[0,0,-70,12,97],[0,0,0,53,43]]
Func _blas_fromArray($aArray, $nMode = 0, Const $sType = "DOUBLE", $iKL = 0, $iKU = 0)
	If $sType <> "DOUBLE" And $sType <> "FLOAT" Then Return SetError(1, 0, Null)

	; you can define the array directly as a string
	If IsString($aArray) Then $aArray = __blas_ArrayFromString($aArray)
	If @error Then Return SetError(3, @error, Null)

	; only 1D (Vector) or 2D (Matrix) arrays allowed
	If UBound($aArray, 0) > 2 Or UBound($aArray, 0) < 1 Then Return SetError(2, UBound($aArray, 0), Null) 

	; variables used in this function
	Local $tStruct, _
	      $mRet[], _             ; the return object/map structure
	      $iRows, $iCols, $iN, $iLDA, _
		  $iX = 1, _                ; current position in struct
		  $iMin, $iMax

	; Vector
	If UBound($aArray, 0) = 1 Then
		$iN = UBound($aArray, 1)
		$tStruct = DllStructCreate(StringFormat("%s[%d]", $sType, $iN))
		If @error Then Return SetError(4, @error, Null)

		For $i = 0 To $iN - 1
			DllStructSetData($tStruct, 1, $aArray[$i], $i + 1)
		Next

		$mRet.elements    = $iN
		$mRet.size        = $iN
		$mRet.rows        = $iN
		$mRet.cols        = 1
		$mRet.storageType = 0
		$mRet.datatype    = $sType
		$mRet.struct      = $tStruct
		$mRet.ptr         = DllStructGetPtr($tStruct)

		Return SetExtended($iN, $mRet)
	EndIf

	$iRows = UBound($aArray, 1)
	$iCols = UBound($aArray, 2)

	$nMode += $__g_BLAS_STYPE_MATRIX
	$mRet.elements    = $iRows * $iCols
	$mRet.rows        = $iRows
	$mRet.cols        = $iCols
	$mRet.datatype    = $sType

	If BitAND($nMode, $__g_BLAS_STYPE_PACKED) Then
	; packed matrices
		Select
			Case BitAND($nMode, $__g_BLAS_STYPE_BAND)
				If BitAND($nMode, $__g_BLAS_STYPE_SYMMETRIC) Then
				; upper: https://www.ibm.com/docs/en/essl/6.3?topic=representation-upper-band-packed-storage-mode
				; lower: https://www.ibm.com/docs/en/essl/6.3?topic=representation-lower-band-packed-storage-mode
					; no extra case for symmetric band matrix because it`s the same algorithm as for lower/upper

					; set upper as default for symmetric band matrix if user not set lower/upper
					If Not (BitAND($nMode, $__g_BLAS_STYPE_UPPER) Or BitAND($nMode, $__g_BLAS_STYPE_LOWER)) Then $nMode = BitOR($nMode, $__g_BLAS_STYPE_UPPER)
				EndIf

				If BitAND($nMode, $__g_BLAS_STYPE_UPPER) Then
				; upper packed band matrix (sub diagonals = $iKL)
				; https://www.ibm.com/docs/en/essl/6.3?topic=representation-upper-triangular-band-packed-storage-mode
					If $iKU = 0 Then $iKU = $iKL
					$iKL = 0
					; we use the same algorithm for general band and lower/upper packed band matrices

				ElseIf BitAND($nMode, $__g_BLAS_STYPE_LOWER) Then
				; lower packed band matrix (sub diagonals = $iKL)
				; https://www.ibm.com/docs/en/essl/6.3?topic=representation-lower-triangular-band-packed-storage-mode

					$iKU = 0
					; we use the same algorithm for general band and lower/upper packed band matrices

				EndIf

				; general packed band matrix (sub diagonals = $iKL, super diagonals = $iKL)
				; https://www.ibm.com/docs/en/essl/6.3?topic=representation-blas-general-band-storage-mode

				$iLDA      = $iKL + $iKU + 1	; size of rows in packed band matrix array ( A(m,n) --> ASB(LDA,n) )
				$mRet.size = $iLDA * $iCols
				$mRet.kl = $iKL
				$mRet.ku = $iKU

				$tStruct = DllStructCreate(StringFormat("%s[%d]", $sType, $mRet.size))
				If @error Then Return SetError(4, @error, Null)

				For $j = 0 To $iCols - 1
					$iMin = ($iRows - 1) > ($j + $iKL) ? $j + $iKL : $iRows - 1
					$iMax = 0 > ($j - $iKU) ? 0 : $j - $iKU
				    For $i = $iMax To $iMin
						DllStructSetData($tStruct, 1, $aArray[$i][$j], ($iKU - $j + $i) + $j * $iLDA + 1)
				    Next
				Next

			Case BitAND($nMode, $__g_BLAS_STYPE_TRIDIAGONAL)
			; not implemented yet
				Return SetError(6,0, Null)

			Case BitAND($nMode, $__g_BLAS_STYPE_TRIANGLE + $__g_BLAS_STYPE_SYMMETRIC) ; triangle and symmetric matrices are stored the same way - only the flag differs

				If BitAND($nMode, $__g_BLAS_STYPE_LOWER) Then
				; triangular: https://www.ibm.com/docs/en/essl/6.3?topic=representation-lower-triangular-packed-storage-mode
				; symmetric:  https://www.ibm.com/docs/en/essl/6.3?topic=representation-lower-packed-storage-mode

					$iN = $iRows <= $iCols ? $iRows * ($iRows + 1) / 2 : $iRows * $iCols - ($iCols * ($iCols-1) ) / 2
					$mRet.size = $iN
					$tStruct   = DllStructCreate(StringFormat("%s[%d]", $sType, $iN))
					If @error Then Return SetError(4, @error, Null)

					For $j = 0 To $iCols - 1
						For $i = $j To $iRows - 1
							DllStructSetData($tStruct, 1, $aArray[$i][$j], $iX)
							$iX += 1
						Next
					Next

				Else ; = $__g_BLAS_STYPE_UPPER (or default)
				; triangular: https://www.ibm.com/docs/en/essl/6.3?topic=representation-upper-triangular-packed-storage-mode
				; symmetric:  https://www.ibm.com/docs/en/essl/6.3?topic=representation-upper-packed-storage-mode

					$iN = $iCols <= $iRows ? $iCols * ($iCols + 1) / 2 : $iCols * $iRows - ($iRows * ($iRows-1) ) / 2
					$mRet.size = $iN
					$tStruct   = DllStructCreate(StringFormat("%s[%d]", $sType, $iN))
					If @error Then Return SetError(4, @error, Null)

					For $j = 0 To $iCols - 1
						For $i = 0 To ($j >= $iRows ? $iRows - 1 : $j)
							DllStructSetData($tStruct, 1, $aArray[$i][$j], $iX)
							$iX += 1
						Next
					Next

				EndIf

			Case BitAND($nMode, $__g_BLAS_STYPE_DIAGONAL) ; diagonal vector as matrix
				$iN = $iRows < $iCols ? $iRows : $iCols
				$tStruct = DllStructCreate(StringFormat("%s[%d]", $sType, $iN))
				For $i = 0 To $iN - 1
					DllStructSetData($tStruct, 1, $aArray[$i][$i], $i + 1)
				Next

			Case Else
			; there is no storage representation for packed matrices which are not symmetric, triangular, tridiagonal or banded
			Return SetError(5, $nMode, Null)

		EndSelect

	Else ; unpacked matrices
	; https://www.ibm.com/docs/en/essl/6.3?topic=matrices-matrix-storage-representation
		$iN        = $iRows * $iCols
		$mRet.size = $iN
		$tStruct   = DllStructCreate(StringFormat("%s[%d]", $sType, $iN))
		If @error Then Return SetError(4, @error, Null)

		Select
			Case BitAND($nMode, $__g_BLAS_STYPE_BAND)
				If BitAND($nMode, $__g_BLAS_STYPE_UPPER) Then
					If $iKU = 0 Then $iKU = $iKL
					$iKL = 0
				EndIf
				If BitAND($nMode, $__g_BLAS_STYPE_LOWER) Then $iKU = 0
				$mRet.kl = BitAND($nMode, $__g_BLAS_STYPE_UPPER) ? 0 : $iKL
				$mRet.ku = $iKU

				; iterate over relevant elements only
				For $j = 0 To $iCols - 1
					$iMin = ($iRows - 1) > ($j + $iKL) ? $j + $iKL : $iRows - 1
					$iMax = 0 > ($j - $iKU) ? 0 : $j - $iKU
					For $i = $iMax To $iMin
						DllStructSetData($tStruct, 1, $aArray[$i][$j], $iX + $i)
					Next
					$iX += $iRows
				Next

			Case BitAND($nMode, $__g_BLAS_STYPE_TRIDIAGONAL)
			; not implemented yet
				Return SetError(6,0,Null)

			Case BitAND($nMode, $__g_BLAS_STYPE_TRIANGLE + $__g_BLAS_STYPE_SYMMETRIC) ; triangle and symmetric matrices are stored the same way - only the flag differs
				; symmetric:  https://www.ibm.com/docs/en/essl/6.3?topic=representation-lower-storage-mode
				; triangular: https://www.ibm.com/docs/en/essl/6.3?topic=representation-lower-triangular-storage-mode

				If BitAND($nMode, $__g_BLAS_STYPE_LOWER) Then

					; iterate over relevant elements only
					If BitAND($nMode, $__g_BLAS_STYPE_NONUNIT) Then
						; only the lower left elements without the main diagonal are processed
						For $j = 0 To $iCols - 1
							$iX += $j + 1
							For $i = $j + 1 To $iRows - 1
								DllStructSetData($tStruct, 1, $aArray[$i][$j], $iX)
								$iX += 1
							Next
						Next

					Else
						; only the lower left elements including the main diagonal are processed
						For $j = 0 To $iCols - 1
							$iX += $j
							For $i = $j To $iRows - 1
								DllStructSetData($tStruct, 1, $aArray[$i][$j], $iX)
								$iX += 1
							Next
						Next
					EndIf

				ElseIf BitAND($nMode, $__g_BLAS_STYPE_UPPER) Then
					; symmetric:  https://www.ibm.com/docs/en/essl/6.3?topic=representation-upper-storage-mode#am5gr_utsm
					; triangular: https://www.ibm.com/docs/en/essl/6.3?topic=representation-upper-triangular-storage-mode

					If BitAND($nMode, $__g_BLAS_STYPE_NONUNIT) Then
						; only the upper right elements without the main diagonal are processed
						For $j = 0 To $iCols - 1
							For $i = 0 To (($j-1) >= $iRows ? $iRows - 1 : ($j - 1))
								DllStructSetData($tStruct, 1, $aArray[$i][$j], $iX)
								$iX += 1
							Next
							$iX += $iRows - $j
						Next

					Else
						; only the upper right elements including the main diagonal are processed
						For $j = 0 To $iCols - 1
							For $i = 0 To ($j >= $iRows ? $iRows - 1 : $j)
								DllStructSetData($tStruct, 1, $aArray[$i][$j], $iX)
								$iX += 1
							Next
							$iX += $iRows - $j - 1
						Next
					EndIf

				Else
					ContinueCase ; handle as full matrix if no upper/lower defined

				EndIf

			Case Else
				; full unpacked representation of a matrix - process all elements
				; https://www.ibm.com/docs/en/essl/6.3?topic=matrices-matrix-storage-representation

				; store in column major order(Fortran order)
				For $j = 0 To $iCols - 1
					For $i = 0 To $iRows - 1
						DllStructSetData($tStruct, 1, $aArray[$i][$j], $iX)
						$iX += 1
					Next
				Next

		EndSelect
	EndIf

	$mRet.storageType = $nMode
	$mRet.struct      = $tStruct
	$mRet.ptr         = DllStructGetPtr($tStruct)

	Return SetExtended($iN, $mRet)
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_fromStruct()
; Description ...: creates a matrix/vector map from a DllStruct as used here in the UDF.
; Syntax ........: _blas_fromStruct($tStruct, $iRows, [$iCols = 0, [$sDatatype = "DOUBLE", [$nMode = $__g_BLAS_STYPE_MATRIX, [$iKL = Default, [$iKU = Default]]]]])
; Parameters ....: tStruct   - [DllStruct] the DllStruct variable with the payload data of the matrix/vector
;                  iRows     - [Int] Number of rows in the matrix/vector
;                  iCols     - [Int] (Default: 0)
;                            ↳ Number of columns in the matrix/vector
;                              If 0: The result will be a vector
;                  sDatatype - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
;                  nMode     - (Default: $__g_BLAS_STYPE_MATRIX)
;                            ↳ BLAS memory layout for the matrix (see $__g_BLAS_STYPE_XXX flags)
;                  iKL       - [Int] (Default: Default)
;                            ↳ if nMode contains $__g_BLAS_STYPE_BAND: number of lower band diagonals
;                  iKU       - [Int] (Default: Default)
;                            ↳ if nMode contains $__g_BLAS_STYPE_BAND: number of upper band diagonals
; Return value ..: [Map] BLAS Matrix Map as it is structured in this UDF
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: Yes
;                  Global $tStruct = DllStructCreate("DOUBLE[10]")
;                  For $i = 1 To 10
;                     DllStructSetData($tStruct, 1, $i * 2, $i)
;                  Next
;                  Global $mBandPacked = _blas_fromStruct($tStruct, 10, 10, "DOUBLE", $__g_BLAS_STYPE_BAND + $__g_BLAS_STYPE_PACKED, 0, 0)
;                  _blas_display($mBandPacked)
; ===============================================================================================================================
Func _blas_fromStruct(ByRef $tStruct, $iRows, $iCols = 0, $sDatatype = "DOUBLE", $nMode = $__g_BLAS_STYPE_MATRIX, $iKL = Default, $iKU = Default)
	Local $mRet[]

	$mRet.rows        = $iRows
	$mRet.cols        = ($iCols <= 1 ? 1 : $iCols)
	$mRet.elements    = $iRows * $mRet.cols
	$mRet.size        = DllStructGetSize($tStruct) / ($sDataType = "DOUBLE" ? $iBLAS_SIZE_DOUBLE : $iBLAS_SIZE_FLOAT)
	$mRet.storageType = ($iCols <= 1 ? 0 : $nMode)
	$mRet.datatype    = $sDatatype
	$mRet.struct      = $tStruct
	$mRet.ptr         = DllStructGetPtr($tStruct)

	; if type is packed band matrix
	If IsKeyword($iKL) <> 1 Then $mRet.kl = $iKL
	If IsKeyword($iKU) <> 1 Then $mRet.ku = $iKU

	Return $mRet
EndFunc   ;==>_la_fromStruct


; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_duplicate()
; Description ...: creates an independent copy of a matrix/vector map
; Syntax ........: _blas_duplicate($mMatrix)
; Parameters ....: mMatrix - [Map] a matrix/vector as a map, as structured here in the UDF
; Return value ..: Success: [Map] the copy
;                  Failure: Null and set @error to:
;                           | 1: $mMatrix is not a map
;                           | 2: error during RtlCopyMemory (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://learn.microsoft.com/en-us/windows-hardware/drivers/ddi/wdm/nf-wdm-rtlcopymemory
; Example .......: Yes
;                  Global $mMatrix = _blas_fromArray("[[4,-1,1],[-1,4,-2],[1,-2,5]]", "")
;                  Global $mDouble = _blas_duplicate($mMatrix)
;                  _blas_display($mDouble)
; ===============================================================================================================================
Func _blas_duplicate(Const $mMatrix)
	If Not IsMap($mMatrix) Then Return SetError(1, 0, Null)

	; create destination memory
	Local $tStruct = DllStructCreate(StringFormat("%s[%d]", $mMatrix.datatype, $mMatrix.size))
	Local $hPtr = DllStructGetPtr($tStruct)

	; copy content
	DllCall('msvcrt.dll', 'NONE:cdecl', 'memcpy', 'PTR', $hPtr, 'PTR', $mMatrix.ptr, 'ULONG_PTR', DllStructGetSize($tStruct))
	If @error Then Return SetError(2, @error, Null)

	Local $mRet = $mMatrix
	$mRet.struct = $tStruct
	$mRet.ptr = $hPtr

	Return $mRet
EndFunc   ;==>_blas_duplicate

; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_display()
; Description ...: displays a matrix/vector map, similar to _ArrayDisplay
; Syntax ........: _blas_display($mData, [$sTitle = "", [iDecimalPlaces = 5, [$iFlags = 64]]])
; Parameters ....: mData          - [Map] a matrix/vector as a map, as structured here in the UDF
;                  sTitle         - [String] (Default: "")
;                                 ↳ the window title to be displayed
;                  iDecimalPlaces - [UInt] (Default: 5)
;                                 ↳ number of decimal places to which the figures are to be rounded
;                  iFlags         - [Int] (Default: 64)
;                                 ↳ display options - see $iFlags for _ArrayDisplay()
; Return value ..: Success: return value of _ArrayDisplay
;                  Failure: False and set @error to:
;                           | 1: error during _blas_toArray() (@extended: @error from _blas_toArray)
;                           | 2: @error during _ArrayDisplay (@extended: @error from _ArrayDisplay )
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......: _ArrayDisplay()
; Link ..........:
; Example .......: Yes
;                  _blas_display(_blas_fromArray("[[1,2,3],[4,5,6],[7,8,9]]"))
; ===============================================================================================================================
Func _blas_display($mData, $sTitle = "", $iDecimalPlaces = 5, $iFlags = 64)
	Local $aArray = _blas_toArray($mData)
	If @error Then Return SetError(1, @error, False)

	If $iDecimalPlaces > 0 Then
		If UBound($aArray, 0) = 1 Then
			For $i = 0 To UBound($aArray, 1) - 1
				If IsNumber($aArray[$i]) Then $aArray[$i] = StringFormat("%." & $iDecimalPlaces & "g", $aArray[$i])
			Next
		Else
			For $i = 0 To UBound($aArray, 1) - 1
				For $j = 0 To UBound($aArray, 2) - 1
					If IsNumber($aArray[$i][$j]) Then $aArray[$i][$j] = StringFormat("%." & $iDecimalPlaces & "g", $aArray[$i][$j])
				Next
			Next
		EndIf
	EndIf

	Local $iRet = _ArrayDisplay($aArray, $sTitle, "", $iFlags)
	Return @error ? SetError(2, @error, $iRet) : $iRet
EndFunc

#EndRegion

#Region BLAS Level 1 (Vector Operations)


; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_axpy()
; Description ...: calculate a * X + Y  and store it in Y (if a = 1 then X + Y)
; Syntax ........: _blas_axpy($mVecX, $mVecY, [$fScalar = 1, [$iStartX = 0, [$iStartY = 0, [$iIncX = 1, [$iIncY = 1, [$iN = Default, [$sDataType = "DOUBLE"]]]]]]])
; Parameters ....: mVecX     - [Map] matrix/vector x as a map, DllStruct or pointer
;                  mVecY     - [Map] matrix/vector y as a map, DllStruct or pointer (will be overwritten)
;                  fScalar   - [Float] (Default: 1)
;                            ↳ scalar value a with which x is scaled
;                  iStartX   - [Int] (Default: 0)
;                            ↳ start element index in X (0-based)
;                  iStartY   - [Int] (Default: 0)
;                            ↳ start element index in Y (0-based)
;                  iIncX     - [Int] (Default: 1)
;                            ↳ storage spacing between elements of X (can be used to handle parts of a matrix as a vector)
;                  iIncY     - [Int] (Default: 1)
;                            ↳ storage spacing between elements of Y (can be used to handle parts of a matrix as a vector)
;                  iN        - [Int] (Default: Default)
;                            ↳ number of elements in input vector(s)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of axpy (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d5/d4b/group__axpy_gadb136e14634fe4b772a0b034f52f3939.html
; Example .......: Yes
;                  Global $mY = _blas_fromArray("[10,20,30,40,50]")
;                  _blas_axpy(_blas_fromArray("[1,2,3,4,5]"), $mY, 2)
;                  _blas_display($mY, "a*X + Y")
; ===============================================================================================================================
Func _blas_axpy($mVecX, $mVecY, $fScalar = 1, $iStartX = 0, $iStartY = 0, $iIncX = 1, $iIncY = 1, $iN = Default, $sDataType = "DOUBLE")
	Local $pX ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mVecX)
			$sDataType = $mVecX.datatype
			If IsKeyword($iN) = 1 Then $iN = $mVecX.elements
			$pX = $mVecX.ptr
		Case IsPtr($mVecX)
			$pX = $mVecX
		Case IsDllStruct($mVecX)
			$pX = DllStructGetPtr($mVecX)
	EndSelect
	Local $pY = IsMap($mVecY) ? $mVecY.ptr : (IsPtr($mVecY) ? $mVecY : DllStructGetPtr($mVecY))

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d", _
	            $dSize   = ($sDataType = "DOUBLE") ? $iBLAS_SIZE_DOUBLE : $iBLAS_SIZE_FLOAT

	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "axpy", _
		"INT*", $iN, _                     ; number of elements in vectors x and y
		$sDataType & "*", $fScalar, _      ; the scalar a
		"PTR",  $pX + $dSize * $iStartX, _ ; start ptr of x
		"INT*", $iIncX, _                  ; increment in x
		"PTR",  $pY + $dSize * $iStartY, _ ; start ptr of y
		"INT*", $iIncY _                   ; increment in y
	)
	Return @error ? SetError(1, @error, False) : True
EndFunc



; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_dot()
; Description ...: calculate the dot product of two vectors X and Y
; Syntax ........: _blas_dot($mVecX, $mVecY, [$iStartX = 0, [$iStartY = 0, [$iIncX = 1, [$iIncY = 1, [$iN = Default, [$sDataType = "DOUBLE"]]]]]])
; Parameters ....: mVecX     - [Map] matrix/vector x as a map, DllStruct or pointer
;                  mVecY     - [Map] matrix/vector y as a map, DllStruct or pointer
;                  iStartX   - [Int] (Default: 0)
;                            ↳ start element index in X (0-based)
;                  iStartY   - [Int] (Default: 0)
;                            ↳ start element index in Y (0-based)
;                  iIncX     - [Int] (Default: 1)
;                            ↳ storage spacing between elements of X (can be used to handle parts of a matrix as a vector)
;                  iIncY     - [Int] (Default: 1)
;                            ↳ storage spacing between elements of Y (can be used to handle parts of a matrix as a vector)
;                  iN        - [Int] (Default: Default)
;                            ↳ number of elements in input vector(s)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: [Number] the dot product
;                  Failure: Null and set @error to:
;                           | 1: error during DllCall of dot (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d1/dcc/group__dot_ga2a42ecc597403b22ad786715c739196b.html#ga2a42ecc597403b22ad786715c739196b
; Example .......: Yes
;                  Global $fDot = _blas_dot(_blas_fromArray("[1,2,-4,4,5]"), _blas_fromArray("[9,8,7,-6,5]")) ; --> -2
;                  ConsoleWrite("X · Y = " & $fDot & @CRLF)
; ===============================================================================================================================
Func _blas_dot($mVecX, $mVecY, $iStartX = 0, $iStartY = 0, $iIncX = 1, $iIncY = 1, $iN = Default, $sDataType = "DOUBLE")
	Local $pX ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mVecX)
			$sDataType = $mVecX.datatype
			If IsKeyword($iN) = 1 Then $iN = $mVecX.elements
			$pX = $mVecX.ptr
		Case IsPtr($mVecX)
			$pX = $mVecX
		Case IsDllStruct($mVecX)
			$pX = DllStructGetPtr($mVecX)
	EndSelect
	Local $pY = IsMap($mVecY) ? $mVecY.ptr : (IsPtr($mVecY) ? $mVecY : DllStructGetPtr($mVecY))

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d", _
	            $dSize   = ($sDataType = "DOUBLE") ? $iBLAS_SIZE_DOUBLE : $iBLAS_SIZE_FLOAT

	Local $aDLL = DllCall($__g_hBLAS_DLL, $sDataType & ":cdecl", $cPrefix & "dot", _
		"INT*", $iN, _                     ; number of elements in vectors x and y
		"PTR",  $pX + $dSize * $iStartX, _ ; start ptr of x
		"INT*", $iIncX, _                  ; increment in x
		"PTR",  $pY + $dSize * $iStartY, _ ; start ptr of y
		"INT*", $iIncY _                   ; increment in y
	)
	Return @error ? SetError(1, @error, Null) : $aDLL[0]
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_rot()
; Description ...: applies a plane rotation to coordinate-pairs
; Syntax ........: _blas_rot($mVecX, $mVecY, $fAlpha, [$iStartX = 0, [$iStartY = 0, [$iIncX = 1, [$iIncY = 1, [$iN = Default, [$sDataType = "DOUBLE"]]]]]])
; Parameters ....: mVecX     - [Map] vector with x-coordinates as a map, as structured here in the UDF (will be overwritten)
;                  mVecY     - [Map] vector with y-coordinates as a map, as structured here in the UDF (will be overwritten)
;                  fAlpha    - [Float] rotation angle [radian]
;                  iStartX   - [Int] (Default: 0)
;                            ↳ start element index in X (0-based)
;                  iStartY   - [Int] (Default: 0)
;                            ↳ start element index in Y (0-based)
;                  iIncX     - [Int] (Default: 1)
;                            ↳ storage spacing between elements of X (can be used to handle parts of a matrix as a vector)
;                  iIncY     - [Int] (Default: 1)
;                            ↳ storage spacing between elements of Y (can be used to handle parts of a matrix as a vector)
;                  iN        - [Int] (Default: Default)
;                            ↳ number of elements in input vector(s)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of rot (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d1/d45/group__rot_gae48ef017306866ac2d5a8c5a52617858.html#gae48ef017306866ac2d5a8c5a52617858
; Example .......: Yes
;                  Global $mX = _blas_fromArray("[1,2,3,4,5]")
;                  Global $mY = _blas_fromArray("[-1,-2,-3,-4,-5]")
;                  _blas_rot($mX, $mY, 1.0471975512)
;                  _blas_display($mX, "rotated x coordinates")
;                  _blas_display($mY, "rotated y coordinates")
; ===============================================================================================================================
Func _blas_rot(ByRef $mVecX, $mVecY, $fAlpha, $iStartX = 0, $iStartY = 0, $iIncX = 1, $iIncY = 1, $iN = Default, $sDataType = "DOUBLE")
	Local $pX ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mVecX)
			$sDataType = $mVecX.datatype
			If IsKeyword($iN) = 1 Then $iN = $mVecX.elements
			$pX = $mVecX.ptr
		Case IsPtr($mVecX)
			$pX = $mVecX
		Case IsDllStruct($mVecX)
			$pX = DllStructGetPtr($mVecX)
	EndSelect
	Local $pY = IsMap($mVecY) ? $mVecY.ptr : (IsPtr($mVecY) ? $mVecY : DllStructGetPtr($mVecY))

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d", _
	            $sType   = $sDataType & "*", _
	            $dSize   = ($sDataType = "DOUBLE") ? $iBLAS_SIZE_DOUBLE : $iBLAS_SIZE_FLOAT

	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "rot", _
		"INT*", $iN, _                     ; number of elements in vectors x and y
		"PTR",  $pX + $dSize * $iStartX, _ ; start ptr of x
		"INT*", $iIncX, _                  ; increment in x
		"PTR",  $pY + $dSize * $iStartY, _ ; start ptr of y
		"INT*", $iIncY, _                  ; increment in y
		$sType, Cos($fAlpha), _            ; rotation angle
		$sType, Sin($fAlpha) _             ; rotation angle
	)
	Return @error ? SetError(1, @error, False) : True
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_rotg()
; Description ...: constructs a plane rotation
; Syntax ........: _blas_rotg($a, $b, [$sType = "DOUBLE"])
; Parameters ....: a     - [Number] the scalar a
;                  b     - [Number] the scalar b
;                  sType -[String] (Default: "DOUBLE")
;                        ↳ data type - either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: array[4] with the values [r,z,c,s]
;                  Failure: Null and set @error to:
;                           | 1: error during DllCall of rotg (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d7/dc5/group__rotg_gaafa91c51f75df6c3f2182032a221c2db.html#gaafa91c51f75df6c3f2182032a221c2db
; Example .......: Yes
;                  Global $aRZCS = _blas_rotg(11, 3)
;                  _ArrayDisplay($aRZCS, "R Z C S")
; ===============================================================================================================================
Func _blas_rotg($a, $b, $sType = "DOUBLE")
	Local Const $cPrefix = ($sType = "FLOAT") ? "s" : "d"
	$sType &= "*"

	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "rotg", $sType, $a, $sType, $b, $sType, 0, $sType, 0)
	If @error Then Return SetError(1, @error, Null)
	Local $aRet[4] = [$aDLL[1], $aDLL[2], $aDLL[3], $aDLL[4]]
	Return $aRet
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_swap()
; Description ...: interchanges two vectors (useful for matrix manipulations)
; Syntax ........: _blas_swap($mMatrixA, $mMatrixB, [$iStartA = 0, [$iStartB = 0, [$iIncA = 1, [$iIncB = 1, [$iN = Default, [$sDataType = "DOUBLE"]]]]]])
; Parameters ....: mMatrixA  - [Map] matrix/vector A as a map, DllStruct or pointer
;                  mMatrixB  - [Map] matrix/vector B as a map, DllStruct or pointer
;                  iStartA   - [Int] (Default: 0)
;                            ↳ start element index in A (0-based)
;                  iStartB   - [Int] (Default: 0)
;                            ↳ start element index in B (0-based)
;                  iIncA     - [Int] (Default: 1)
;                            ↳ storage spacing between elements of A (can be used to handle parts of a matrix as a vector)
;                  iIncB     - [Int] (Default: 1)
;                            ↳ storage spacing between elements of B (can be used to handle parts of a matrix as a vector)
;                  iN        - [Int] (Default: Default)
;                            ↳ number of elements in input vector(s)
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of swap (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d7/d51/group__swap_ga780475990528dce288cf4f7bba36c90f.html#ga780475990528dce288cf4f7bba36c90f
; Example .......: Yes
;                  Global $mMatrix = _blas_fromArray("[[1,2,3,4,5],[6,7,8,9,10],[11,12,13,14,15],[16,17,18,19,20],[21,22,23,24,25]]")
;                  _blas_swap($mMatrix, $mMatrix, 0, 15, 1, 1, 5)
;                  _blas_display($mMatrix, "swapped columns 0 & 3")
; ===============================================================================================================================
Func _blas_swap($mMatrixA, $mMatrixB, $iStartA = 0, $iStartB = 0, $iIncA = 1, $iIncB = 1, $iN = Default, $sDataType = "DOUBLE")
	Local $pA ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrixA)
			$sDataType = $mMatrixA.datatype
			; calculate number of elements if not defined
			If IsKeyword($iN) = 1 Then $iN = Floor(($mMatrixA.elements - $iStartA) / $iIncA)
			$pA = $mMatrixA.ptr
		Case IsPtr($mMatrixA)
			$pA = $mMatrixA
		Case IsDllStruct($mMatrixA)
			$pA = DllStructGetPtr($mMatrixA)
	EndSelect
	Local $pB = IsMap($mMatrixB) ? $mMatrixB.ptr : (IsPtr($mMatrixB) ? $mMatrixB : DllStructGetPtr($mMatrixB))

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d", _
	            $dSize   = ($sDataType = "DOUBLE") ? $iBLAS_SIZE_DOUBLE : $iBLAS_SIZE_FLOAT

	; call the function
	DllCall($__g_hBLAS_DLL, $sDataType & ":cdecl", $cPrefix & "swap", _
		"INT*", $iN, _                     ; number of elements to check
		"PTR",  $pA + $dSize * $iStartA, _ ; start ptr to read
		"INT*", $iIncA, _                  ; increment
		"PTR",  $pB + $dSize * $iStartB, _ ; start ptr to read
		"INT*", $iIncB _                   ; increment
	)
	Return @error ? SetError(1, @error, False) : True
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_copy()
; Description ...: copies a vector x to a vector y (useful for extracting parts of a matrix to other areas)
; Syntax ........: _blas_copy($mMatrix, [$iStartX = 0, [$iIncX = 1, [$iStartY = 0, [$iIncY = 1, [$iN = Default, [$mTarget = Default, [$bRetMap = True, [$sDataType = "DOUBLE"]]]]]]]])
; Parameters ....: mMatrix   - [Map] matrix/vector x as a map, DllStruct or pointer
;                  iStartX   - [Int] (Default: 0)
;                            ↳ start element index in X (0-based)
;                  iStartY   - [Int] (Default: 0)
;                            ↳ start element index in Y (0-based)
;                  iIncX     - [Int] (Default: 1)
;                            ↳ storage spacing between elements of X (can be used to handle parts of a matrix as a vector)
;                  iIncY     - [Int] (Default: 1)
;                            ↳ storage spacing between elements of Y (can be used to handle parts of a matrix as a vector)
;                  iN        - [Int] (Default: Default)
;                            ↳ number of elements in input vector(s)
;                  mTarget   - [Map] (Default: Default)
;                            ↳ target matrix/vector x as a map, DllStruct or pointer
;                              if Default: a new matrix/vector is returned
;                  bRetMap   - [Bool] (Default: True)
;                            ↳ if true: a new map/vector is created if mTarget is a DllStruct or a point only
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: If $bRetMap: True result matrix/vector - else: True
;                  Failure: Null and set @error to:
;                           | 1: error during DllCall of copy (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d5/d2b/group__copy_gafc29c83942509cba5f3764115f0471d5.html#gafc29c83942509cba5f3764115f0471d5
; Example .......: Yes
;                  Global $mMatrix = _blas_fromArray("[[1,2,3,4,5],[6,7,8,9,10],[11,12,13,14,15],[16,17,18,19,20],[21,22,23,24,25]]")
;                  Global $mRow = _blas_copy($mMatrix, 1, $mMatrix.cols) ; extract row
;                  _blas_display($mRow)
; ===============================================================================================================================
Func _blas_copy(Const $mMatrix, Const $iStartX = 0, Const $iIncX = 1, Const $iStartY = 0, Const $iIncY = 1, $iN = Default, $mTarget = Default, $bRetMap = True, $sDataType = "DOUBLE")
	Local $pA, $pT ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrix)
			$sDataType = $mMatrix.datatype
			; determine the maximum number of elements
			If IsKeyword($iN) = 1 Then $iN = Int(($mMatrix.elements - $iStartX - 1) / $iIncX) + 1
			$pA = $mMatrix.ptr
		Case IsPtr($mMatrix)
			$pA = $mMatrix
		Case IsDllStruct($mMatrix)
			$pA = DllStructGetPtr($mMatrix)
	EndSelect
	Select
		Case IsKeyword($mTarget) = 1
			$mTarget = _blas_createVector($iN, $sDataType)
			$pT = $mTarget.ptr
		Case IsMap($mTarget)
			$pT = $mTarget.ptr
		Case IsPtr($mTarget)
			$pT = $mTarget
			If $bRetMap Then _blas_createVector($iN, $sDataType, $pT)
		Case IsDllStruct($mTarget)
			$pT = DllStructGetPtr($mTarget)
			If $bRetMap Then _blas_createVector($iN, $sDataType, $pT)
	EndSelect

	Local Const $cPrefix = $sDataType = "FLOAT" ? "s" : "d", _
	            $dSize   = ($sDataType = "DOUBLE") ? $iBLAS_SIZE_DOUBLE : $iBLAS_SIZE_FLOAT

	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "copy", _
		"INT*", $iN, _                     ; number of elements to read
		"PTR",  $pA + $dSize * $iStartX, _ ; start ptr to read
		"INT*", $iIncX, _                  ; increment in x direction
		"PTR",  $pT + $dSize * $iStartY, _ ; target matrix/vector
		"INT*", $iIncY _                   ; storage spacing between elements
	)
	Return @error ? SetError(1, @error, Null) : ($bRetMap ? $mTarget : True)
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_scal()
; Description ...: scales a vector with a scalar (X = a * X)
; Syntax ........: _blas_scal($mVector, $fScale, [$iStart = 0, [$iInc = 1, [$iN = Default, [$sDataType = "DOUBLE"]]]])
; Parameters ....: mVector   - [Map] matrix/vector x as a map, DllStruct or pointer (will be overwritten)
;                  fScale    - [Float] the scaling value
;                  iStart    - [Int] (Default: 0)
;                            ↳ start element index in X (0-based)
;                  iInc      - [Int] (Default: 1)
;                            ↳ storage spacing between elements of X
;                  iN        - [Int] (Default: Default)
;                            ↳ number of elements in input vector
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of scal (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d2/de8/group__scal.html
; Example .......: Yes
;                  Global $mVector = _blas_fromArray("[1,2,3,4]")
;                  _blas_scal($mVector, 3)
;                  _blas_display($mVector, "scaled vector")
; ===============================================================================================================================
Func _blas_scal($mVector, $fScale, $iStart = 0, $iInc = 1, $iN = Default, $sDataType = "DOUBLE")
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

	DllCall($__g_hBLAS_DLL, $sDataType & ":cdecl", $cPrefix & "scal", _
		"INT*",           $iN, _                    ; number of elements to check
		$sDataType & "*", $fScale, _                ; scalar
		"PTR",            $pV + $dSize * $iStart, _ ; start ptr to read
		"INT*",           $iInc _                   ; increment
	)
	Return @error ? SetError(1, @error, False) : True
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_nrm2()
; Description ...: calculate the euclidean norm of a vector
; Syntax ........: _blas_nrm2($mVector, [$iStart = 0, [$iInc = 1, [$iN = Default, [$sDataType = "DOUBLE"]]]])
; Parameters ....: mVector   - [Map] matrix/vector as a map, DllStruct or pointer (will be overwritten)
;                  iStart    - [Int] (Default: 0)
;                            ↳ start index
;                  iInc      - [Int] (Default: 1)
;                            ↳ storage spacing between elements
;                  iN        - [Int] (Default: Default)
;                            ↳ number of elements in input vector
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: [Number] the euclidian norm value
;                  Failure: Null and set @error to:
;                           | 1: error during DllCall of nrm2 (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d1/d2a/group__nrm2_gab5393665c8f0e7d5de9bd1dd2ff0d9d0.html#gab5393665c8f0e7d5de9bd1dd2ff0d9d0
; Example .......: Yes
;                  Global $fNorm = _blas_nrm2(_blas_fromArray("[1,-7,-2,-3]")) ; --> 7.9372
;                  ConsoleWrite("Norm(x): " & $fNorm & @CRLF)
; ===============================================================================================================================
Func _blas_nrm2($mVector, $iStart = 0, $iInc = 1, $iN = Default, $sDataType = "DOUBLE")
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

	Local $aDLL = DllCall($__g_hBLAS_DLL, $sDataType & ":cdecl", $cPrefix & "nrm2", _
		"INT*", $iN, _                    ; number of elements to check
		"PTR",  $pV + $dSize * $iStart, _ ; start ptr to read
		"INT*", $iInc _                   ; increment
	)
	Return @error ? SetError(1, @error, Null) : $aDLL[0]
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_asum()
; Description ...: calculate the sum of the absolute(!) values of a matrix/vector
; Syntax ........: _blas_asum($mMatrix, [$iStart = 0, [$iInc = 1, [$iN = Default, [$sDataType = "DOUBLE"]]]])
; Parameters ....: mMatrix   - [Map] matrix/vector as a map, DllStruct or pointer
;                  iStart    - [Int] (Default: 0)
;                            ↳ start index
;                  iInc      - [Int] (Default: 1)
;                            ↳ storage spacing between elements
;                  iN        - [Int] (Default: Default)
;                            ↳ number of elements in input vector
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: [Number] absolute sum of the values
;                  Failure: Null and set @error to:
;                           | 1: error during DllCall of asum (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d5/d72/group__asum_ga829029987b14b622f355aacf54a8e4b9.html#ga829029987b14b622f355aacf54a8e4b9
; Example .......: Yes
;                  Global $fSum = _blas_asum(_blas_fromArray('[1,-7,-2,-3]')) ; --> 13
;                  ConsoleWrite("Sum(|x|): " & $fSum & @CRLF)
; ===============================================================================================================================
Func _blas_asum($mMatrix, $iStart = 0, $iInc = 1, $iN = Default, $sDataType = "DOUBLE")
	Local $pM ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrix)
			$sDataType = $mMatrix.datatype
			; calculate number of elements if not defined
			If IsKeyword($iN) = 1 Then $iN = Floor(($mMatrix.elements - $iStart) / $iInc)
			$pM = $mMatrix.ptr
		Case IsPtr($mMatrix)
			$pM = $mMatrix
		Case IsDllStruct($mMatrix)
			$pM = DllStructGetPtr($mMatrix)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d", _
	            $dSize   = ($sDataType = "DOUBLE") ? $iBLAS_SIZE_DOUBLE : $iBLAS_SIZE_FLOAT

	; call the function
	Local $aDLL = DllCall($__g_hBLAS_DLL, $sDataType & ":cdecl", $cPrefix & "asum", _
		"INT*", $iN, _                    ; number of elements to check
		"PTR",  $pM + $dSize * $iStart, _ ; start ptr to read
		"INT*", $iInc _                   ; increment
	)
	Return @error ? SetError(1, @error, Null) : $aDLL[0]
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_amax()
; Description ...: finds the first element having the maximum absolute(!) value
; Syntax ........: _blas_amax($mMatrix, [$iStart = 0, [$iInc = 1, [$iN = Default, [$sDataType = "DOUBLE"]]]])
; Parameters ....: mMatrix   - [Map] matrix/vector as a map, DllStruct or pointer
;                  iStart    - [Int] (Default: 0)
;                            ↳ start index
;                  iInc      - [Int] (Default: 1)
;                            ↳ storage spacing between elements
;                  iN        - [Int] (Default: Default)
;                            ↳ number of elements in input vector
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: [Number] the maximum absolute value (@extended = index of first element having this value)
;                  Failure: Null and set @error to:
;                           | 1: error during DllCall of amax (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/dd/d52/group__iamax_gacec03c5109f531c06b4fb301cf1a2d7a.html#gacec03c5109f531c06b4fb301cf1a2d7a
; Example .......: Yes
;                  Global $fMax = _blas_amax(_blas_fromArray('[2,14,5,-6,7]')) ; --> 14 (index = 2)
;                  ConsoleWrite("Max: " & $fMax & " (Index: " & @extended & ")" & @CRLF)
; ===============================================================================================================================
Func _blas_amax($mMatrix, $iStart = 0, $iInc = 1, $iN = Default, $sDataType = "DOUBLE")
	Local $pM ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrix)
			$sDataType = $mMatrix.datatype
			; calculate number of elements if not defined
			If IsKeyword($iN) = 1 Then $iN = Floor(($mMatrix.elements - $iStart) / $iInc)
			$pM = $mMatrix.ptr
		Case IsPtr($mMatrix)
			$pM = $mMatrix
		Case IsDllStruct($mMatrix)
			$pM = DllStructGetPtr($mMatrix)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d", _
	            $dSize   = ($sDataType = "DOUBLE") ? $iBLAS_SIZE_DOUBLE : $iBLAS_SIZE_FLOAT

	Local $aDLL = DllCall($__g_hBLAS_DLL, "INT:cdecl", "i" & $cPrefix & "amax", _
		"INT*", $iN, _                    ; number of elements to check
		"PTR",  $pM + $dSize * $iStart, _ ; start ptr to read
		"INT*", $iInc _                   ; increment
	)
	If @error Then Return SetError(1, @error, Null)

	; extract the return value
	Local $fRet
	Select
		Case IsMap($mMatrix)
			$fRet = DllStructGetData($mMatrix.struct, 1, $aDLL[0] + $iStart)
		Case IsPtr($mMatrix)
			Local $tStruct = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iStart + $aDLL[0]), $mMatrix)
			$fRet = DllStructGetData($tStruct, 1, $aDLL[0] + $iStart)
		Case IsDllStruct($mMatrix)
			$fRet = DllStructGetData($mMatrix, 1, $aDLL[0] + $iStart)
	EndSelect

	Return SetExtended($aDLL[0], $fRet)
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_amin()
; Description ...: finds the first element having the minimum absolute(!) value
; Syntax ........: _blas_amin($mMatrix, [$iStart = 0, [$iInc = 1, [$iN = Default, [$sDataType = "DOUBLE"]]]])
; Parameters ....: mMatrix   - [Map] matrix/vector as a map, DllStruct or pointer
;                  iStart    - [Int] (Default: 0)
;                            ↳ start index
;                  iInc      - [Int] (Default: 1)
;                            ↳ storage spacing between elements
;                  iN        - [Int] (Default: Default)
;                            ↳ number of elements in input vector
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: [Number] the minimum absolute value (@extended = index of first element having this value)
;                  Failure: Null and set @error to:
;                           | 1: error during DllCall of amin (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2023-0/i-amin.html
; Example .......: Yes
;                  Global $fMax = _blas_amin(_blas_fromArray('[2,14,5,-1,7]')) ; --> -1 (index = 4)
;                  ConsoleWrite("Min: " & $fMax & " (Index: " & @extended & ")" & @CRLF)
; ===============================================================================================================================
Func _blas_amin($mMatrix, $iStart = 0, $iInc = 1, $iN = Default, $sDataType = "DOUBLE")
	Local $pM ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrix)
			$sDataType = $mMatrix.datatype
			; calculate number of elements if not defined
			If IsKeyword($iN) = 1 Then $iN = Floor(($mMatrix.elements - $iStart) / $iInc)
			$pM = $mMatrix.ptr
		Case IsPtr($mMatrix)
			$pM = $mMatrix
		Case IsDllStruct($mMatrix)
			$pM = DllStructGetPtr($mMatrix)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d", _
	            $dSize   = ($sDataType = "DOUBLE") ? $iBLAS_SIZE_DOUBLE : $iBLAS_SIZE_FLOAT

	Local $aDLL = DllCall($__g_hBLAS_DLL, "INT:cdecl", "i" & $cPrefix & "amin", _
		"INT*", $iN, _                    ; number of elements to check
		"PTR",  $pM + $dSize * $iStart, _ ; start ptr to read
		"INT*", $iInc _                   ; increment
	)
	If @error Then Return SetError(1, @error, Null)

	; extract the return value
	Local $fRet
	Select
		Case IsMap($mMatrix)
			$fRet = DllStructGetData($mMatrix.struct, 1, $aDLL[0] + $iStart)
		Case IsPtr($mMatrix)
			Local $tStruct = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iStart + $aDLL[0]), $mMatrix)
			$fRet = DllStructGetData($tStruct, 1, $aDLL[0] + $iStart)
		Case IsDllStruct($mMatrix)
			$fRet = DllStructGetData($mMatrix, 1, $aDLL[0] + $iStart)
	EndSelect

	Return SetExtended($aDLL[0], $fRet)
EndFunc

#EndRegion



#Region BLAS Level 2 (Matrix-Vector Operations)

; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_gemv()
; Description ...: matrix-vector multiplication:      y := alpha*A*x + beta*y,   or   y := alpha*Aᵀ*x + beta*y
;                  where A is a general unpacked matrix
; Syntax ........: _blas_gemv($mMatrix, $mVecX, $mVecY, [$fAlpha = 1, [$fBeta = 0, [$cTransposed = "N", [$iIncX = 1, [$iIncY = 1, [$iM = Default, [$iN = Default, [$iLDA = Default, [$sDataType = "DOUBLE"]]]]]]]]])
; Parameters ....: mMatrix     - [Map] matrix A as a map, DllStruct or pointer
;                  mVecX       - [Map] vector x as a map, DllStruct or pointer
;                  mVecY       - [Map] vector y as a map, DllStruct or pointer (will be overwritten)
;                  fAlpha      - [Float] (Default: 1)
;                              ↳ scalar value alpha with which x is scaled
;                  fBeta       - [Float] (Default: 0)
;                              ↳ scalar value beta with which y is scaled
;                  cTransposed - [Char] (Default: "N")
;                              ↳ if "T": instead of A the transposed matrix Aᵀ is used
;                  iIncX       - [Int] (Default: 1)
;                              ↳ storage spacing between elements in x
;                  iIncY       - [Int] (Default: 1)
;                              ↳ storage spacing between elements in y
;                  iM          - [Int] (Default: Default)
;                              ↳ number of rows in A
;                  iN          - [Int] (Default: Default)
;                              ↳ number of columns in A
;                  iLDA        - [Int] (Default: Default)
;                              ↳ first("leading") dimension size of A
;                  sDataType   - [String] (Default: "DOUBLE")
;                              ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of gemv (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d7/dda/group__gemv_ga4ac1b675072d18f902db8a310784d802.html#ga4ac1b675072d18f902db8a310784d802
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[8,1,2,3], [4,6,2,3], [2,7,3,7]]")
;                  Global $mX = _blas_fromArray("[3,2,1,2]")
;                  Global $mY = _blas_fromArray("[5,3,2]")
;                  _blas_gemv($mA, $mX, $mY)
;                  _blas_display($mY, "y <- beta * y + alpha * A * x")
; ===============================================================================================================================
Func _blas_gemv($mMatrix, $mVecX, $mVecY, $fAlpha = 1, $fBeta = 0, $cTransposed = "N", Const $iIncX = 1, Const $iIncY = 1, $iM = Default, $iN = Default, $iLDA = Default, $sDataType = "DOUBLE")
	Local $pM ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrix)
			$sDataType = $mMatrix.datatype

			; handle default parameter
			If IsKeyword($iM)   = 1 Then $iM   = $cTransposed = "N" ? $mMatrix.rows : $mMatrix.cols
			If IsKeyword($iN)   = 1 Then $iN   = $cTransposed = "N" ? $mMatrix.cols : $mMatrix.rows
			If IsKeyword($iLDA) = 1 Then $iLDA = $cTransposed = "N" ? $iM : $iN

			$pM = $mMatrix.ptr
		Case IsPtr($mMatrix)
			$pM = $mMatrix
		Case IsDllStruct($mMatrix)
			$pM = DllStructGetPtr($mMatrix)
	EndSelect
	Local $pX = IsMap($mVecX) ? $mVecX.ptr : (IsPtr($mVecX) ? $mVecX : DllStructGetPtr($mVecX))
	Local $pY = IsMap($mVecY) ? $mVecY.ptr : (IsPtr($mVecY) ? $mVecY : DllStructGetPtr($mVecY))

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d", _
	            $sType   = $sDataType & "*"

	; set char buffers
	DllStructSetData($tBLASCHAR1, 1, $cTransposed = "T" ? "T" : "N")

	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "gemv", _
		"PTR",  $pBLASCHAR1, _   ; char *transa  ( If "N" then A is used; if "T" then Aᵀ is used)
		"INT*", $iM, _           ; number of rows in matrix A
		"INT*", $iN, _           ; number of cols in matrix A
		$sType, $fAlpha, _       ; the scaling constant alpha for
		"ptr",  $pM, _           ; A
		"INT*", $iLDA, _         ; lda: length of a column in a (used to calc the element index - can be used to change indexing behavior)
		"ptr",  $pX, _           ; the vector x
		"INT*", $iIncX, _        ; the stride for vector x
		$sType, $fBeta, _        ; the scaling constant beta
		"ptr",  $pY, _           ; the vector y
		"INT*", $iIncY _         ; the stride for vector y
	)
	Return @error ? SetError(1, @error, False) : True
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_trmv()
; Description ...: matrix-vector multiplication:      x := A*x,   or   x := Aᵀ*x
;                  where A is a unpacked upper or lower triangular matrix
; Syntax ........: _blas_trmv($mMatrix, $mVecX, [$cUPLO = "U", [$cDIAG = "N", [$cTRANS = "N", [$iN = Default, [$iIncX = 1, [$iLDA = Default, [$sDataType = "DOUBLE"]]]]]]])
; Parameters ....: mMatrix   - [Map] matrix A as a map, DllStruct or pointer
;                  mVecX     - [Map] vector x as a map, DllStruct or pointer (gets overwritten)
;                  cUPLO     - [Char] (Default: "U")
;                            ↳ "U": A is an upper triangular matrix
;                              "L": A is an lower triangular matrix
;                  cDIAG     - [Char] (Default: "N")
;                            ↳ "U": A is a unit triangular matrix (diagonal = identity = 1)
;                              "N": A is not a unit triangular matrix
;                  cTRANS    - [Char] (Default: "N")
;                            ↳ "N":      x := A * x
;                              "T"/"C":  x := Aᵀ * x
;                  iN        - [Int] (Default: Default)
;                            ↳ order of matrix A
;                  iIncX     - [Int] (Default: 1)
;                            ↳ storage spacing between elements in x
;                  iLDA      - [Int] (Default: Default)
;                            ↳ first("leading") dimension size of A
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of trmv (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d6/d1c/group__trmv_ga73370bd6dca01abe05d54ecd1d91ce9a.html#ga73370bd6dca01abe05d54ecd1d91ce9a
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[1,0,0,0],[1,1,0,0],[2,3,1,0],[3,4,3,1]]")
;                  Global $mX = _blas_fromArray("[1,2,3,4]")
;                  _blas_trmv($mA, $mX, "L", "N")
;                  _blas_display($mX, "A*x")
; ===============================================================================================================================
Func _blas_trmv($mMatrix, $mVecX, $cUPLO = "U", $cDIAG = "N", $cTRANS = "N", $iN = Default, Const $iIncX = 1, $iLDA = Default, $sDataType = "DOUBLE")
	Local $pM ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrix)
			$sDataType = $mMatrix.datatype
			If IsKeyword($iN) = 1 Then $iN = $cTRANS = "N" ? $mMatrix.rows : $mMatrix.cols
			If IsKeyword($iLDA) = 1 Then $iLDA = $iN
			$pM = $mMatrix.ptr
		Case IsPtr($mMatrix)
			$pM = $mMatrix
		Case IsDllStruct($mMatrix)
			$pM = DllStructGetPtr($mMatrix)
	EndSelect
	Local $pV = IsMap($mVecX) ? $mVecX.ptr : (IsPtr($mVecX) ? $mVecX : DllStructGetPtr($mVecX))

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; set char buffers and input healing
	DllStructSetData($tBLASCHAR1, 1, $cUPLO = "U" ? "U" : "L")
	DllStructSetData($tBLASCHAR2, 1, $cTRANS = "N" ? "N" : ($cTRANS = "T" ? "T" :"C"))
	DllStructSetData($tBLASCHAR3, 1, $cDIAG = "N" ? "N" : "U")

	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "trmv", _
		"PTR",  $pBLASCHAR1, _     ; UPLO -> 'U': A = upper triangular, 'L': A = lower triangular
		"PTR",  $pBLASCHAR2, _     ; TRANS -> 'N': x = A*x, 'T': x = Aᵀ*x, 'C': x = Aᵀ*x
		"PTR",  $pBLASCHAR3, _     ; DIAG -> 'U': A = unit triangular, 'N': A = other elements on diagonal
		"INT*", $iN, _             ; N
		"ptr",  $pM, _             ; A
		"INT*", $iLDA, _           ; lda: length of a column/row (leading dimension)
		"ptr",  $pV, _             ; the vector x
		"INT*", $iIncX _           ; the stride for vector x
	)
	Return @error ? SetError(1, @error, False) : True
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_symv()
; Description ...: matrix-vector multiplication:      y := alpha*A*x + beta*y
;                  where A is an unpacked symmetric matrix
; Syntax ........: _blas_symv($mMatrix, $mVecX, $mVecY, [$fAlpha = 1, [$fBeta = 0, [$cUPLO = "U", [$iN = Default, [$iLDA = Default, [$iIncX = 1, [$iIncY = 1, [$sDataType = "DOUBLE"]]]]]]]])
; Parameters ....: mMatrix   - [Map] matrix A as a map, DllStruct or pointer
;                  mVecX     - [Map] vector x as a map, DllStruct or pointer
;                  mVecY     - [Map] vector y as a map, DllStruct or pointer (will be overwritten)
;                  fAlpha    - [Float] (Default: 1)
;                            ↳ scalar value alpha with which x is scaled
;                  fBeta     - [Float] (Default: 0)
;                            ↳ scalar value beta with which y is scaled
;                  cUPLO     - [Char] (Default: "U")
;                            ↳ "U": upper triangular parts of A are used
;                              "L": lower triangular parts of A are used
;                  iN        - [Int] (Default: Default)
;                            ↳ order of matrix A
;                  iLDA      - [Int] (Default: Default)
;                            ↳ first("leading") dimension size of A
;                  iIncX     - [Int] (Default: 1)
;                            ↳ storage spacing between elements in x
;                  iIncY     - [Int] (Default: 1)
;                            ↳ storage spacing between elements in y
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of symv (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/db/d17/group__hemv_ga0b20bcf6e94079dce2f3d035798e9738.html#ga0b20bcf6e94079dce2f3d035798e9738
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[8,1,2,3], [1,5,7,9], [2,7,3,7], [3,9,7,8]]")
;                  Global $mX = _blas_fromArray("[3,2,1,2]")
;                  Global $mY = _blas_fromArray("[5,3,2,5]")
;                  _blas_symv($mA, $mX, $mY, 1, 1)
;                  _blas_display($mY, "y <- beta * y + alpha * A * x")
; ===============================================================================================================================
Func _blas_symv($mMatrix, $mVecX, $mVecY, $fAlpha = 1, $fBeta = 0, $cUPLO = "U", $iN = Default, $iLDA = Default, Const $iIncX = 1, Const $iIncY = 1, $sDataType = "DOUBLE")
	Local $pM ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrix)
			$sDataType = $mMatrix.datatype
			If IsKeyword($iN) = 1 Then $iN = $mMatrix.rows
			If IsKeyword($iLDA) = 1 Then $iLDA = $iN
			$pM = $mMatrix.ptr
		Case IsPtr($mMatrix)
			$pM = $mMatrix
		Case IsDllStruct($mMatrix)
			$pM = DllStructGetPtr($mMatrix)
	EndSelect
	Local $pX = IsMap($mVecX) ? $mVecX.ptr : (IsPtr($mVecX) ? $mVecX : DllStructGetPtr($mVecX))
	Local $pY = IsMap($mVecY) ? $mVecY.ptr : (IsPtr($mVecY) ? $mVecY : DllStructGetPtr($mVecY))

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; set char buffers and input healing
	DllStructSetData($tBLASCHAR1, 1, $cUPLO = "U" ? "U" : "L")

	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "symv", _
		"PTR",            $pBLASCHAR1, _     ; UPLO -> 'U': A = upper triangular, 'L': A = lower triangular
		"INT*",           $iN, _             ; N
		$sDataType & "*", $fAlpha, _         ; ALPHA
		"ptr",            $pM, _             ; A
		"INT*",           $iLDA, _           ; lda: length of a column/row (leading dimension)
		"ptr",            $pX, _             ; the vector x
		"INT*",           $iIncX, _          ; the stride for vector x
		$sDataType & "*", $fBeta, _          ; BETA
		"ptr",            $pY, _             ; the vector y
		"INT*",           $iIncY _           ; the stride for vector y
	)
	Return @error ? SetError(1, @error, False) : True
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_gbmv()
; Description ...: matrix-vector multiplication:      y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y
;                  where A is an packed(!) banded matrix
; Syntax ........: _blas_gbmv($mMatrix, $mVecX, $mVecY, [$fAlpha = 1.0, [$fBeta = 0, [$iKL = 0, [$iKU = 0, [$cTRANS = "N", [$iM = Default, [$iN = Default, [$iIncX = 1, [$iIncY = 1, [$iLDA = Default, [$sDataType = "DOUBLE"]]]]]]]]]]])
; Parameters ....: mMatrix   - [Map] packed(!) matrix A as a map, DllStruct or pointer
;                  mVecX     - [Map] vector x as a map, DllStruct or pointer
;                  mVecY     - [Map] vector y as a map, DllStruct or pointer (will be overwritten)
;                  fAlpha    - [Float] (Default: 1)
;                            ↳ scalar value alpha with which x is scaled
;                  fBeta     - [Float] (Default: 0)
;                            ↳ scalar value beta with which y is scaled
;                  iKL       - [Int] (Default: 0
;                            ↳ number of sub-diagonals (lower) of the matrix A
;                  iKU       - [Int] (Default: 0)
;                            ↳ number of super-diagonals (upper) of the matrix A
;                  cTRANS    - [Char] (Default: "N")
;                            ↳ "N":      x := A * x
;                              "T"/"C":  x := Aᵀ * x
;                  iM        - [Int] (Default: Default)
;                            ↳ number of rows in A
;                  iN        - [Int] (Default: Default)
;                            ↳ number of columns in A
;                  iIncX     - [Int] (Default: 1)
;                            ↳ storage spacing between elements in x
;                  iIncY     - [Int] (Default: 1)
;                            ↳ storage spacing between elements in y
;                  iLDA      - [Int] (Default: Default)
;                            ↳ first("leading") dimension size of A
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of gbmv (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/dd/df4/group__gbmv_ga7001c2a185bcc8a3b6731f5d1ea7093e.html#ga7001c2a185bcc8a3b6731f5d1ea7093e
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[2,0,0,0], [0,3,0,0], [0,0,4,0], [0,0,0,5]]", $__g_BLAS_STYPE_BAND + $__g_BLAS_STYPE_PACKED, "DOUBLE", 0, 0)
;                  Global $mX = _blas_fromArray("[3,2,5,2]")
;                  Global $mY = _blas_fromArray("[0,0,0]")
;                  _blas_gbmv($mA, $mX, $mY, 1, 0, 0, 0)
;                  _blas_display($mY, "y <- beta * y + alpha * A * x")
; ===============================================================================================================================
Func _blas_gbmv($mMatrix, $mVecX, $mVecY, $fAlpha = 1.0, $fBeta = 0, $iKL = 0, $iKU = 0, $cTRANS = "N", $iM = Default, $iN = Default, Const $iIncX = 1, Const $iIncY = 1, $iLDA = Default, $sDataType = "DOUBLE")
	Local $pM ; pointer to the data in memory

	; set char buffers and input healing
	$cTRANS = $cTRANS = "N" ? "N" : ($cTRANS = "T" ? "T" : "C")
	DllStructSetData($tBLASCHAR1, 1, $cTRANS)

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrix)
			$sDataType = $mMatrix.datatype
			If IsKeyword($iM)   = 1 Then $iM   = $cTRANS = "N" ? $mMatrix.rows : $mMatrix.cols
			If IsKeyword($iN)   = 1 Then $iN   = $cTRANS = "N" ? $mMatrix.cols : $mMatrix.rows
			If IsKeyword($iLDA) = 1 Then $iLDA = $iKL + $iKU + 1
			$pM = $mMatrix.ptr
		Case IsPtr($mMatrix)
			$pM = $mMatrix
		Case IsDllStruct($mMatrix)
			$pM = DllStructGetPtr($mMatrix)
	EndSelect
	Local $pX = IsMap($mVecX) ? $mVecX.ptr : (IsPtr($mVecX) ? $mVecX : DllStructGetPtr($mVecX)), _
	      $pY = IsMap($mVecY) ? $mVecY.ptr : (IsPtr($mVecY) ? $mVecY : DllStructGetPtr($mVecY))

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d", _
                $sType   = $sDataType & "*"


	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "gbmv", _
		"PTR",  $pBLASCHAR1, _    ; TRANS -> 'N': y := alpha*A*x + beta*y, 'T'/'C': y := alpha*A**T*x + beta*y
		"INT*", $iM, _            ; M
		"INT*", $iN, _            ; N
		"INT*", $iKL, _           ; KL
		"INT*", $iKU, _           ; KU
		$sType, $fAlpha, _        ; Alpha
		"PTR",  $pM, _            ; A
		"INT*", $iLDA, _          ; lda: length of a column/row (leading dimension)
		"PTR",  $pX, _            ; X
		"INT*", $iIncX, _         ; the stride for vector x
		$sType, $fBeta, _         ; BETA
		"PTR",  $pY, _            ; Y
		"INT*", $iIncY _          ; the stride for vector y
	)
	Return @error ? SetError(1, @error, False) : True
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_sbmv()
; Description ...: matrix-vector multiplication:      alpha*A*x + beta*y --> y
;                  where A is an packed(!) symmetric band matrix
; Syntax ........: _blas_sbmv($mMatrix, $mVecX, $mVecY, [$fAlpha = 1.0, [$fBeta = 1.0, [$iK = 0, [$cUPLO = "U", [$iN = Default, [$iIncX = 1, [$iIncY = 1, [$iLDA = Default, [$sDataType = "DOUBLE"]]]]]]]]])
; Parameters ....: mMatrix   - [Map] packed(!) matrix A as a map, DllStruct or pointer
;                  mVecX     - [Map] vector x as a map, DllStruct or pointer
;                  mVecY     - [Map] vector y as a map, DllStruct or pointer (will be overwritten)
;                  fAlpha    - [Float] (Default: 1)
;                            ↳ scalar value alpha with which x is scaled
;                  fBeta     - [Float] (Default: 0)
;                            ↳ scalar value beta with which y is scaled
;                  iK        - [Int] (Default: 0)
;                            ↳ number of super-diagonals of the matrix A
;                  cUPLO     - [Char] (Default: "U")
;                            ↳ "U": upper triangular parts of A are used
;                              "L": lower triangular parts of A are used
;                  iN        - [Int] (Default: Default)
;                            ↳ order of matrix A
;                  iIncX     - [Int] (Default: 1)
;                            ↳ storage spacing between elements in x
;                  iIncY     - [Int] (Default: 1)
;                            ↳ storage spacing between elements in y
;                  iLDA      - [Int] (Default: Default)
;                            ↳ first("leading") dimension size of A
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of sbmv (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/da/dd4/group__hbmv_ga3ddb7bcc544c5881a9365cb78d642228.html#ga3ddb7bcc544c5881a9365cb78d642228
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[2,3,0,0],[3,4,5,0],[0,5,6,7],[0,0,7,8]]", $__g_BLAS_STYPE_BAND + $__g_BLAS_STYPE_PACKED + $__g_BLAS_STYPE_SYMMETRIC + $__g_BLAS_STYPE_UPPER, "DOUBLE", 1) ; packed symmetric upper band matrix
;                  Global $mX = _blas_fromArray("[2,3,4,5]")
;                  Global $mY = _blas_fromArray("[6,7,8,9]")
;                  _blas_sbmv($mA, $mX, $mY, 1, 1, 1)
;                  _blas_display($mY, "y <- beta * y + alpha * A * x")
; ===============================================================================================================================
Func _blas_sbmv($mMatrix, $mVecX, $mVecY, $fAlpha = 1.0, $fBeta = 1.0, $iK = 0, $cUPLO = "U", $iN = Default, Const $iIncX = 1, Const $iIncY = 1, $iLDA = Default, $sDataType = "DOUBLE")
	Local $pM ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrix)
			$sDataType = $mMatrix.datatype
			If IsKeyword($iN) = 1 Then $iN = $mMatrix.rows
			$pM = $mMatrix.ptr
		Case IsPtr($mMatrix)
			$pM = $mMatrix
		Case IsDllStruct($mMatrix)
			$pM = DllStructGetPtr($mMatrix)
	EndSelect
	Local $pX = IsMap($mVecX) ? $mVecX.ptr : (IsPtr($mVecX) ? $mVecX : DllStructGetPtr($mVecX)), _
	      $pY = IsMap($mVecY) ? $mVecY.ptr : (IsPtr($mVecY) ? $mVecY : DllStructGetPtr($mVecY))

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d", _
                $sType   = $sDataType & "*"

	If IsKeyword($iLDA) = 1 Then $iLDA = $iK + 1

	; set char buffers and input healing
	DllStructSetData($tBLASCHAR1, 1, $cUPLO = "U" ? "U" : "L")

	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "sbmv", _
		"PTR",  $pBLASCHAR1, _    ; UPLO -> 'U': A = upper triangular, 'L': A = lower triangular
		"INT*", $iN, _            ; N
		"INT*", $iK, _            ; K
		$sType, $fAlpha, _        ; Alpha
		"PTR",  $pM, _            ; A
		"INT*", $iLDA, _          ; lda: length of a column/row (leading dimension)
		"PTR",  $pX, _            ; X
		"INT*", $iIncX, _         ; the stride for vector x
		$sType, $fBeta, _         ; BETA
		"PTR",  $pY, _            ; Y
		"INT*", $iIncY _          ; the stride for vector y
	)
	Return @error ? SetError(1, @error, False) : True
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_tbsv()
; Description ...: solves one of the systems of equations:    A * x = b,   or   Aᵀ * x = b
;                  Where A is an packed(!) lower or upper triangular band matrix
; Syntax ........: _blas_tbsv($mMatrix, $mVecX, [$iK = 0, [$cUPLO = "U", [$cDIAG = "N", [$cTRANS = "N", [$iN = Default, [$iIncX = 1, [$iLDA = Default, [$sDataType = "DOUBLE"]]]]]]]])
; Parameters ....: mMatrix   - [Map] packed(!) matrix A as a map, DllStruct or pointer
;                  mVecX     - [Map] vector b as a map, DllStruct or pointer (will be overwritten)
;                  iK        - [Int] (Default: 0)
;                            ↳ number of super-diagonals of the matrix A
;                  cUPLO     - [Char] (Default: "U")
;                            ↳ "U": A is an upper triangular band matrix
;                              "L": A is an upper triangular band matrix
;                  cDIAG     - [Char] (Default: "N")
;                            ↳ "U": A is a unit triangular matrix (diagonal = identity = 1)
;                              "N": A is not a unit triangular matrix
;                  cTRANS    - [Char] (Default: "N")
;                            ↳ "N":      A   * x = b
;                              "T"/"C":  Aᵀ * x = b
;                  iN        - [Int] (Default: Default)
;                            ↳ order of matrix A
;                  iIncX     - [Int] (Default: 1)
;                            ↳ storage spacing between elements in x
;                  iLDA      - [Int] (Default: Default)
;                            ↳ first("leading") dimension size of A
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of tbsv (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d4/dcd/group__tbsv_gaf2fd78c89c68cdbc9f697ea37c075ece.html#gaf2fd78c89c68cdbc9f697ea37c075ece
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[4,2,1],[0,3,5],[0,0,2]]", $__g_BLAS_STYPE_BAND + $__g_BLAS_STYPE_PACKED + $__g_BLAS_STYPE_UPPER, "DOUBLE", 2, 0) ; must be in special compressed upper triangular band matrix format
;                  Global $mX = _blas_fromArray("[11,16,4]")
;                  _blas_tbsv($mA, $mX, 2)
;                  _blas_display($mX, "solve(A*x = b)")
; ===============================================================================================================================
Func _blas_tbsv($mMatrix, $mVecX, $iK = 0, $cUPLO = "U", $cDIAG = "N", $cTRANS = "N", $iN = Default, Const $iIncX = 1, $iLDA = Default, $sDataType = "DOUBLE")
	Local $pM ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrix)
			$sDataType = $mMatrix.datatype
			If IsKeyword($iN)   = 1 Then $iN   = $cTRANS = "N" ? $mMatrix.rows : $mMatrix.cols
			If IsKeyword($iLDA) = 1 Then $iLDA = $cTRANS = "N" ? $mMatrix.rows : $mMatrix.cols
			$pM = $mMatrix.ptr
		Case IsPtr($mMatrix)
			$pM = $mMatrix
		Case IsDllStruct($mMatrix)
			$pM = DllStructGetPtr($mMatrix)
	EndSelect
	Local $pX = IsMap($mVecX) ? $mVecX.ptr : (IsPtr($mVecX) ? $mVecX : DllStructGetPtr($mVecX))

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; set char buffers and input healing
	DllStructSetData($tBLASCHAR1, 1, $cUPLO  = "U" ? "U" : "L")
	DllStructSetData($tBLASCHAR2, 1, $cTRANS = "N" ? "N" : ($cTRANS = "T" ? "T" :"C"))
	DllStructSetData($tBLASCHAR3, 1, $cDIAG  = "N" ? "N" : "U")

	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "tbsv", _
		"PTR",  $pBLASCHAR1, _   ; UPLO -> 'U': A = upper triangular, 'L': A = lower triangular
		"PTR",  $pBLASCHAR2, _   ; TRANS -> 'N': x = A*x, 'T': x = Aᵀ*x, 'C': x = Aᵀ*x
		"PTR",  $pBLASCHAR3, _   ; DIAG -> 'U': A = unit triangular, 'N': A = other elements on diagonal
		"INT*", $iN, _           ; N
		"INT*", $iK, _           ; Super-diagonals (if UPLO = U) or sub-diagonals (if UPLO = L)
		"ptr",  $pM, _           ; A
		"INT*", $iLDA, _         ; lda: length of a column/row (leading dimension)
		"ptr",  $pX, _           ; the vector x
		"INT*", $iIncX _         ; the stride for vector x
	)
	Return @error ? SetError(1, @error, False) : True
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_trsv()
; Description ...: solves one of the systems of equations:    A * x = b,   or   Aᵀ * x = b
;                  Where A is an unpacked lower or upper triangular matrix
; Syntax ........: _blas_trsv($mMatrix, $mVecX, [$cUPLO = "U", [$cDIAG = "N", [$cTRANS = "N", [$iN = Default, [$iIncX = 1, [$iLDA = Default, [$sDataType = "DOUBLE"]]]]]]])
; Parameters ....: mMatrix   - [Map] matrix A as a map, DllStruct or pointer
;                  mVecX     - [Map] vector b as a map, DllStruct or pointer (will be overwritten)
;                  cUPLO     - [Char] (Default: "U")
;                            ↳ "U": A is an upper triangular band matrix
;                              "L": A is an upper triangular band matrix
;                  cDIAG     - [Char] (Default: "N")
;                            ↳ "U": A is a unit triangular matrix (diagonal = identity = 1)
;                              "N": A is not a unit triangular matrix
;                  cTRANS    - [Char] (Default: "N")
;                            ↳ "N":      A   * x = b
;                              "T"/"C":  Aᵀ * x = b
;                  iN        - [Int] (Default: Default)
;                            ↳ order of matrix A
;                  iIncX     - [Int] (Default: 1)
;                            ↳ storage spacing between elements in x
;                  iLDA      - [Int] (Default: Default)
;                            ↳ first("leading") dimension size of A
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of trsv (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/dd/dc3/group__trsv_ga7a7dcbb8745b4776ce13063ab031141f.html#ga7a7dcbb8745b4776ce13063ab031141f
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[4,2,1],[0,3,5],[0,0,2]]", $__g_BLAS_STYPE_TRIANGLE + $__g_BLAS_STYPE_UPPER) ; "upper triangle" not necessary - just little speedup
;                  Global $mX = _blas_fromArray("[11,16,4]")
;                  _blas_trsv($mA, $mX)
;                  _blas_display($mX, "solve(A*x = b)") ; --> [1.25, 2, 2]
; ===============================================================================================================================
Func _blas_trsv($mMatrix, $mVecX, $cUPLO = "U", $cDIAG = "N", $cTRANS = "N", $iN = Default, Const $iIncX = 1, $iLDA = Default, $sDataType = "DOUBLE")
	Local $pM ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrix)
			$sDataType = $mMatrix.datatype
			If IsKeyword($iN)   = 1 Then $iN   = $cTRANS = "N" ? $mMatrix.rows : $mMatrix.cols
			If IsKeyword($iLDA) = 1 Then $iLDA = $cTRANS = "N" ? $mMatrix.rows : $mMatrix.cols
			$pM = $mMatrix.ptr
		Case IsPtr($mMatrix)
			$pM = $mMatrix
		Case IsDllStruct($mMatrix)
			$pM = DllStructGetPtr($mMatrix)
	EndSelect
	Local $pX = IsMap($mVecX) ? $mVecX.ptr : (IsPtr($mVecX) ? $mVecX : DllStructGetPtr($mVecX))

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; set char buffers and input healing
	DllStructSetData($tBLASCHAR1, 1, $cUPLO  = "U" ? "U" : "L")
	DllStructSetData($tBLASCHAR2, 1, $cTRANS = "N" ? "N" : ($cTRANS = "T" ? "T" :"C"))
	DllStructSetData($tBLASCHAR3, 1, $cDIAG  = "N" ? "N" : "U")

	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "trsv", _
		"PTR",  $pBLASCHAR1, _     ; UPLO -> 'U': A = upper triangular, 'L': A = lower triangular
		"PTR",  $pBLASCHAR2, _     ; TRANS -> 'N': x = A*x, 'T': x = Aᵀ*x, 'C': x = Aᵀ*x
		"PTR",  $pBLASCHAR3, _     ; DIAG -> 'U': A = unit triangular, 'N': A = other elements on diagonal
		"INT*", $iN, _             ; N
		"ptr",  $pM, _             ; A
		"INT*", $iLDA, _           ; lda: length of a column/row (leading dimension)
		"ptr",  $pX, _             ; the vector x
		"INT*", $iIncX _           ; the stride for vector x
	)
	Return @error ? SetError(1, @error, False) : True
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_ger()
; Description ...: calculate the rank 1 operation: A := alpha * x * yᵀ + A
; Syntax ........: _blas_ger($mA, $mVecX, $mVecY, [$fAlpha = 1.0, [$iM = Default, [$iN = Default, [$incX = 1, [$IncY = 1, [$iLDA = Default, [$sDataType = "DOUBLE"]]]]]]])
; Parameters ....: mA        - [Map] matrix A as a map, DllStruct or pointer (will be overwritten)
;                  mVecX     - [Map] vector x as a map, DllStruct or pointer
;                  mVecY     - [Map] vector x as a map, DllStruct or pointer
;                  fAlpha    - [Float] (Default: 1.0)
;                            ↳ scalar value to scale x with
;                  iM        - [Int] (Default: Default)
;                            ↳ number of rows in A
;                  iN        - [Int] (Default: Default)
;                            ↳ number of columns in A
;                  incX      - [Int] (Default: 1)
;                            ↳ storage spacing between elements in x
;                  IncY      - [Int] (Default: 1)
;                            ↳ storage spacing between elements in y
;                  iLDA      - [Int] (Default: Default)
;                            ↳ first("leading") dimension size of A
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of ger (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d8/d75/group__ger_gaef5d248da0fdfb62bccb259725935cb8.html#gaef5d248da0fdfb62bccb259725935cb8
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[1,2,3], [2,2,4], [3,2,2], [4,2,1]]")
;                  Global $mX = _blas_fromArray("[3,2,1,4]")
;                  Global $mY = _blas_fromArray("[1,2,3]")
;                  _blas_ger($mA, $mX, $mY)
;                  _blas_display($mA, "A <-- alpha*x*y'+ A")
; ===============================================================================================================================
Func _blas_ger($mA, $mVecX, $mVecY, Const $fAlpha = 1.0, $iM = Default, $iN = Default, Const $incX = 1, Const $IncY = 1, $iLDA = Default, $sDataType = "DOUBLE")
	Local $pA ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype

			; handle default parameter
			If IsKeyword($iM)   = 1 Then $iM   = $mA.rows
			If IsKeyword($iN)   = 1 Then $iN   = $mA.cols
			If IsKeyword($iLDA) = 1 Then $iLDA = $mA.rows

			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect
	Local $pX = IsMap($mVecX) ? $mVecX.ptr : (IsPtr($mVecX) ? $mVecX : DllStructGetPtr($mVecX))
	Local $pY = IsMap($mVecY) ? $mVecY.ptr : (IsPtr($mVecY) ? $mVecY : DllStructGetPtr($mVecY))

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "ger", _
		"INT*",           $iM, _      ; number of rows in matrix A
		"INT*",           $iN, _      ; number of columns in matrix A
		$sDataType & "*", $fAlpha, _  ; the scaling constant alpha for
		"ptr",            $pX, _      ; the vector x
		"INT*",           $incX, _    ; the stride for vector x
		"ptr",            $pY, _      ; the vector y
		"INT*",           $IncY, _    ; the stride for vector y
		"ptr",            $pA, _      ; A
		"INT*",           $iLDA _     ; lda: length of a column in a (used to calc the element index - can be used to change indexing behavior)
	)
	Return @error ? SetError(1, @error, False) : True
EndFunc

#EndRegion

#Region BLAS Level 3 (Matrix-Matrix Operations)

; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_gemm()
; Description ...: matrix multiplication:           C := alpha * op(A) * op(B) + beta * C   ; op(X) is X or Xᵀ
;                  where A, B and C are general m × n matrices
; Syntax ........: _blas_gemm($mMatrixA, $mMatrixB, $mMatrixC, [$fAlpha = 1, [$fBeta = 0, [$cTRANSA = "N", [$cTransposedB = "N", [$iM = Default, [$iN = Default, [$iK = Default, [$iLDA = Default, [$iLDB = Default, [$iLDC = Default, [$sDataType = "DOUBLE"]]]]]]]]]]])
; Parameters ....: mMatrixA     - [Map] matrix A as a map, DllStruct or pointer
;                  mMatrixB     - [Map] matrix B as a map, DllStruct or pointer
;                  mMatrixC     - [Map] matrix C as a map, DllStruct or pointer (will be overwritten)
;                  fAlpha       - [Float] (Default: 1)
;                               ↳ scalar value to scale A with
;                  fBeta        - [Float] (Default: 0)
;                               ↳ scalar value to scale C with
;                  cTRANSA - [Char] (Default: "N")
;                               ↳ "N":      op(A) = A
;                                 "T"/"C":  op(A) = Aᵀ
;                  cTransposedB - [Char] (Default: "N")
;                               ↳ "N":      op(B) = B
;                                 "T"/"C":  op(B) = Bᵀ
;                  iM           - [Int] (Default: Default)
;                               ↳ number of rows of op(A) and C
;                  iN           - [Int] (Default: Default)
;                               ↳ number columns of op(B) and C
;                  iK           - [Int] (Default: Default)
;                               ↳ number of columns of op(A) and rows of op(B)
;                  iLDA         - [Int] (Default: Default)
;                               ↳ first("leading") dimension size of A
;                  iLDB         - [Int] (Default: Default)
;                               ↳ first("leading") dimension size of B
;                  iLDC         - [Int] (Default: Default)
;                               ↳ first("leading") dimension size of C
;                  sDataType    - [String] (Default: "DOUBLE")
;                               ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of gemm (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/dd/d09/group__gemm_ga1e899f8453bcbfde78e91a86a2dab984.html#ga1e899f8453bcbfde78e91a86a2dab984
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[1, -3], [2, 4], [1, -1]]")
;                  Global $mB = _blas_fromArray("[[1, -3], [2, 4], [1, -1]]")
;                  Global $mC = _blas_fromArray("[[0.5, 0.5, 0.5], [0.5, 0.5, 0.5], [0.5, 0.5, 0.5]]")
;                  _blas_gemm($mA, $mB, $mC, 1, 2, "N", "T")
;                  _blas_display($mC, "C <- alpha * A * B + beta * C")
; ===============================================================================================================================
Func _blas_gemm($mMatrixA, $mMatrixB, $mMatrixC, $fAlpha = 1, $fBeta = 0, $cTRANSA = "N", $cTransposedB = "N", $iM = Default, $iN = Default, $iK = Default, $iLDA = Default, $iLDB = Default, $iLDC = Default, $sDataType = "DOUBLE")
	Local $pA, $pB, $pC ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrixA)
			$sDataType = $mMatrixA.datatype
			If IsKeyword($iM)   = 1 Then $iM   = $cTRANSA = "N" ? $mMatrixA.rows : $mMatrixA.cols
			If IsKeyword($iK)   = 1 Then $iK   = $cTRANSA = "N" ? $mMatrixA.cols : $mMatrixA.rows
			If IsKeyword($iLDA) = 1 Then $iLDA = $cTRANSA = "N" ? $iM : $iK
			$pA = $mMatrixA.ptr
		Case IsPtr($mMatrixA)
			$pA = $mMatrixA
		Case IsDllStruct($mMatrixA)
			$pA = DllStructGetPtr($mMatrixA)
	EndSelect
	Select
		Case IsMap($mMatrixB)
			If IsKeyword($iN)   = 1 Then $iN   = $cTransposedB = "N" ? $mMatrixB.cols : $mMatrixB.rows
			If IsKeyword($iK)   = 1 Then $iK   = $cTransposedB = "N" ? $mMatrixB.rows : $mMatrixB.cols
			If IsKeyword($iLDB) = 1 Then $iLDB = $cTransposedB = "N" ? $iK : $iN
			$pB = $mMatrixB.ptr
		Case IsPtr($mMatrixB)
			$pB = $mMatrixB
		Case IsDllStruct($mMatrixB)
			$pB = DllStructGetPtr($mMatrixB)
	EndSelect
	Select
		Case IsMap($mMatrixC)
			If IsKeyword($iN)   = 1 Then $iN   = $mMatrixC.cols
			If IsKeyword($iK)   = 1 Then $iK   = $mMatrixC.rows
			If IsKeyword($iLDC) = 1 Then $iLDC = $iM
			$pC = $mMatrixC.ptr
		Case IsPtr($mMatrixC)
			$pC = $mMatrixC
		Case IsDllStruct($mMatrixC)
			$pC = DllStructGetPtr($mMatrixC)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d", _
	            $sType   = $sDataType & "*"

	; set char buffers
	DllStructSetData($tBLASCHAR1, 1, $cTRANSA = "T" ? "T" : "N")
	DllStructSetData($tBLASCHAR2, 1, $cTransposedB = "T" ? "T" : "N")

	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "gemm", _
		"PTR",  $pBLASCHAR1, _   ; char *transa  ( If "N" then A is used; if "T" then Aᵀ is used)
		"PTR",  $pBLASCHAR2, _   ; char *transb  ( If "N" then B is used; if "T" then Bᵀ is used)
		"INT*", $iM, _           ; M - number of rows in matrix A and matrix C
		"INT*", $iN, _           ; N - number of columns in matrix B and matrix C
		"INT*", $iK, _           ; K - number of columns in matrix A and number of rows in matrix C
		$sType, $fAlpha, _       ; the scaling constant alpha
		"ptr",  $pA, _           ; A
		"INT*", $iLDA, _         ; lda
		"ptr",  $pB, _           ; B
		"INT*", $iLDB, _         ; ldb
		$sType, $fBeta, _        ; the scaling constant beta
		"ptr",  $pC, _           ; C
		"INT*", $iLDC _          ; ldc
	)
	Return @error ? SetError(1, @error, False) : True
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_symm()
; Description ...: matrix multiplication:           C := alpha * A * B + beta * C   or C := alpha * B * A + beta * C
;                  where A is an unpacked symmetric matrix, B and C are general m × n matrices
; Syntax ........: _blas_symm($mMatrixA, $mMatrixB, $mMatrixC, [$fAlpha = 1, [$fBeta = 0, [$cUPLO = "U", [$cSIDE = "L", [$iM = Default, [$iN = Default, [$iLDA = Default, [$iLDB = Default, [$iLDC = Default, [$sDataType = "DOUBLE"]]]]]]]]]])
; Parameters ....: mMatrixA  - [Map] unpacked symmatric matrix A as a map, DllStruct or pointer
;                  mMatrixB  - [Map] matrix B as a map, DllStruct or pointer
;                  mMatrixC  - [Map] matrix A as a map, DllStruct or pointer (will be overwritten)
;                  fAlpha    - [Float] (Default: 1)
;                            ↳ scalar value to scale A/B with
;                  fBeta     - [Float] (Default: 0)
;                            ↳ scalar value to scale C with
;                  cUPLO     - [Char] (Default: "U")
;                            ↳ "U": values are stored in upper triangular part of A
;                              "L": values are stored in lower triangular part of A
;                  cSIDE     - [Char] (Default: "L")
;                            ↳ "L": C := alpha * A * B + beta * C
;                              "R": C := alpha * B * A + beta * C
;                  iM        - [Int] (Default: Default)
;                            ↳ numer of rows of C
;                  iN        - [Int] (Default: Default)
;                            ↳ number of columns of C
;                  iLDA      - [Int] (Default: Default)
;                            ↳ first("leading") dimension size of A
;                  iLDB      - [Int] (Default: Default)
;                            ↳ first("leading") dimension size of B
;                  iLDC      - [Int] (Default: Default)
;                            ↳ first("leading") dimension size of C
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: @error ? SetError(1, @error, False) : True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of symm (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d0/d16/group__hemm_ga83617b92d007d90eda87d135bf08175b.html#ga83617b92d007d90eda87d135bf08175b
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[8,1,2,3], [0,5,7,9], [0,0,3,7], [0,0,0,8]]")
;                  Global $mB = _blas_fromArray("[[1,2,3],[4,5,6],[7,8,9],[10,11,12]]")
;                  Global $mC = _blas_createMatrix(4, 3, $mA.datatype) ; empty because beta = 0: means only A*B
;                  _blas_symm($mA, $mB, $mC, 1, 0, "U", "L")
;                  _blas_display($mC, "C <- alpha * A * B    (A = symmetric)")
; ===============================================================================================================================
Func _blas_symm($mMatrixA, $mMatrixB, $mMatrixC, $fAlpha = 1, $fBeta = 0, $cUPLO = "U", $cSIDE = "L", $iM = Default, $iN = Default, $iLDA = Default, $iLDB = Default, $iLDC = Default, $sDataType = "DOUBLE")
	Local $pA, $pB, $pC ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrixC)
			If IsKeyword($iM)   = 1 Then $iM   = $mMatrixC.rows
			If IsKeyword($iN)   = 1 Then $iN   = $mMatrixC.cols
			If IsKeyword($iLDC) = 1 Then $iLDC = $iM
			$pC = $mMatrixC.ptr
		Case IsPtr($mMatrixC)
			$pC = $mMatrixC
		Case IsDllStruct($mMatrixC)
			$pC = DllStructGetPtr($mMatrixC)
	EndSelect
	Select
		Case IsMap($mMatrixB)
			If IsKeyword($iM)   = 1 Then $iM   = $mMatrixB.rows
			If IsKeyword($iN)   = 1 Then $iN   = $mMatrixB.cols
			If IsKeyword($iLDB) = 1 Then $iLDB = $iM
			$pB = $mMatrixB.ptr
		Case IsPtr($mMatrixB)
			$pB = $mMatrixB
		Case IsDllStruct($mMatrixB)
			$pB = DllStructGetPtr($mMatrixB)
	EndSelect
	Select
		Case IsMap($mMatrixA)
			$sDataType = $mMatrixA.datatype
			If IsKeyword($iM)   = 1 Then $iM   = $mMatrixA.rows
			If IsKeyword($iN)   = 1 Then $iN   = $mMatrixA.cols
			If IsKeyword($iLDA) = 1 Then $iLDA = $cSIDE = "L" ? $iM : $iN
			$pA = $mMatrixA.ptr
		Case IsPtr($mMatrixA)
			$pA = $mMatrixA
		Case IsDllStruct($mMatrixA)
			$pA = DllStructGetPtr($mMatrixA)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d", _
	            $sType   = $sDataType & "*"

	; set char buffers
	DllStructSetData($tBLASCHAR1, 1, $cSIDE = "L" ? "L" : "R")
	DllStructSetData($tBLASCHAR2, 1, $cUPLO = "U" ? "U" : "L")

	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "symm", _
		"PTR",  $pBLASCHAR1, _   ; SIDE
		"PTR",  $pBLASCHAR2, _   ; UPLO
		"INT*", $iM, _           ; M - number of rows in matrix A and matrix C
		"INT*", $iN, _           ; N - number of columns in matrix B and matrix C
		$sType, $fAlpha, _       ; the scaling constant alpha
		"ptr",  $pA, _           ; the general matrix A
		"INT*", $iLDA, _         ; lda
		"ptr",  $pB, _           ; the general matrix B
		"INT*", $iLDB, _         ; ldb
		$sType, $fBeta, _        ; the scaling constant beta
		"ptr",  $pC, _           ; the general matrix C
		"INT*", $iLDC _          ; ldc
	)
	Return @error ? SetError(1, @error, False) : True
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_trmm()
; Description ...: matrix multiplication:   B := alpha * op(A) * B   or   B := alpha * B * op(A)      where op(A) = A or Aᵀ
;                  where A is an unpacked lower or upper triangular matrix, B is a general m × n matrix
; Syntax ........: _blas_trmm($mA, $mB, [$fAlpha = 1, [$cSide = "L", [$cUplo = "U", [$cDiag = "N", [$cTRANSA = "N", [$iM = Default, [$iN = Default, [$iLDA = Default, [$iLDB = Default, [$sDataType = "DOUBLE"]]]]]]]]]])
; Parameters ....: mA           - [Map] unpacked triangular matrix A as a map, DllStruct or pointer
;                  mB           - [Map] matrix B as a map, DllStruct or pointer (will be overwritten)
;                  fAlpha       - [Float] (Default: 1)
;                               ↳ scalar value to scale A/B with
;                  cSide        - [Char] (Default: "L")
;                               ↳ "L": B := alpha * op(A) * B
;                                 "R": B := alpha * B * op(A)
;                  cUplo        - [Char] (Default: "U")
;                               ↳ "U": A is an upper triangular matrix
;                                 "L": A is an lower triangular matrix
;                  cDiag        - [Char] (Default: "N")
;                               ↳ "U": A is a unit triangular matrix (diagonal = identity = 1)
;                                 "N": A is not a unit triangular matrix
;                  cTRANSA      - [Char] (Default: "N")
;                               ↳ "N":      op(A) = A
;                                 "T"/"C":  op(A) = Aᵀ
;                  iM           - [Int] (Default: Default)
;                               ↳ number of rows of B
;                  iN           - [Int] (Default: Default)
;                               ↳ number of columns of B
;                  iLDA         - [Int] (Default: Default)
;                               ↳ first("leading") dimension size of A
;                  iLDB         - [Int] (Default: Default)
;                               ↳ first("leading") dimension size of B
;                  sDataType    - [String] (Default: "DOUBLE")
;                               ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of trmm (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/dd/dab/group__trmm_ga4d2f76d6726f53c69031a2fe7f999add.html#ga4d2f76d6726f53c69031a2fe7f999add
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[1,2,3],[0,4,5],[0,0,6]]")
;                  Global $mB = _blas_fromArray("[[1,2,3],[4,5,6],[7,8,9]]")
;                  _blas_trmm($mA, $mB, 1)
;                  _blas_display($mB, "B <-- A * B with A = triangular")
; ===============================================================================================================================
Func _blas_trmm($mA, $mB, $fAlpha = 1, $cSide = "L", $cUplo = "U", $cDiag = "N", $cTRANSA = "N", $iM = Default, $iN = Default, $iLDA = Default, $iLDB = Default, $sDataType = "DOUBLE")
	Local $pA, $pB ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iLDA) = 1 Then $iLDA = $cTRANSA = "T" ? $mA.cols : $mA.rows
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect
	Select
		Case IsMap($mB)
			If IsKeyword($iM)   = 1 Then $iM   = $mB.rows
			If IsKeyword($iN)   = 1 Then $iN   = $mB.cols
			$pB = $mB.ptr
		Case IsPtr($mB)
			$pB = $mB
		Case IsDllStruct($mB)
			$pB = DllStructGetPtr($mB)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	If IsKeyword($iLDA) = 1 Then $iLDA = $cTRANSA ? $iM : $iN
	If IsKeyword($iLDB) = 1 Then $iLDB = $iM

	; set char buffers
	DllStructSetData($tBLASCHAR1, 1, $cSide = "R" ? "R" : "L")
	DllStructSetData($tBLASCHAR2, 1, $cUplo = "L" ? "L" : "U")
	DllStructSetData($tBLASCHAR3, 1, $cTRANSA = "T" ? "T" : "N")
	DllStructSetData($tBLASCHAR4, 1, $cDiag = "U" ? "U" : "N")

	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "trmm", _
		"PTR",            $pBLASCHAR1, _    ; char *side (if "R" formula changes to B <- alpha * B . A[ᵀ])
		"PTR",            $pBLASCHAR2, _    ; char *uplo (indicates if A is a upper or lower triangular matrix)
		"PTR",            $pBLASCHAR3, _    ; if "T" the formula is used:      alpha*Aᵀ.A + beta * C   if "N": alpha*A.Aᵀ + beta * C
		"PTR",            $pBLASCHAR4, _    ; "U": A is a unit triangular matrix   "N": A is not a unit triangular matrix
		"INT*",           $iM, _            ; m: number of rows in matrix B
		"INT*",           $iN, _            ; n: number of columns in matrix B
		$sDataType & "*", $fAlpha, _        ; the scaling constant alpha
		"ptr",            $pA, _            ; the triangular (unpacked) matrix A
		"INT*",           $iLDA, _          ; lda
		"ptr",            $pB, _            ; the general matrix B
		"INT*",           $iLDB _           ; ldb
	)
	Return @error ? SetError(1, @error, False) : True
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_trsm()
; Description ...: solves one of the matrix equations: op(A) * X = alpha * B,   or   X * op(A) = alpha * B     where op(A) = A or Aᵀ
;                  where A is an unpacked lower or upper triangular matrix, B is a general m × n matrix
; Syntax ........: _blas_trsm($mA, $mB, [$fAlpha = 1.0, [$cSIDE = "L", [$cUPLO = "U", [$cDIAG = "N", [$cTRANSA = "N", [$iM = Default, [$iN = Default, [$iLDA = Default, [$iLDB = Default, [$sDataType = "DOUBLE"]]]]]]]]]])
; Parameters ....: mA        - [Map] unpacked triangular matrix A as a map, DllStruct or pointer
;                  mB        - [Map] unpacked triangular matrix/vector B as a map, DllStruct or pointer (will be overwritten)
;                  fAlpha    - [Float] (Default: 1.0)
;                            ↳ scalar value to scale B with
;                  cSIDE     - [Char] (Default: "L")
;                            ↳ "L": op(A) * X = alpha * B
;                              "R": X * op(A) = alpha * B
;                  cUPLO     - [Char] (Default: "U")
;                            ↳ "U": A is an upper triangular matrix
;                              "L": A is an lower triangular matrix
;                  cDIAG     - [Char] (Default: "N")
;                            ↳ "U": A is a unit triangular matrix (diagonal = identity = 1)
;                              "N": A is not a unit triangular matrix
;                  cTRANSA   - [Char] (Default: "N")
;                            ↳ "N":      op(A) = A
;                              "T"/"C":  op(A) = Aᵀ
;                  iM        - [Int] (Default: Default)
;                            ↳ number of rows of B
;                  iN        - [Int] (Default: Default)
;                            ↳ number of columns of B
;                  iLDA      - [Int] (Default: Default)
;                            ↳ first("leading") dimension size of A
;                  iLDB      - [Int] (Default: Default)
;                            ↳ first("leading") dimension size of B
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of trsm (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d9/de5/group__trsm_ga7120d931d7b1a15e12d50d328799df8a.html#ga7120d931d7b1a15e12d50d328799df8a
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[4,2,1],[0,3,5],[0,0,2]]", $__g_BLAS_STYPE_TRIANGLE + $__g_BLAS_STYPE_UPPER) ; "upper triangle" not necessary - just little speedup
;                  Global $mX = _blas_fromArray("[[11,16,4],[5,9,12],[19,12,14]]")
;                  _blas_trsm($mA, $mX)
;                  _blas_display($mX)
; ===============================================================================================================================
Func _blas_trsm($mA, $mB, $fAlpha = 1.0, $cSIDE = "L", $cUPLO = "U", $cDIAG = "N", $cTRANSA = "N", $iM = Default, $iN = Default, $iLDA = Default, $iLDB = Default, $sDataType = "DOUBLE")
	Local $pA, $pB ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mB)
			$sDataType = $mB.datatype
			If IsKeyword($iM)   = 1 Then $iM   = $mB.rows
			If IsKeyword($iN)   = 1 Then $iN   = $mB.cols
			If IsKeyword($iLDB) = 1 Then $iLDB = $iM
			$pB = $mB.ptr
		Case IsPtr($mB)
			$pB = $mB
		Case IsDllStruct($mB)
			$pB = DllStructGetPtr($mB)
	EndSelect
	Select
		Case IsMap($mA)
			$sDataType = $mA.datatype
			If IsKeyword($iLDA) = 1 Then $iLDA = $cSIDE = "L" ? $iM : $iN
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; set char buffers and input healing
	DllStructSetData($tBLASCHAR1, 1, $cSIDE  = "L" ? "L" : "R")
	DllStructSetData($tBLASCHAR2, 1, $cUPLO  = "U" ? "U" : "L")
	DllStructSetData($tBLASCHAR3, 1, $cTRANSA = "N" ? "N" : ($cTRANSA = "T" ? "T" :"C"))
	DllStructSetData($tBLASCHAR4, 1, $cDIAG  = "N" ? "N" : "U")

	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "trsm", _
		"PTR",            $pBLASCHAR1, _     ; SIDE   -> 'L': A*X = B,   'R': X*A = B
		"PTR",            $pBLASCHAR2, _     ; UPLO   -> 'U': A = upper triangular, 'L': A = lower triangular
		"PTR",            $pBLASCHAR3, _     ; TRANSA -> 'N': x = A*x, 'T': x = Aᵀ*x, 'C': x = Aᵀ*x
		"PTR",            $pBLASCHAR4, _     ; DIAG   -> 'U': A = unit triangular, 'N': A = other elements on diagonal
		"INT*",           $iM, _             ; M
		"INT*",           $iN, _             ; N
		$sDataType & "*", $fAlpha, _         ; ALPHA
		"ptr",            $pA, _             ; A
		"INT*",           $iLDA, _           ; LDB: length of a column/row (leading dimension)
		"ptr",            $pB, _             ; B
		"INT*",           $iLDB _            ; LDB: length of a column/row (leading dimension)
	)
	Return @error ? SetError(1, @error, False) : True
EndFunc


; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_syrk()
; Description ...: symmetric rank k operations: C := alpha * A * Aᵀ + beta * C    or    alpha * Aᵀ * A + beta * C
;                  where A is an unpacked symmetric matrix, C is a general m × n matrix
; Syntax ........: _blas_syrk($mMatrixA, $mMatrixC, [$fAlpha = 1, [$fBeta = 0, [$cUPLO = "U", [$cTRANS = "N", [$iN = Default, [$iK = Default, [$iLDA = Default, [$iLDC = Default, [$sDataType = "DOUBLE"]]]]]]]]])
; Parameters ....: mMatrixA  - [Map] unpacked symmetric matrix A as a map, DllStruct or pointer
;                  mMatrixC  - [Map] matrix C as a map, DllStruct or pointer (will be overwritten)
;                  fAlpha    - [Float] (Default: 1)
;                            ↳ scalar value to scale A with
;                  fBeta     - [Float] (Default: 0)
;                            ↳ scalar value to scale C with
;                  cUPLO     - [Char] (Default: "U")
;                            ↳ "U": values are stored in upper triangular part of A
;                              "L": values are stored in lower triangular part of A
;                  cTRANS    - [Char] (Default: "N")
;                            ↳ "N":      op(A) = A
;                              "T"/"C":  op(A) = Aᵀ
;                  iN        - [Int] (Default: Default)
;                            ↳ order of matrix C
;                  iK        - [Int] (Default: Default)
;                            ↳ number of columns in A
;                  iLDA      - [Int] (Default: Default)
;                            ↳ first("leading") dimension size of A
;                  iLDC      - [Int] (Default: Default)
;                            ↳ first("leading") dimension size of C
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of syrk (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d4/d6e/group__herk_ga09ec411a6b845c47efb0bc8bd3283b03.html#ga09ec411a6b845c47efb0bc8bd3283b03
; Example .......: Yes
;                  Global $mA = _blas_fromArray("[[8,1,2,3], [1,5,7,9], [2,7,3,7], [3,9,7,8]]")
;                  Global $mC = _blas_createMatrix(4, 4, $mA.datatype)
;                  _blas_syrk($mA, $mC)
;                  $mC.storageType = BitOR($mC.storageType, $__g_BLAS_STYPE_SYMMETRIC + $__g_BLAS_STYPE_UPPER) ; C is symmetric and only upper part is stored here
;                  _blas_display($mC, "C <- A * Aᵀ (result is symmetric)")
; ===============================================================================================================================
Func _blas_syrk($mMatrixA, $mMatrixC, $fAlpha = 1, $fBeta = 0, $cUPLO = "U", $cTRANS = "N", $iN = Default, $iK = Default, $iLDA = Default, $iLDC = Default, $sDataType = "DOUBLE")
	Local $pA, $pC ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrixA)
			$sDataType = $mMatrixA.datatype
			If IsKeyword($iN)   = 1 Then $iN   = $cTRANS = "N" ? $mMatrixA.rows : $mMatrixA.cols
			If IsKeyword($iK)   = 1 Then $iK   = $cTRANS = "N" ? $mMatrixA.cols : $mMatrixA.rows
			If IsKeyword($iLDA) = 1 Then $iLDA = $cTRANS = "N" ? $mMatrixA.rows : $mMatrixA.cols
			$pA = $mMatrixA.ptr
		Case IsPtr($mMatrixA)
			$pA = $mMatrixA
		Case IsDllStruct($mMatrixA)
			$pA = DllStructGetPtr($mMatrixA)
	EndSelect
	Select
		Case IsMap($mMatrixC)
			If IsKeyword($iLDC) = 1 Then $iLDC = $iN
			$pC = $mMatrixC.ptr
		Case IsPtr($mMatrixC)
			$pC = $mMatrixC
		Case IsDllStruct($mMatrixC)
			$pC = DllStructGetPtr($mMatrixC)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d", _
	            $sType   = $sDataType & "*"

	; set char buffers
	DllStructSetData($tBLASCHAR1, 1, $cUPLO = "U" ? "U" : "L")
	DllStructSetData($tBLASCHAR2, 1, $cTRANS = "N" ? "N" : ($cTRANS = "T" ? "T" : "C"))

	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "syrk", _
		"PTR",  $pBLASCHAR1, _   ; UPLO
		"PTR",  $pBLASCHAR2, _   ; TRANS
		"INT*", $iN, _           ; N 
		"INT*", $iK, _           ; K
		$sType, $fAlpha, _       ; the scaling constant alpha
		"ptr",  $pA, _           ; A
		"INT*", $iLDA, _         ; lda
		$sType, $fBeta, _        ; the scaling constant beta
		"ptr",  $pC, _           ; C
		"INT*", $iLDC _          ; ldc
	)
	Return @error ? SetError(1, @error, False) : True
EndFunc

#EndRegion


#Region additional functions (depends on library)

; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_tpttr()
; Description ...: copies a triangular matrix from the standard packed format (TP) to the standard full format (TR)
; Syntax ........: _blas_tpttr($mAP, [$cUPLO = Default, [$mA = Default, [$iN = Default, [$iLDA = Default, [$bRetMap = True, [$sDataType = "DOUBLE"]]]]]])
; Parameters ....: mAP       - [Map] packed upper or lower triangular matrix as a map, DllStruct or pointer (will be overwritten)
;                  cUPLO     - [Char] (Default: Default)
;                            ↳ "U": A is an upper triangular matrix
;                              "L": A is an lower triangular matrix
;                  mA        - [Map] (Default: Default)
;                            ↳ target map where the values should be written
;                              if Default: function returns a new matrix
;                  iN        - [Int] (Default: Default)
;                            ↳ order of the matrix A (number of rows)
;                  iLDA      - [Int] (Default: Default)
;                            ↳ first("leading") dimenson size of A
;                  bRetMap   - [Bool] (Default: True)
;                            ↳ if true: function returns a unpacked upper or lower triangular matrix as a map
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: $bRetMap = True:  unpacked upper or lower triangular matrix as a map
;                           $bRetMap = False: True
;                  Failure: $bRetMap ? Null : False and set @error to:
;                           | 1: error during DllCall of tpttr (@extended: @error from DllCall)
;                           | 2: error inside tpttr (@extended: INFO-value from tpttr)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/d2/d04/group__tpttr_gad8388ffa6869aea9194cc0e2070c731b.html#gad8388ffa6869aea9194cc0e2070c731b
; Example .......: Yes
;                  Global $mTP = _blas_fromArray('[[611,197,-192,407,-8,-52,-49,29],[197,899,113,-192,-71,-43,-8,-44],[-192,113,899,196,61,49,8,52],[407,-192,196,611,8,44,59,-23],[-8,-71,61,8,411,-599,208,208],[-52,-43,49,44,-599,411,208,208],[-49,-8,8,59,208,208,99,-911],[29,-44,52,-23,208,208,-911,99]]', $__g_BLAS_STYPE_TRIANGLE + $__g_BLAS_STYPE_PACKED + $__g_BLAS_STYPE_UPPER)
;                  _blas_display($mTP, "packed (" & $mTP.size & " elements)")
;                  Global $mU = _blas_tpttr($mTP)
;                  _blas_display($mU, "unpacked (" & $mU.size & " elements)")
; ===============================================================================================================================
Func _blas_tpttr($mAP, $cUPLO = Default, $mA = Default, $iN = Default, $iLDA = Default, $bRetMap = True, $sDataType = "DOUBLE")
	Local $pAP, $pA ; pointer to the data in memory
	Local $bSpecialCase = False

	If IsKeyword($cUPLO) = 1 Then $cUPLO = BitAND($mAP.storageType, $__g_BLAS_STYPE_LOWER) ? "L" : "U"

	; Set parameters depending on the input type
	Select
		Case IsMap($mAP)
			$sDataType = $mAP.datatype
			If IsKeyword($iN)   = 1 Then $iN   = $mAP.rows
			If IsKeyword($iLDA) = 1 Then $iLDA = $iN
			$pAP = $mAP.ptr
		Case IsPtr($mAP)
			$pAP = $mAP
		Case IsDllStruct($mAP)
			$pAP = DllStructGetPtr($mAP)
	EndSelect
	Select
		Case IsKeyword($mA) = 1
			If $bRetMap Then
				$mA = _blas_createMatrix($iLDA, $iN, $sDataType, BitOR($__g_BLAS_STYPE_MATRIX, $__g_BLAS_STYPE_TRIANGLE, $cUPLO = "L" ? $__g_BLAS_STYPE_LOWER : $__g_BLAS_STYPE_UPPER))
				$pA = $mA.ptr
			Else
				Local $tStruct = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iLDA * $iN))
				$pA = DllStructGetPtr($tStruct)
				$bSpecialCase = True
			EndIf
		Case IsMap($mA)
			$pA = $mA.ptr
		Case IsPtr($mA)
			$pA = $mA
		Case IsDllStruct($mA)
			$pA = DllStructGetPtr($mA)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; set char buffers
	DllStructSetData($tBLASCHAR1, 1, $cUPLO)

	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "tpttr", _
		"PTR",  $pBLASCHAR1, _ ; UPLO
		"INT*", $iN, _         ; N
		"PTR",  $pAP, _        ; AP
		"PTR",  $pA, _         ; A
		"INT*", $iLDA, _       ; lda
		"INT*", 0 _            ; INFO
	)
	If @error        Then Return SetError(1, @error, $bRetMap ? Null : False)
	If $aDLL[6] <> 0 Then Return SetError(2, $aDLL[6], $bRetMap ? Null : False)

	If $bRetMap Then Return $bSpecialCase ? _blas_fromStruct($tStruct, $iLDA, $iN, $sDataType, BitOR($__g_BLAS_STYPE_MATRIX, $__g_BLAS_STYPE_TRIANGLE, $cUPLO = "L" ? $__g_BLAS_STYPE_LOWER : $__g_BLAS_STYPE_UPPER)) : $mA
	Return $bSpecialCase ? $tStruct : True
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_trttp()
; Description ...: copies a triangular matrix from the standard full format (TR) to the standard packed format (TP)
; Syntax ........: _blas_trttp($mMatrix, [$cUPLO = "U", [$mAP = Default, [$iN = Default, [$iLDA = Default, [$bRetMap = True, [$sDataType = "DOUBLE"]]]]]])
; Parameters ....: mMatrix   - [Map] unpacked upper or lower triangular matrix
;                  cUPLO     - [Char] (Default: "U")
;                            ↳ "U": matrix is an upper triangular matrix
;                              "L": matrix is an lower triangular matrix
;                  mAP       - [Map] (Default: Default)
;                            ↳ target map where the values should be written
;                              if Default: function returns a new matrix
;                  iN        - [Int] (Default: Default)
;                            ↳ order of the matrix A (number of rows)
;                  iLDA      - [Int] (Default: Default)
;                            ↳ first("leading") dimension size of A
;                  bRetMap   - [Bool] (Default: True)
;                            ↳ if true: function returns an packed upper or lower triangular matrix as a map
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: $bRetMap = True:  packed upper or lower triangular matrix as a map
;                           $bRetMap = False: True
;                  Failure: $bRetMap ? Null : False and set @error to:
;                           | 1: error during DllCall of trttp (@extended: @error from DllCall)
;                           | 2: error inside trttp (@extended: INFO-value from trttp)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/df/d41/group__trttp_ga5eb049cce8be54e0f7c75b7d998a8679.html#ga5eb049cce8be54e0f7c75b7d998a8679
; Example .......: Yes
;                  Global $mA = _blas_fromArray('[[611,197,-192,407,-8,-52,-49,29],[197,899,113,-192,-71,-43,-8,-44],[-192,113,899,196,61,49,8,52],[407,-192,196,611,8,44,59,-23],[-8,-71,61,8,411,-599,208,208],[-52,-43,49,44,-599,411,208,208],[-49,-8,8,59,208,208,99,-911],[29,-44,52,-23,208,208,-911,99]]')
;                  _blas_display($mA, "unpacked matrix")
;                  Global $mT = _blas_trttp($mA)
;                  _blas_display($mT, "packed triangular matrix")
; ===============================================================================================================================
Func _blas_trttp($mMatrix, $cUPLO = "U", $mAP = Default, $iN = Default, $iLDA = Default, $bRetMap = True, $sDataType = "DOUBLE")
	Local $pAP, $pA, _ ; pointer to the data in memory
	      $bSpecialCase = False

	 If IsKeyword($cUPLO) = 1 Then $cUPLO = BitAND($mMatrix.storageType, $__g_BLAS_STYPE_LOWER) ? "L" : "U"

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrix)
			$sDataType = $mMatrix.datatype
			If IsKeyword($iN)   = 1 Then $iN   = $mMatrix.cols
			If IsKeyword($iLDA) = 1 Then $iLDA = $mMatrix.rows
			$pA = $mMatrix.ptr
		Case IsPtr($mMatrix)
			$pA = $mMatrix
		Case IsDllStruct($mMatrix)
			$pA = DllStructGetPtr($mMatrix)
	EndSelect
	Local $iNNew = $iN * ($iN + 1) / 2
	Select
		Case IsKeyword($mAP) = 1
			If $bRetMap Then
				$mAP = _blas_createVector($iNNew, $sDataType)
				; adjust matrix object settings
				$mAP.rows        = $mMatrix.rows
				$mAP.cols        = $mMatrix.cols
				$mAP.storageType = BitOR($__g_BLAS_STYPE_MATRIX, $__g_BLAS_STYPE_PACKED, $__g_BLAS_STYPE_TRIANGLE, $cUPLO = "L" ? $__g_BLAS_STYPE_LOWER : $__g_BLAS_STYPE_UPPER)  
				$mAP.elements    = $iN
				$mAP.size        = $iNNew

				$pAP = $mAP.ptr
			Else
				Local $tStruct = DllStructCreate(StringFormat("%s[%d]", $sDataType, $iNNew))
				$pAP = DllStructGetPtr($tStruct)
				$bSpecialCase = True
			EndIf
		Case IsMap($mAP)
			$pAP = $mAP.ptr
		Case IsPtr($mAP)
			$pAP = $mAP
		Case IsDllStruct($mAP)
			$pAP = DllStructGetPtr($mAP)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; set char buffers and input healing
	DllStructSetData($tBLASCHAR1, 1, $cUPLO = "U" ? "U" : "L")

	Local $aDLL = DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "trttp", _
		"PTR",  $pBLASCHAR1, _ ; UPLO
		"INT*", $iN, _         ; N
		"PTR",  $pA, _         ; A
		"INT*", $iLDA, _       ; lda
		"PTR",  $pAP, _        ; AP
		"INT*", 0 _            ; INFO
	)
	If @error        Then Return SetError(1, @error, $bRetMap ? Null : False)
	If $aDLL[6] <> 0 Then Return SetError(2, $aDLL[6], $bRetMap ? Null : False)

	If $bRetMap Then Return SetExtended($iNNew, $bSpecialCase ? _blas_fromStruct($tStruct, $iLDA, $iN, $sDataType, BitOR($__g_BLAS_STYPE_MATRIX, $__g_BLAS_STYPE_PACKED, $__g_BLAS_STYPE_TRIANGLE, $cUPLO = "L" ? $__g_BLAS_STYPE_LOWER : $__g_BLAS_STYPE_UPPER)) : $mAP)
	Return SetExtended($iNNew, $bSpecialCase ? $tStruct : True)
EndFunc



; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_imatcopy()
; Description ...: performs scaling and in-place transposition/copying of matrices
; Syntax ........: _blas_imatcopy($mMatrix, [$fAlpha = 1, [$cTransposed = "M", [$iM = Default, [$iN = Default, [$iLDA = Default, [$iLDB = Default, [$sDataType = "DOUBLE"]]]]]]])
; Parameters ....: mMatrix     - [Map] matrix A as a map, DllStruct or pointer (gets overwritten)
;                  fAlpha      - [Float] (Default: 1)
;                              ↳ scalar value to scale the matrix with
;                  cTransposed - [Char] (Default: "T")
;                              ↳ if "T": the matrix is transposed during copying
;                                if "N": the matrix is copied as it is
;                  iM          - [Int] (Default: Default)
;                              ↳ number of rows of the matrix
;                  iN          - [Int] (Default: Default)
;                              ↳ number of columns of the matrix
;                  iLDA        - [Int] (Default: Default)
;                              ↳ first("leading") dimension size of the source matrix
;                  iLDB        - [Int] (Default: Default)
;                              ↳ first("leading") dimension size of the destination matrix
;                  sDataType   - [String] (Default: "DOUBLE")
;                              ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of imatcopy (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......: Not included in all BLAS implementations (at least included in Intel MKL and OpenBLAS)
; Related .......:
; Link ..........: https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-0/mkl-imatcopy.html
; Example .......: Yes
;                  Global $mMatrix = _blas_fromArray('[[1,2,3,4,5],[6,7,8,9,10],[11,12,13,14,15],[16,17,18,19,20]]')
;                  _blas_display($mMatrix, "before")
;                  _blas_imatcopy($mMatrix)
;                  Global $iTmp = $mMatrix.rows   ; swap .rows and .cols
;                  $mMatrix.rows = $mMatrix.cols
;                  $mMatrix.cols = $iTmp
;                  _blas_display($mMatrix, "transposed")
; ===============================================================================================================================
Func _blas_imatcopy($mMatrix, $fAlpha = 1, $cTransposed = "T", $iM = Default, $iN = Default, $iLDA = Default, $iLDB = Default, $sDataType = "DOUBLE")
	Local $pA ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrix)
			$sDataType = $mMatrix.datatype
			; calculate number of elements if not defined
			If IsKeyword($iM)   = 1 Then $iM   = $mMatrix.rows
			If IsKeyword($iN)   = 1 Then $iN   = $mMatrix.cols
			If IsKeyword($iLDA) = 1 Then $iLDA = $iM
			If IsKeyword($iLDB) = 1 Then $iLDB = $iN

			$pA = $mMatrix.ptr
		Case IsPtr($mMatrix)
			$pA = $mMatrix
		Case IsDllStruct($mMatrix)
			$pA = DllStructGetPtr($mMatrix)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; set char buffers
	DllStructSetData($tBLASCHAR1, 1,"C")
	DllStructSetData($tBLASCHAR2, 1, $cTransposed = "T" ? "T" : "N")

	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "imatcopy", _
		"PTR",            $pBLASCHAR1, _   ; order (here column major)
		"PTR",            $pBLASCHAR2, _   ; char *transa  ( If "N" then A is used; if "T" then Aᵀ is used)
		"INT*",           $iM, _           ; number of rows in matrix A
		"INT*",           $iN, _           ; number of columns in matrix A
		$sDataType & "*", $fAlpha, _       ; the scaling constant alpha for
		"ptr",            $pA, _           ; the matrix
		"INT*",           $iLDA, _         ; lda: length of a column in a (used to calc the element index - can be used to change indexing behavior)
		"INT*",           $iLDB _          ; ldb: length of a column in target matrix (used to calc the element index - can be used to change indexing behavior)
	)
	Return @error ? SetError(1, @error, False) : True
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _blas_omatcopy()
; Description ...: performs scaling and out-place transposition/copying of matrices.
; Syntax ........: _blas_omatcopy($mMatrix, [$fAlpha = 1, [$cTransposed = "N", [$mRet = Default, [$iM = Default, [$iN = Default, [$iLDA = Default, [$iLDB = Default, [$sDataType = "DOUBLE"]]]]]]]])
; Parameters ....: mMatrix     - [Map] matrix A as a map, DllStruct or pointer
;                  fAlpha      - [Float] (Default: 1)
;                              ↳ scalar value to scale the matrix with
;                  cTransposed - [Char] (Default: "N")
;                              ↳ if "T": the matrix is transposed during copying
;                                if "N": the matrix is copied as it is
;                  mRet        - [Map] (Default: Default)
;                              ↳ target map with appropriate dimensioning where the values should be written
;                                if Default: function returns a new matrix
;                  iM          - [Int] (Default: Default)
;                              ↳ number of rows of the matrix
;                  iN          - [Int] (Default: Default)
;                              ↳ number of columns of the matrix
;                  iLDA        - [Int] (Default: Default)
;                              ↳ first("leading") dimension size of the source matrix
;                  iLDB        - [Int] (Default: Default)
;                              ↳ first("leading") dimension size of the destination matrix
;                  sDataType   - [String] (Default: "DOUBLE")
;                              ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: $mRet = Default: matrix as a map
;                           Else: True
;                  Failure: $mRet = Default ? Null : False and set @error to:
;                           | 1: error during DllCall of omatcopy (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......: Not included in all BLAS implementations (at least included in Intel MKL and OpenBLAS)
; Link ..........: https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-0/mkl-omatcopy.html
; Example .......: Yes
;                  Global $mMatrix = _blas_fromArray('[[1,2,3,4,5],[6,7,8,9,10],[11,12,13,14,15],[16,17,18,19,20]]')
;                  _blas_display($mMatrix, "before")
;                  Global $mTransposed = _blas_omatcopy($mMatrix, 1, "T")
;                  _blas_display($mTransposed, "after (transposed)")
; ===============================================================================================================================
Func _blas_omatcopy($mMatrix, $fAlpha = 1, $cTransposed = "N", $mRet = Default, $iM = Default, $iN = Default, $iLDA = Default, $iLDB = Default, $sDataType = "DOUBLE")
	Local $pA, $pR ; pointer to the data in memory
	Local $bRetMap = False

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrix)
			$sDataType = $mMatrix.datatype
			; calculate number of elements if not defined
			If IsKeyword($iM)   = 1 Then $iM   = $mMatrix.rows
			If IsKeyword($iN)   = 1 Then $iN   = $mMatrix.cols
			If IsKeyword($iLDA) = 1 Then $iLDA = $iM
			If IsKeyword($iLDB) = 1 Then $iLDB = $iN

			$pA = $mMatrix.ptr
		Case IsPtr($mMatrix)
			$pA = $mMatrix
		Case IsDllStruct($mMatrix)
			$pA = DllStructGetPtr($mMatrix)
	EndSelect
	Select
		Case IsKeyWord($mRet) = 1
			$mRet = _blas_createMatrix($iN, $iM, $sDataType)
			$pR = $mRet.ptr
			$bRetMap = True
		Case IsMap($mRet)
			$pR = $mRet.ptr
		Case IsPtr($mRet)
			$pR = $mRet
		Case IsDllStruct($mRet)
			$pR = DllStructGetPtr($mRet)
	EndSelect

	Local Const $cPrefix = ($sDataType = "FLOAT") ? "s" : "d"

	; set char buffers
	DllStructSetData($tBLASCHAR1, 1,"C")
	DllStructSetData($tBLASCHAR2, 1, $cTransposed = "T" ? "T" : "N")

	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "omatcopy", _
		"PTR",            $pBLASCHAR1, _ ; order (here column major)
		"PTR",            $pBLASCHAR2, _ ; char *transa  ( If "N" then A is used; if "T" then Aᵀ is used)
		"INT*",           $iM, _         ; number of rows in matrix A
		"INT*",           $iN, _         ; number of columns in matrix A
		$sDataType & "*", $fAlpha, _     ; the scaling constant alpha for
		"ptr",            $pA, _         ; the matrix
		"INT*",           $iLDA, _       ; lda: length of a column in a (used to calc the element index - can be used to change indexing behavior)
		"ptr",            $pR, _         ; target matrix
		"INT*",           $iLDB _        ; ldb: length of a column in target matrix (used to calc the element index - can be used to change indexing behavior)
	)
	Return @error ? SetError(1, @error, $bRetMap ? Null : False) : $bRetMap ? $mRet : True
EndFunc

#EndRegion


#Region helper functions

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __blas_ArrayFromString()
; Description ...: creates an array from an array definition in AutoIt syntax, which is passed as a string
; Syntax ........: __blas_ArrayFromString($sArrayDef)
; Parameters ....: sArrayDef - [String] array definition in AutoIt syntax as a string
; Return value ..: Success: SetExtended(UBound($aRet), $aRet)
;                  Failure: Null and set @error to:
;                           | 1: Invalid form of the string (@extended: @error from StringRegExp())
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Link ..........: https://www.autoitscript.com/autoit3/docs/intro/lang_variables.htm#ArrayMaps
; Example .......: Yes
;                  Global $aMatrix = __blas_ArrayFromString("[[2,3,4],[2^2, 3^2, 3^4],[sqrt(2), sqrt(3), sqrt(4)]]")
;                  _ArrayDisplay($aMatrix)
; ===============================================================================================================================
Func __blas_ArrayFromString($sArrayDef)
	$sArrayDef = StringRegExpReplace($sArrayDef, '(^\s*\[|\]\s*$)', '')
	Local $aVals = StringRegExp($sArrayDef, '(?sx)(?(DEFINE)   (?<string> ''(?>[^'']+|'''')*'' | "(?>[^"]+|"")*")   (?<value> \s*(?>\g<string> | [^,\[\]]+)\s* )   (?<subarray> \s*\K\[ \g<value>(?:, \g<value>)* \])   (?<outervalue> \g<subarray> | \g<value> ))\g<outervalue>', 3)
	If @error Then Return SetError(1, @error, Null)

	Local $n2Dim = 0

	For $i = 0 To UBound($aVals) - 1
		If StringRegExp($aVals[$i], '^\s*\[') Then
			; count elements of second dimension:
			StringRegExpReplace($aVals[$i], '(?sx)(?(DEFINE)   (?<string> ''(?>[^'']+|'''')*'' | "(?>[^"]+|"")*")   (?<value> \s*(?>\g<string> | [^,\[\]]+)\s* ))\g<value>', '')
			If @extended > $n2Dim Then $n2Dim = @extended
		EndIf
	Next

	If $n2Dim > 0 Then ; 2D-array
		Local $aRet[UBound($aVals)][$n2Dim], $aSubVals

		For $i = 0 To UBound($aVals) - 1
			$aSubVals = StringRegExp($aVals[$i], '(?sx)(?(DEFINE)   (?<string> ''(?>[^'']+|'''')*'' | "(?>[^"]+|"")*")   (?<value> \s*(?>\g<string> | [^,\[\]]+)\s* ))\g<value>', 3)
			For $j = 0 To UBound($aSubVals) - 1
				$aRet[$i][$j] = Execute($aSubVals[$j])
			Next
		Next

		Return SetExtended(UBound($aRet), $aRet)
	Else ; 1D-Array
		Local $aRet[UBound($aVals)]

		For $i = 0 To UBound($aVals) - 1
			$aRet[$i] = Execute($aVals[$i])
		Next
		Return SetExtended(UBound($aRet), $aRet)
	EndIf
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __blas_fillWithScalar()
; Description ...: fills a matrix, a vector or parts thereof with a specific value
; Syntax ........: __blas_fillWithScalar($mMatrix, [$fScalar = 1.0, [$iStart = 0, [$iInc = 1, [$iN = Default, [$sDataType = "DOUBLE"]]]]])
; Parameters ....: mMatrix   - [Map] matrix A as a map, DllStruct or pointer (gets overwritten)
;                  fScalar   - [Float] (Default: 1.0)
;                            ↳ the scalar value with which the elements are to be overwritten
;                  iStart    - [Int] (Default: 0)
;                            ↳ start index
;                  iInc      - [Int] (Default: 1)
;                            ↳ storage spacing between elements
;                  iN        - [Int] (Default: Default)
;                            ↳ number of elements in input vector
;                  sDataType - [String] (Default: "DOUBLE")
;                            ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: True (@extended = number of elements overwritten)
;                  Failure: False and set @error to:
;                           | 1: error during DllCall of copy (@extended: @error from DllCall)
; Author ........: AspirinJunkie
; Modified.......: 2024-08-27
; Remarks .......:
; Related .......:
; Link ..........: https://www.netlib.org/lapack/explore-html/da/d55/group__gelss_gac6159de3953ae0386c2799294745ac90.html
; Example .......: Yes
;                  Global $mMatrix = _blas_createMatrix(10, 10)
;                  __blas_fillWithScalar($mMatrix, 7.3, 0, $mMatrix.rows + 1)
;                  _blas_display($mMatrix)
; ===============================================================================================================================
Func __blas_fillWithScalar(ByRef $mMatrix, $fScalar = 1.0, $iStart = 0, $iInc = 1, $iN = Default, $sDataType = "DOUBLE")
	Local $pM ; pointer to the data in memory

	; Set parameters depending on the input type
	Select
		Case IsMap($mMatrix)
			$sDataType = $mMatrix.datatype
			If IsKeyword($iN) = 1 Then $iN = Int(($mMatrix.size - $iStart - 1) / $iInc) + 1
			$pM = $mMatrix.ptr
		Case IsPtr($mMatrix)
			$pM = $mMatrix
		Case IsDllStruct($mMatrix)
			$pM = DllStructGetPtr($mMatrix)
	EndSelect

	Local Const $cPrefix = $sDataType = "FLOAT" ? "s" : "d", _
	            $dSize   = ($sDataType = "DOUBLE") ? $iBLAS_SIZE_DOUBLE : $iBLAS_SIZE_FLOAT

	; run blas copy to fill the target
	DllCall($__g_hBLAS_DLL, "NONE:cdecl", $cPrefix & "copy", _
		"INT*",           $iN, _                    ; N - number of elements to fill
		$sDataType & "*", $fScalar, _               ; DX - pointer to source value (here direct because we only have 1 element)
		"INT*",           0, _                      ; INCX - no increment in source (because we only have one value)
		"PTR",            $pM + $dSize * $iStart, _ ; DY - pointer to target matrix/vector
		"INT*",           $iInc _                  ; INCY - element step in the target matrix/vector
	)
	Return @error ? SetError(1, @error, False) : SetExtended($iN, True)
EndFunc

; #INTERNAL_USE_ONLY# ===========================================================================================================
; Name ..........: __blas_GBSfromArray()
; Description ...: converts an AutoIt array into a banded matrix, which is stored in "General-Band Storage Mode"
; Syntax ........: __blas_GBSfromArray($aArray, [$kl = 0, [$ku = 0, [$sType = "DOUBLE"]]])
; Parameters ....: aArray - [Array] Array holding the matrix
;                  kl     - (Default: 0)
;                         ↳ number of subdiagonals
;                  ku     - (Default: 0)
;                         ↳ number of super-diagonals
;                  sType  - [String] (Default: "DOUBLE")
;                         ↳ data type of the individual elements of the matrix/vector. Either "DOUBLE" or "FLOAT" possible.
; Return value ..: Success: [Map] BLAS Matrix Map as it is structured in this UDF
;                  Failure: Null and set @error to:
;                           | 1: matrix must be quadratic (@extended: number of cols of aArray)
;                           | 2: error during DllStructCreate (@extended: @error from DllStructCreate())
; Author ........: AspirinJunkie
; Modified.......: 2024-10-01
; Remarks .......: The GBS-Mode is NOT the "BLAS-General-Band Storage Mode" (see link for differences)
; Related .......:
; Link ..........: https://www.ibm.com/docs/en/essl/6.3?topic=representation-general-band-storage-mode
; Example .......: Yes
;                  Global $aArray[][] = [[19,-80,-55,0,0],[29,-92,-67,-94,0],[0,-29,77,-7,40],[0,0,-70,12,97],[0,0,0,53,43]]
;                  Global $mBanded = __blas_GBSfromArray($aArray, 1, 2)
;                  _blas_display($mBanded)
; ===============================================================================================================================
Func __blas_GBSfromArray($aArray, $kl = 0, $ku = 0, Const $sType = "DOUBLE")
	Local $iLDAB = 2 * $kl + $ku + 1
	Local $iN = UBound($aArray, 1)

	If UBound($aArray, 2) <> $iN Then Return SetError(1, UBound($aArray, 2), Null)

	Local $tStruct = DllStructCreate(StringFormat("%s[%d]", $sType, $iLDAB * $iN))
	If @error Then Return SetError(2, @error, Null)

	Local $iIndex
    For $iJ = 1 To $iN
        For $iI = (1 >= $iJ - $ku ? 1 : $iJ - $ku) To $iN <= $iJ + $kl ? $iN : $iJ + $kl
            $iIndex = ($iI - $iJ + $kl + $ku + 1) + ($iJ - 1) * $iLDAB
			DllStructSetData($tStruct, 1, $aArray[$iI - 1][$iJ - 1], $iIndex)
        Next
    Next

	Local $mRet[]
	$mRet.struct      = $tStruct
	$mRet.ptr         = DllStructGetPtr($tStruct)
	$mRet.elements    = $iLDAB * $iN
	$mRet.size        = $iLDAB * $iN
	$mRet.rows        = $iLDAB
	$mRet.cols        = $iN
	$mRet.kl          = $kl
	$mRet.ku          = $ku
	$mRet.storageType = $__g_BLAS_STYPE_MATRIX
	$mRet.datatype    = $sType

	Return $mRet
EndFunc

#EndRegion