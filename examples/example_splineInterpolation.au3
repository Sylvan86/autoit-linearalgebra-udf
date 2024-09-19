; load the BLAS/LAPACK-Dll (MUST done BEFORE #include "LinearAlgebra"!)
Global $__g_hBLAS_DLL = DllOpen(FileGetLongName("../libopenblas.dll"))

#include "../LinearAlgebra.au3"

; ------- Task -----------------
; A cubic spline is to be calculated for a list of X-Y value pairs in order to interpolate between the interpolation points.

; the input values
Global $aPoints[][2] = [ _
	[1,	2.3], _
	[2,	1.3], _
	[3,	0.7], _
	[4,	1.8], _
 	[5,	1.3]  _
]

; calculate the polynomials of the cubic spline
Global $aPolynoms = _spline_cubic($aPoints)

; print polynoms:
For $i = 0 To UBound($aPolynoms) - 1
	ConsoleWrite(StringFormat("polynom %d [%g,%g]: %s\n", $i+1, $aPolynoms[$i].from, $aPolynoms[$i].to, $aPolynoms[$i].Func))
Next
ConsoleWrite(@CRLF)

; interpolate between the grid points
Global $fX = 2.5, $fY = _calcSplinePoint($aPolynoms, $fX)
ConsoleWrite("f(" & $fx & "): " & $fY & @CRLF & @CRLF)





; function to interpolate in a spline
Func _calcSplinePoint($aPolynoms, $fX)
	Local $aPoly
	For $aPoly In $aPolynoms
		If $aPoly.from > $fX or $aPoly.to < $fX Then ContinueLoop
		Return $aPoly.d + $aPoly.c * ($fx - $aPoly.from) + $aPoly.b * ($fx - $aPoly.from)^2 + $aPoly.a * ($fx - $aPoly.from)^3
	Next
	Return SetError(1, 0, Null)
EndFunc

; #FUNCTION# ====================================================================================================================
; Name ..........: _spline_cubic()
; Description ...: calculates the polynomials of a cubic spline for a list of X-Y pairs
; Syntax ........: _spline_cubic($aPoints, [$fDisplayDecimals = 2, [$fM0 = 0.0, [$fMn = 0.0]]])
; Parameters ....: aPoints          - [2D-Array] X-Y value pairs of the sampling points
;                  fDisplayDecimals - [Float] (Default: 2)
;                                   ↳ number of relevant decimal places in the output of the function string
;                  fM0              - [Float] (Default: 0.0)
;                                   ↳ given 2nd derivation at the first grid point
;                  fMn              - [Float] (Default: 0.0)
;                                   ↳ given 2nd derivation at the last grid point
; Return value ..: Success: [Array of Maps]: [
;                              {"a": param a, "b": param b, "c": param c, "d": param d, "Func": polynom as a string, "from": Lower definition range limit, "to": upper definition range limit}, ...
;                           ]
;                  Failure: Null and set @error to:
;                           |1X: error X during _lp_gtsv() (@extended: @extended from _lp_gtsv())
; Author ........: AspirinJunkie
; Modified.......: 2024-09-19
; Remarks .......:
; Related .......:
; Link ..........:
; Example .......: No
; ===============================================================================================================================
Func _spline_cubic($aPoints, $fDisplayDecimals = 2, $fM0 = 0.0, $fMn = 0.0)
	Local $iN = UBound($aPoints)

	; determine the increments between the interpolation points
	Local $aSteps[$iN - 1]
	For $i = 0 To $iN - 2
		$aSteps[$i] = $aPoints[$i + 1][0] - $aPoints[$i][0]
	Next

	; create system of equations
	; The function matrix A is a tridiagonal matrix. For reasons of efficiency, we therefore want to solve the system of equations using _lp_gtsv().
	; To do this, we consider the 3 diagonals of the matrix A as separate vectors.
	Local $tD  = DllStructCreate("DOUBLE[" & ($iN - 2) & "]")
	Local $tDU = DllStructCreate("DOUBLE[" & ($iN - 3) & "]")
	Local $tDL = DllStructCreate("DOUBLE[" & ($iN - 3) & "]")
	Local $tB  = DllStructCreate("DOUBLE[" & ($iN - 2) & "]")
	For $i = 0 To $iN - 3
		; main diagonal of A
		DllStructSetData($tD, 1,  2 * ($aSteps[$i] + $aSteps[$i + 1]), $i + 1)
		; supper diagonal of A
		If $i > 0 Then DllStructSetData($tDL, 1, $aSteps[$i], $i)
		; super diagonal of A
		If $i < $iN - 3 Then DllStructSetData($tDU, 1, $aSteps[$i + 1], $i + 1)
		; fill the observation vector b
		DllStructSetData($tB, 1, 6 * ( (($aPoints[$i + 2][1] - $aPoints[$i + 1][1]) / $aSteps[$i + 1]) - (($aPoints[$i + 1][1] - $aPoints[$i][1]) / $aSteps[$i]) ), $i + 1)
	Next

	; calculate the 2nd derivatives at the interpolation points and stores them in $aM
	_lp_gtsv($tD, $tB, $tDU, $tDL, 1, $iN - 2)
	If @error Then Return SetError(10 + @error, @extended, Null)
	Local $aM[$iN] = [$fM0]
	$aM[$iN - 1] = $fMn
	For $i = 1 To $iN - 2
		; the 2nd derivates are stored in the result vector
		$aM[$i] = DllStructGetData($tB, 1, $i)
	Next

	; calculate the polynom parameters and prepare for return
	Local $aReturn[$iN - 1], $sFunc
	For $i = 1 To $iN - 1
		Local $mPolynom[]

		; the limits of the scope of the polynomial
		$mPolynom.from = $aPoints[$i - 1][0]
		$mPolynom.to = $aPoints[$i][0]

		; the parameters of the cubic polynomial
		$mPolynom.a = ($aM[$i] - $aM[$i - 1]) / 6 * $aSteps[$i - 1]
		$mPolynom.b = $aM[$i - 1] / 2
		$mPolynom.c = ($aPoints[$i][1] - $aPoints[$i - 1][1]) / $aSteps[$i - 1] - $aSteps[$i - 1] / 6   * (2 * $aM[$i - 1] + $aM[$i])
		$mPolynom.d = $aPoints[$i - 1][1]

		$mPolynom.Func = ($mPolynom.d = 0 ? "" : StringFormat("%." & $fDisplayDecimals & "g", $mPolynom.d)) & _
		         ($mPolynom.c = 0 ? "" : ($mPolynom.c >= 0 ? " + " : " - ") & StringFormat("%." & $fDisplayDecimals & "g * (x ", Abs($mPolynom.c)) & ($aPoints[$i - 1][0] >= 0 ? "- " : "+ ") & StringFormat("%." & $fDisplayDecimals & "g)", Abs($aPoints[$i - 1][0]))) & _
				 ($mPolynom.b = 0 ? "" : ($mPolynom.b >= 0 ? " + " : " - ") & StringFormat("%." & $fDisplayDecimals & "g * (x ", Abs($mPolynom.b)) & ($aPoints[$i - 1][0] >= 0 ? "- " : "+ ") & StringFormat("%." & $fDisplayDecimals & "g)^2", Abs($aPoints[$i - 1][0]))) & _
				 ($mPolynom.a = 0 ? "" : ($mPolynom.a >= 0 ? " + " : " - ") & StringFormat("%." & $fDisplayDecimals & "g * (x ", Abs($mPolynom.a)) & ($aPoints[$i - 1][0] >= 0 ? "- " : "+ ") & StringFormat("%." & $fDisplayDecimals & "g)^3", Abs($aPoints[$i - 1][0])))

		$aReturn[$i - 1] = $mPolynom
	Next

	Return $aReturn
EndFunc

