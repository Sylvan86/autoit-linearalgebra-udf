; load the BLAS/LAPACK-Dll (MUST done BEFORE #include "LinearAlgebra"!)
$__g_hBLAS_DLL = DllOpen(FileGetLongName("../libopenblas.dll"))

#include "../LinearAlgebra.au3"

; ------- Task -----------------
; The coordinates of 4 points were determined, all of which we know lie on a circle.
; The coordinates of these points should now be balanced so that this condition is met.
; The expected standard deviation of the coordinate components is 0.3 m.

Global $iFlags = 2 * $__LA_LSTSQ_RANK - 1 ; = all

Global $aPoints[4][2] = [ _
	[ 46.0, 78.3], _    ; P1
	[103.2, 66.5], _    ; P2
	[112.7, 22.9], _    ; P3
	[ 20.5, 23.5] _     ; P4
]

; add the observations (functional relationship, observation value and a priori standard deviation of the observation)
Global $sFunc, $mObs
For $i = 0 To UBound($aPoints) - 1
	; add the point coordinates as observations for the target parameters
	_la_adj_addObservation($mObs, StringFormat("X_%i", $i + 1), $aPoints[$i][0], 0.3)
	_la_adj_addObservation($mObs, StringFormat("Y_%i", $i + 1), $aPoints[$i][1], 0.3)

	; formulate the condition that the points lie on a circle as an equation.
	$sFunc = StringFormat("sqrt((X_m - X_%i)^2 + (Y_m - Y_%i)^2) - R", $i + 1, $i + 1)
	_la_adj_addObservation($mObs, $sFunc, 0, 0.01)
Next

; create parameter list and add approximate values for the parameters.
Global $mParams = __la_adj_getParamList($mObs)
__la_adj_setApproxValue($mParams, "X_M", 30)
__la_adj_setApproxValue($mParams, "Y_M", 20)
__la_adj_setApproxValue($mParams, "R", 1)
For $i = 0 To 3
	__la_adj_setApproxValue($mParams, StringFormat("X_%i", $i + 1), $aPoints[$i][0])
	__la_adj_setApproxValue($mParams, StringFormat("Y_%i", $i + 1), $aPoints[$i][1])
Next

; do the adjustment and show the results (since there is a linear functional relationship, we do not need approximate values for the parameters)
Global $mLstSq = _la_adjustment($mObs, $mParams, "GM", "QR", $iFlags)
If @error Then Exit MsgBox(16, "error", "error during _la_adjustment()" & @CRLF & "@error: " & @error & @CRLF & "@extended: " & @extended)

ConsoleWrite("s0: " & $mLstSq.s0 & @CRLF)
ConsoleWrite("r^TPr: " & $mLstSq.r2sum & @CRLF)
_la_adj_showResult($mLstSq, "x")
;~ _la_adj_showResult($mLstSq, "Qx")
;~ _la_adj_showResult($mLstSq, "r")
;~ _la_adj_showResult($mLstSq, "sdY")