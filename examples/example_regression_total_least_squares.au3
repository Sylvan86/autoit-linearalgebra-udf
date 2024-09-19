; load the BLAS/LAPACK-Dll (MUST done BEFORE #include "LinearAlgebra"!)
Global $__g_hBLAS_DLL = DllOpen(FileGetLongName("../libopenblas.dll"))

#include "../LinearAlgebra.au3"

; ------- Task -----------------
; In a normal regression, you have x-y pairs of measured values and the linear equation to be determined y = a * x + b
; If the parameters are determined from this using a normal regression, it is implicitly assumed that only the y-values are subject to error.
; The x-values are simply entered into the function equation as given.
; They then initially act as constants, although they are probably also subject to errors.
; The alternative to this is "orthogonal regression" - or for more general cases "Total Least Squares" (TLS)
; This always occurs when observation values appear in the observation equations for other observations.
; The problem: The calculation requires specially adapted calculation methods, such as the Gauss-Helmert model of adjustment or special TLS models.
; However, these are difficult to cast into generally usable functions, which is actually the aim of this UDF.
; Instead, these problems can also be calculated using the Gauss-Markov model used here by extending the observation equation and adding pseudo-observations:
;
; Let us assume the observation equation y = a * x + b.
; In the Gauss-Markov model, this implicitly becomes y + vy = a * x + b
; However, since we now say that x is also subject to errors, we extend the equation to include residuals for x:
; y + vy = a * (x + vx) + b
; We simply convert this to 0:
; 0 = a * (x + vx) + b - (y + vy)
; The number of parameters is therefore increased by 2n and the Jacobian matrix is therefore significantly larger.
; The observation vector here consists of nothing but zeros, while the x and y values are entered directly into the Jacobian matrix.
; To make the system definite again, we add pseudo-observations for each vx_i and vy_i with the estimated value 0
; and the highest possible standard deviation so that these observations do not influence the solution in any relevant way.
;
; In principle, this basic concept can be applied to all TLS problems.
; However, the disadvantage is that the Jacobian matrix quickly becomes very large as the number of parameters increases significantly.
; As long as the computing power is available, this approach is a good alternative to explicit TLS approaches.

Global $iFlags = $__LA_LSTSQ_SDX

Global $aPoints[][2] = [ _
	[5,		80], _
	[23,	78], _
	[25,	60], _
	[48,	53], _
	[17,	85], _
	[8,		84], _
	[4,		73], _
	[26,	79], _
	[11,	81], _
	[19,	75], _
	[14,	68], _
	[35,	72], _
	[29,	58], _
	[4,		92], _
	[23,	65] _
]

; add the observations (direct observations + pseudo observations)
Global $sFunc, $mObs ; the variable holding the observations
For $i = 0 To UBound($aPoints) - 1
	; add the observation (but with the observed value 0, as we have converted the equation to 0 and included the y-value in the equation itself)
	$sFunc = StringFormat("a * (%.16g + vx_%d) + b - %.16g - vy_%d ", $aPoints[$i][0], $i, $aPoints[$i][1], $i)
	_la_adjustment_addObservation($mObs, $sFunc, 0, 1)

	; add pseudo observations for vy_i & vx_i = 0 (with a high standard deviation so that it does not have a relevant influence on the result)
	_la_adjustment_addObservation($mObs, StringFormat("vy_%d", $i), 0, 100)
	_la_adjustment_addObservation($mObs, StringFormat("vx_%d", $i), 0, 100)
Next

; do the regression and show the results
Global $mLstSq = _la_adjustment($mObs, Default, "GN", "Cholesky", $iFlags, 1e-4)
If @error Then Exit MsgBox(16, "error", "error during _la_adjustment()" & @CRLF & "@error: " & @error & @CRLF & "@extended: " & @extended)
_la_adj_showResult($mLstSq, "x", "orthogonal regression")

; display the resulting function
Global $aX = _la_toArray($mLstSq.x)
Global $fA = $aX[_ArraySearch($mLstSq.params, "A")]
Global $fB = $aX[_ArraySearch($mLstSq.params, "B")]
ConsoleWrite(StringFormat("% 22s: y = %.2g * x + %.2g\n", "orthogonal regression", $fA, $fB))



; for comparison: do a normal regression with the given values
Global $aX[UBound($aPoints)], $aY[UBound($aPoints)]
For $i = 0 To UBound($aPoints) - 1
	$aX[$i] = $aPoints[$i][0]
	$aY[$i] = $aPoints[$i][1]
Next
; calculate the regression
Global $mRegression = _la_regression("a * $x + b", $aX, $aY, $iFlags)
ConsoleWrite(StringFormat("% 22s: y = %s\n\n", "normal regression", $mRegression.Func))
_la_adj_showResult($mRegression, "x", "normal regression")
