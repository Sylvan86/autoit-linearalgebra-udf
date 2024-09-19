; load the BLAS/LAPACK-Dll (MUST done BEFORE #include "LinearAlgebra"!)
$__g_hBLAS_DLL = DllOpen(FileGetLongName("../libopenblas.dll"))

#include "../LinearAlgebra.au3"

Global $iFlags = 2 * $__LA_LSTSQ_RANK - 1 ; = all


Global $sFunc = "a * $x + b * $y + c"
Global $aX[] = [ 0,       0,     0,      0,       1,       1,      1,      1,       2,      2,      2,      2,      3,      3,      3,      3]
Global $aY[] = [ 0,       1,     2,      3,       0,       1,      2,      3,       0,      1,      2,      3,      0,      1,      2,      3]
Global $aZ[] = [-1.0015, -0.67, -0.3309, 0.0001, -0.5023, -0.1669, 0.1606, 0.5033, -0.0029, 0.3313, 0.6707, 0.9973, 0.4988, 0.8318, 1.1715, 1.5002]
Global $mVars[]
$mVars.X = $aX
$mVars.Y = $aY


Global $mRegression = _la_regression($sFunc, $mVars, $aZ, $iFlags)

ConsoleWrite(StringFormat("\n% 11s: %4.6f\n% 11s: %4.6f\n% 11s: %4.6f\n\n", "R²", $mRegression.R2, "R² adjusted", $mRegression.R2_corr, "Pearson", $mRegression.pearson))
ConsoleWrite("result function: y = " & $mRegression.Func & @CRLF & @CRLF)
_la_adj_showResult($mRegression, "x")
_la_adj_showResult($mRegression, "r")
_la_display($mRegression.Qx, "Qx")
