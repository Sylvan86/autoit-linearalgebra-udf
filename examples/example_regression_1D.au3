; load the BLAS/LAPACK-Dll (MUST done BEFORE #include "LinearAlgebra"!)
$__g_hBLAS_DLL = DllOpen(FileGetLongName("../libopenblas.dll"))

#include "../LinearAlgebra.au3"

Global $iFlags = 2 * $__LA_LSTSQ_RANK - 1 ; = all

; regression function (a, b = parameters to be estimated; $x = placeholder for the x-values)
Global $sFunc = "a * $x + b"
Global $aX[] = [ 0,    2,    4,    6,    8,    10,    12,    15,    20,    25]   ; x-values
Global $aY[] = [17.2, 34.1, 50.9, 68.5, 86.2, 103.1, 121.0, 143.8, 185.3, 228.7] ; y-values (1D-Array means 1D-regression)

; calculate the regression
Global $mRegression = _la_regression($sFunc, $aX, $aY, $iFlags)
ConsoleWrite(StringFormat("\n% 11s: %4.5f\n% 11s: %4.5f\n% 11s: %4.5f\n\n", "R²", $mRegression.R2, "R² adjusted", $mRegression.R2_corr, "Pearson", $mRegression.pearson))
ConsoleWrite("result function: y = " & $mRegression.Func & @CRLF & @CRLF)
_la_adj_showResult($mRegression, "x")
_la_adj_showResult($mRegression, "r")
