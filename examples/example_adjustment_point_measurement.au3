; load the BLAS/LAPACK-Dll (MUST done BEFORE #include "LinearAlgebra"!)
$__g_hBLAS_DLL = DllOpen(FileGetLongName("../libopenblas.dll"))

#include "../LinearAlgebra.au3"

; ------- Task -----------------
; Distances to a new point were measured from 2 known points and the angles between this point and the other point.
; The coordinates of the new point are to be determined from this.
;
;              N
;             ╱ ╲
;            ╱   ╲
;           ╱     ╲
;	     D1╱       ╲D2
;         ╱         ╲
;        ╱           ╲
;       ╱╮α1       α2╭╲
;      P1 ――――――――――― P2


Global Const $fPI = ACos(-1)
Global $iFlags = 2 * $__LA_LSTSQ_RANK - 1 ; = all

; coordinates of the fixed points
Global $aFixedPoints[4][2] = [ _  ; [m]
	[414.06, 589.27], _ ; P1
	[495.07, 821.01], _ ; P2
	[263.77, 824.34], _ ; P3
	[227.28, 601.58] _  ; P4
]
; distance measurements    -->P1    -->P2     -->P3     -->P4
Global $aDistances[4] =  [173.518, 210.891,  107.586,  147.535]  ; [m]
; direction ("angles")     -->P1    -->P2     -->P3     -->P4
Global $aDirections[4] = [65.8211, 153.2858, 249.7261, 387.1499] ; [gon]

; add the observations (distances and directions)
Global $sFunc, $mObs ; the variable holding the observations
For $i = 0 To 3
	; add distance observation:  sqrt((y_n - y_i)^2 + (x_n - x_i)^2)
	$sFunc = StringFormat("sqrt((x_n - %.16g)^2 + (y_n - %.16g)^2)", $aFixedPoints[$i][1], $aFixedPoints[$i][0]) ; functional relationship between observation an target parameters
	_la_adj_addObservation($mObs, $sFunc, $aDistances[$i], 0.01, "sD") ; "sD" means: observation belongs to the observation group “sD”. This is used to examine the accuracy of the groups independently (variance component estimation)
	; 0.01 means: expected  standard deviation of the observation - estimate how accurate the observation is. Specified in [m]

	; add direction observation: 200 / Pi * arctan( (x_n - x_i) / (y_n - y_i)) - r0
	; r0 = initial direction on the viewpoint
	$sFunc = StringFormat('_calcDirection_gon( x_n, y_n, %.16g, %.16g, r0)', $aFixedPoints[$i][1], $aFixedPoints[$i][0] ) ; functional relationship between observation an target parameters
	_la_adj_addObservation($mObs, $sFunc, $aDirections[$i], 0.001, "sR") ; "sR" means: observation belongs to the observation group “sD”. This is used to examine the accuracy of the groups independently (variance component estimation)
	; 0.001 means: expected  standard deviation of the observation - estimate how accurate the observation is. Specified in [gon]
Next

; Approximate values for the new point coordinates and the direction (e.g. can be derived from the triangle P4-P3-N)
Global $mParams = __la_adj_getParamList($mObs)
__la_adj_setApproxValue($mParams, "x_n", 725.8)
__la_adj_setApproxValue($mParams, "y_n", 306.9)
__la_adj_setApproxValue($mParams, "r0", 277)

; do the adjustment and show the results
Global $mLstSq = _la_adjustment($mObs, $mParams, "GN", "QR", $iFlags)
If @error Then Exit MsgBox(16, "error", "error during _la_adjustment()" & @CRLF & "@error: " & @error & @CRLF & "@extended: " & @extended)
_la_adj_showResult($mLstSq, "x")
_la_adj_showResult($mLstSq, "Qx")
_la_adj_showResult($mLstSq, "r")




#Region helper functions
; calc direction ("angle") in [gon] from one point to another by knowing their coordinates and the initial direction
Func _calcDirection_gon($fXfrom, $fYfrom, $fxTo, $fYTo, $fR0 = 0)
	Local Static $fRho = 200 / $fPI
	Return Mod($fRho * atan2($fxTo - $fXfrom, $fYTo - $fYfrom) - $fR0 + 800, 400)
EndFunc

Func atan2($y, $x)
	Return (2 * ATan($y / ($x + Sqrt($x * $x + $y * $y))))
EndFunc
#EndRegion


