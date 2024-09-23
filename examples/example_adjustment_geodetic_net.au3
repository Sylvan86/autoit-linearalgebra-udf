; load the BLAS/LAPACK-Dll (MUST done BEFORE #include "LinearAlgebra"!)
$__g_hBLAS_DLL = DllOpen(FileGetLongName("../libopenblas.dll"))

#include "../LinearAlgebra.au3"

; ------- Task -----------------
; The distance and direction to a total of 4 points with known coordinates were measured on a new point.
; The coordinates of the position of this new point and the orientation of the instrument are to be calculated from this.
; We estimate that the distances were measured to approximately ± 0.01 m and the directions to ± 0.001 gon.
; However, as we are not entirely sure about this, the values are to be derived using a variance component estimate.
;
;                                   P1
;                                    ╲
;                                     ╲d1
;                                      ╲
;                                       ╲ ╮    d2
;                                        PN ―――――――――P2
;                                      ╱ │
;                                   d4╱  │
;                                    ╱   │d3
;                                   P4   │
;                                        │
;                                        P3

Global Const $fPI = ACos(-1)
Global $iFlags = 2 * $__LA_LSTSQ_RANK - 1 ; = all

; coordinates of the fixed points
Global $aFixedPoints[4][2] = [ _  ; [m]
	[314.06, 489.27], _ ; P1
	[395.07, 721.01], _ ; P2
	[163.77, 724.34], _ ; P3
	[127.28, 501.58] _  ; P4
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
; Important: Lower case letters must be converted to upper case letters for the parameter names
Global $mParams = __la_adj_getParamList($mObs)
__la_adj_setApproxValue($mParams, "x_n", 625.8)
__la_adj_setApproxValue($mParams, "y_n", 206.9)
__la_adj_setApproxValue($mParams, "r0", 277)

; do the adjustment and show the results
Global $mLstSq = _la_adjustment($mObs, $mParams, "GN", "QR", $iFlags, 1e-9, "Central", 1, 50, "DOUBLE")
If @error Then Exit MsgBox(16, "error", "error during _la_adjustment()" & @CRLF & "@error: " & @error & @CRLF & "@extended: " & @extended)
_la_adj_showResult($mLstSq, "x")



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


