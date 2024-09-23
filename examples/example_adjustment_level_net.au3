; load the BLAS/LAPACK-Dll (MUST done BEFORE #include "LinearAlgebra"!)
$__g_hBLAS_DLL = DllOpen(FileGetLongName("../libopenblas.dll"))

#include "../LinearAlgebra.au3"

; ------- Task -----------------
; A net was measured between a total of 6 points during leveling.
; The height of point no. 6 is given as 67.228 m.
; The absolute heights of the other points are now to be determined from the height differences.
; It is assumed that the greater the distance of the leveling, the lower the accuracy.
; The weight of the measurement is therefore assumed to be inversely proportional to the square root of the distance.
; In addition, we would like to introduce the concept of pseudo-observations to formulate conditions.
;
;              24     45
;           2――――――4――――――5
;          ╱ ╲  II │ III ╱ ╲
;         ╱   ╲    │    ╱   ╲
;      12╱   23╲   │34 ╱     ╲56
;       ╱       ╲  │  ╱35     ╲
;      ╱    I    ╲ │ ╱    IV   ╲
;     ╱           ╲│╱           ╲
;    1―――――――――――――3―――――――――――――6
;           13           36

Global $iFlags = 2 * $__LA_LSTSQ_RANK - 1 ; = all

Global $aMeasurements[9][4] = [ _
 _; from to   Δh[m]  D[km]
	[1,  2,  -8.206, 0.62], _
	[1,  3,  -5.734, 1.2], _
	[2,  3,   2.481, 0.45], _
	[2,  4,  -4.433, 0.8], _
	[3,  4,  -6.909, 1.0], _
	[3,  5, -18.872, 1.1], _
	[3,  6,   4.035, 0.44], _
	[4,  5, -11.962, 0.72], _
	[5,  6,  22.904, 0.83] _
]

; add the observations (functional relationship, observation value and a priori standard deviation of the observation)
Global $sFunc, $mObs
For $i = 0 To 8
	$sFunc = StringFormat("H%d - H%d", $aMeasurements[$i][1], $aMeasurements[$i][0])
	_la_adjustment_addObservation($mObs, $sFunc, $aMeasurements[$i][2], sqrt($aMeasurements[$i][3]))
Next

; The height of point 6 is known and is considered fixed. Therefore, a pseudo-observation with this height and a very small standard deviation is added for H6.
; With this the absolute height of the other points can be directly calculated.
_la_adjustment_addObservation($mObs, "H6", 67.228, 1e-5)

; the network structure results in additional conditions (loops must sum up to 0), which further increase the quality of the solution.
; These are formulated and added as pseudo-observations. (for seeing effect try once with the pseudo-observations and once without)
_la_adjustment_addObservation($mObs, "(H2 - H1) + (H3 - H2) - (H3 - H1)", 0, 1e-5)                                     ; I
_la_adjustment_addObservation($mObs, "(H4 - H2) - (H4 - H3) - (H3 - H2)", 0, 1e-5)                                     ; II
_la_adjustment_addObservation($mObs, "(H4 - H3) + (H5 - H4) - (H5 - H3)", 0, 1e-5)                                     ; III
_la_adjustment_addObservation($mObs, "(H5 - H3) + (H6 - H5) - (H6 - H3)", 0, 1e-5)                                     ; IV
_la_adjustment_addObservation($mObs, "(H3 - H1) + (H4 - H3) - (H2 - H1) - (H4 - H2)", 0, 1e-5)                         ; I, II
_la_adjustment_addObservation($mObs, "(H3 - H1) + (H5 - H3) - (H5 - H4) - (H4 - H2) - (H2 - H1)", 0, 1e-5)             ; I, II, III
_la_adjustment_addObservation($mObs, "(H3 - H1) + (H6 - H3) - (H6 - H5) - (H5 - H4) - (H4 - H2) - (H2 - H1)", 0, 1e-5) ; I,II,III,IV
_la_adjustment_addObservation($mObs, "(H5 - H3) - (H5 - H4) - (H4 - H2) + (H3 - H2)", 0, 1e-5)                         ; II, III
_la_adjustment_addObservation($mObs, "(H6 - H3) - (H6 - H5) - (H5 - H4) - (H4 - H2) + (H3 - H2)", 0, 1e-5)             ; II, III, IV
_la_adjustment_addObservation($mObs, "(H6 - H3) - (H6 - H5) - (H5 - H4) - (H4 - H3)", 0, 1e-5)                         ; III, IV

; do the adjustment and show the results (since there is a linear functional relationship, we do not need approximate values for the parameters)
Global $mLstSq = _la_adjustment($mObs, Default, "GN", "QR", $iFlags)
If @error Then Exit MsgBox(16, "error", "error during _la_adjustment()" & @CRLF & "@error: " & @error & @CRLF & "@extended: " & @extended)

ConsoleWrite("s0: " & $mLstSq.s0 & @CRLF)
ConsoleWrite("r^TPr: " & $mLstSq.r2sum & @CRLF)
_la_adj_showResult($mLstSq, "x")
;~ _la_adj_showResult($mLstSq, "Qx")
;~ _la_adj_showResult($mLstSq, "Qy")
;~ _la_adj_showResult($mLstSq, "Qr")
;~ _la_adj_showResult($mLstSq, "r")
;~ _la_adj_showResult($mLstSq, "sdY")
;~ _la_adj_showResult($mLstSq, "sdX")
;~ _la_adj_showResult($mLstSq, "sdR")