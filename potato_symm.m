(*Define constants for the implicit curve*)a = 1025/1000;
b = 9/100;
a2 = a^2;

(*Implicit curve equation F(x,y)=0-as numerical functions*)

F[x_?NumericQ, y_?NumericQ] := y^2 - 2*b*x*y - a2*x + a2*x^3;

(*Partial derivatives-computed numerically*)

Fx[x_?NumericQ, y_?NumericQ] := -2*b*y - a2 + 3*a2*x^2;
Fy[x_?NumericQ, y_?NumericQ] := 2*y - 2*b*x;
Fxx[x_?NumericQ, y_?NumericQ] := 6*a2*x;
Fyy[x_?NumericQ, y_?NumericQ] := 2.0;
Fxy[x_?NumericQ, y_?NumericQ] := -2.0*b;

(*Solve for y given x*)

solveYForX[xval_?NumericQ] := 
  Module[{discriminant, sqrtD, y1, y2}, 
   discriminant = (2*b*xval)^2 - 4*1*(-a2*xval + a2*xval^3);
   If[discriminant < 1*^-9, discriminant = 0];
   sqrtD = Sqrt[discriminant];
   y1 = (2*b*xval + sqrtD)/2.0;
   y2 = (2*b*xval - sqrtD)/2.0;
   {y1, y2}];

(*Adaptive sampling function*)

adaptiveXSampling[xMax_, baseResolution_ : 800] := 
  Module[{region1, region2, region3, region4, region5, nDense, 
    nSparse}, 
   nDense = 
    Floor[baseResolution*1.5];(*Dense sampling for special regions*)
   nSparse = 
    Floor[baseResolution];(*Sparse sampling for other \
regions*)(*Dense sampling near 0:[-0.02,0.02]-
   adjust start if negative*)
   region1 = Range[Max[0, -0.02], 0.02, 0.04/nDense];
   (*Sparse sampling:(0.02,0.98)*)
   region2 = Range[0.02 + 0.96/nSparse, 0.98, 0.96/nSparse];
   (*Dense sampling near 1:[0.98,1.02]*)
   region3 = 
    Range[0.98, Min[1.02, xMax], Min[0.04, 1.02 - 0.98]/nDense];
   (*Sparse sampling from 1.02 to xMax if xMax>1.02*)
   If[xMax > 1.02, 
    region4 = 
     Range[1.02 + (xMax - 1.02)/nSparse, 
      xMax, (xMax - 1.02)/nSparse];
    region5 = {}, region4 = {};
    region5 = 
     If[xMax > 1.02, {}, 
      Range[Max[region3] + (xMax - Max[region3])/nSparse, 
       xMax, (xMax - Max[region3])/nSparse]]];
   (*Combine and remove duplicates*)
   Union[Join[region1, region2, region3, region4, region5]]];

(*Trace the curve loop with adaptive sampling*)

traceCurveLoop[resolution_ : 800] := 
  Module[{roots, xMax, xValsForward, upperBranch, lowerBranch, y1, 
    y2},(*Find the positive root of-a^2*x^2+b^2*x+a^2=0*)
   roots = x /. Solve[-a2*x^2 + b^2*x + a2 == 0, x, Reals];
   xMax = Max[roots];
   (*Use adaptive sampling*)
   xValsForward = adaptiveXSampling[xMax, resolution];
   (*Trace upper branch*)
   upperBranch = Table[{y1, y2} = solveYForX[xval];
     {xval, y1}, {xval, xValsForward}];
   (*Trace lower branch in reverse*)
   lowerBranch = Table[{y1, y2} = solveYForX[xval];
     {xval, y2}, {xval, Reverse[xValsForward]}];
   (*Combine,removing duplicate endpoint*)
   Join[upperBranch, Rest[lowerBranch]]];

(*Compute tangent and normal for implicit curve*)

computeTangentNormalImplicit[p_] := 
  Module[{grad, normGrad, normal, tangent}, 
   grad = {Fx[p[[1]], p[[2]]], Fy[p[[1]], p[[2]]]};
   normGrad = Norm[grad];
   If[normGrad < 1*^-10, Return[{Null, Null}]];
   normal = grad/normGrad;
   tangent = {-normal[[2]], normal[[1]]};
   {tangent, normal}];

(*Compute curvature for implicit curve*)

computeCurvatureImplicit[p_] := 
  Module[{x, y, fx, fy, fxx, fyy, fxy, numerator, denominator}, 
   x = p[[1]]; y = p[[2]];
   fx = Fx[x, y]; fy = Fy[x, y];
   fxx = Fxx[x, y]; fyy = Fyy[x, y]; fxy = Fxy[x, y];
   numerator = fxx*fy^2 - 2*fxy*fx*fy + fyy*fx^2;
   denominator = (fx^2 + fy^2)^1.5;
   If[Abs[denominator] < 1*^-10, Return[0.0]];
   -numerator/denominator];

(*Compute focal point for implicit curve*)

computeFocalPointImplicit[p_, RMax_ : 10.0] := 
  Module[{T, N, kappa, R}, {T, N} = computeTangentNormalImplicit[p];
   If[N === Null, Return[Null]];
   kappa = computeCurvatureImplicit[p];
   If[Abs[kappa] < 1*^-10, Return[Null]];
   R = 1.0/kappa;
   If[Abs[R] > RMax, Return[Null]];
   p + N*R];

(*Compute focal set*)

computeFocalSetImplicit[points_, RMax_ : 3.0] := 
  Module[{focalPoints = {}, fp}, 
   Do[fp = computeFocalPointImplicit[p, RMax];
    If[fp =!= Null && VectorQ[fp, NumericQ], 
     AppendTo[focalPoints, fp]];, {p, points}];
   focalPoints];

(*Compute symmetry set*)

computeSymmetrySetImplicit[points_, lambdaMax_ : 3.0] := 
  Module[{centers = {}, nPoints, thresholdIdx, i, j, p1, p2, T1, N1, 
    T2, N2, prevCheck1, prevCheck2, prevP2, diff, Tsum, Tdiff, 
    eqCheck1, eqCheck2, alpha, pInterp, TInterp, NInterp, Nvec, denom,
     lam}, nPoints = Length[points];
   If[nPoints < 2, Return[{}]];
   thresholdIdx = Floor[0.05*nPoints];
   Do[p1 = points[[i]];
    {T1, N1} = computeTangentNormalImplicit[p1];
    If[T1 === Null, Continue[]];
    prevCheck1 = Null;
    prevCheck2 = Null;
    prevP2 = Null;
    Do[p2 = points[[j]];
     {T2, N2} = computeTangentNormalImplicit[p2];
     If[T2 === Null, Continue[]];
     (*Skip if tangents are parallel or anti-parallel*)
     If[Norm[T1 - T2] < 1*^-5 || Norm[T1 + T2] < 1*^-5, 
      prevCheck1 = Null;
      prevCheck2 = Null;
      prevP2 = p2;
      Continue[]];
     (*Compute tangency conditions*)diff = p1 - p2;
     Tsum = T1 + T2;
     Tdiff = T1 - T2;
     eqCheck1 = diff.Tsum;
     eqCheck2 = diff.Tdiff;
     (*Check for sign changes*)
     If[prevP2 =!= Null,(*Check eq1 sign change*)
      If[prevCheck1 =!= Null && prevCheck1*eqCheck1 < 0, 
       alpha = Abs[prevCheck1]/(Abs[prevCheck1] + Abs[eqCheck1]);
       pInterp = prevP2 + alpha*(p2 - prevP2);
       {TInterp, NInterp} = computeTangentNormalImplicit[pInterp];
       If[NInterp =!= Null, Nvec = N1 + NInterp;
        denom = Nvec.Nvec;
        If[Abs[denom] > 1*^-8, lam = (pInterp - p1).Nvec/denom;
         If[Abs[lam] < lambdaMax, AppendTo[centers, p1 + lam*N1]]]]];
      (*Check eq2 sign change*)
      If[prevCheck2 =!= Null && prevCheck2*eqCheck2 < 0, 
       alpha = Abs[prevCheck2]/(Abs[prevCheck2] + Abs[eqCheck2]);
       pInterp = prevP2 + alpha*(p2 - prevP2);
       {TInterp, NInterp} = computeTangentNormalImplicit[pInterp];
       If[NInterp =!= Null, Nvec = N1 - NInterp;
        denom = Nvec.Nvec;
        If[Abs[denom] > 1*^-8, lam = (pInterp - p1).Nvec/denom;
         If[Abs[lam] < lambdaMax, 
          AppendTo[centers, p1 + lam*N1]]]]]];
     (*Update previous values*)prevCheck1 = eqCheck1;
     prevCheck2 = eqCheck2;
     prevP2 = p2;, {j, i + thresholdIdx, nPoints}];, {i, 1, 
     nPoints}];
   centers];

(*Main execution*)
Print["1. Tracing the implicit curve..."];
resolution = 5000;
curvePoints = traceCurveLoop[resolution];
Print["Found ", Length[curvePoints], " points on the curve."];

Print["2. Computing focal set..."];
RMax = 10.0;
focalPoints = computeFocalSetImplicit[curvePoints, RMax];
Print["Found ", Length[focalPoints], " points in the focal set."];

Print["3. Computing symmetry set..."];
lambdaMax = \[Infinity];
symmetryPoints = computeSymmetrySetImplicit[curvePoints, lambdaMax];
Print["Found ", Length[symmetryPoints], 
  " points in the symmetry set."];

formattedCurve = Insert[#, "curve", 1] & /@ curvePoints;
formattedFocal = Insert[#, "focal", 1] & /@ focalPoints;
formattedSymmetry = Insert[#, "symmetry", 1] & /@ symmetryPoints;

(*Combine all data into one big list*)

allData = Join[formattedCurve, formattedFocal, formattedSymmetry];

(*Add a header row for clarity in the CSV*)

header = {"Type", "X", "Y"};
exportData = Prepend[allData, header];
Export["potato_plot_new.csv", exportData]