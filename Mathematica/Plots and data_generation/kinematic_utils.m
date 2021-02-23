(* ::Package:: *)

BeginPackage["`kinematicutils`"]


(*rotx::usage = " rotation matrix corresponding to x  axis "
roty::usage = " rotation matrix corresponding to y  axis ";
rotz::usage = " rotation matrix corresponding to z  axis ";*)
axis2skew::usage = " converts the vector c1, c2, c3 to a skew symmetric matrix ";
cayley2r::usage = "takes as input, the elements of the skew symmetric matrix c1, c2 , c3 and returns the corresponding rotation matrix";
bezoutmatrix::usage = "bezoutmatrix[poly1,poly2,x] returns the bezout matrix formed from the coefficients of polynomial1 and polynomial2 expressed in the variable x";
polymult::usage = "polymult[poly1,poly2] multiplies 2 univariate polynomials, used only in Bezout module"; (* To be replaced later *)
polymultuv::usage = "polymultuv[poly1,poly2,x] multiplies 2 univariate polynomials";
polyadduv::usage = "polyadduv[poly1,poly2,x] adds 2 univariate polynomials";
linvelocityjacobian::usage = " returns the linear velocity jacobian for a particular manipulator. ";
angvelocityjacobian::usage = " returns J\[Omega]s or J\[Omega]b depending on the choice given. J\[Omega]s ->1 or J\[Omega]b -> 2 ";
mpaux1::usage = " Called by monopoly ";
monopoly::usage = " Splits a polynonomial into its monomials and power products ";
simpolytrig::usage = " monopoly followed by simplification ";
parsimpolytrig::usage = " parallel simpolytrig ";



size[x_] :=  Print["Size of the expression is:", N[ByteCount[x]/1024],  " KB."];

tic[]:=Block[{},$mytime1=SessionTime[];$mytime2=TimeUsed[];];
toc[]:=Block[{t1,t2,hrs,mins,secs},
t1=SessionTime[]-$mytime1;
t2=TimeUsed[]-$mytime2;
If[t1>=3600,
hrs=Quotient[t1,3600];
mins=Quotient[Mod[t1,3600],60];
secs=Mod[t1,60];
Print["Actual time consumed: ",hrs," hours, ",mins," minutes, ",secs," seconds."];];
If[3600>=t1>=60,
mins=Quotient[t1,60];
secs=Mod[t1,60];
Print["Actual time consumed: ",mins," minutes, ",secs," seconds."];];
Print["Actual time consumed: ",t1," seconds."];

If[t2>=3600,
hrs=Quotient[t2,3600];
mins=Quotient[Mod[t2,3600],60];
secs=Mod[t2,60];
Print["CPU time consumed: ",hrs," hours, ",mins," minutes, ",secs," seconds."];];
If[3600>=t2>=60,
mins=Quotient[t2,60];
secs=Mod[t2,60];
Print["CPU time consumed: ",mins," minutes, ",secs," seconds."];];
Print["CPU time consumed: ",t2," seconds."];
];


Begin["`Private`"]


(*rotx[\[Alpha]_]:=Module[{temp},
temp={{1,0,0},{0,Cos[\[Alpha]],-Sin[\[Alpha]]},{0,Sin[\[Alpha]],Cos[\[Alpha]]}};
Return[temp];
];*)


(*roty[\[Beta]_]:=Module[{temp},
temp={{Cos[\[Beta]],0,Sin[\[Beta]]},{0,1,0},{-Sin[\[Beta]],0,Cos[\[Beta]]}};
Return[temp];
];*)


(*rotz[\[Gamma]_]:=Module[{temp},
temp={{Cos[\[Gamma]],-Sin[\[Gamma]],0},{Sin[\[Gamma]],Cos[\[Gamma]],0},{0,0,1}};
Return[temp];
];*)


axis2skew[x_] := Module[{},
      
      (*First check to make sure we have a skew symmetric matrix*)
      
      If[ !VectorQ[x] || Dimensions[x] != 3,
         Message[utils::wronginput, "axis2skew", x, 
         "is not a 3-vector"];
        Return[];
        ];
      
      (*Now extract the appropriate component*)
      
      Return[{{0, -x[[3]], x[[2]]},
              {x[[3]], 0, -x[[1]]},
              {-x[[2]], x[[1]], 0}}];       
      ];


linvelocityjacobian[p_, q_] := Module[{temp, x},
   x = Length[p];
   temp = Table[D[p[[i]], {q}], {i, 1, x}];
   Return[temp];
   ];


angvelocityjacobian[R_, q_, x_] := 
  Module[{temp, \[CapitalOmega]s, \[CapitalOmega]b, \[Omega]s, \[Omega]b, J\[Omega]s, J\[Omega]b},
   
  If[R == IdentityMatrix[3],
    temp = ConstantArray[0, {3, Length[q]}];
    Return[temp];
    ];
   \[CapitalOmega]s = Simplify[D[R,t].Transpose[R]];
   \[CapitalOmega]b = Simplify[Transpose[R].D[R,t]];
   \[Omega]s = {\[CapitalOmega]s[[3,2]], \[CapitalOmega]s[[1,3]], \[CapitalOmega]s[[2,1]]};
   \[Omega]b = {\[CapitalOmega]b[[3,2]], \[CapitalOmega]b[[1,3]], \[CapitalOmega]b[[2, 1]]};

   Print[R];
   Print[D[R,t]];
   
   J\[Omega]s = Normal[CoefficientArrays[{\[Omega]s[[1]], \[Omega]s[[2]], \[Omega]s[[3]]}, q]][[2]];
   J\[Omega]b = Normal[CoefficientArrays[{\[Omega]b[[1]], \[Omega]b[[2]], \[Omega]b[[3]]}, q]][[2]];
   If[x == 1, temp = J\[Omega]s;];
   If[x == 2, temp = J\[Omega]b;];
   Return[temp];
   ];



cayley2r[{c1_, c2_, c3_}] := 
    Simplify[(IdentityMatrix[3] + axis2skew[{c1, c2, c3}]).Inverse[IdentityMatrix[3] - axis2skew[{c1,c2,c3}]]];


mpaux1[coeffs_, powers_, x_, simplify_:Simplify] :=
    
    Module[{n, coeffnew, pownew, cx, m, coeff, pow, ii},
      (* 
        inputs : 
          coeffs : list of coefficients of the input power products 
              powers : power products 
        *)
      
      n = Length[powers];
      
      coeffnew = {};
      pownew = {};
      
      For[ii = 1, ii <= n,
        cx = CoefficientList[coeffs[[ii]], x];
        m = Length[cx];
        pow = Complement[Range[m], Flatten[Position[cx, 0]]] - 1;
        coeff = cx[[pow + 1]];
        pownew = Flatten[Append[pownew, powers[[ii]]*x^pow]];
        coeffnew = Flatten[Append[coeffnew, coeff]];
        ii++;];
      
      Return[{simplify[coeffnew], pownew}];
      ];

(* the main function *)




monopoly[poly_, vars_, simplify_:Identity] :=
    
    Module[{x, cx, n, nonzerocoeff, pow, nzpowcoeffx, powx, coeffs, powers,jj},
      If[poly == 0, Return[0]];
      
        If[vars ==  {},
            Print["monopoly: variable list is empty"];
            Return[]
            ];
        
      If[And @@ Map[FreeQ[poly, #] &, vars] == True,
	(*        
	Print["monopoly: input 1 is not a polynomial in input variable(s)"];
	*)  
      Return[poly];
        ];
      
      (* take out the first variable *)
      
      x = vars[[1]];
      
      (* find the coefficientlist *)
      
      cx = CoefficientList[poly, x];
      
      n = Length[cx];
      
      (* find out the nonzero coefficeints of powers of x *) 
      
      (* set of powers of x *)
      
      pow = Complement[Range[n], Flatten[Position[cx, 0]]] - 1;
      
      (* corresponding set of coefficients, non zero coeffs. of x *)
      
      nzcoeffx = cx[[pow + 1]];
      
      powx = x^pow;
      
      (* Return	at this stage if univariate *)

      If[Length[vars] == 1, 
	Return[{nzcoeffx, powx}];
	];
          
      (* pass it onto the next level for expanding wrt the next variable *)
  
    
      varlength = Length[vars];
      
      {coeffs, powers} = {nzcoeffx, powx};
      
      For[jj = 2, jj <= varlength,
        
        {coeffs, powers} = mpaux1[coeffs, powers, vars[[jj]], simplify];
        
             jj++;
        ];
      Return[{simplify[coeffs], powers}];
      ];



simpolytrig[expr_, vars_List,  simplify_:Simplify, expand_:Identity] :=
Module[{}, 
(* multivariate case *)
 
 If[And @@ Map[FreeQ[expr, #] &, vars] == True,
	(*        
	Print["monopoly: input 1 is not a polynomial in input variable(s)"];
	*)  
      Return[expr];
        ];

Return[Dot@@monopoly[expr,vars,simplify]]; 
];


parsimpolytrig[expr_,vars_List,N_,expand_: Identity]:=Module[{mp,mp1},If[And@@Map[FreeQ[expr,#]&,vars]==True,Return[expr];];
If[N==0,Return[Dot@@monopoly[expr,vars,Simplify]]];
If[N!=0,mp=monopoly[expr,vars];
LaunchKernels[N];
SetSharedVariable[mp1];
mp1=ConstantArray[0,Length[mp[[1]]]];
ParallelDo[mp1[[i]]=Simplify[mp[[1]][[i]]],{i,Length[mp[[1]]]}];
CloseKernels[];
(*Clear[mp[[1]]];*)Return[Dot[mp1,mp[[2]]]];];];


(* ::Input:: *)
(**)


(*simpolytrig[expr_, var_,  simplify_:Simplify, expand_:Identity] :=
Module[{}, 
(* univariate case *)

Return[Collect[expand[expr], var, simplify]];
];*)


polyadduv[poly1_, poly2_, x_] := 
  Module[{coeffvec1, coeffvec2, i, j, sum, temp, temp1, polysum, n1, 
    n2, p, len, len1}, 
   coeffvec1 = Reverse[CoefficientList[poly1, x]];
   coeffvec2 = Reverse[CoefficientList[poly2, x]];
   n1 = Length[coeffvec1];
   n2 = Length[coeffvec2];
   len = Max[n1, n2];
   sum = ConstantArray[0, len];
   If[n1 == n2, sum = coeffvec1 + coeffvec2;];
   If[n1 > n2, temp1 = ConstantArray[0, n1 - n2];
    sum = coeffvec1 + RotateRight[Join[coeffvec2, temp1], n1 - n2];];
   If[n2 > n1, temp1 = ConstantArray[0, n2 - n1];
    sum = RotateRight[Join[coeffvec1, temp1], n2 - n1] + coeffvec2;];
   temp = Reverse[RotateRight[Append[x^Range[len - 1], 1], 1]];
   polysum = sum.temp;
   Return[polysum];
   ];


(*Module to multiply two polynomials*)

polymultuv[poly1_, poly2_, x_] := 
  Module[{coeffvec1, coeffvec2, coeffvec3, i, j, sum, temp, temp1, 
    polyprod, n1, n2, p, p1, p2, len, len1, x1},
   coeffvec1 = Reverse[CoefficientList[poly1, x]];
   coeffvec2 = Reverse[CoefficientList[poly2, x]];
   n1 = Length[coeffvec1];
   n2 = Length[coeffvec2];
   len1 = n1 + n2 - 1;
   If[(n1 == n2) || (n1 > n2), (p = 
       Outer[Times, coeffvec1, coeffvec2]) && (x1 = len1 - n2), (p = 
       Outer[Times, coeffvec2, coeffvec1]) && (x1 = len1 - n1)];
   sum = ConstantArray[0, len1];
   len = Length[p];
   temp1 = ConstantArray[0, x1];
   For[i = 1, i <= len, i++, 
    sum = sum + RotateRight[Join[p[[i]], temp1], i - 1];];
   temp = Reverse[RotateRight[Append[x^Range[n1 + n2 - 2], 1], 1]];
   polyprod = sum.temp;
   Return[polyprod];
   ];


(* Module to multiply two univariate polynomials by manipulation of the coefficient vectors alone *)
polymult[poly1_,poly2_]:= Module[{coeffvec3,i,j,coeffvec,temp,temp1,augment,productcoeffvec,n1,n2,p,p1,p2,len,nproduct,x1},

n1=Length[poly1];

n2=Length[poly2];

(* Multiply the 2 lists using outer , define the length of the number of zeros to be appended by the variable x1 *)
nproduct =n1+n2-1;
If[(n1==n2)||(n1>n2),
(p=Outer[Times,poly1,poly2])&&(x1=nproduct-n2),
(p=Outer[Times,poly2,poly1])&&(x1=nproduct-n1)];

(* productcoeffvec defines the final product coefficient vector, the length of the product vector is given by nproduct *)
productcoeffvec=ConstantArray[0,nproduct];
len=Length[p];
(* Each p, individual product vector will be augmented by x1 zeros, the list being defined by augment *)
augment=ConstantArray[0,x1]; 
(* Loop to shift and add *)
(* Take each product given by p[i], augment with zeros to bring it to the length of the final product vector, shift right by (i-1)th instance and add to the final coefficient vector each time *)
For[i=1,i<=len,i++,
productcoeffvec=productcoeffvec+RotateRight[Join[p[[i]],augment],i-1];
];
(* productcoeffvec givesthe final product coefficient list *)

Return[productcoeffvec];
]


(*Module to calculate Bezout's Matrix from 2 given polynomials*)
bezoutmatrix[p1_,p2_,x_]:=Module[{a,b,m,n,elim1,elim2,res,rem,i,c,p3,coeff,temp,bigger,smaller,c1,p,matrix,deg1,deg2,ap,poly1,poly2,poly3,poly4,poly5,poly6},
If[Exponent[p1,x]>Exponent[p2,x],
bigger=p1;
smaller=p2;,
bigger=p2;
smaller=p1;
];
a=Reverse[CoefficientList[bigger,x]];
b=Reverse[CoefficientList[smaller,x]];
If[Exponent[p1,x]==Exponent[p2,x],
a=Reverse[CoefficientList[p1,x]];
b=Reverse[CoefficientList[p2,x]];
];
m=Length[a];
deg1=m-1;(*degree of polynomial 1*)
n=Length[b];
deg2=n-1;(*degree of polynomial 2*)

(* Case 1: When degree of two polynomials are equal,deg1=deg2, then n-1 linearly independent equations are generated by eliminating each monomial*)
If[m==n,
For[i=1,i<=n-1,i++,
c=i;
poly1=a[[1;;i]];
poly2=b[[i+1;;n]];
poly3=polymult[poly1,poly2];
poly4=b[[1;;i]];
poly5=a[[i+1;;n]];
poly6=polymult[poly4,poly5];
coeff[i]=Reverse[poly3-poly6];
];
];
(* Case 2: When degree of one polynomial is greater than the other one, deg1>deg2, then at first m-n linearly independent equations are generated *)
If[m>n,
p3=Expand[smaller*x^(m-n)];
For[i=1,i<=n-1,i++,
c=i;
elim1=a[[1;;i]].Reverse[RotateRight[x^Append[Range[i-1],0]]];
elim2=b[[1;;i]].Reverse[RotateRight[x^Append[Range[i-1],0]]];
res=CoefficientList[Collect[(p3*elim1-bigger*elim2),x],x];
coeff[i]=res;
];
(*the remaining equations are generated by multiplying the lower degree equation by powers of m-n-1,m-n-2....0 of the univariate, x in this case*)
ap=ConstantArray[0,m-n-1];
rem=Join[Reverse[b],ap];
For[i=1,i<=m-n,i++,
If[i==1,coeff[c+1]=rem;,
coeff[c+1]=RotateRight[rem];
];
rem=coeff[c+1];
(*Print[coeff[c+1]];*)
c=c+1;
];
];
temp=coeff[1];
For[i=2,i<=c,i++,
temp=Join[temp,coeff[i]];
];
matrix=Partition[temp,c];
Return[matrix];
];


End[]
EndPackage[]
