(* ::Package:: *)

SetDirectory[NotebookDirectory[]]
<<"kinematic_utils.m"
(*BeginPackage["`sylvester`"]*)


dumeq::usage = " dummification of equations "
sylvmat::usage = " generates sylvester matrix for three equations in three unknowns "



dumeq[eqn_,vars_,dumvar_]:=Module[{eqmono,len,dumlist,i,varmap,dumexp},
eqmono=monopoly[eqn,vars];
(*Print[eqmono];*)
len=Length[eqmono[[2]]];
dumlist=Table[ToExpression[dumvar<>ToString[i]],{i,1,len}];
varmap=Inner[Rule,dumlist,eqmono[[1]],List];
dumexp=dumlist.eqmono[[2]];
Return[{varmap,dumexp}];
];


(*creating the partition list*)
sylvmat[syseqns_,deg_,totvar_,var_]:=Module[{deglist,partlist,lenlist,i,j,k,part1,partnew,newlist,evenc,oddc,binarr,bincoeff,lenarr,eqlist,eqlisttemp,leneqlistcheck,leneq,lenpart,detpart,p,temp,lentemp,detstage,parttemp,dettemp,eqnstemp,coeffarr,count,expvar,templist,coeff,subslist,eqmod,rem,sum,ptemp,deteq,detp,lenp,lendeteq,powprod,toteqns,leneqns,eq,powlist,match,augcoeff,list,str,strcoeff,l,matchtemp,matchcoeff,leneqlist,powprod1,a,dimmat,eq4dum,eq5dum,eq6dum,subsd,subse,subsf,eqns,subsaugcoeff1,subsaugcoeff1a,subsaugcoeff2,subsaugcoeff2a,subsaugcoeff3,subsaugcoeff3a},
(*deg=6;
totvar={c1,c2,c3};
var={c2,c3};*)
eq4dum=dumeq[syseqns[[1]],{var[[1]],var[[2]]},"d"];
eq5dum=dumeq[syseqns[[2]],{var[[1]],var[[2]]},"e"];
eq6dum=dumeq[syseqns[[3]],{var[[1]],var[[2]]},"f"];
(*Print[eq4dum];*)
subsd=eq4dum[[1]];
subse=eq5dum[[1]];
subsf=eq6dum[[1]];
eqns={eq4dum[[2]],eq5dum[[2]],eq6dum[[2]]};
deglist=Reverse[Table[i,{i,2,deg+1}]];
partlist=ConstantArray[0,Length[deglist]];
lenlist=Length[deglist];
For[i=1,i<=lenlist,i++,
part1=IntegerPartitions[deglist[[i]],{2}];
partnew=ConstantArray[0,Length[part1]];
newlist=ConstantArray[0,2*Length[part1]];
partnew=Table[Reverse[part1[[j]]],{j,1,Length[part1]}];
evenc=1;
oddc=1;
For[j=1,j<=Length[newlist],j++,
If[OddQ[j],newlist[[j]]=part1[[evenc++]],newlist[[j]]=partnew[[oddc++]]];
];
newlist=DeleteDuplicates[newlist];
partlist[[i]]=newlist;
];

(*generate power product terms of a polynomial of degree: n-2*)
binarr=ConstantArray[0,(deg-2)];
bincoeff=monopoly[Expand[(1+var[[1]]+var[[2]])^(deg-2)],var];
binarr=bincoeff[[2]];
binarr=Flatten[binarr];
lenarr=Length[binarr];
(*multiply each equation by each term of a polynomial degree n-2 to generate (3/2)n(n-1) independent equations*)
Print["multiplying each equation by each term of a polynomial degree ", deg-2," to generate ", (3/2)*deg*(deg-1)," independent equations"];
Print["....."];
eqlist=ConstantArray[0,Length[eqns]];
For[i=1,i<=Length[eqns],i++,
eqlisttemp=ConstantArray[0,lenarr];
For[j=1,j<=lenarr,j++,
eqlisttemp[[j]]=Expand[eqns[[i]]*binarr[[j]]];
];
eqlist[[i]]=eqlisttemp;
];
eqlist=Flatten[eqlist];
leneqlistcheck=Length[eqlist];
Print["Generated ",leneqlistcheck," equations"];

(*group the terms according to the degree partition and form the determinants to find (1/2)n(n+1) independent equations*)
(*Print["Now grouping the terms according to the degree partition and form the determinants to find ",(1/2)*deg*(deg+1)," independent equations"];
Print["...."];*)
leneq=Length[eqns];
lenpart=Length[partlist];
detpart=ConstantArray[0,lenpart];
For[p=1,p<=lenpart,p++,
temp=partlist[[p]];
lentemp=Length[temp];
detstage=ConstantArray[0,lentemp];
For[l=1,l<=lentemp,l++,
parttemp=temp[[l]];
(*Print[" grouping : ",parttemp];*)
dettemp=ConstantArray[0,3];
For[k=1,k<=leneq,k++,
eqnstemp=eqns[[k]];
(*Print[" Now grouping equation ",k,"...."];*)
coeffarr=ConstantArray[0,3];
count=1;
For[i=1,i<=2,i++,
(*Print["exponent of ",var[[i]],"  required : ",parttemp[[i]]];*)
expvar=Exponent[eqnstemp,var[[i]]];
templist=Select[Table[j-parttemp[[i]],{j,parttemp[[i]],expvar}],#>= 0&];
(*Print[templist];*)
If[Length[templist]==0,
(*Print["coefficient of the exponent does not exist"];*)
coeff=0,
subslist=Table[var[[i]]^(parttemp[[i]]+templist[[j]])->var[[i]]^templist[[j]]*a,{j,1,Length[templist]}];
(*Print["substitution list "];*)
(*Print[subslist];*)
eqmod=eqnstemp/.subslist;
(*Print[" equation after substituting ..... \n",eqmod];*)
coeff=Coefficient[eqmod,a];
];
(*Print[" coefficient extracted ..... \n",coeff];
Print[ "end of iteration ",i,"....."];*)
rem=Simplify[eqnstemp-Expand[(coeff*a)/.{a->var[[i]]^parttemp[[i]]}]];
eqnstemp=Expand[rem];
coeffarr[[count++]]=coeff;
];
(*Print[" Remainder term after grouping .... "];
Print[Expand[rem]];*)
coeffarr[[count]]=Expand[rem];
dettemp[[k]]=coeffarr;
];
detstage[[l]]=dettemp;
(*Print[" Determinant generated : "];
Print[dettemp];
Print["................................................................................."];*)
];
detpart[[p]]=detstage;
];

sum=0;
For[i=1,i<=lenpart,i++,
ptemp=partlist[[i]];
sum=sum+Length[ptemp];
];

deteq=ConstantArray[0,sum];
count=1;
For[i=1,i<=lenpart,i++,
detp=detpart[[i]];
lenp=Length[detp];
For[j=1,j<=lenp,j++,
deteq[[count++]]=Expand[Det[detp[[j]]]];
];
];
(*Print[" Generated ",count-1," equations "];
Print[" ************** "];
Print[" Generating list of power products of the ",(3/2)*deg*(deg-1)," equations ....."];*)
leneqlist=Length[eqlist];
For[i=1,i<=leneqlist,i++,
powprod1=monopoly[eqlist[[i]],{var[[1]],var[[2]]}][[2]];
(*Print[powprod1];
Print[Exponent[eqlist[[i]],{c2,c3}]];*)
];
Print[" ****************** "];
Print[" Now generating list of power products of the ",(1/2)*deg*(deg+1)," equations ....."];
lendeteq=Length[deteq];
For[i=1,i<=lendeteq,i++,
powprod=monopoly[deteq[[i]],{var[[1]],var[[2]]}][[2]];
(*Print[powprod];
Print[Exponent[deteq[[i]],{c2,c3}]];*)
];
(* Join the equations list, to form the complete set of the equations*)
toteqns=Join[eqlist,deteq];
leneqns=Length[toteqns];
(*Print[" Total number of equations generated: ",leneqns];*)
eq=Table[monopoly[toteqns[[i]],{var[[1]],var[[2]]}],{i,1,leneqns}];
powlist=Reverse[monopoly[Expand[(1+var[[1]]+var[[2]])^((2*deg)-2)],{var[[1]],var[[2]]}][[2]]];
Length[powlist];
(*Print[" ********** "];
Print[" Forming the augmented matrix ...... "];*)
match=ConstantArray[0,leneqns];
augcoeff=ConstantArray[0,leneqns];
list=Table[i,{i,1,leneqns}];
For[i=1,i<=leneqns,i++,
str=eq[[i]][[2]];
strcoeff=eq[[i]][[1]];
l=Length[str];
matchtemp=ConstantArray[0,leneqns];
matchcoeff=ConstantArray[0,leneqns];
For[j=1,j<=Length[powlist],j++,
For[k=1,k<=l,k++,
If[powlist[[j]]===str[[k]],matchtemp[[j]]=j;matchcoeff[[j]]=strcoeff[[k]];];
];
];
(*Print[matchtemp];*)
match[[i]]=Complement[list,DeleteCases[matchtemp,0]];
augcoeff[[i]]=matchcoeff;
(*Print[augcoeff[[i]]];*)
];
dimmat=Dimensions[augcoeff];
Print[" Dimension of the Sylvester matrix generated: ",dimmat];
subsaugcoeff1=augcoeff/.subsd;
subsaugcoeff1a=Table[Expand[subsaugcoeff1[[i]]],{i,1,dimmat[[1]]}];
subsaugcoeff2=subsaugcoeff1a/.subse;
subsaugcoeff2a=Table[Expand[subsaugcoeff2[[i]]],{i,1,dimmat[[1]]}];
subsaugcoeff3=subsaugcoeff2a/.subsf;
subsaugcoeff3a=Table[Expand[subsaugcoeff3[[i]]],{i,1,dimmat[[1]]}];
Return[{subsaugcoeff3a,powlist}];
];


(*End[]*)
(*EndPackage[]*)
