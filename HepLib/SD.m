(* ::Package:: *)

(* ::Section:: *)
(*M2C*)


Unprotect[pow];
Clear[pow];
pow::usage="pow in GiNaC";
Protect[pow];


Clear[M2C];
M2C[exp_]:=Module[{ret},
ret=exp/.Rule->List;
ret=ret/.{EulerGamma->Euler,Zeta->zeta,Gamma->tgamma,Power->pow,Log->log,Sqrt->sqrt};
ret=ToString[InputForm[ret],PageWidth->Infinity];
ret=StringReplace[ret,{"["->"(", "]"->")"}];
ret
];


(* ::Section:: *)
(*C2M*)


Clear[C2M]
C2M[cex_String]:=Module[{ret},
ret=StringReplace[cex,{
"VE"->"VE","iEpsilon"->"iEpsilon","Euler"->"EulerGamma","zeta"->"Zeta","tgamma"->"Gamma","pow"->"Power","log"->"Log","sqrt"->"Sqrt","E"->"*10^","+-"->"+0*"
}];
ret=ToExpression[ret,TraditionalForm];
ret=ret/.{psi[ex_]:>PolyGamma[1,ex]};
Return[ret];
];


RE[cex_String]:=C2M[cex]/.Complex[r_,_]:>r


(* ::Section::Closed:: *)
(*VE*)


Unprotect[VE];
Clear[VE];
VE::usage="Value and Error Wrapper";
VE/:VE[0,0]:=0;
VE/:VE[Complex[vr_,vi_],Complex[er_,ei_]]:=VE[vr,er]+I VE[vi,ei];
VE/:VE[v_?NumericQ,e_]/;Element[v,Reals]&&v<0:=-VE[-v,e];
VE/:MakeBoxes[VE[v_?NumericQ,0],TraditionalForm]:=StyleBox[MakeBoxes[v,TraditionalForm],FontColor->Blue];
VE/:MakeBoxes[VE[v_?NumericQ,e_/;Abs[e]<=0],TraditionalForm]:=StyleBox[MakeBoxes[v,TraditionalForm],FontColor->Blue];


Clear[CLog];
CLog[ex_]:=Ceiling[Log[10,Abs[ex]]];
CLog[ex_]/;Abs[Rationalize[ex,0]]<Power[10,-100]:=-100;
VEWidth=1;


VE/:MakeBoxes[VE[v0_/;Element[v0,Reals],e0_?NumericQ],TraditionalForm]:=Module[{v,e,n,v1,e1,a,en},
n=Max[CLog[v0],CLog[e0]];
{v,e}={Abs[v0],Abs[e0]}/Power[10,n];
e1=N[Rationalize[e,0],VEWidth];
e1=AccountingForm[e1,VEWidth]//ToString;
If[StringStartsQ[e1,"1."],
e1=StringReplace[e1,"1."->"0.1"];
n=n+1;
v=v/10;
];
a=StringLength[e1]-2;
v=N[Round[v,Power[10,-a]],a+2];
v1=AccountingForm[v]//ToString;
If[v1==="0",v1="0."];
On[Assert];Assert[Or[StringStartsQ[v1,"1."],StringStartsQ[v1,"0."]]];
If[StringLength[v1]<a+2,v1=v1<>StringJoin[Table["0",a+2]]];
v1=StringTake[v1,a+2];
e1=StringTake[e1,-VEWidth];
If[n=!=0,RowBox[{v1,"(",e1,")","\[Times]",SuperscriptBox[10,n]}],RowBox[{v1,"(",e1,")"}]]
];


Protect[VE];


Clear[ChopVE];
ChopVE[exp_,diff_:Power[10,-6]]:=exp/.VE[v_/;Abs[v]<diff,e_/;Abs[e]<diff]->0;


Clear[SimplifyVE];
SimplifyVE[exp_]:=Module[{tmp,VEH,VF,IH},
tmp=Cases[exp,s_Symbol/;Not[NumericQ[s]],{0,Infinity}]//Union;
If[Length[tmp]=!=0,Print["Extra Symbols Found in SimplifyVE: ",tmp];Return[tmp]];
tmp=exp/.VE->VEH;
tmp=tmp/.VEH[Complex[r_,i_],e_]:>VEH[r,e]+I VEH[i,e];
tmp=Collect[tmp,_VEH,VF];
tmp=tmp/.VF[0]->0;
tmp=tmp//.Dispatch[{
VEH[v1_,e1_]^2:>VEH[v1^2,Abs[2v1 e1]],
VEH[v1_,e1_] VEH[v2_,e2_]:>VEH[v1 v2,Sqrt[Abs[v1 e2]^2+Abs[v2 e1]^2]]
}];
tmp=tmp//.Dispatch[{
VF[Complex[r_,i_]]:>(VF[r]+IH VF[i]),
VF[-c_]:>-VF[c]
}];
tmp=Collect[tmp,{IH,_VEH}];
tmp=tmp/.VF[0]->0;
tmp=tmp/.VF[c_] VEH[v_,e_]:>VEH[c v,Abs[c e]];
tmp=tmp/.VF[c_]:>c;
tmp=tmp//.Dispatch[{
VEH[v1_,e1_]+VEH[v2_,e2_]:>VEH[v1+v2,Sqrt[e1^2+e2^2]],
r_+VEH[v_,e_]/;Element[r,Reals]:>VEH[r+v,e]
}];
tmp/.{VEH->VE,IH->I}
];


Clear[VESimplify];
VESimplify[exp_,ri_:Null]:=Module[{tmp,VF,VE2,IM,resR,resI},
On[Assert];Assert[SubsetQ[{True},Union[Map[NumericQ,Union[Cases[exp,_Symbol,{0,Infinity}]]]]]];
tmp=Expand[exp];
tmp=tmp/.Complex[a_,b_]:>a+b IM;
tmp=Expand[tmp,_VE|IM];
tmp=Distribute[VF[tmp]];
tmp=tmp/.VF[0]->0;
tmp=tmp/.VF[ex_]/;FreeQ[ex,_VE]:>VF[VE[ex,0]];
tmp=tmp/.{VF[c_. VE[e_,v_]]/;FreeQ[c,IM]:>VE2[c e,Power[Abs[c] v,2]]};
tmp=tmp/.{VF[c_. IM VE[e_,v_]]:>IM VE2[c e,Power[Abs[c] v,2]]};
On[Assert];Assert[FreeQ[tmp,_VF]];
resR=tmp/.IM->0/.VE2->List;
If[resR=!=0,
On[Assert];Assert[Length[resR]===2];
resR[[2]]=Sqrt[resR[[2]]];
];
resR=VE@@resR;
resI=Coefficient[tmp,IM]/.VE2->List;
If[resI=!=0,
On[Assert];Assert[Length[resI]===2];
resI[[2]]=Sqrt[resI[[2]]];
];
resI=VE@@resI;
tmp=Switch[ri
,Null,resR+resI I
,Re,resR
,Im,I resI];
Return[tmp];
];


(* ::Section:: *)
(*UF*)


Clear[UF];
UF[ls_,tls_,ps_,ns_,lr_:{},tlr_:{},nr_:{},isQuasi_:False]:=Module[{exe,proc,err,is,os,es,res},
exe="/usr/local/feng/bin/UF";
proc=StartProcess[exe,ProcessEnvironment-><|"DYLD_LIBRARY_PATH"->"/usr/local/feng/lib"|>];
is=ProcessConnection[proc,"StandardInput"];
os=ProcessConnection[proc,"StandardOutput"];
es=ProcessConnection[proc,"StandardError"];
err=ReadString[es,EndOfBuffer];
If[StringLength[StringTrim[err]]>0,Print[err];Abort[];Return[]];
WriteLine[is,M2C[If[isQuasi,1,0]]];
WriteLine[is,".end"];
WriteLine[is,M2C[ls]];
WriteLine[is,".end"];
WriteLine[is,M2C[tls]];
WriteLine[is,".end"];
WriteLine[is,M2C[ps]];
WriteLine[is,".end"];
WriteLine[is,M2C[ns]];
WriteLine[is,".end"];
WriteLine[is,M2C[lr]];
WriteLine[is,".end"];
WriteLine[is,M2C[tlr]];
WriteLine[is,".end"];
WriteLine[is,M2C[nr]];
WriteLine[is,".end"];
res=ReadString[os,EndOfFile];
If[res===EndOfFile||StringLength[StringTrim[res]]==0,Print[ReadString[es,EndOfBuffer]];Abort[]];
KillProcess[proc];
res
];
