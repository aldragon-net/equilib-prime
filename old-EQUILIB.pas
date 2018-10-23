PROGRAM EQUILIB;      { Расчет равновесных параметров газа за ударной волной }
uses CRT,DOS,EQLEXTRA,CONSTANT;{ Все величины записаны в системе СИ          }
label  3,4,5;         { Атмосфеpы физические                                 }
const                 { Alfa:=Alfa*Pi/180; u1:=u1*Sin(Alfa)                  }
     ns = 100;        { Delta:=Alfa-atan(sin(Alfa)/cos(Alfa)/r2r1)           }
     ms = 15;         { Delta:=Delta*180/Pi                                  }
{----------------------------------------------------------------------------}
type
    TInfoSubst = record FName :string[12]; SPart :single
                 end;
    TFile = array[1..ms] of TInfoSubst;
    TInfoSpec = record Name :string[12]; Part,Weight :single; Cp,Entl :array[1..ns] of single
                end;
    TArrayRec = array[1..ms] of TInfoSpec;
    F3V = Function(x1, x2 :ABE; x3 :extended) :extended;
{----------------------------------------------------------------------------}
var
   F                            :text;
   Subst1,Subst4                :TFile;
   Mix1,Mix4                    :TArrayRec;
{   StrT1,StrT4,StrP1,StrP4      :string[2];}
   UnitU                        :string[5];
   UnitP                        :string[6];
   UnitT1,UnitT4                :string[3];
   StFOut,StFDat                       :string[12];
   StPath                              :string;
   Temp                                :array[1..ns] of single;
   EntMix1,EntMix4,Gmix1,Gmix2,Gmix4,T :ABE;
   WeightMix1,WeightMix4,H,G1,G2:extended;
   T1,T2,T4,T5,r2r1,r5r1,x,base :extended;
   P1,P4,P2P1,P4P1,P5P1,n1      :extended;
   u1,M,Mr,a1,a2,a4,Alfa,Delta  :extended;
   base1,base2,base3,G4         :extended;
   i,ModeCalc,is1,is4,jmax,mm   :byte;
   Number                       :word;
{----------------------------------------------------------------------------}
                                  {$F+}
 Function IntLin(y,x:ABE; s:extended):extended;
 Var i :byte;
 Begin {--------------------------IntLin-------------------------------------}
   i:=1;
   While (x[i]<s)and(i<61) do
        Inc(i);
   IntLin:=(y[i]-y[i-1])*(s-x[i-1])/(x[i]-x[i-1])+y[i-1]
 End; {---------------------------IntLin-------------------------------------}
                                  {$F-}
{----------------------------------------------------------------------------}
 Procedure OutRes( var F:text);
 Var T:ABE;
 Begin {---------------------------OutRes------------------------------------}
     Window(38,5,78,21); ClrScr;
     Writeln(F,'    Температура       T2 = ',T2:6:0,' [K]');
     Writeln(F,'      Отношение             ');
     Writeln(F,'     плотностей  ro2/ro1 = ',r2r1:6:3);
     Writeln(F,'       Давление       Р2 = ',P2P1*P1/1.01325E5:6:3,' [атм]');
     Writeln(F,'  Концентpация [см-3] n2 =',P2P1*P1/(Kb*T2)*1E-6:12);
     Writeln(F,'       Скорость       u1 = ',u1:6:0,' [м/c]');
     Writeln(F,'     Число Маха        M = ',M:6:3);
     Writeln(F);
     Writeln(F,'    Температура       T5 = ',T5:6:0,' [K]');
     Writeln(F,'      Отношение            ');
     Writeln(F,'     плотностей  ro5/ro1 = ',r5r1:6:3);
     Writeln(F,'       Давление       P5 = ',P5P1*P1/1.01325E5:6:3,' [атм]');
     Writeln(F,'  Концентpация [см-3] n5 =',P5P1*P1/(Kb*T5)*1E-6:12);
     Writeln(F,'     Число Маха       Mr = ',Mr:6:3);
     Writeln(F)
 end; {---------------------------OutRes-------------------------------------}

 Procedure IterTemp5(IntLin:F3V; Sign:F2V; var T5,r5r1,P5P1:extended);
 Var
{    T,H : ABE; }
    x,T5a,T5b,H_2,H5,s0,s1,GT2 :extended;
    js                         :shortint;
 Begin {--------------------------IterTemp5----------------------------------}
   x:=20;
   T5a:=T1*(2*(G1-1)*M*M+3-G1)*(M*M*(3*G1-1)-2*G1+2)/SQR((G1+1)*M);
   js:=0;
   H_2:=IntLin(EntMix1,T,T2);
   Repeat H5:=IntLin(EntMix1,T,T5a);
          s0:=5E-4*SQR(u1*(1-T2/T1/P2P1));
          s1:=(H5-H_2)/WeightMix1;
          T5b:=(s1-s0)*(T2/(s1+s0)+1E3*WeightMix1/R);
          If js*Sign(T5a,T5b)<=0 then x:=0.1*x;
          js:=Sign(T5a,T5b);
          If T5a>T5b then T5a:=T5a+x
                     else T5a:=T5a-x
   Until Abs(T5a-T5b)<1E-8;
   T5:=0.5*(T5a+T5b);
   r5r1:=r2r1*(s1+s0)/(s1-s0);
   P5P1:=r5r1*T5/T1;{ ?   }
{   Mr:=SQRT((0.5+0.5/GT2)*P5P1/P2P1+0.5-0.5/GT2);}

 end; {---------------------------IterTemp5----------------------------------}
 Procedure CalcVel(var M,P4P1 :extended);
 Var s,x,P        :extended;
     js           :shortint;
     i            :byte;
 Begin {--------------------------CalcVel------------------------------------}
   M:=1.5; x:=1; i:=0;
   P4P1:=P4{*1.01325E5}/P1; js:=0; P:=0;
   While (Abs(P4P1-P)>1E-9)and(i<250) do begin
         Inc(i);
         P:=(2*M*M*G1-G1+1)/(G1+1)*
             Exp(-2*G4/(G4-1)*Ln(1-(M-1/M)*(G4-1)/(G1+1)*a1/a4));
         If js*Sign(P4P1,P)<=0 then x:=0.1*x;
         js:=Sign(P4P1,P);
         If P4P1<P then M:=M-x
                   else M:=M+x
                            end { While }
 end; {---------------------------CalcVel------------------------------------}
 Procedure RevCalc_u1(var M:extended);
  Label 1;
  Var T2a,T2b,Wm4_1,x       :extended;
      js                    :shortint;
      i,i1,j                :word;
  Begin {-------------------------RevCalc_u1------------------------------------}
{  1: Window(3,e+5,35,23); TextBackGround(1); TextColor(14);  ClrScr;}
     Write(' Температура     T2 =  ');
     Read(T2a);
{     If T2a<1.2*T1 then Goto 1;}
     M:=1.5; x:=1; T2b:=100; js:=0;
     While Abs(T2a-T2b)>1E-6 do begin
          GetTime(i1, i1, i, i1);
          If i-j<>0 then
            Write('');
          j:=i;
{          IterPress(IntLin,Sign,u1,T2b,r2r1,P2P1);}
          If js*Sign(T2a,T2b)<=0 then x:=0.1*x;
          js:=Sign(T2a,T2b);
          If T2a<T2b then M:=M-x
                     else M:=M+x
                                end; { While
    Write(' Давление  P4 [атм] =  ');
   Read(P4c); P4c:=P4c*1.01325E5;
    Wm4_1:=GT4/GT1*Sqr((GT1+1)/(GT4-1)*(1/(M-1/M)))*
            Sqr(1-Exp((GT4-1)/(2*GT4)*Ln(P1/P4c*(2*M*M/(1+1/GT1)-(GT1-1)/(GT1+1)))));
     If jrap=1 then
       x:=(Wm4_1*Wmix*1E3-Wh2)/(28.97-Wh2)
               else
       x:=(Wm4_1*Wmix*1E3-Whe)/(28.97-Whe);
     Writeln(' Молекулярный вес      [г/моль]   толкающего газа W4 = ',Wm4_1*Wmix*1E3:5:2);
     Writeln('      Доля воздуха  = ',x*100:6:2,'%')}
  end; {--------------------------RevCalc_u1------------------------------------}
 Procedure InpPres(var P1,n1:extended);
 Begin {--------------------------InpPres------------------------------------}
{   Window(3,e+3,35,23);}
   TextBackground(1);Textcolor(14);
   ClrScr;
   Write('  Начальное                           давление  P1 =  ','    [торр]');
   GotoXY(21,2);
   Readln(P1);
   P1:=P1*133.322;
   n1:=P1/(Kb*T1)*1E-6;
   Write('  Hачальная                        концентpация n1 =',n1:12,
                                         '                         [см-3]');
 end; {---------------------------InpPres------------------------------------}
 Procedure OutFres(var P4P1:extended);
 Var i          :word;
     x          :extended;
 Begin {--------------------------OutFres------------------------------------}
     Write('  Hомеp экспеpимента?  N=  ');Read(i);
     Writeln(F,'  N  ',i);
{     x:=(GT1+1)/(GT4-1)*a4/a1;}
     If M>=0.5*(x+Sqrt(x*x+4)) then begin
     Writeln(' Hепpавильный выбоp состава      и/или темпеpатуpы газа T4 в КВД ');
     P4P1:=0                        end
                               else
{     P4P1:=(2*M*M*GT1-GT1+1)/(GT1+1)*
           Exp(-2*GT4/(GT4-1)*Ln(1-(M-1/M)/x));}
     Writeln(F,i:6,'  ',T1:5:1,' ',T2:5:0,' ',T5:6:0,' ',u1:5:0,' ',M:6:3,
                  ' ',P1/133.322:5:0,' ',P4:5:1,' ',P4*1.01325E5/P1:10:3,' ',
                  P4P1:10:3,' ',P2P1*P1/1.01325E5:8:3,' ',
                  P5P1*P1/1.01325E5:8:3,' ',r2r1:6:3,'  ',r5r1:7:3,' ',n1:12);
 end; {---------------------------OutFres------------------------------------
 Procedure RevCalc_T2(var M,T2,T5:extended);
 Var T5b,x                          :extended;
     js                             :shortint;
     i,j,i1,i0                      :word;
 Begin --------------------------RevCalc_T2---------------------------------
   Window(3,e+7,35,23); TextBackGround(1); TextColor(14);  ClrScr;
   Write(' Температура    T5 =  ');
   Read(T5);
   If T5<400 then Exit;
   Window(17,e+8,23,e+8); TextBackGround(4); TextColor(14);  ClrScr;
   GotoXY(4,1);
   M:=1.5; x:=1; T5b:=100; js:=0; i0:=0;
   While Abs(T5-T5b)>1E-6 do begin
        GetTime(i1, i1, i, i1);
        If i-j<>0 then begin
          Write(i0);
          If i0<10 then GotoXY(WhereX-1,WhereY)
                   else GotoXY(WhereX-2,WhereY);
          Inc(i0)            end;
        j:=i;
        IterPress(IntLin,Sign,u1,T2,r2r1,P2P1);
        IterTemp5(IntLin,Sign,a2,u1,T5b,r5r1,P5P1,Mr);
        If js*Sign(T5,T5b)<=0 then x:=0.1*x;
        js:=Sign(T5,T5b);
        If T5<T5b then M:=M-x
                  else M:=M+x
                              end; { While
   Window(3,e+8,35,23);
   TextColor(11); TextBackGround(1); ClrScr
 end; ---------------------------RevCalc_T2---------------------------------}
{----------------------------------------------------------------------------}
 Procedure ReadFileInp(var F :text; var is :byte; var Subst :TFile);
 Begin {--------------------------ReadFileInp--------------------------------}
   is:=0;
   Repeat
         Inc(is);
         Readln(F, Subst[is].FName);
         if Subst[is].FName<>'END' then
           Readln(F, Subst[is].SPart)
   Until Subst[is].FName='END';
   Dec(is)
 End; {---------------------------ReadFileInp--------------------------------}
{----------------------------------------------------------------------------}
 Procedure InputDate;
 Var F           :text;
     ModeCalcStr :string;
 Begin {--------------------------InputDate----------------------------------}
    Assign(F,'equilib.inp');
    Reset(F);
    Readln(F, StPath);
    Readln(F, T1);
    Readln(F, UnitT1);
    If UnitT1='[C]' then T1:=T1+273.15;
    Readln(F, P1);
    Readln(F, UnitP);
    If UnitP='[torr]' then P1:=P1/760;
    Readln(F, u1);
    Readln(F, UnitU,base1,base2,base3);
{    Readln(F, Alfa);}
    READFILEINP(F, is1, Subst1);
    Readln(F, StFOut);
    Readln(F, StFDat);
    READFILEINP(F, is4, Subst4);
    Readln(F, T4);
    Readln(F, UnitT4);
{    If UnitT4='[C]' then T4:=T4+273.15;}
    If T4=0 then T4:=T1;
    Readln(F, P4); { атм }
    Readln(F, Number);{ Number-номеp эксп.? }
    Readln(F, ModeCalcStr);
    If ModeCalcStr='u1'    then ModeCalc:=1;
    If ModeCalcStr='T2'    then ModeCalc:=2;
    If ModeCalcStr='T5'    then ModeCalc:=3;
    If ModeCalcStr='P4/P1' then ModeCalc:=4;
    Close(F);
 End; {---------------------------InputDate----------------------------------}
{----------------------------------------------------------------------------}
 Procedure WriteFileInp;
 Var F           :text;
     i           :byte;
     ModeCalcStr :string;
 Begin {--------------------------WriteFileInp-------------------------------}
   Assign(F,'equilib.inp');
   Rewrite(F);
   Writeln(F, StPath);
   Writeln(F, T1:6:2);
   Writeln(F, UnitT1);
   Writeln(F, P1:7:2);
   Writeln(F, UnitP);
   Writeln(F, u1:8:3);
   Writeln(F, UnitU,base1,base2,base3);
   Writeln(F, Alfa);
   For i:=1 to is1 do begin
      Writeln(F, Subst1[i].FName);
      Writeln(F, Subst1[i].SPart:10:5)
                      end;
   Writeln(F, 'END');
   Writeln(F, StFOut);
   Writeln(F, StFDat);
   For i:=1 to is4 do begin
      Writeln(F, Subst4[i].FName);
      Writeln(F, Subst4[i].SPart:10:5)
                      end;
   Writeln(F, 'END');
   Writeln(F, T4:6:2);
   Writeln(F, UnitT4);
   Writeln(F, P4:6:2);
   Writeln(F, Number);
   If ModeCalc=1 then ModeCalcStr:='u1';
   If ModeCalc=2 then ModeCalcStr:='T2';
   If ModeCalc=3 then ModeCalcStr:='T5';
   If ModeCalc=4 then ModeCalcStr:='P4/P1';
   Writeln(F, ModeCalcStr);
   Close(F)
 End; {---------------------------WriteFileInp-------------------------------}
{----------------------------------------------------------------------------}
 Procedure ReadFilEql(is :byte; Subst :TFile; var Mix :TArrayRec;
                       var jmax :byte);
 Var i :byte;
     F :text;
 Begin {--------------------------ReadFileql---------------------------------}
   For i:=1 to is do begin
      Mix[i].Part:=Subst[i].SPart;
      Assign(F, StPath+Subst[i].FName);
      Reset(F);
      Readln(F, Mix[i].Name);
      Readln(F, Mix[i].Weight);
      Readln(F);
      jmax:=0;
      Repeat
            Inc(jmax);
            Readln(F, T[jmax], Mix[i].Cp[jmax], Mix[i].Entl[jmax])
      Until SeekEof(F);{               end;}
      Close(F)       end {For}
 End; {--------------------------ReadFileql----------------------------------}
{----------------------------------------------------------------------------}
 Procedure CalcMixDate(var H,Gamma :ABE; Mix :TArrayRec; is :byte;
                       jmax :byte; var WeightMix :extended);
 Var     i,j :byte;
     SumPart :extended;
 Begin {-------------------------CalcMixDate---------------------------------}
   WeightMix:=0;
   SumPart:=0;
   For i:=1 to is do
      SumPart:=SumPart+Mix[i].Part;
   For i:=1 to is do
      WeightMix:=WeightMix+Mix[i].Weight*(Mix[i].Part/SumPart);
   WeightMix:=1E-3*WeightMix;
   For j:=1 to jmax do begin
      H[j]:=0;
      Gamma[j]:=0;
      For i:=1 to is do begin
         H[j]:=H[j]+(Mix[i].Entl[j]-Mix[i].Entl[3])*
               (Mix[i].Part/SumPart);
         Gamma[j]:=Gamma[j]+Mix[i].Cp[j]*(Mix[i].Part/SumPart)
                        end;
      Gamma[j]:=Gamma[j]/(Gamma[j]-R)
                       end
 End; {--------------------------CalcMixDate---------------------------------}
{----------------------------------------------------------------------------}
 Procedure IterPress(IntLin:F3V; Sign:F2V;
                      var T2,r2r1,P2P1:extended);
 Var x,Pa,Pb,H0                 :extended;
     js                         :shortint;
 Begin {------------------------IterPress------------------------------------}
{   u1:=M*SQRT(GMix1*T1*R/WeightMix1);}
   r2r1:=(G1+1)/(G1-1+2/SQR(M));
   Pa:=1; Pb:=2; js:=0; x:=0.1*r2r1;
   While Abs(Pa-Pb)>=1E-9 do begin
         H0:=(0.5*u1*u1*(1-SQR(1/r2r1))+IntLin(Entmix1,T,T1))*1E-3*Weightmix1;
         T2:=IntLin(T,Entmix1,H0);
         Pa:=r2r1*(T2/T1);
         Pb:=1+u1*u1*Weightmix1/T1/R*(1-1/r2r1);
         If js*Sign(Pa,Pb)<=0 then x:=0.1*x;
         js:=Sign(Pa,Pb);
         If Pa<Pb then r2r1:=r2r1+x
                  else r2r1:=r2r1-x
                             end;
   P2P1:=0.5*(Pa+Pb)
  end; {--------------------------IterPress----------------------------------}
{============================================================================}
BEGIN {---------------------------main---------------------------------------}
    INPUTDATE;
    TextBackGround(black);
    ClrScr;
    TextBackGround(white);
    TextColor(red);
    InsLine;
    GotoXY(27,1);
    Writeln('Conditions behind Shock Wave');
    TextColor(yellow);
    TextBackGround(blue);
    GotoXY(1,3);
    Writeln('Mixture:');
    Window(15,3,80,11);
    ClrScr;
    Writeln(' Initial   conditions           Incident SW         Reflected SW');
    Writeln(' ───────────────────────────────────────────────────────────────');
    Writeln(' Temperature T1=       [K]      T2=                 T5=');
    Writeln(' Pressure    P1=       [atm]    P2=                 P5=');
    Writeln(' Pressure    P4=       [atm]                           ');
    Writeln(' Ratio of densities        ro2/ro1=            ro5/ro1=');
    Writeln(' Concentration [cm-3]           n2=                 n5= ');
    Writeln(' Velocity SW    [m/s]           u1=                 u2=');
    Write('                  [M]            M=                 Mr=');
    Window(1,1,80,25);
    GotoXY(15,13);
    Writeln(' Sound velocity [m/s] a1=       a2=           ');
    If FSearch(StFOut,'')='' then
      begin Assign(F,StFOut); Rewrite(F); Close(F) end;
    If FSearch(StFDat,'')='' then
      begin Assign(F,StFDat); Rewrite(F); Close(F) end;
    READFILEQL(is1, Subst1, Mix1, jmax);
    READFILEQL(is4, Subst4, Mix4, jmax);
    CALCMIXDATE(EntMix1, GMix1, Mix1, is1, jmax, WeightMix1);
    CALCMIXDATE(EntMix4, GMix4, Mix4, is4, jmax, WeightMix4);
    G1:=IntLin(GMix1,T,T1);
    a1:=SQRT(G1*T1*R/WeightMix1);
    G4:=IntLin(GMix4,T,T4);
    a4:=SQRT(G4*T4*R/WeightMix4);
    If UnitU='[M]'    then u1:=M*a1;
    If UnitU='[mks]'  then begin u1:=base1/u1*1E3; M:=u1/a1 end;
    If UnitU='[m/s]'  then M:=u1/a1;
    Case ModeCalc of
        1: ; { u1 }
        2: ; { T2 }
        3: ; { T5 }
        4:   { P4/P1 }
    End;
    ITERPRESS(IntLin, Sign, T2, r2r1, P2P1);
    G2:=IntLin(GMix1,T,T2);
    a2:=SQRT(G2*T2*R/WeightMix1);
    x:=a1/a2*(G2+1)/(G1+1)*(M-1/M);
    Mr:=0.5*(x+SQRT(x*x+4));
    GotoXY(31,6);
    Writeln(P1:6:3);
    GotoXY(31,5);
    Writeln(T1:4:0);
    GotoXY(50,5);
    Writeln(T2:5:0);
    GotoXY(51,6);
    Writeln(P2P1*P1:5:3);
    GotoXY(50,8);
    Writeln(r2r1:6:3);
    GotoXY(50,9);
    Writeln(P2P1*P1/(Kb*T2)*1.01325E-1:6);
    GotoXY(50,10);
    Writeln(u1:5:0);
    GotoXY(51,11);
    Writeln(M:5:3);
    GotoXY(41,13);
    Writeln(a1:4:0);
    ITERTEMP5(IntLin, Sign, T5, r5r1, P5P1);
    GotoXY(70,5);
    Writeln(T5:5:0);
    GotoXY(71,6);
    Writeln(P5P1*P1:5:3);
    GotoXY(71,8);
    Writeln(r5r1:6:3);
    GotoXY(70,9);
    Writeln(P5P1*P1/(Kb*T2)*1.01325E-1:6);
    GotoXY(71,10);
    Writeln(a2*Mr:5:0);
    GotoXY(71,11);
    Writeln(Mr:5:3);
    GotoXY(50,13);
    Writeln(a2:4:0);
    For i:=1 to is1 do begin
       GotoXY(1,4+i);
       Writeln(Mix1[i].Name,' = ',Mix1[i].Part:3:0)
                       end;
    CALCVEL(M, P4P1);
    GotoXY(32,7);
    Writeln(P4:5:3);
    GotoXY(15,15);
    Writeln(' Calculating velocity SW [m/s]  u1=',M*a1:5:0,'      ');
    GotoXY(15,16);
    Writeln('                                 M=',M:6:3,'     ');
    Readln;
END.
