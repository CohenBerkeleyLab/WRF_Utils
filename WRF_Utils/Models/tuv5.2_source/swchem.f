      SUBROUTINE swchem(nw,wl,nz,tlev,airden,
     $     j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Load various "weighting functions" (products of cross section and        =*
*=  quantum yield at each altitude and each wavelength).  The altitude       =*
*=  dependence is necessary to ensure the consideration of pressure and      =*
*=  temperature dependence of the cross sections or quantum yields.          =*
*=  The actual reading, evaluation and interpolation is done in separate     =*
*=  subroutines for ease of management and manipulation.  Please refer to    =*
*=  the inline documentation of the specific subroutines for detail          =*
*=  information.                                                             =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section * quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'params'

* input

      INTEGER nw
      REAL wl(kw)
      
      INTEGER nz
      REAL tlev(kz)
      REAL airden(kz)

* weighting functions

      CHARACTER*50 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* local:
      REAL wc(kw)
      INTEGER iw
*_______________________________________________________________________

* complete wavelength grid

      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
 5    CONTINUE

*____________________________________________________________________________


******** Ox Photochemistry

* A1.  O2 + hv -> O + O
* reserve first position.  Cross section parameterization in Schumman-Runge and 
* Lyman-alpha regions are zenith-angle dependent, will be written in 
* subroutine seto2.f.
 
      j = 1
      jlabel(j) = 'O2 -> O + O'

*A2.  O3 + hv ->  (both channels)
      CALL r01(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

******** HOx Photochemistry

*B1. HO2 + hv -> OH + O
      CALL r39(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*B3. H2O2 + hv -> 2 OH
      CALL r08(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

******** NOx Photochemistry

*C1.  NO2 + hv -> NO + O(3P)
      CALL r02(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*C2.  NO3 + hv ->  (both channels)
      CALL r03(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)
     
*C3.  N2O + hv -> N2 + O(1D)
      CALL r44(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*C5.  N2O5 + hv -> (both channels)
      CALL r04(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*C6.  HNO2 + hv -> OH + NO
      CALL r05(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*C7.  HNO3 + hv -> OH + NO2
      CALL r06(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*C8.  HNO4 + hv -> HO2 + NO2
      CALL r07(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* NO3-(aq) + hv -> NO2 + O-     (for snow)
* NO3-(aq) + hv -> NO2- + O(3P) (for snow)
      CALL r118(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

******** Organic Photochemistry

*D1.  CH2O + hv -> (both channels)
c      CALL r10(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)
      CALL pxCH2O(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*D2.  CH3CHO + hv -> (all three channels)
      CALL r11(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*D3.  C2H5CHO + hv -> C2H5 + HCO
      CALL r12(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*D9.  CH3(OOH) + hv -> CH3O + OH
      CALL r16(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*D10. HOCH2OOH -> HOCH2O. + OH
      CALL r121(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*D12. CH3(ONO2) + hv -> CH3O + NO2
      CALL r17(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*D13. CH3(OONO2) -> CH3(OO) + NO2
      call r134(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*     CH3CH2(ONO2) -> CH3CH2O + NO2
      CALL r106(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*     C2H5(ONO2) -> C2H5O + NO2
      call r141(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*     n-C3H7ONO2 -> n-C3H7O + NO2
      call r142(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*     1-C4H9ONO2 -> 1-C4H9O + NO2
      call r143(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*     2-C4H9ONO2 -> 2-C4H9O + NO2
      call r144(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*     CH3CH(ONO2)CH3 -> CH3CHOCH3 + NO2
      CALL r107(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*     CH2(OH)CH2(ONO2) -> CH2(OH)CH2(O.) + NO2
      CALL r108(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*     CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2
      CALL r109(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*     C(CH3)3(ONO2) -> C(CH3)3(O.) + NO2
      CALL r110(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*     C(CH3)3(ONO) -> C(CH3)3(O) + NO
      call r135(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*D14. PAN + hv -> CH3CO(OO) + NO2
*     PAN + hv -> CH3CO(O) + NO3
      CALL r18(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*D15. CH3CH2COO2NO2 -> CH3CH2CO(OO) + NO2
*     CH3CH2COO2NO2 -> CH2CH2CO(O) + NO3
      CALL r120(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*D16. CH2=CHCHO + hv -> Products
      CALL r122(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*D17. CH2=C(CH3)CHO + hv -> Products
      CALL r104(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*D18. CH3COCH=CH2 + hv -> Products
      CALL r103(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*D19. CH2(OH)CHO + hv -> Products
      CALL r101(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*D20. CH3COCH3 + hv -> Products
      CALL r15(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*     CH3COCH2CH3 -> CH3CO + CH2CH3
      CALL r119(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*D21. CH2(OH)COCH3 -> CH3CO + CH2(OH)
*     CH2(OH)COCH3 -> CH2(OH)CO + CH3
      CALL r112(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*D22. CHOCHO + hv -> Products
      CALL r13(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*D23. CH3COCHO + hv -> Products
      CALL r14(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*     CH3COCOCH3 + hv -> Products
      CALL r102(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*D25. CH3CO(OH) + hv -> Products
      CALL r138(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*D26. CH3CO(OOH) + hv -> Products
      CALL r123(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*D28. CH3COCO(OH) + hv -> Products
      CALL r105(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

* (CH3)2NNO -> products
      call r124(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

******** FOx Photochemistry

*E12. CF2O + hv -> Products
      CALL r22(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

******** ClOx Photochemistry

*F1.  Cl2 + hv -> Cl + Cl

      CALL r47(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F2.  ClO -> Cl + O
      call r125(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F3.  ClOO + hv -> Products
      CALL r31(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F4.  OCLO -> Products
      call r132(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F7.  ClOOCl -> Cl + ClOO
      CALL r111(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F13. HCl -> H + Cl
      CALL r137(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F14. HOCl -> HO + Cl
      call r130(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F15. NOCl -> NO + Cl
      call r131(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F16. ClNO2 -> Cl + NO2
      call r126(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F17. ClONO -> Cl + NO2
      call r136(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)      

*F18. ClONO2 + hv -> Products
      CALL r45(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F19. CCl4 + hv -> Products
      CALL r20(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F20. CH3OCl + hv -> Cl + CH3O
      CALL r139(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F21. CHCl3 -> Products
      CALL r140(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F23. CH3Cl + hv -> Products
      CALL r30(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F24. CH3CCl3 + hv -> Products
      CALL r29(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F30. CCl2O + hv -> Products
      CALL r19(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F32. CClFO + hv -> Products
      CALL r21(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F33. CCl3F (CFC-11) + hv -> Products
      CALL r26(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F34. CCl2F2 (CFC-12) + hv -> Products
      CALL r27(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F36. CF2ClCFCl2 (CFC-113) + hv -> Products
      CALL r23(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F37. CF2ClCF2Cl (CFC-114) + hv -> Products
      CALL r24(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F38. CF3CF2Cl (CFC-115) + hv -> Products
      CALL r25(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F40. CHClF2 (HCFC-22) + hv -> Products
      CALL r38(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F42. CF3CHCl2 (HCFC-123) + hv -> Products
      CALL r32(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F43. CF3CHFCl (HCFC-124) + hv -> Products
      CALL r33(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F45. CH3CFCl2 (HCFC-141b) + hv -> Products
      CALL r34(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F46. CH3CF2Cl (HCFC-142b) + hv -> Products
      CALL r35(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F56. CF3CF2CHCl2 (HCFC-225ca) + hv -> Products
      CALL r36(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*F57. CF2ClCF2CHFCl (HCFC-225cb) + hv -> Products
      CALL r37(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

**** BrOx Photochemistry

*G1.  Br2 -> Br + Br
      CALL r115(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*G3.  BrO -> Br + O
      CALL r114(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*G6.  HOBr -> OH + Br
      CALL r113(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*G7.  BrNO -> Br + NO
      call r127(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*G8.  BrONO -> Br + NO2
*     BrONO -> BrO + NO
      call r129(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*G9.  BrNO2 -> Br + NO2
      call r128(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*G10. BrONO2 + hv -> Products
      CALL r46(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*G11. BrCl -> Br + Cl
      call r133(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*G13. CH3Br + hv -> Products
      CALL r28(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*G15. CHBr3 + hv -> Products
      CALL r09(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*G25. CF2Br2 (Halon-1202) + hv -> Products
      CALL r40(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*G26. CF2BrCl (Halon-1211) + hv -> Products
      CALL r41(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*G27. CF3Br (Halon-1301) + hv -> Products
      CALL r42(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*G35. CF2BrCF2Br (Halon-2402) + hv -> Products
      CALL r43(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

*H24  perfluoro n-iodo propane -> products
      CALL r145(nw,wl,wc,nz,tlev,airden,j,sq,jlabel)

****************************************************************

      IF (j .GT. kj) STOP '1002'
      RETURN
      END
