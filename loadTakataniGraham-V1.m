% loadTakataniGraham.m
% Returs nmTG(), muaoxyTG(), muadeoxyTG()
% Tabulated Molar Extinction Coefficient for Hemoglobin in Water
%These values for the molar extinction coefficient e in [cm-1/(moles/liter)] were compiled by Scott Prahl (prahl@ece.ogi.edu) using data from 
% S. Takatani and M. D. Graham, "Theoretical analysis of diffuse reflectance from a two-layer tissue model," IEEE Trans. Biomed. Eng., BME-26, 656--664, (1987). 
%To convert this data to absorbance A, multiply by the molar concentration and the pathlength. For example, if x is the number of grams per liter and a 1 cm cuvette is being used, then the absorbance is given by 
%        (e) [(1/cm)/(moles/liter)] (x) [g/liter] (1) [cm]
%  A =  ---------------------------------------------------
%                          64,500 [g/mole]
%using 64,500 as the gram molecular weight of hemoglobin. 
%To convert this data to absorption coefficient in (cm-1), multiply by the molar concentration and 2.303, 
% mua = log(10)*extTG*150/64500
%where x is the number of grams per liter. A typical value of x for whole blood is x=150 g Hb/liter. 
%lambda	HbO2	Hb	
%nm	cm-1/M	cm-1/M	
extTG = [
200 50      50
300 50      50
400 50      50
450	68000	58000
460	45040	20600
480	27360	13360
500	20200	16360
507	19240	19240
510	19040	20000
520	23520	25080
522	25680	25680
540	57080	41120
542	57480	44000
549	49840	49840
555	36000	52160
560	33880	50160
569	45080	45080
577	61480	36800
579	54920	35440
586	28920	28920
600	3200	14600
605	1860	9496
615	1152	5776
625	732	4400
635	488	3796
645	396	3436
655	340	3244
665	292	3156
675	288	3028
685	272	2796
695	280	2424
705	300	1988
715	328	1628
725	368	1464
735	412	1464
745	480	1616
755	556	1756
765	616	1640
775	684	1340
785	736	1040
795	776	964
805	880	896
815	880	880
825	952	832
835	996	820
845	1048	820
855	1068	820
865	1116	820
875	1140	848
885	1168	832
895	1188	884
905	1208	896
915	1220	924
925	1228	860
935	1216	848
945	1212	756
955	1196	704
965	1176	616
975	1148	552
985	1108	424
995	1052	372
1100 100    100
1200 50     50
1700 50  50];

nmTG    = extTG(:,1);
muaoxyTG = log(10)*extTG(:,2)*150/64500;
muadeoxyTG = log(10)*extTG(:,3)*150/64500;
clear extTG