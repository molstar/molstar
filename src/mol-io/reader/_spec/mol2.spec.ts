
import { parseMol2 } from '../mol2/parser';

const Mol2String = `@<TRIPOS>MOLECULE
5816
 26 26 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 O           1.7394   -2.1169   -1.0894 O.3     1  LIG1       -0.3859
      2 O          -2.2941    1.0781   -1.7979 O.3     1  LIG1       -0.5033
      3 O          -3.6584    0.5842    0.5722 O.3     1  LIG1       -0.5033
      4 N           2.6359    1.0243    0.7030 N.3     1  LIG1       -0.3162
      5 C           1.6787   -1.1447   -0.0373 C.3     1  LIG1        0.0927
      6 C           0.2684   -0.6866    0.1208 C.ar    1  LIG1       -0.0143
      7 C           2.6376    0.0193   -0.3576 C.3     1  LIG1        0.0258
      8 C          -0.3658   -0.0099   -0.9212 C.ar    1  LIG1       -0.0109
      9 C          -0.4164   -0.9343    1.3105 C.ar    1  LIG1       -0.0524
     10 C          -1.6849    0.4191   -0.7732 C.ar    1  LIG1        0.1586
     11 C          -1.7353   -0.5053    1.4585 C.ar    1  LIG1       -0.0162
     12 C          -2.3696    0.1713    0.4166 C.ar    1  LIG1        0.1582
     13 C           3.5645    2.1013    0.3950 C.3     1  LIG1       -0.0157
     14 H           2.0210   -1.6511    0.8741 H       1  LIG1        0.0656
     15 H           2.3808    0.4742   -1.3225 H       1  LIG1        0.0453
     16 H           3.6478   -0.3931   -0.4831 H       1  LIG1        0.0453
     17 H           0.1501    0.1801   -1.8589 H       1  LIG1        0.0659
     18 H           0.0640   -1.4598    2.1315 H       1  LIG1        0.0622
     19 H           2.9013    0.5888    1.5858 H       1  LIG1        0.1217
     20 H          -2.2571   -0.7050    2.3907 H       1  LIG1        0.0655
     21 H           2.6646   -2.4067   -1.1652 H       1  LIG1        0.2103
     22 H           3.2862    2.6124   -0.5325 H       1  LIG1        0.0388
     23 H           4.5925    1.7346    0.3078 H       1  LIG1        0.0388
     24 H           3.5401    2.8441    1.1985 H       1  LIG1        0.0388
     25 H          -3.2008    1.2997   -1.5231 H       1  LIG1        0.2923
     26 H          -3.9690    0.3259    1.4570 H       1  LIG1        0.2923
@<TRIPOS>BOND
     1     1     5    1
     2     1    21    1
     3     2    10    1
     4     2    25    1
     5     3    12    1
     6     3    26    1
     7     4     7    1
     8     4    13    1
     9     4    19    1
    10     5     6    1
    11     5     7    1
    12     5    14    1
    13     6     8   ar
    14     6     9   ar
    15     7    15    1
    16     7    16    1
    17     8    10   ar
    18     8    17    1
    19     9    11   ar
    20     9    18    1
    21    10    12   ar
    22    11    12   ar
    23    11    20    1
    24    13    22    1
    25    13    23    1
    26    13    24    1`;

const Mol2StringMultiBlocks = `@<TRIPOS>MOLECULE
5816
 26 26 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 O           1.7394   -2.1169   -1.0894 O.3     1  LIG1       -0.3859
      2 O          -2.2941    1.0781   -1.7979 O.3     1  LIG1       -0.5033
      3 O          -3.6584    0.5842    0.5722 O.3     1  LIG1       -0.5033
      4 N           2.6359    1.0243    0.7030 N.3     1  LIG1       -0.3162
      5 C           1.6787   -1.1447   -0.0373 C.3     1  LIG1        0.0927
      6 C           0.2684   -0.6866    0.1208 C.ar    1  LIG1       -0.0143
      7 C           2.6376    0.0193   -0.3576 C.3     1  LIG1        0.0258
      8 C          -0.3658   -0.0099   -0.9212 C.ar    1  LIG1       -0.0109
      9 C          -0.4164   -0.9343    1.3105 C.ar    1  LIG1       -0.0524
     10 C          -1.6849    0.4191   -0.7732 C.ar    1  LIG1        0.1586
     11 C          -1.7353   -0.5053    1.4585 C.ar    1  LIG1       -0.0162
     12 C          -2.3696    0.1713    0.4166 C.ar    1  LIG1        0.1582
     13 C           3.5645    2.1013    0.3950 C.3     1  LIG1       -0.0157
     14 H           2.0210   -1.6511    0.8741 H       1  LIG1        0.0656
     15 H           2.3808    0.4742   -1.3225 H       1  LIG1        0.0453
     16 H           3.6478   -0.3931   -0.4831 H       1  LIG1        0.0453
     17 H           0.1501    0.1801   -1.8589 H       1  LIG1        0.0659
     18 H           0.0640   -1.4598    2.1315 H       1  LIG1        0.0622
     19 H           2.9013    0.5888    1.5858 H       1  LIG1        0.1217
     20 H          -2.2571   -0.7050    2.3907 H       1  LIG1        0.0655
     21 H           2.6646   -2.4067   -1.1652 H       1  LIG1        0.2103
     22 H           3.2862    2.6124   -0.5325 H       1  LIG1        0.0388
     23 H           4.5925    1.7346    0.3078 H       1  LIG1        0.0388
     24 H           3.5401    2.8441    1.1985 H       1  LIG1        0.0388
     25 H          -3.2008    1.2997   -1.5231 H       1  LIG1        0.2923
     26 H          -3.9690    0.3259    1.4570 H       1  LIG1        0.2923
@<TRIPOS>BOND
     1     1     5    1
     2     1    21    1
     3     2    10    1
     4     2    25    1
     5     3    12    1
     6     3    26    1
     7     4     7    1
     8     4    13    1
     9     4    19    1
    10     5     6    1
    11     5     7    1
    12     5    14    1
    13     6     8   ar
    14     6     9   ar
    15     7    15    1
    16     7    16    1
    17     8    10   ar
    18     8    17    1
    19     9    11   ar
    20     9    18    1
    21    10    12   ar
    22    11    12   ar
    23    11    20    1
    24    13    22    1
    25    13    23    1
    26    13    24    1
@<TRIPOS>MOLECULE
5816
 26 26 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 O           1.7394   -2.1169   -1.0894 O.3     1  LIG1       -0.3859
      2 O          -2.2941    1.0781   -1.7979 O.3     1  LIG1       -0.5033
      3 O          -3.6584    0.5842    0.5722 O.3     1  LIG1       -0.5033
      4 N           2.6359    1.0243    0.7030 N.3     1  LIG1       -0.3162
      5 C           1.6787   -1.1447   -0.0373 C.3     1  LIG1        0.0927
      6 C           0.2684   -0.6866    0.1208 C.ar    1  LIG1       -0.0143
      7 C           2.6376    0.0193   -0.3576 C.3     1  LIG1        0.0258
      8 C          -0.3658   -0.0099   -0.9212 C.ar    1  LIG1       -0.0109
      9 C          -0.4164   -0.9343    1.3105 C.ar    1  LIG1       -0.0524
     10 C          -1.6849    0.4191   -0.7732 C.ar    1  LIG1        0.1586
     11 C          -1.7353   -0.5053    1.4585 C.ar    1  LIG1       -0.0162
     12 C          -2.3696    0.1713    0.4166 C.ar    1  LIG1        0.1582
     13 C           3.5645    2.1013    0.3950 C.3     1  LIG1       -0.0157
     14 H           2.0210   -1.6511    0.8741 H       1  LIG1        0.0656
     15 H           2.3808    0.4742   -1.3225 H       1  LIG1        0.0453
     16 H           3.6478   -0.3931   -0.4831 H       1  LIG1        0.0453
     17 H           0.1501    0.1801   -1.8589 H       1  LIG1        0.0659
     18 H           0.0640   -1.4598    2.1315 H       1  LIG1        0.0622
     19 H           2.9013    0.5888    1.5858 H       1  LIG1        0.1217
     20 H          -2.2571   -0.7050    2.3907 H       1  LIG1        0.0655
     21 H           2.6646   -2.4067   -1.1652 H       1  LIG1        0.2103
     22 H           3.2862    2.6124   -0.5325 H       1  LIG1        0.0388
     23 H           4.5925    1.7346    0.3078 H       1  LIG1        0.0388
     24 H           3.5401    2.8441    1.1985 H       1  LIG1        0.0388
     25 H          -3.2008    1.2997   -1.5231 H       1  LIG1        0.2923
     26 H          -3.9690    0.3259    1.4570 H       1  LIG1        0.2923
@<TRIPOS>BOND
     1     1     5    1
     2     1    21    1
     3     2    10    1
     4     2    25    1
     5     3    12    1
     6     3    26    1
     7     4     7    1
     8     4    13    1
     9     4    19    1
    10     5     6    1
    11     5     7    1
    12     5    14    1
    13     6     8   ar
    14     6     9   ar
    15     7    15    1
    16     7    16    1
    17     8    10   ar
    18     8    17    1
    19     9    11   ar
    20     9    18    1
    21    10    12   ar
    22    11    12   ar
    23    11    20    1
    24    13    22    1
    25    13    23    1
    26    13    24    1`;

const Mol2StringMinimal = `@<TRIPOS>MOLECULE
5816
 26 26 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 O           1.7394   -2.1169   -1.0894 O.3
      2 O          -2.2941    1.0781   -1.7979 O.3
      3 O          -3.6584    0.5842    0.5722 O.3
      4 N           2.6359    1.0243    0.7030 N.3
      5 C           1.6787   -1.1447   -0.0373 C.3
      6 C           0.2684   -0.6866    0.1208 C.ar
      7 C           2.6376    0.0193   -0.3576 C.3
      8 C          -0.3658   -0.0099   -0.9212 C.ar
      9 C          -0.4164   -0.9343    1.3105 C.ar
     10 C          -1.6849    0.4191   -0.7732 C.ar
     11 C          -1.7353   -0.5053    1.4585 C.ar
     12 C          -2.3696    0.1713    0.4166 C.ar
     13 C           3.5645    2.1013    0.3950 C.3
     14 H           2.0210   -1.6511    0.8741 H
     15 H           2.3808    0.4742   -1.3225 H
     16 H           3.6478   -0.3931   -0.4831 H
     17 H           0.1501    0.1801   -1.8589 H
     18 H           0.0640   -1.4598    2.1315 H
     19 H           2.9013    0.5888    1.5858 H
     20 H          -2.2571   -0.7050    2.3907 H
     21 H           2.6646   -2.4067   -1.1652 H
     22 H           3.2862    2.6124   -0.5325 H
     23 H           4.5925    1.7346    0.3078 H
     24 H           3.5401    2.8441    1.1985 H
     25 H          -3.2008    1.2997   -1.5231 H
     26 H          -3.9690    0.3259    1.4570 H
@<TRIPOS>BOND
     1     1     5    1
     2     1    21    1
     3     2    10    1
     4     2    25    1
     5     3    12    1
     6     3    26    1
     7     4     7    1
     8     4    13    1
     9     4    19    1
    10     5     6    1
    11     5     7    1
    12     5    14    1
    13     6     8   ar
    14     6     9   ar
    15     7    15    1
    16     7    16    1
    17     8    10   ar
    18     8    17    1
    19     9    11   ar
    20     9    18    1
    21    10    12   ar
    22    11    12   ar
    23    11    20    1
    24    13    22    1
    25    13    23    1
    26    13    24    1`;

const Mol2StringCrysin = `@<TRIPOS>MOLECULE
1144204
    12    11     2     0     0
SMALL
USER_CHARGES
****
Generated from the CSD

@<TRIPOS>ATOM
     1 Cl1      0.0925   3.6184   1.9845   Cl        1 RES1  -1.0000
     2 C1      -4.7391   0.3350   0.4215   C.ar      2 RES2   0.0000
     3 C2      -3.4121   0.2604   0.9351   C.ar      2 RES2   0.0000
     4 C3      -2.9169   1.2555   1.7726   C.ar      2 RES2   0.0000
     5 C4      -3.7118   2.3440   2.1099   C.ar      2 RES2   0.0000
     6 C5      -5.0314   2.4052   1.6209   C.ar      2 RES2   0.0000
     7 C6      -5.5372   1.4057   0.7962   C.ar      2 RES2   0.0000
     8 C7      -6.9925   1.4547   0.3334   C.3       2 RES2   0.0000
     9 C8      -7.8537   0.5554   1.1859   C.3       2 RES2   0.0000
    10 N1      -9.3089   0.7134   0.8192   N.3       2 RES2   1.0000
    11 O1      -2.6613  -0.8147   0.5707   O.3       2 RES2   0.0000
    12 O2      -1.6204   1.0919   2.2584   O.3       2 RES2   0.0000
@<TRIPOS>BOND
     1     2     3   ar
     2     3     4   ar
     3     4     5   ar
     4     5     6   ar
     5     6     7   ar
     6     7     2   ar
     7     8     7    1
     8     9     8    1
     9    10     9    1
    10    11     3    1
    11    12     4    1
@<TRIPOS>SUBSTRUCTURE
     1 RES1        1 GROUP             0 ****  ****    0
     2 RES2        2 GROUP             0 ****  ****    0
@<TRIPOS>CRYSIN
   10.5150   11.1300    7.9380   90.0000   90.0000   90.0000    29     5
@<TRIPOS>MOLECULE
   1144204
    12    11     2     0     0
SMALL
USER_CHARGES
****
Generated from the CSD

@<TRIPOS>ATOM
     1 Cl1      0.0925   3.6184   1.9845   Cl        1 RES1  -1.0000
     2 C1      -4.7391   0.3350   0.4215   C.ar      2 RES2   0.0000
     3 C2      -3.4121   0.2604   0.9351   C.ar      2 RES2   0.0000
     4 C3      -2.9169   1.2555   1.7726   C.ar      2 RES2   0.0000
     5 C4      -3.7118   2.3440   2.1099   C.ar      2 RES2   0.0000
     6 C5      -5.0314   2.4052   1.6209   C.ar      2 RES2   0.0000
     7 C6      -5.5372   1.4057   0.7962   C.ar      2 RES2   0.0000
     8 C7      -6.9925   1.4547   0.3334   C.3       2 RES2   0.0000
     9 C8      -7.8537   0.5554   1.1859   C.3       2 RES2   0.0000
    10 N1      -9.3089   0.7134   0.8192   N.3       2 RES2   1.0000
    11 O1      -2.6613  -0.8147   0.5707   O.3       2 RES2   0.0000
    12 O2      -1.6204   1.0919   2.2584   O.3       2 RES2   0.0000
@<TRIPOS>BOND
     1     2     3   ar
     2     3     4   ar
     3     4     5   ar
     4     5     6   ar
     5     6     7   ar
     6     7     2   ar
     7     8     7    1
     8     9     8    1
     9    10     9    1
    10    11     3    1
    11    12     4    1
@<TRIPOS>SUBSTRUCTURE
     1 RES1        1 GROUP             0 ****  ****    0
     2 RES2        2 GROUP             0 ****  ****    0
@<TRIPOS>CRYSIN
   10.5150   11.1300    7.9380   90.0000   90.0000   90.0000    29     5
`;
const Mol2AtomWithStatusBit = `@<TRIPOS>MOLECULE
LIG
   46    49     1     0     0
SMALL
No Charge or Current Charge


@<TRIPOS>ATOM
      1 O           14.1460   158.8920    36.9950 O          0 LIG      -0.603100  
      2 C           13.1110   158.9860    37.6350 C          0 LIG       0.663700   BACKBONE|DICT|DIRECT
      3 N           12.4620   160.2570    37.7450 N          0 LIG      -0.561900   
      4 C1          13.0120   161.4260    37.0760 C          0 LIG       0.076300   
      5 C2          12.5590   157.7730    38.3130 C          0 LIG      -0.124600   
      6 C3          12.2660   157.8050    39.6760 C          0 LIG      -0.072000   
      7 C4          11.7670   156.6730    40.3080 C          0 LIG      -0.141000   
      8 C5          11.5580   155.5010    39.5740 C          0 LIG      -0.101000   
      9 C6          12.3520   156.5920    37.5790 C          0 LIG       0.003900   
     10 C7          11.8560   155.4710    38.2080 C          0 LIG      -0.142000   
     11 S           12.7010   156.5180    35.8190 S          0 LIG      -0.113800   
     12 C8          11.4490   155.4060    35.1740 C          0 LIG      -0.047100   
     13 C9          10.0940   155.6050    35.5100 C          0 LIG      -0.123000   
     14 C10          9.1060   154.7310    35.0210 C          0 LIG      -0.081000   
     15 C11         11.8290   154.3550    34.3550 C          0 LIG      -0.031000   
     16 C12         10.8200   153.4660    33.8530 C          0 LIG      -0.096000   
     17 N1          10.9130   152.3630    33.0240 N          0 LIG      -0.249400   
     18 N2           9.6460   151.8660    32.8370 N          0 LIG       0.063800   
     19 C13          8.7500   152.6100    33.5180 C          0 LIG      -0.082100   
     20 C14          9.4900   153.6500    34.1780 C          0 LIG      -0.099000   
     21 C15          7.2670   152.3500    33.5450 C          0 LIG      -0.076000   
     22 C16          6.3800   153.1480    34.2220 C          0 LIG      -0.215800   
     23 C17          4.9250   152.7650    34.2910 C          0 LIG       0.449100   
     24 N3           4.4220   151.8900    33.3980 N          0 LIG      -0.668000   
     25 C18          3.1180   151.5270    33.4310 C          0 LIG       0.401200   
     26 C19          2.2730   152.0660    34.3980 C          0 LIG      -0.242300   
     27 C20          2.7710   152.9720    35.3270 C          0 LIG      -0.092000   
     28 C21          4.1180   153.3270    35.2680 C          0 LIG      -0.233300   
     29 H           13.9080   161.5960    37.4060 H          0 LIG       0.044033   
     30 H1          13.0450   161.2660    36.1200 H          0 LIG       0.044033   
     31 H2          12.4500   162.1960    37.2560 H          0 LIG       0.044033   
     32 H3          11.7420   160.3300    38.2100 H          0 LIG       0.320500   
     33 H4          12.4110   158.6190    40.1840 H          0 LIG       0.156000   
     34 H5          11.5640   156.6950    41.2560 H          0 LIG       0.139000   
     35 H6          11.2070   154.7100    40.0120 H          0 LIG       0.136000   
     36 H7          11.7130   154.6570    37.7010 H          0 LIG       0.146000   
     37 H8           6.6900   153.9580    34.6560 H          0 LIG       0.134000   
     38 H9           6.9230   151.5830    33.0610 H          0 LIG       0.152000   
     39 H10          9.4520   151.1800    32.3560 H          0 LIG       0.328700   
     40 H11         12.7620   154.2220    34.1250 H          0 LIG       0.161000   
     41 H12          8.1740   154.8630    35.2550 H          0 LIG       0.137000   
     42 H13          9.8420   156.3480    36.0810 H          0 LIG       0.141000   
     43 H14          2.1910   153.3510    36.0060 H          0 LIG       0.141000   
     44 H15          1.3380   151.8090    34.4230 H          0 LIG       0.146000   
     45 H16          2.7730   150.8930    32.7830 H          0 LIG       0.024100   
     46 H17          4.4860   153.9600    35.9040 H          0 LIG       0.143000   
@<TRIPOS>BOND
     1     1     2 2   
     2     2     3 1   
     3     2     5 1   
     4     3     4 1   
     5     3    32 1   
     6     4    29 1   
     7     4    30 1   
     8     4    31 1   
     9     5     6 2  
    10     5     9 1  
    11     6     7 1  
    12     6    33 1   
    13     7     8 2  
    14     7    34 1   
    15     8    10 1  
    16     8    35 1   
    17     9    10 2  
    18     9    11 1   
    19    10    36 1   
    20    11    12 1   
    21    12    13 1  
    22    12    15 2
    23    13    14 2  
    24    13    42 1   
    25    14    20 1  
    26    14    41 1   
    27    15    16 1  
    28    15    40 1   
    29    16    17 2   
    30    16    20 1  
    31    17    18 1   
    32    18    19 1   
    33    18    39 1   
    34    19    20 2   
    35    19    21 1   
    36    21    22 2   
    37    21    38 1   
    38    22    23 1   
    39    22    37 1   
    40    23    24 1  
    41    23    28 2  
    42    24    25 2  
    43    25    26 1  
    44    25    45 1   
    45    26    27 2  
    46    26    44 1   
    47    27    28 1  
    48    27    43 1   
    49    28    46 1   
@<TRIPOS>SUBSTRUCTURE
     1 LIG         1 TEMP              0 ****  ****    0 ROOT
`;

const Mol2BondWithStatusBit = `@<TRIPOS>MOLECULE
CHEMBL3264998.sdf
63   67    1
SMALL
USER_CHARGES


@<TRIPOS>ATOM
      1 O1        -28.2626   19.0003    5.7733 O.3       1 UNK        -0.3497
      2 O2        -29.5974   17.6450    3.7070 O.3       1 UNK        -0.3546
      3 O3        -21.8499   21.1925   11.8327 O.2       1 UNK        -0.5085
      4 O4        -22.2187   26.1390   17.0059 O.co2     1 UNK        -0.7591
      5 O5        -21.4758   25.2834   18.9072 O.co2     1 UNK        -0.7591
      6 N6        -21.1110   17.6090    6.6310 N.pl3     1 UNK         0.4919
      7 N7        -23.8100   18.2190    6.3820 N.2       1 UNK        -0.3081
      8 N8        -24.1350   16.7140    4.6410 N.pl3     1 UNK        -0.6273
      9 N9        -21.2130   16.1100    4.9760 N.2       1 UNK        -0.4516
     10 N10       -23.8663   22.0782   12.4561 N.am      1 UNK        -0.3464
     11 C11       -21.8890   16.9780    5.6880 C.ar      1 UNK         0.1777
     12 C12       -22.9770   18.8390    7.3180 C.ar      1 UNK         0.1139
     13 C13       -23.2980   17.3130    5.5550 C.ar      1 UNK         0.2664
     14 C14       -23.5990   19.8300    8.2100 C.ar      1 UNK        -0.0403
     15 C15       -21.6650   18.5530    7.4790 C.ar      1 UNK        -0.4461
     16 C16       -22.9920   20.1240    9.4510 C.ar      1 UNK        -0.0690
     17 C17       -19.8460   17.0980    6.5270 C.ar      1 UNK        -0.5699
     18 C18       -23.6360   21.0620   10.2600 C.ar      1 UNK        -0.1521
     19 C19       -25.4970   16.9560    4.4550 C.ar      1 UNK         0.2316
     20 C20       -24.8050   20.4470    7.8230 C.ar      1 UNK        -0.0994
     21 C21       -19.9390   16.1850    5.4940 C.ar      1 UNK         0.0647
     22 C22       -24.8320   21.6830    9.9010 C.ar      1 UNK        -0.0876
     23 C23       -25.4380   21.3760    8.6720 C.ar      1 UNK        -0.1418
     24 C24       -23.0307   21.4535   11.5935 C.2       1 UNK         0.5396
     25 C25       -26.2230   17.8610    5.2300 C.ar      1 UNK        -0.2280
     26 C26       -26.1560   16.2880    3.4260 C.ar      1 UNK        -0.1397
     27 C27       -23.6286   22.5325   13.7671 C.ar      1 UNK        -0.0415
     28 C28       -27.5820   18.0940    4.9860 C.ar      1 UNK         0.1575
     29 C29       -23.3495   23.5327   16.3634 C.ar      1 UNK        -0.0012
     30 C30       -28.2580   17.4240    3.9550 C.ar      1 UNK         0.1058
     31 C31       -27.5070   16.5190    3.1870 C.ar      1 UNK        -0.1481
     32 C32       -24.5744   23.4206   14.2849 C.ar      1 UNK        -0.1252
     33 C33       -22.5445   22.1320   14.5645 C.ar      1 UNK        -0.1252
     34 C34       -24.4354   23.9214   15.5761 C.ar      1 UNK        -0.1256
     35 C35       -22.4087   22.6336   15.8562 C.ar      1 UNK        -0.1256
     36 C36       -23.2024   24.0723   17.7882 C.3       1 UNK        -0.0609
     37 C37       -22.2122   25.2558   17.8954 C.2       1 UNK         0.6125
     38 C38       -27.5050   19.6025    6.7571 C.3       1 UNK         0.1379
     39 C39       -30.1253   16.9147    2.6618 C.3       1 UNK         0.1389
     40 H40       -21.0096   19.0039    8.2097 H         1 UNK         0.2436
     41 H41       -23.6625   16.0546    4.0394 H         1 UNK         0.3281
     42 H42       -22.0806   19.6089    9.7198 H         1 UNK         0.1477
     43 H43       -19.0305   17.4102    7.1624 H         1 UNK         0.2701
     44 H44       -25.2551   20.2179    6.8677 H         1 UNK         0.1309
     45 H45       -19.1399   15.5750    5.0983 H         1 UNK         0.1637
     46 H46       -25.3072   22.4040   10.5504 H         1 UNK         0.1428
     47 H47       -26.3672   21.8362    8.3698 H         1 UNK         0.1332
     48 H48       -25.7364   18.3777    6.0394 H         1 UNK         0.1344
     49 H49       -25.6216   15.5882    2.7995 H         1 UNK         0.1646
     50 H50       -24.7927   22.2801   12.1104 H         1 UNK         0.2846
     51 H51       -27.9781   15.9788    2.3801 H         1 UNK         0.1434
     52 H52       -25.4203   23.7337   13.6893 H         1 UNK         0.1175
     53 H53       -21.8054   21.4274   14.2187 H         1 UNK         0.1175
     54 H54       -25.1696   24.6135   15.9622 H         1 UNK         0.1418
     55 H55       -21.5717   22.3223   16.4642 H         1 UNK         0.1418
     56 H56       -24.1746   24.3788   18.1778 H         1 UNK         0.0601
     57 H57       -22.8929   23.2438   18.4298 H         1 UNK         0.0601
     58 H58       -28.1279   20.2966    7.3214 H         1 UNK         0.0380
     59 H59       -26.6792   20.1461    6.2982 H         1 UNK         0.0380
     60 H60       -27.1097   18.8405    7.4288 H         1 UNK         0.0380
     61 H61       -31.1815   17.1576    2.5458 H         1 UNK         0.0378
     62 H62       -30.0181   15.8504    2.8716 H         1 UNK         0.0378
     63 H63       -29.5931   17.1591    1.7425 H         1 UNK         0.0378
@<TRIPOS>BOND
     1    1   28 1
     2    1   38 1
     3    2   30 1
     4    2   39 1
     5    3   24 2
     6    4   37 1
     7    5   37 2
     8    6   11 ar
     9    6   15 ar
    10    6   17 ar
    11    7   12 ar
    12    7   13 ar
    13    8   13 1
    14    8   19 1
    15    8   41 1
    16    9   11 ar
    17    9   21 ar
    18   10   24 am   BACKBONE|DICT|INTERRES
    19   10   27 1
    20   10   50 1
    21   11   13 ar
    22   12   14 1
    23   12   15 ar
    24   14   16 ar
    25   14   20 ar
    26   15   40 1
    27   16   18 ar
    28   16   42 1
    29   17   21 ar
    30   17   43 1
    31   18   22 ar
    32   18   24 1
    33   19   25 ar
    34   19   26 ar
    35   20   23 ar
    36   20   44 1
    37   21   45 1
    38   22   23 ar
    39   22   46 1
    40   23   47 1
    41   25   28 ar
    42   25   48 1
    43   26   31 ar
    44   26   49 1
    45   27   32 ar
    46   27   33 ar
    47   28   30 ar
    48   29   34 ar
    49   29   35 ar
    50   29   36 1
    51   30   31 ar
    52   31   51 1
    53   32   34 ar
    54   32   52 1
    55   33   35 ar
    56   33   53 1
    57   34   54 1
    58   35   55 1
    59   36   37 1
    60   36   56 1
    61   36   57 1
    62   38   58 1
    63   38   59 1
    64   38   60 1
    65   39   61 1
    66   39   62 1
    67   39   63 1
@<TRIPOS>SUBSTRUCTURE
     1 UNK         1 GROUP             0       ****    0 ROOT
`;

describe('mol2 reader', () => {
    it('basic', async () => {
        const parsed = await parseMol2(Mol2String, '').run();
        if (parsed.isError) {
            throw new Error(parsed.message);
        }
        const mol2File = parsed.result;

        // number of structures
        expect(mol2File.structures.length).toBe(1);

        const data = mol2File.structures[0];
        const { molecule, atoms, bonds } = data;

        // molecule fields
        expect(molecule.mol_name).toBe('5816');
        expect(molecule.num_atoms).toBe(26);
        expect(molecule.num_bonds).toBe(26);
        expect(molecule.num_subst).toBe(0);
        expect(molecule.num_feat).toBe(0);
        expect(molecule.num_sets).toBe(0);
        expect(molecule.mol_type).toBe('SMALL');
        expect(molecule.charge_type).toBe('GASTEIGER');
        expect(molecule.status_bits).toBe('');
        expect(molecule.mol_comment).toBe('');

        // required atom fields
        expect(atoms.count).toBe(26);
        expect(atoms.atom_id.value(0)).toBe(1);
        expect(atoms.atom_name.value(0)).toBe('O');
        expect(atoms.x.value(0)).toBeCloseTo(1.7394, 0.001);
        expect(atoms.y.value(0)).toBeCloseTo(-2.1169, 0.0001);
        expect(atoms.z.value(0)).toBeCloseTo(-1.0893, 0.0001);
        expect(atoms.atom_type.value(0)).toBe('O.3');

        // optional atom fields
        expect(atoms.subst_id.value(0)).toBe(1);
        expect(atoms.subst_name.value(0)).toBe('LIG1');
        expect(atoms.charge.value(0)).toBeCloseTo(-0.3859);
        expect(atoms.status_bit.value(0)).toBe('');

        // required bond fields
        expect(bonds.count).toBe(26);
        expect(bonds.bond_id.value(0)).toBe(1);
        expect(bonds.origin_atom_id.value(0)).toBe(1);
        expect(bonds.target_atom_id.value(0)).toBe(5);
        expect(bonds.bond_type.value(0)).toBe('1');

        // optional bond fields
        expect(bonds.status_bits.value(0)).toBe('');
    });

    it('multiblocks', async () => {
        const parsed = await parseMol2(Mol2StringMultiBlocks, '').run();
        if (parsed.isError) {
            throw new Error(parsed.message);
        }
        const mol2File = parsed.result;

        // number of structures
        expect(mol2File.structures.length).toBe(2);

        const data = mol2File.structures[1];
        const { molecule, atoms, bonds } = data;

        // molecule fields
        expect(molecule.mol_name).toBe('5816');
        expect(molecule.num_atoms).toBe(26);
        expect(molecule.num_bonds).toBe(26);
        expect(molecule.num_subst).toBe(0);
        expect(molecule.num_feat).toBe(0);
        expect(molecule.num_sets).toBe(0);
        expect(molecule.mol_type).toBe('SMALL');
        expect(molecule.charge_type).toBe('GASTEIGER');
        expect(molecule.status_bits).toBe('');
        expect(molecule.mol_comment).toBe('');

        // required atom fields
        expect(atoms.count).toBe(26);
        expect(atoms.atom_id.value(0)).toBe(1);
        expect(atoms.atom_name.value(0)).toBe('O');
        expect(atoms.x.value(0)).toBeCloseTo(1.7394, 0.001);
        expect(atoms.y.value(0)).toBeCloseTo(-2.1169, 0.0001);
        expect(atoms.z.value(0)).toBeCloseTo(-1.0893, 0.0001);
        expect(atoms.atom_type.value(0)).toBe('O.3');

        // optional atom fields
        expect(atoms.subst_id.value(0)).toBe(1);
        expect(atoms.subst_name.value(0)).toBe('LIG1');
        expect(atoms.charge.value(0)).toBeCloseTo(-0.3859);
        expect(atoms.status_bit.value(0)).toBe('');

        // required bond fields
        expect(bonds.count).toBe(26);
        expect(bonds.bond_id.value(0)).toBe(1);
        expect(bonds.origin_atom_id.value(0)).toBe(1);
        expect(bonds.target_atom_id.value(0)).toBe(5);
        expect(bonds.bond_type.value(0)).toBe('1');

        // optional bond fields
        expect(bonds.status_bits.value(0)).toBe('');
    });

    it('minimal', async () => {
        const parsed = await parseMol2(Mol2StringMinimal, '').run();
        if (parsed.isError) {
            throw new Error(parsed.message);
        }
        const mol2File = parsed.result;

        // number of structures
        expect(mol2File.structures.length).toBe(1);

        const data = mol2File.structures[0];
        const { molecule, atoms, bonds } = data;

        // molecule fields
        expect(molecule.mol_name).toBe('5816');
        expect(molecule.num_atoms).toBe(26);
        expect(molecule.num_bonds).toBe(26);
        expect(molecule.num_subst).toBe(0);
        expect(molecule.num_feat).toBe(0);
        expect(molecule.num_sets).toBe(0);
        expect(molecule.mol_type).toBe('SMALL');
        expect(molecule.charge_type).toBe('GASTEIGER');
        expect(molecule.status_bits).toBe('');
        expect(molecule.mol_comment).toBe('');

        // required atom fields
        expect(atoms.count).toBe(26);
        expect(atoms.atom_id.value(0)).toBe(1);
        expect(atoms.atom_name.value(0)).toBe('O');
        expect(atoms.x.value(0)).toBeCloseTo(1.7394, 0.001);
        expect(atoms.y.value(0)).toBeCloseTo(-2.1169, 0.0001);
        expect(atoms.z.value(0)).toBeCloseTo(-1.0893, 0.0001);
        expect(atoms.atom_type.value(0)).toBe('O.3');

        // optional atom fields
        expect(atoms.subst_id.value(0)).toBe(0);
        expect(atoms.subst_name.value(0)).toBe('');
        expect(atoms.charge.value(0)).toBeCloseTo(0);
        expect(atoms.status_bit.value(0)).toBe('');

        // required bond fields
        expect(bonds.count).toBe(26);
        expect(bonds.bond_id.value(0)).toBe(1);
        expect(bonds.origin_atom_id.value(0)).toBe(1);
        expect(bonds.target_atom_id.value(0)).toBe(5);
        expect(bonds.bond_type.value(0)).toBe('1');

        // optional bond fields
        expect(bonds.status_bits.value(0)).toBe('');
    });

    it('atom status_bit', async () => {
        const parsed = await parseMol2(Mol2AtomWithStatusBit, '').run();
        if (parsed.isError) {
            throw new Error(parsed.message);
        }
        const mol2File = parsed.result;
        const data = mol2File.structures[0];
        const statusBits = data.atoms.status_bit.toArray();
        expect(statusBits.length).toEqual(data.atoms.count);
        expect(statusBits[1]).toEqual('BACKBONE|DICT|DIRECT');
        for (let i = 0; i < data.atoms.count; i++) {
            if (i !== 1) expect(statusBits[i]).toEqual('');
        }
    });

    it('bond status_bit', async () => {
        const parsed = await parseMol2(Mol2BondWithStatusBit, '').run();
        if (parsed.isError) {
            throw new Error(parsed.message);
        }
        const mol2File = parsed.result;
        const data = mol2File.structures[0];
        const statusBits = data.bonds.status_bits.toArray();
        expect(statusBits.length).toEqual(data.bonds.count);
        expect(statusBits[17]).toEqual('BACKBONE|DICT|INTERRES');
        for (let i = 0; i < data.bonds.count; i++) {
            if (i !== 17) expect(statusBits[i]).toEqual('');
        }
    });

    it('crysin', async () => {
        const parsed = await parseMol2(Mol2StringCrysin, '').run();
        if (parsed.isError) {
            throw new Error(parsed.message);
        }
        const mol2File = parsed.result;

        // number of structures
        expect(mol2File.structures.length).toBe(2);

        // crysin fields
        for (const data of mol2File.structures) {
            expect(data.crysin).toEqual({
                a: 10.5150,
                b: 11.1300,
                c: 7.9380,
                alpha: 90.0,
                beta: 90.0,
                gamma: 90.0,
                spaceGroup: 29,
                setting: 5
            });
        }
    });
});
