/*
 This software is part of iowahills_dsp, a set of DSP routines under MIT License.
 2016 By Daniel Klostermann, Iowa Hills Software, LLC  IowaHills.com
 Copyright (c) 2021  Hayati Ayguen <h_ayguen@web.de>
 All rights reserved.


 This is a C translation of the Parks McClellan algorithm originally done in Fortran.

 The original fortran code came from the Parks McClellan article on Wikipedia.
 http://en.wikipedia.org/wiki/Parks%E2%80%93McClellan_filter_design_algorithm

 This code is quite different from the original. The original code had 69 goto statements,
 which made it nearly impossible to follow. And of course, it was Fortran code, so many changes
 had to be made regardless of style.

 Apparently, Fortran doesn't use global variables. Instead, is uses something called
 common memory space. e.g. COMMON PI2,AD,DEV,X,Y,GRID,DES,WT,ALPHA,IEXT,NFCNS,NGRID
 I simply converted these to globals. It isn't pretty, but for this purpose, who cares?

 The first step was to get a C version of the code working with as few changes as possible.
 That version is also available on:  http://www.iowahills.com/A7ExampleCodePage.html
 Then, in our desire to see if the code could be made more understandable, we decided to
 remove as many goto statements as possible. We checked our work by comparing the coefficients
 between this code and our original translation on more than 1000 filters while varying all the parameters.

 Ultimately, we were able to reduce the goto count from 69 to 7, all of which are in the Remez
 function. Of the 7 remaining, 3 of these are at the very bottom of the function, and go
 back to the very top of the function. These could have been removed, but our goal was to
 clarify the code, not restyle it, and since they are clear, we let them be.

 The other 4 goto statements are intertwined in a rather nasty way. We recommend you print out
 the Remez code, tape the sheets end to end, and trace out the goto's. It wasn't apparent to
 us that they can be removed without an extensive study of the code.

 For better or worse, we also removed any code that was obviously related to Hilbert transforms
 and Differentiators. We did this because we aren't interested in these, and we also don't
 believe this algorithm does a very good job with them (far too much ripple).

 We added the functions CalcCoefficients() and ErrTest() as a way to simplify things a bit.

 We also found 3 sections of code that never executed. Two of the sections were just a few lines
 that the goto's always went around. The third section represented nearly half of the CalcCoefficients()
 function. This statement always tested the same, which never allowed the code to execute.
 if(GRID[1] < 0.01 && GRID[NGRID] > 0.49) KKK = 1;
 This may be due to the 0.01 minimum width limit we set for the bands.

 Note our use of MIN_TEST_VAL. The original code wasn't guarding against division by zero.
 Limiting the return values as we have also helped the algorithm's convergence behavior.

 In an effort to improve readability, we made a large number of variable name changes and also
 deleted a large number of variables. We left many variable names in tact, in part as an aid when
 comparing to the original code, and in part because a better name wasn't obvious.

 This code is essentially straight c, and should compile with few, if any changes. Note the error
 message in CalcParkCoeff2. It warns of the possibility of convergence failure, but you will
 find that the iteration count NITER, isn't always an indicator of convergence problems when
 it is less than 3, as stated in the original Fortran code comments.

*/


#include <iowahills/parks_mcclellan.h>

#include <cmath>

//---------------------------------------------------------------------------

#define ITRMAX 50               // Max Number of Iterations. Some Notch and BPF are running ~ 43
#define MIN_TEST_VAL 1.0E-6     // Min value used in LeGrangeInterp and GEE

//---------------------------------------------------------------------------

#define M_2PI  6.28318530717958647692

static double LeGrangeInterp2(const ParksGlobals *PtrGlobals, int K, int N, int M);
static double GEE2(const ParksGlobals *PtrGlobals, int K, int N);
static int Remez2(ParksGlobals *PtrGlobals, int GridIndex);
static bool ErrTest(const ParksGlobals *PtrGlobals, int k, int Nut, double Comp, double *Err);
static void CalcCoefficients(ParksGlobals *PtrGlobals);
static int CalcParkCoeff2(ParksGlobals *PtrGlobals, int NBANDS, int NFILT, double *FirCoeff);

//---------------------------------------------------------------------------

int ParksMcClellanEx(ParksGlobals *PtrGlobals, double *FirCoeff, int NumTaps, TFIRPassTypes PassType, double OmegaC, double BW, double ParksWidth)
{
 ParksGlobals &glb = *PtrGlobals;
 if(NumTaps > PARKS_MAX_NUM_TAPS)
   return 1;
 int j, NumBands = 0;

 // Note: There is no feedback to the caller if ParksWidth or NumTaps are modified here.
 if(PassType == firBPF || PassType == firNOTCH || NumTaps > 70)
  {
   if(ParksWidth > 0.15)
     ParksWidth = 0.15; // Else we have covergence problems.
  }

 if(PassType == firNOTCH || PassType == firHPF)
  {
   if(NumTaps % 2 == 0)
     NumTaps++;  // High pass and notch filters must have odd tap counts.
  }

 if(NumTaps > PARKS_MAX_NUM_TAPS)
   NumTaps = PARKS_MAX_NUM_TAPS;

 // It helps the algorithm a great deal if each band is at least 0.01 wide.
 // The weights used here came from the orig PM code.
 if(PassType == firLPF)
  {
   NumBands = 2;
   glb.Edge[1] = 0.0;                       // Omega = 0
   glb.Edge[2] = OmegaC;                    // Pass band edge
   if(glb.Edge[2] < 0.01)   glb.Edge[2] = 0.01;
   if(glb.Edge[2] > 0.98)   glb.Edge[2] = 0.98;
   glb.Edge[3] = glb.Edge[2] + ParksWidth;      // Stop band edge
   if(glb.Edge[3] > 0.99)   glb.Edge[3] = 0.99;
   glb.Edge[4] = 1.0;                       // Omega = Pi
   glb.BandMag[1] = 1.0;
   glb.BandMag[2] = 0.0;
   glb.InitWeight[1] = 1.0;
   glb.InitWeight[2] = 10.0;
  }

 if(PassType == firHPF)
  {
   NumBands = 2;
   glb.Edge[1] = 0.0;                       // Omega = 0
   glb.Edge[3] = OmegaC;                    // Pass band edge
   if(glb.Edge[3] > 0.99)   glb.Edge[3] = 0.99;
   if(glb.Edge[3] < 0.02)   glb.Edge[3] = 0.02;
   glb.Edge[2] = glb.Edge[3] - ParksWidth;      // Stop band edge
   if(glb.Edge[2] < 0.01)   glb.Edge[2] = 0.01;
   glb.Edge[4] = 1.0;                       // Omega = Pi
   glb.BandMag[1] = 0.0;
   glb.BandMag[2] = 1.0;
   glb.InitWeight[1] = 10.0;
   glb.InitWeight[2] = 1.0;
  }

 if(PassType == firBPF)
  {
   NumBands = 3;
   glb.Edge[1] = 0.0;                       // Omega = 0
   glb.Edge[3] = OmegaC - BW/2.0;           // Left pass band edge.
   if(glb.Edge[3] < 0.02)   glb.Edge[3] = 0.02;
   glb.Edge[2] = glb.Edge[3] - ParksWidth;      // Left stop band edge
   if(glb.Edge[2] < 0.01)   glb.Edge[2] = 0.01;
   glb.Edge[4] = OmegaC + BW/2.0;           // Right pass band edge
   if(glb.Edge[4] > 0.98)   glb.Edge[4] = 0.98;
   glb.Edge[5] = glb.Edge[4] + ParksWidth;      // Right stop band edge
   if(glb.Edge[5] > 0.99)   glb.Edge[5] = 0.99;
   glb.Edge[6] = 1.0;                       // Omega = Pi

   glb.BandMag[1] = 0.0;
   glb.BandMag[2] = 1.0;
   glb.BandMag[3] = 0.0;
   glb.InitWeight[1] = 10.0;
   glb.InitWeight[2] = 1.0;
   glb.InitWeight[3] = 10.0;
  }

 // This algorithm tends to have problems with narrow band notch filters.
 if(PassType == firNOTCH)
  {
   NumBands = 3;
   glb.Edge[1] = 0.0;                        // Omega = 0
   glb.Edge[3] = OmegaC - BW/2.0;            // Left stop band edge.
   if(glb.Edge[3] < 0.02)   glb.Edge[3] = 0.02;
   glb.Edge[2] = glb.Edge[3] - ParksWidth;       // Left pass band edge
   if(glb.Edge[2] < 0.01)   glb.Edge[2] = 0.01;
   glb.Edge[4] = OmegaC + BW/2.0;            // Right stop band edge
   if(glb.Edge[4] > 0.98)   glb.Edge[4] = 0.98;
   glb.Edge[5] = glb.Edge[4] + ParksWidth;       // Right pass band edge
   if(glb.Edge[5] > 0.99)   glb.Edge[5] = 0.99;
   glb.Edge[6] = 1.0;                        // Omega = Pi

   glb.BandMag[1] = 1.0;
   glb.BandMag[2] = 0.0;
   glb.BandMag[3] = 1.0;
   glb.InitWeight[1] = 1.0;
   glb.InitWeight[2] = 10.0;
   glb.InitWeight[3] = 1.0;
  }

 // Parks McClellan's edges are based on 2Pi, we are based on Pi.
 for(j=1; j<=2*NumBands; j++)
   glb.Edge[j] /= 2.0;

 return CalcParkCoeff2(&glb, NumBands, NumTaps, FirCoeff);
}
//---------------------------------------------------------------------------

int ParksMcClellan(double *FirCoeff, int NumTaps, TFIRPassTypes PassType, double OmegaC, double BW, double ParksWidth)
{
 if(!NumTaps)
   return 1;
 ParksGlobals *PtrGlobals = new ParksGlobals();
 if (!PtrGlobals)
   return 1;
 int ret = ParksMcClellanEx(PtrGlobals, FirCoeff, NumTaps, PassType, OmegaC, BW, ParksWidth);
 delete PtrGlobals;
 return ret;
}

//---------------------------------------------------------------------------

static int CalcParkCoeff2(ParksGlobals *PtrGlobals, int NumBands, int TapCount, double *FirCoeff)
{
 ParksGlobals &glb = *PtrGlobals;
 int j, k, GridCount, GridIndex, BandIndex, NumIterations;
 double LowFreqEdge, UpperFreq, TempVar, Change;
 bool OddNumTaps;
 GridCount = 16;               // Grid Density

 if(TapCount % 2)OddNumTaps = true;
 else OddNumTaps = false;

 glb.HalfTapCount = TapCount/2;
 if(OddNumTaps) glb.HalfTapCount++;

 glb.Grid[1] = glb.Edge[1];
 LowFreqEdge = GridCount * glb.HalfTapCount;
 LowFreqEdge = 0.5 / LowFreqEdge;
 j = 1;
 k = 1;
 BandIndex = 1;
 while(BandIndex <= NumBands)
  {
   UpperFreq = glb.Edge[k+1];
   while(glb.Grid[j] <= UpperFreq)
	{
     TempVar = glb.Grid[j];
     glb.DesiredMag[j] = glb.BandMag[BandIndex];
     glb.Weight[j] = glb.InitWeight[BandIndex];
	 j++;;
     glb.Grid[j] = TempVar + LowFreqEdge;
	}

   glb.Grid[j-1] = UpperFreq;
   glb.DesiredMag[j-1] = glb.BandMag[BandIndex];
   glb.Weight[j-1] = glb.InitWeight[BandIndex];
   k+=2;
   BandIndex++;
   if(BandIndex <= NumBands)    glb.Grid[j] = glb.Edge[k];
  }

 GridIndex = j-1;
 if(!OddNumTaps && glb.Grid[GridIndex] > (0.5-LowFreqEdge)) GridIndex--;

 if(!OddNumTaps)
  {
   for(j=1; j<=GridIndex; j++)
	{
     Change = cos(M_PI * glb.Grid[j] );
     glb.DesiredMag[j] = glb.DesiredMag[j] / Change;
     glb.Weight[j] = glb.Weight[j] * Change;
	}
  }

 TempVar = (double)(GridIndex-1)/(double)glb.HalfTapCount;
 for(j=1; j<=glb.HalfTapCount; j++)
  {
   glb.ExchangeIndex[j] = (int)( (double)(j-1) * TempVar + 1.0 );
  }
 glb.ExchangeIndex[glb.HalfTapCount+1] = GridIndex;

 NumIterations = Remez2(&glb, GridIndex);
 CalcCoefficients(&glb);

 // Calculate the impulse response.
 if(OddNumTaps)
  {
   for(j=1; j<=glb.HalfTapCount-1; j++)
	{
     glb.Coeff[j] = 0.5 * glb.Alpha[glb.HalfTapCount+1-j];
	}
   glb.Coeff[glb.HalfTapCount] = glb.Alpha[1];
  }
 else
  {
   glb.Coeff[1] = 0.25 * glb.Alpha[glb.HalfTapCount];
   for(j=2; j<=glb.HalfTapCount-1; j++)
	{
     glb.Coeff[j] = 0.25 * (glb.Alpha[glb.HalfTapCount+1-j] + glb.Alpha[glb.HalfTapCount+2-j]);
	}
   glb.Coeff[glb.HalfTapCount] = 0.5 * glb.Alpha[1] + 0.25 * glb.Alpha[2];
  }


 // Output section.
 for(j=1; j<=glb.HalfTapCount; j++)FirCoeff[j-1] = glb.Coeff[j];
 if(OddNumTaps)
  for(j=1; j<glb.HalfTapCount; j++)FirCoeff[glb.HalfTapCount+j-1] = glb.Coeff[glb.HalfTapCount-j];
 else
  for(j=1; j<=glb.HalfTapCount; j++)FirCoeff[glb.HalfTapCount+j-1] = glb.Coeff[glb.HalfTapCount-j+1];

 // Display the iteration count.
 if(NumIterations < 3)
  {
   // ShowMessage("Parks McClellan unable to coverge");
   return 2;
  }

 return 0;
}

//---------------------------------------------------------------------------------------
int Remez2(ParksGlobals *PtrGlobals, int GridIndex)
{
 ParksGlobals &glb = *PtrGlobals;
 int j, JET, K, k, NU, JCHNGE, K1, KNZ, KLOW, NUT, KUP;
 int NUT1 = 0, LUCK, KN, NITER;
 double Deviation, DNUM, DDEN, TempVar;
 double DEVL, COMP = 0, YNZ = 0, Y1 = 0, ERR;

 LUCK = 0;
 DEVL = -1.0;
 NITER = 1; // Init this to 1 to be consistent with the orig code.

 TOP_LINE:  // We come back to here from 3 places at the bottom.
 glb.ExchangeIndex[glb.HalfTapCount+2] = GridIndex + 1;

 for(j=1; j<=glb.HalfTapCount+1; j++)
  {
   TempVar = glb.Grid[ glb.ExchangeIndex[j] ];
   glb.CosOfGrid[j] = cos(TempVar * M_2PI);
  }

 JET = (glb.HalfTapCount-1)/15 + 1;
 for(j=1; j<=glb.HalfTapCount+1; j++)
  {
   glb.LeGrangeD[j] = LeGrangeInterp2(&glb, j,glb.HalfTapCount+1,JET);
  }

 DNUM = 0.0;
 DDEN = 0.0;
 K = 1;
 for(j=1; j<=glb.HalfTapCount+1; j++)
  {
   k = glb.ExchangeIndex[j];
   DNUM += glb.LeGrangeD[j] * glb.DesiredMag[k];
   DDEN += (double)K * glb.LeGrangeD[j]/glb.Weight[k];
   K = -K;
  }
 Deviation = DNUM / DDEN;

 NU = 1;
 if(Deviation > 0.0) NU = -1;
 Deviation = -(double)NU * Deviation;
 K = NU;
 for(j=1; j<=glb.HalfTapCount+1; j++)
  {
   k = glb.ExchangeIndex[j];
   TempVar = (double)K * Deviation/glb.Weight[k];
   glb.DesPlus[j] = glb.DesiredMag[k] + TempVar;
   K = -K;
  }

 if(Deviation <= DEVL)return(NITER); // Ouch

 DEVL = Deviation;
 JCHNGE = 0;
 K1 = glb.ExchangeIndex[1];
 KNZ = glb.ExchangeIndex[glb.HalfTapCount+1];
 KLOW = 0;
 NUT = -NU;

 //Search for the extremal frequencies of the best approximation.

 j=1;
 while(j<glb.HalfTapCount+2)
  {
   KUP = glb.ExchangeIndex[j+1];
   k = glb.ExchangeIndex[j] + 1;
   NUT = -NUT;
   if(j == 2) Y1 = COMP;
   COMP = Deviation;

   if(k < KUP && !ErrTest(&glb, k, NUT, COMP, &ERR))
	{
	 L210:
	 COMP = (double)NUT * ERR;
	 for(k++; k<KUP; k++)
	  {
       if( ErrTest(&glb, k, NUT, COMP, &ERR) )break; // for loop
	   COMP = (double)NUT * ERR;
	  }

     glb.ExchangeIndex[j] = k-1;
	 j++;
	 KLOW = k - 1;
	 JCHNGE++;
	 continue;  // while loop
	}

   k--;

   L225: k--;
   if(k <= KLOW)
	{
     k = glb.ExchangeIndex[j] + 1;
	 if(JCHNGE > 0)
	  {
       glb.ExchangeIndex[j] = k-1;
	   j++;
	   KLOW = k - 1;
	   JCHNGE++;
	   continue;  // while loop
	  }
	 else  // JCHNGE <= 0
	  {
	   for(k++; k<KUP; k++)
		{
         if( ErrTest(&glb, k, NUT, COMP, &ERR) )continue; // for loop
		 goto L210;
		}

        KLOW = glb.ExchangeIndex[j];
		j++;
		continue; // while loop
	   }
	}
   // Can't use a do while loop here, it would outdent the two continue statements.
   if( ErrTest(&glb, k, NUT, COMP, &ERR) && JCHNGE <= 0)goto L225;

   if( ErrTest(&glb, k, NUT, COMP, &ERR) )
	{
     KLOW = glb.ExchangeIndex[j];
	 j++;
	 continue; // while loop
	}

   COMP = (double)NUT * ERR;

   L235:
   for(k--; k>KLOW; k--)
	{
     if( ErrTest(&glb, k, NUT, COMP, &ERR) ) break; // for loop
	 COMP = (double)NUT * ERR;
	}

   KLOW = glb.ExchangeIndex[j];
   glb.ExchangeIndex[j] = k + 1;
   j++;
   JCHNGE++;
 }  // end while(j<HalfTapCount

 if(j == glb.HalfTapCount+2) YNZ = COMP;

 while(j <= glb.HalfTapCount+2)
  {
   if(K1 > glb.ExchangeIndex[1]) K1 = glb.ExchangeIndex[1];
   if(KNZ < glb.ExchangeIndex[glb.HalfTapCount+1]) KNZ = glb.ExchangeIndex[glb.HalfTapCount+1];
   NUT1 = NUT;
   NUT = -NU;
   k = 0 ;
   KUP = K1;
   COMP = YNZ * 1.00001;
   LUCK = 1;

   for(k++; k<KUP; k++)
	{
     if( ErrTest(&glb, k, NUT, COMP, &ERR) ) continue; // for loop
     j = glb.HalfTapCount+2;
	 goto L210;
	}
   LUCK = 2;
   break;  // break while(j <= HalfTapCount+2) loop
  } // end while(j <= HalfTapCount+2)

 if(LUCK == 1 || LUCK == 2)
  {
   if(LUCK == 1)
	{
	 if(COMP > Y1) Y1 = COMP;
     K1 = glb.ExchangeIndex[glb.HalfTapCount+2];
	}

   k = GridIndex + 1;
   KLOW = KNZ;
   NUT = -NUT1;
   COMP = Y1 * 1.00001;

   for(k--; k>KLOW; k--)
	{
     if( ErrTest(&glb, k, NUT, COMP, &ERR) )continue;  // for loop
     j = glb.HalfTapCount+2;
	 COMP = (double)NUT * ERR;
	 LUCK = 3;   // last time in this if(LUCK == 1 || LUCK == 2)
	 goto L235;
	}

   if(LUCK == 2)
	{
	 if(JCHNGE > 0 && NITER++ < ITRMAX) goto TOP_LINE;
	 else return(NITER);
	}

   for(j=1; j<=glb.HalfTapCount; j++)
	{
     glb.ExchangeIndex[glb.HalfTapCount+2 - j] = glb.ExchangeIndex[glb.HalfTapCount+1 - j];
	}
   glb.ExchangeIndex[1] = K1;
   if(NITER++ < ITRMAX) goto TOP_LINE;
  }  // end if(LUCK == 1 || LUCK == 2)



 KN = glb.ExchangeIndex[glb.HalfTapCount+2];
 for(j=1; j<=glb.HalfTapCount; j++)
  {
   glb.ExchangeIndex[j] = glb.ExchangeIndex[j+1];
  }
 glb.ExchangeIndex[glb.HalfTapCount+1] = KN;
 if(NITER++ < ITRMAX) goto TOP_LINE;


 return(NITER);

}

//-----------------------------------------------------------------------
// Function to calculate the lagrange interpolation coefficients for use in the function gee.
static double LeGrangeInterp2(const ParksGlobals *PtrGlobals, int K, int N, int M) // D
{
 const ParksGlobals &glb = *PtrGlobals;
 int j, k;
 double Dee, Q;
 Dee = 1.0;
 Q = glb.CosOfGrid[K];
 for(k=1; k<=M; k++)
 for(j=k; j<=N; j+=M)
  {
   if(j != K)Dee = 2.0 * Dee * (Q - glb.CosOfGrid[j]);
  }
 if(fabs(Dee) < MIN_TEST_VAL )
  {
   if(Dee < 0.0)Dee = -MIN_TEST_VAL;
   else         Dee =  MIN_TEST_VAL;
  }
 return(1.0/Dee);
}

//-----------------------------------------------------------------------
// Function to evaluate the frequency response using the Lagrange interpolation
// formula in the barycentric form.
static double GEE2(const ParksGlobals *PtrGlobals, int K, int N)
{
 const ParksGlobals &glb = *PtrGlobals;
 int j;
 double P,C,Dee,XF;
 P = 0.0;
 XF = glb.Grid[K];
 XF = cos(M_2PI * XF);
 Dee = 0.0;
 for(j=1; j<=N; j++)
  {
   C = XF - glb.CosOfGrid[j];
   if(fabs(C) < MIN_TEST_VAL )
	{
	 if(C < 0.0)C = -MIN_TEST_VAL;
	 else       C =  MIN_TEST_VAL;
	}
   C = glb.LeGrangeD[j] / C;
   Dee = Dee + C;
   P = P + C*glb.DesPlus[j];
  }
 if(fabs(Dee) < MIN_TEST_VAL )
  {
   if(Dee < 0.0)Dee = -MIN_TEST_VAL;
   else         Dee =  MIN_TEST_VAL;
  }
 return(P/Dee);
}

//-----------------------------------------------------------------------

static bool ErrTest(const ParksGlobals *PtrGlobals, int k, int Nut, double Comp, double *Err)
{
 const ParksGlobals &glb = *PtrGlobals;
 *Err = GEE2(&glb, k, glb.HalfTapCount+1);
 *Err = (*Err - glb.DesiredMag[k]) * glb.Weight[k];
 if((double)Nut * *Err - Comp <= 0.0) return(true);
 else return(false);
}

//-----------------------------------------------------------------------

// Calculation of the coefficients of the best approximation using the inverse discrete fourier transform.
static void CalcCoefficients(ParksGlobals *PtrGlobals)
{
 ParksGlobals &glb = *PtrGlobals;
 int j, k, n;
 double GTempVar, OneOverNumTaps;
 double Omega, TempVar, FreqN, TempX,  GridCos;
 double GeeArray[PARKS_SMALL];

 GTempVar = glb.Grid[1];
 glb.CosOfGrid[glb.HalfTapCount+2] = -2.0;
 OneOverNumTaps = 1.0 /(double)(2*glb.HalfTapCount-1);
 k = 1;

 for(j=1; j<=glb.HalfTapCount; j++)
  {
   FreqN = (double)(j-1) * OneOverNumTaps;
   TempX = cos(M_2PI * FreqN);

   GridCos = glb.CosOfGrid[k];
   if(TempX <= GridCos)
	{
	 while(TempX <= GridCos && (GridCos-TempX) >= MIN_TEST_VAL) // MIN_TEST_VAL = 1.0E-6
	  {
	   k++;;
       GridCos = glb.CosOfGrid[k];
	  }
	}
   if(TempX <= GridCos || (TempX-GridCos) < MIN_TEST_VAL)
	{
     GeeArray[j] = glb.DesPlus[k]; // Desired Response
	}
   else
	{
     glb.Grid[1] = FreqN;
     GeeArray[j] = GEE2(&glb, 1, glb.HalfTapCount+1);
	}
   if(k > 1) k--;
  }

 glb.Grid[1] = GTempVar;
 for(j=1; j<=glb.HalfTapCount; j++)
  {
   TempVar = 0.0;
   Omega = (double)(j-1) * M_2PI * OneOverNumTaps;
   for(n=1; n<=glb.HalfTapCount-1; n++)
	{
	 TempVar += GeeArray[n+1] * cos(Omega * (double)n);
	}
   TempVar = 2.0 * TempVar + GeeArray[1];
   glb.Alpha[j] = TempVar;
  }

 glb.Alpha[1] = glb.Alpha[1] * OneOverNumTaps;
 for(j=2; j<=glb.HalfTapCount; j++)
  {
   glb.Alpha[j] = 2.0 * glb.Alpha[j] * OneOverNumTaps;
  }

}

//-----------------------------------------------------------------------

