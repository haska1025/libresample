/**********************************************************************

  resamplesubs.c

  Real-time library interface by Dominic Mazzoni

  Based on resample-1.7:
    http://www-ccrma.stanford.edu/~jos/resample/

  Dual-licensed as LGPL and BSD; see README.md and LICENSE* files.

  This file provides the routines that do sample-rate conversion
  on small arrays, calling routines from filterkit.

**********************************************************************/

/* Definitions */
#include "resample_defs.h"

#include "filterkit.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/* Sampling rate up-conversion only subroutine;
 * Slightly faster than down-conversion;
 */
int lrsSrcUp(float X[],
             float Y[],
             double factor,
             double *TimePtr,
             UWORD Nx,
             UWORD Nwing,
             float LpScl,
             float Imp[],
             float ImpD[],
             BOOL Interp)
{
    float *Xp, *Ystart;
    float v;
    
    double CurrentTime = *TimePtr;
    double dt;                 /* Step through input signal */ 
    double endTime;            /* When Time reaches EndTime, return to user */
    
    /**
     * 这一段逻辑主要是进行插值。
     * fo(frequence old)：重采样前，信号 x[n] 的采样频率
     * to(time old)：重采样前，采样间隔。to = 1/fo
     * fn(frequence new): 重采样后，信号 x[n] 的采样频率
     * tn(time new): 重采样后，采样间隔。tn = 1/fn
     * L = fn/fo，这里的 factor 就是 fn/fo，factor = (1/tn) / (1/to) = to/tn, tn = to/factor，可以得知重采样的时间间隔是采样前的 1/factor
     */
    dt = 1.0/factor;           /* Output sampling period */
    
    Ystart = Y;
    endTime = CurrentTime + Nx;
    while (CurrentTime < endTime)
    {
        double LeftPhase = CurrentTime-floor(CurrentTime);
        double RightPhase = 1.0 - LeftPhase;

        // 插值，连续插入 factor 个 x[(int)currentTime]
        Xp = &X[(int)CurrentTime]; /* Ptr to current input sample */
        /* Perform left-wing inner product */
        v = lrsFilterUp(Imp, ImpD, Nwing, Interp, Xp,
                        LeftPhase, -1);
        /* Perform right-wing inner product */
        v += lrsFilterUp(Imp, ImpD, Nwing, Interp, Xp+1, 
                         RightPhase, 1);

        v *= LpScl;   /* Normalize for unity filter gain */

        *Y++ = v;               /* Deposit output */
        CurrentTime += dt;      /* Move to next sample by time increment */
    }

    *TimePtr = CurrentTime;
    return (Y - Ystart);        /* Return the number of output samples */
}

/* Sampling rate conversion subroutine */

int lrsSrcUD(float X[],
             float Y[],
             double factor,
             double *TimePtr,
             UWORD Nx,
             UWORD Nwing,
             float LpScl,
             float Imp[],
             float ImpD[],
             BOOL Interp)
{
    float *Xp, *Ystart;
    float v;

    double CurrentTime = (*TimePtr);
    double dh;                 /* Step through filter impulse response */
    double dt;                 /* Step through input signal */
    double endTime;            /* When Time reaches EndTime, return to user */
    
    dt = 1.0/factor;            /* Output sampling period */
    
    dh = MIN(Npc, factor*Npc);  /* Filter sampling period */
    
    Ystart = Y;
    endTime = CurrentTime + Nx;
    while (CurrentTime < endTime)
    {
        double LeftPhase = CurrentTime-floor(CurrentTime);
        double RightPhase = 1.0 - LeftPhase;

        Xp = &X[(int)CurrentTime];     /* Ptr to current input sample */
        /* Perform left-wing inner product */
        v = lrsFilterUD(Imp, ImpD, Nwing, Interp, Xp,
                        LeftPhase, -1, dh);
        /* Perform right-wing inner product */
        v += lrsFilterUD(Imp, ImpD, Nwing, Interp, Xp+1, 
                         RightPhase, 1, dh);

        v *= LpScl;   /* Normalize for unity filter gain */
        *Y++ = v;               /* Deposit output */
        
        CurrentTime += dt;      /* Move to next sample by time increment */
    }

    *TimePtr = CurrentTime;
    return (Y - Ystart);        /* Return the number of output samples */
}
