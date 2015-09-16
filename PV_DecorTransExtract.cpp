/**
 * Transient extraction plugin for supercollider
 * By Bryan Dalle Molle (bryandallemolle@gmail.com)
 * University of Illinois at Chicago 2015
 * 
 * Based on approach detailed by Ross Penniman in 
 * "A High-Quality General-Purpose
 * Decorrelator With Transient Fidelity"
 *
 */

#include "SC_PlugIn.h"
#include "FFT_UGens.h"

InterfaceTable *ft;

#define MIN_VAL 0.000000001		// this should be better chosen values...
#define MAX_VAL 9.999999999		// as should this...

struct PV_DecorTransExtract : Unit 
{
	// arguments
	bool retTrans;
	float alphaVal;
	float betaVal;
	unsigned int iVal;					
	unsigned int jVal;
	unsigned int dVal;
	unsigned int lowFreqCutVal; 

	// other variables
	unsigned int prevBinsIdx;
	unsigned int littleOmegaIdx;
	unsigned int bigOmegaIdx;
	unsigned int numFreqBins;
	bool initFirstCalc;
	bool initDebugFlag;

	// dynamic memory
	bool* transFlagBins; 
	float* releaseBins;
	float** prevBins;
	float** littleOmegaBins;
	float* bigOmegaBins;
};

static void PV_DecorTransExtract_Ctor(PV_DecorTransExtract *unit);
static void PV_DecorTransExtract_next(PV_DecorTransExtract *unit, int inNumSamples);
static void PV_DecorTransExtract_Flanger_Dtor(PV_DecorTransExtract *unit);
static void firstCalc(PV_DecorTransExtract *unit);
static float calcAverageBins(PV_DecorTransExtract *unit, int freqBinIdx);
static float calcMinBin(PV_DecorTransExtract *unit, int freqBinIdx);
static float calcReleaseBin(PV_DecorTransExtract *unit, int freqBinIdx);
static float min(float a, float b);
static float max(float a, float b);

void PV_DecorTransExtract_Ctor(PV_DecorTransExtract *unit)
{
	// get input variables
	unit->retTrans = IN0(1);
	unit->iVal = IN0(2);
	unit->jVal = IN0(3);
	unit->alphaVal = IN0(4);
	unit->betaVal = IN0(5);
	unit->dVal = IN0(6);
	unit->lowFreqCutVal = IN0(7);

	// set other variables
	unit->numFreqBins = 0;
	unit->prevBinsIdx = 0;
	unit->littleOmegaIdx = 0;
	unit->bigOmegaIdx = 0;
	unit->initDebugFlag = true;
	unit->initFirstCalc = true;

	// set output
	ZOUT0(0) = ZIN0(0);

	// set calculation function (next)
	SETCALC(PV_DecorTransExtract_next);
}

void firstCalc(PV_DecorTransExtract *unit)
{
	// set the number of bins from first calculation
	if (unit->initDebugFlag)
	{
		printf("initial next function call\n");
		printf("ARGS:\n");
		printf(" - numFreqBins = %d\n", unit->numFreqBins);
		printf(" - iVal = %u\n", unit->iVal);
		printf(" - jVal = %u\n", unit->jVal);
		printf(" - alphaVal = %f\n", unit->alphaVal);
		printf(" - betaVal = %f\n", unit->betaVal);
		printf(" - dVal = %d\n", unit->dVal);
		printf(" - lowFreqCutVal = %d\n", unit->lowFreqCutVal);
	}

	// allocate some memory
	unit->prevBins = (float**)RTAlloc(unit->mWorld, unit->iVal * sizeof(float*));
	for (int i = 0; i < unit->iVal; ++i)
		unit->prevBins[i] = (float*)RTAlloc(unit->mWorld, unit->numFreqBins * sizeof(float));

	unit->littleOmegaBins = (float**)RTAlloc(unit->mWorld, unit->jVal * sizeof(float*));
	for (int j = 0; j < unit->jVal; ++j)
		unit->littleOmegaBins[j] = (float*)RTAlloc(unit->mWorld, unit->numFreqBins * sizeof(float));

	unit->bigOmegaBins = (float*)RTAlloc(unit->mWorld, unit->numFreqBins * sizeof(float));
	unit->transFlagBins = (bool*)RTAlloc(unit->mWorld, unit->numFreqBins * sizeof(bool));
	unit->releaseBins = (float*)RTAlloc(unit->mWorld, unit->numFreqBins * sizeof(float));

	// zero out all memory
	for (int i = 0; i < unit->iVal; ++i)
		memset(unit->prevBins[i], 0, unit->numFreqBins * sizeof(float));

	for (int j = 0; j < unit->jVal; ++j)
		memset(unit->littleOmegaBins[j], 0, unit->numFreqBins * sizeof(float));

	memset(unit->bigOmegaBins, 0, unit->numFreqBins * sizeof(float));
	memset(unit->releaseBins, 1.0, unit->numFreqBins * sizeof(float));

	if (unit->initDebugFlag)
	{
		printf("successful init!\n");
		unit->initDebugFlag = false;
	}

}

float min(float a, float b)
{
	if (a < b) return a;
	else return b;
}

float max(float a, float b)
{
	if (a > b) return a;
	else return b;
}

float calcAverageBins(PV_DecorTransExtract *unit, int k)
{
	float tempSum = 0.0;
	float average = 0.0;

	for (int i = 0; i < unit->iVal; ++i)
		tempSum += unit->prevBins[i][k];

	average = tempSum / unit->iVal;

	return average;
}

// k = frequency bin index
// s = time slice of bin
float calcMinBin(PV_DecorTransExtract *unit, int k)
{
	float minBin = MAX_VAL;
	int currBin = unit->littleOmegaIdx;

	for (int j = 0; j < unit->jVal; ++j)
		// ignore the current bin, it can't be the minimum...
		if (j != currBin && unit->littleOmegaBins[j][k] < minBin)
			minBin = unit->littleOmegaBins[j][k];

	return minBin;
}

float calcReleaseBin(PV_DecorTransExtract *unit, int k)
{
	float weightedRelease = unit->betaVal * unit->releaseBins[k];

	// get the [s - 1] little omega bin value
	int prevS = unit->littleOmegaIdx;
	if (prevS < 0) prevS = unit->jVal - 1;

	// if transient bin...
	if (unit->transFlagBins[k])
		return min(weightedRelease, unit->littleOmegaBins[prevS][k]);
	// if not a transient bin...
	else return max(weightedRelease, MIN_VAL);
}

void PV_DecorTransExtract_next(PV_DecorTransExtract *unit, int inNumSamples)
{
	unsigned int numTransBins = 0;
	unsigned int diffBins = 0;
	float kVal = 0.0;
	float outputVal = 0.0;
	float gain = 1.0;

	// this macro gets an FFT chain if it has "fired"
	// otherwise, returns -1 and we go on our way...
	PV_GET_BUF

	// primes the ugen upon first calculation
	if (unit->initFirstCalc)
	{
		unit->numFreqBins = numbins;
		firstCalc(unit);
		unit->initFirstCalc = false;
	}
	
	// check that number of bins is static! this can't change!
	if (unit->numFreqBins != numbins)
	{
		printf("PV_DecorTransExtract: num bins mismatch!\n");
		printf(" - fft chain must be of constant size with PV_DecoreTransExtract()\n");
		printf(" - userbins = %d, numbins = %d\n", unit->numFreqBins, numbins);
		return;
	}
	else if (unit->lowFreqCutVal > unit->numFreqBins)
	{
		printf("PV_DecorTransExtract: low cut bin > number of bins\n");
		printf("setting low cut to 0\n");
		unit->lowFreqCutVal = 0;
		return;
	}

	SCPolarBuf *p = ToPolarApx(buf);

	// store bin values
	for (int i = unit->lowFreqCutVal; i < numbins; ++i)
		unit->prevBins[unit->prevBinsIdx][i] = p->bin[i].mag;

	// calculate average freq history for each bin
	for (int j = unit->lowFreqCutVal; j < numbins; ++j)
		unit->littleOmegaBins[unit->littleOmegaIdx][j] = calcAverageBins(unit, j);

	// reset transient flag
	memset(unit->transFlagBins, 0, unit->numFreqBins * sizeof(bool));

	// printf("mark transients. little omega index = %d\n", unit->littleOmegaIdx);

	for (int i = unit->lowFreqCutVal; i < numbins; ++i)
	{
		unit->bigOmegaBins[i] = calcMinBin(unit, i);

		if (p->bin[i].mag > unit->alphaVal * unit->bigOmegaBins[i])
			numTransBins++;
		else 
		{
			if (numTransBins >= unit->dVal) 
			{
				diffBins = (i - 1) - numTransBins;
				for (int j = (i - 1); j >= diffBins; --j)
					unit->transFlagBins[j] = true;	
			}
			numTransBins = 0;
		}
	}

	// printf("release stage. little omega index = %d\n", unit->littleOmegaIdx);

	// this is where the lower freq bins are cut out
	if (unit->retTrans)
		for (int i = 0; i < unit->lowFreqCutVal; ++i)
		 	p->bin[i].mag = 0.0;

	for (int i = unit->lowFreqCutVal; i < numbins; ++i)
	{
		unit->releaseBins[i] = calcReleaseBin(unit, i);

		kVal = min(unit->releaseBins[i], p->bin[i].mag);
		outputVal = kVal / p->bin[i].mag;

		if (unit->retTrans)
			gain = 1.0 - outputVal;
		else gain = outputVal;

		// apply gain to output
		p->bin[i].mag = gain * p->bin[i].mag;
	}

	// increment circular queue of bins
	if (unit->prevBinsIdx < unit->iVal - 1)
		unit->prevBinsIdx++;
	else unit->prevBinsIdx = 0;

	if (unit->littleOmegaIdx < unit->jVal - 1)
		unit->littleOmegaIdx++;
	else unit->littleOmegaIdx = 0;
}

void PV_DecorTransExtract_Dtor(PV_DecorTransExtract *unit)
{
	RTFree(unit->mWorld, unit->bigOmegaBins);

	for (int i = 0; i < unit->iVal; ++i)
		RTFree(unit->mWorld, unit->littleOmegaBins[i]);
	RTFree(unit->mWorld, unit->littleOmegaBins);

	for (int j = 0; j < unit->jVal; ++j)
		RTFree(unit->mWorld, unit->prevBins[j]);
	RTFree(unit->mWorld, unit->prevBins);

	RTFree(unit->mWorld, unit->transFlagBins);
	RTFree(unit->mWorld, unit->releaseBins);
}

PluginLoad(PV_DecorTransExtract)
{
	ft = inTable;
	DefineDtorUnit(PV_DecorTransExtract);
}