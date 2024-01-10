#include <jni.h>
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>

#include "SavGolJNI.h"

#define UPPER true
#define LOWER false
#define MIN_RADIUS 12

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

struct gh_values // global host values
{
	float *d_values;
	float *h_values;
	float *precalcPow;
};

int blockSizePowers;
int blockSizeSmooth;

/* =============================================================================
* Checks if a CUDA capable device is present
*
* @param env:       JNI environment
* @param thisObj:   reference to "this" Java object
*
* @return:          true if there is a cuda capable device
*/
JNIEXPORT jboolean JNICALL Java_at_tugraz_genome_lda_quantification_SavGolJNI_cudaCapableDeviceNative
(JNIEnv *env, jobject thisObj)
{

	int deviceCount = 0;
	cudaGetDeviceCount(&deviceCount);

	if (deviceCount == 0)
		return false;
	else
		return true;

}

/* =============================================================================
* Checks if an error accured on the device
*
* @param ans: a cuda function call
*/
// http://stackoverflow.com/a/14038590/6914637
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void
gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
	if (code != cudaSuccess)
	{
		printf("GPUassert: %s %s %d\n", cudaGetErrorString(code), file,
			line);
		if (abort)
			exit(code);
	}
}

/* =============================================================================
* Finds the lowest or highest index within the range
*
* @param dtIndx:        The point to extend the range
* @param range:         The time which the index should be extendet to
* @param posDirection:  Gives the information to move in the array up or down
* @param values:        Contains time, raw data and smoothed data
* @param numberOfScans: Length of values
*
* @return the lowest or highest index
*/
int calcBoundIndex(int dtIndx, float range, bool posDirection, float *values,
	int numberOfScans)
{
	int boundIndex = dtIndx;
	if (posDirection == UPPER)
	{
		while (boundIndex < (numberOfScans - 1)
			&& (values[boundIndex * 4] - values[dtIndx * 4]) < range)
		{
			++boundIndex;
		}
	}
	else
	{
		while (boundIndex > 0
			&& (values[dtIndx * 4] - values[boundIndex * 4]) < range)
		{
			--boundIndex;
		}
	}
	return boundIndex;
}

/* =============================================================================
* Pre calculates the 4th root of the raw/smoothed values
*
* @param powers:        holds the result of the calculations
* @param values:        Contains time, raw data and smoothed data
* @param numberOfScans: Length of values
* @param copyDirection: location of the smoothed data
* 						    ( 0->d_values[][2], 1->d_values[][3] )
*/
__global__
void
precalculatePowers(float *powers, float *values, int numberOfScans,
	int copyDirection)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	while (idx < numberOfScans)
	{
		if (values[4 * idx + 2 + copyDirection] > 1)
			powers[idx] = rsqrt(rsqrt(values[4 * idx + 2 + copyDirection]));
		else
			powers[idx] = 1.0f;
		idx += blockDim.x * gridDim.x;
	}
}

/* =============================================================================
* Finds the lowest or highest index within the range
*
* @param boundIndex:    The calculated lowest or highest index
* @param idx:           The point to extend the range
* @param posDirection:  Gives the information to move in the array up or down
* @param numberOfScans: Length of values
* @param values:        Contains time, raw data and smoothed data
* @param range:	        The time which the index should be extendet to
*
* @return the lowest or highest index
*/
__device__
void
deviceCalcBoundIndex(int *boundIndex, int idx, bool posDirection,
	int numberOfScans, float *values, float range)
{
	*boundIndex = idx;
	if (posDirection == UPPER)
	{
		while (*boundIndex < (numberOfScans - 1)
			&& (values[4 * (*boundIndex) + 0] - values[4 * idx + 0]) < range)
		{
			++*boundIndex;
		}
	}
	else
	{
		while (*boundIndex > 0
			&& (values[4 * idx + 0] - values[4 * (*boundIndex) + 0]) < range)
		{
			--*boundIndex;
		}
	}
}

/* =============================================================================
* performs LU decomposition
*
* @param mtrx
* @param order
* @param indx
* @param vv
*/
__device__
void
myLUDcmp(float mtrx[][5], int order, int *indx, float *vv)
{
	float big, dum, sum, temp;
	int i, j, k;
	int imax = 0;
	float d = 1.0f;

	for (i = 0; i < order; i++)
	{
		big = 0.0f;
		for (j = 0; j < order; j++)
			if ((temp = fabs(mtrx[i][j])) > big)
				big = temp;
		vv[i] = 1.0F / big;
	}

	for (j = 0; j < order; j++)
	{
		for (i = 0; i < j; i++)
		{
			sum = mtrx[i][j];
			for (k = 0; k < i; k++)
				sum -= mtrx[i][k] * mtrx[k][j];
			mtrx[i][j] = sum;
		}
		big = 0.0f;
		for (i = j; i < order; i++)
		{
			sum = mtrx[i][j];
			for (k = 0; k < j; k++)
				sum -= mtrx[i][k] * mtrx[k][j];
			mtrx[i][j] = sum;
			dum = vv[i] * fabs(sum);
			if (dum >= big)
			{
				big = dum;
				imax = i;
			}
		}
		if (j != imax)
		{
			for (k = 0; k < order; k++)
			{
				dum = mtrx[imax][k];
				mtrx[imax][k] = mtrx[j][k];
				mtrx[j][k] = dum;
			}
			d *= -1;
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (mtrx[j][j] == 0.0)
			mtrx[j][j] = 1E-20f;
		if (j != order)
		{
			dum = 1.0F / mtrx[j][j];
			for (i = j + 1; i < order; i++)
				mtrx[i][j] *= dum;
		}
	}
}

/* =============================================================================
* performs LU backsubstition
*
* @param mtrx
* @param order
* @param indx
* @param vec
*/
__device__
void
myLUBksb(float mtrx[][5], int order, int *indx, float *vec)
{
	float sum;
	int ii, ip;
	int i, j;

	ii = -1;
	for (i = 0; i < order; i++)
	{
		ip = indx[i];
		sum = vec[ip];
		vec[ip] = vec[i];
		if (ii != -1)
			for (j = ii; j < i; j++)
				sum -= mtrx[i][j] * vec[j];
		else if (sum != 0)
			ii = i;
		vec[i] = sum;
	}
	for (i = order - 1; i >= 0; i--)
	{
		sum = vec[i];
		for (j = i + 1; j < order; j++)
			sum -= mtrx[i][j] * vec[j];
		vec[i] = sum / mtrx[i][i];
	}
}

/* =============================================================================
* Does a Savitzky Golay Filter around a given point
*
* @param val:           The smoothed value at idx
* @param g_values:      Stored on the global memory.
*                           Contains time, raw data and smoothed data
* @param g_powers:      Precalculated powers
* @param idx:           Index of point
* @param lower:         Lower range border
* @param upper:         Upper range border
* @param order:         Order of polynome
* @param copyDirection: Location of the smoothed data
*                           ( 0->g_values[][2], 1->g_values[][3] )
* @param numberOfScans: Length of g_values
* @param radius:        Maximum number of points within the range
* @param shared_memory: Dynamically allocated shared memory
* @param blockSize:     Block size
*/
__device__
void
SavGolFilter(float *val, float *g_values, float *g_powers, int idx, int lower,
	int upper, int order, int copyDirection, int numberOfScans,
	int radius, float *shared_memory, int blockSize)
{
	// -------------------------------------------------------------------------
	// copy values[] and powers[] from the global to the shared memory

	// setting new indices
	int s_idx = threadIdx.x + radius;
	int offset = blockSize * (idx / blockSize) - radius;
	lower = lower - offset;
	upper = upper - offset;

	// init shared variables
	float* values = &shared_memory[0];
	float* powers = &shared_memory[4 * blockSize + 2 * 4 * radius];

	// copy the data
	for (int i = 0; i < 4; i++)
		values[4 * s_idx + i] = g_values[4 * idx + i];
	powers[s_idx] = g_powers[idx];

	// copy the data outside of the radius
	if (threadIdx.x < radius)
	{
		if (idx - radius >= 0)
		{
			for (int i = 0; i < 4; i++)
			{
				values[4 * (s_idx - radius) + i] = g_values[4 * (idx - radius) + i];
			}
			powers[s_idx - radius] = g_powers[idx - radius];
		}
		if (idx + blockSize < numberOfScans)
		{
			for (int i = 0; i < 4; i++)
			{
				values[4 * (s_idx + blockSize) + i] = g_values[4 * (idx + blockSize) + i];
			}
			powers[s_idx + blockSize] = g_powers[idx + blockSize];
		}
	}

	// wait until every thread within a block is done with copying to shared
	// memory
	__syncthreads();
	
	float x = values[4 * s_idx + 0];
	float sum = 0.0f;
	float adding = 0.0f;
	float mtrx[5][5];
	float vec[5];
	int indx[5];
	float vv[4];

	// -------------------------------------------------------------------------
	// get "mtrx"

	float sums[9];

	// calculate the sums for mtrx
	for (int m = 0; m < 2 * order + 1; m++)
	{
		sums[m] = 0.0f;

		for (int k = lower; k <= upper; k++)
		{
			adding = 1.0f;

			for (int l = 0; l < m; l++)
			{
				adding = adding * (values[4 * k + 0] - x);
			}
			if (m == 0)
				adding = 1.0f;

			adding = adding * powers[k];
			sums[m] += adding;

		}
	}

	// assign the sums to the mtrx
	for (int i = 0; i <= order; i++)
	{
		for (int j = i; j <= order; j++)
		{
			mtrx[i][j] = sums[i + j];
			mtrx[j][i] = sums[i + j];
		}
	}

	// -------------------------------------------------------------------------
	// LU decomposition

	myLUDcmp(mtrx, order + 1, indx, vv);

	// -------------------------------------------------------------------------
	// get "vec"

	for (int i = 0; i <= order; i++)
	{
		sum = 0.0f;
		for (int k = lower; k <= upper; k++)
		{
			sum += powf(values[4 * k + 0] - x, i)
				* values[4 * k + 2 + copyDirection] * powers[k];
		}
		vec[i] = sum;
	}

	// -------------------------------------------------------------------------
	// LU backsubstition

	myLUBksb(mtrx, order + 1, indx, vec);

	*val = vec[0];
}

extern __shared__ float shared_memory[];
/* =============================================================================
* smoothes raw spectrum at a specific point.
*
* @param range:         number of seconds around given points
* @param threshold:     minimum value to pass back
* @param powers:        precalculated powers
* @param values:        Contains time, raw data and smoothed data
* @param numberOfScans: Length of g_values
* @param copyDirection: location of the smoothed data
*                           ( 0->g_values[][2], 1->g_values[][3] )
* @param radius:        Maximum number of points within the range
* @param blockSize:     Block size
* @param startOffset:   Difference between precalculation start and smooth start
* @param stopOffset:    Difference between precalculation stop and smooth stop
*/
__global__
void
SmoothDataPoint(float range, float threshold, float *powers, float *values,
	int numberOfScans, int copyDirection, int radius,
	int blockSize, int startOffset, int stopOffset)
{
	int idx = threadIdx.x + blockIdx.x * blockDim.x;
	int order = 4;
	float val = 0;

	while (idx < numberOfScans)
	{
		// ---------------------------------------------------------------------
		// get boundary of indices which should be included into the calculation

		int lower;
		deviceCalcBoundIndex(&lower, idx, LOWER, numberOfScans, values, range);
		int upper;
		deviceCalcBoundIndex(&upper, idx, UPPER, numberOfScans, values, range);

		while (upper - lower < 10)
		{
			if (lower > 0)
				--lower;
			if (upper < numberOfScans - 1)
				++upper;
			if (lower <= 0 && upper >= numberOfScans - 1)
				break;
		}

		// ---------------------------------------------------------------------
		// Get the smoothed value

		order = order < upper - lower - 1 ? order : upper - lower - 1;

		if (order < 1)
		{
			val = values[4 * idx + 2 + copyDirection];
		}
		else
		{
			SavGolFilter(&val, values, powers, idx, lower, upper, order,
				copyDirection, numberOfScans, radius, shared_memory,
				blockSize);
			if (val < threshold)
			{
				val = threshold;
			}
		}
		
		if (idx < startOffset || idx >= (numberOfScans - stopOffset))
		{
			val = values[4 * idx + 2 + copyDirection];
		}

		values[4 * idx + 3 - copyDirection] = val;

		idx += blockDim.x * gridDim.x;
	}
}

/* =============================================================================
* Smooths the raw spectrum
*
* @param values:            Contains time, raw data and smoothed data
* @param numberOfScans:     Length of g_values
* @param range:             number of seconds around given points
* @param repeats:           how often the spectrum should be smoothed
* @param startSmoothScan:   start index of the spectrum
* @param stopSmoothScan:    end index of the spectrum
* @param address:           Adress of values struct for each thread
*/
void prepareSmooth(float *values, int numberOfScans,
	float range, int repeats, int startSmoothScan,
	int stopSmoothScan, long address)
{
	struct gh_values *ptr = (struct gh_values *)address;
	// -------------------------------------------------------------------------
	// Calculate the minimum intensity of the raw data and
	// store it into threshold:

	float threshold = values[1];
	for (int i = 0; i < numberOfScans; i++)
	{
		if (values[i * 4 + 1] < threshold)
		{
			threshold = values[i * 4 + 1];
		}
	}

	// -------------------------------------------------------------------------
	// Set the start and stop point of the smoothing range

	int startScan = 0;
	if (startSmoothScan > -1)
		startScan = startSmoothScan;
	int preCalcStart =
		calcBoundIndex(startScan, range, LOWER, values, numberOfScans) - 10 ;
	if (preCalcStart < 0)
		preCalcStart = 0;

	int stopScan = numberOfScans;
	if (stopSmoothScan > -1)
		stopScan = stopSmoothScan;
	int preCalcStop =
		calcBoundIndex(stopScan, range, UPPER, values, numberOfScans) + 10;
	if (preCalcStop > numberOfScans)
		preCalcStop = numberOfScans;

	// -------------------------------------------------------------------------
	// Adjust the length of values

	int smooth_length = preCalcStop - preCalcStart;
	float *smooth_values = values + preCalcStart * 4;

	// -------------------------------------------------------------------------
	// Setting the size of the block and grid dimension

	int gridSizePowers = (smooth_length + blockSizePowers - 1) / blockSizePowers;
	int gridSizeSmooth = (smooth_length + blockSizeSmooth - 1) / blockSizeSmooth;

	// -------------------------------------------------------------------------
	// Copy the data from the host to the device

	gpuErrchk(
		cudaMemcpy(ptr->d_values, smooth_values, 4 * smooth_length * sizeof(float),
			cudaMemcpyHostToDevice));

	// -------------------------------------------------------------------------
	// Calculate the shared memory allocation size
    
	float min_time_diff = smooth_values[(smooth_length - 1) * 4] - smooth_values[0];
	for (int i = 0; i < smooth_length - 1; i++)
	{
		float time_diff = smooth_values[(i + 1) * 4] - smooth_values[i * 4];
		min_time_diff = MIN(min_time_diff, time_diff);
	}
	int radius = MAX((range / min_time_diff + 1) * 2, MIN_RADIUS);

	int shared_space_values = 4 * blockSizeSmooth + 2 * 4 * radius;
	int shared_space_powers = blockSizeSmooth + 2 * radius;
	int shared_space = shared_space_values + shared_space_powers;

	// -------------------------------------------------------------------------
	// Smooth the raw data

	for (int i_rep = 0; i_rep < repeats; i_rep++)
	{
		precalculatePowers <<<gridSizePowers, blockSizePowers >>>(ptr->precalcPow,
			ptr->d_values, smooth_length, i_rep % 2);
		SmoothDataPoint <<<gridSizeSmooth, blockSizeSmooth,
			shared_space * sizeof(float) >>>(range, threshold, ptr->precalcPow,
			ptr->d_values, smooth_length, i_rep % 2, radius, blockSizeSmooth,
			startScan-preCalcStart, preCalcStop-stopScan);
	}
	gpuErrchk(cudaPeekAtLastError());
	gpuErrchk(cudaDeviceSynchronize());

	// -------------------------------------------------------------------------
	// copying back the data from the device

	gpuErrchk(
		cudaMemcpy(smooth_values, ptr->d_values, 4 * smooth_length * sizeof(float),
			cudaMemcpyDeviceToHost));

	if (repeats % 2 == 1)
	{
		for (int i = 2; i < 4 * smooth_length; i += 4)
		{
			values[i + preCalcStart * 4] = smooth_values[i + 1];
		}
	}
	else
	{
		for (int i = 2; i < 4 * smooth_length; i += 4)
		{
			values[i + preCalcStart * 4] = smooth_values[i];
		}
	}
}

/* =============================================================================
* Allocates memory on the graphicscard
*
* @param env:           JNI environment
* @param thisObj:       reference to "this" Java object
* @param mallocSize:    the size of allocated memory
*/
JNIEXPORT jlong JNICALL Java_at_tugraz_genome_lda_quantification_SavGolJNI_initMallocNative(JNIEnv *env,
	jobject thisObj,
	jint mallocSize)
{
	struct gh_values *ptr = (struct gh_values *)malloc(sizeof(*ptr));

	ptr->h_values = (float *)malloc(4 * (int)mallocSize * sizeof(float));
	gpuErrchk(cudaMalloc(&ptr->d_values, 4 * (int)mallocSize * sizeof(float)));
	gpuErrchk(cudaMalloc(&ptr->precalcPow, (int)mallocSize * sizeof(float)));

	int minGridSize;
	gpuErrchk(cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSizePowers,
		precalculatePowers, 0, 0));
	gpuErrchk(cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSizeSmooth,
		SmoothDataPoint, 0, 0));

	return (jlong)ptr;
}

/* =============================================================================
* Receives variables from the Java class SavGolJNI
*
* @param env:               JNI environment
* @param thisObj:           reference to "this" Java object
* @param j_values:          Contains time, raw data and smoothed data
* @param numberOfScans:	    Length of g_values
* @param range:             number of seconds around given points
* @param repeats:           how often the spectrum should be smoothed
* @param startSmoothScan:   start index of the spectrum
* @param stopSmoothScan:    end index of the spectrum
* @param address:           address of the device variables
*
*/
JNIEXPORT jfloatArray JNICALL Java_at_tugraz_genome_lda_quantification_SavGolJNI_SmoothNative(JNIEnv *env,
	jobject thisObj,
	jobjectArray j_value, jint numberOfScans, jfloat range, jint repeats,
	jint startSmoothScan, jint stopSmoothScan, jlong address)
{ 
	struct gh_values *ptr = (struct gh_values *)address;
	// -------------------------------------------------------------------------
	// Copies the values from Java to C

	// input variables
	int i, j;
	int len2 = 4;
	float* values;
	values = ptr->h_values; 

	// output variables
	jfloatArray j_smoothed = env->NewFloatArray(numberOfScans);
	jfloat *smoothed = env->GetFloatArrayElements(j_smoothed, NULL);

	// copying the values
	for (i = 0; i<numberOfScans; ++i)
	{
		jfloatArray oneDim = (jfloatArray)env->GetObjectArrayElement(j_value, i);
		jfloat *element = env->GetFloatArrayElements(oneDim, 0);

		for (j = 0; j<len2; ++j)
		{
			values[i*len2 + j] = element[j];
		}

		env->ReleaseFloatArrayElements(oneDim, element, JNI_ABORT);
		env->DeleteLocalRef(oneDim);
	}

	// -------------------------------------------------------------------------
	// Smoothing of the raw data

	prepareSmooth(values, (int)numberOfScans, (float)range, (int)repeats,
		(int)startSmoothScan, (int)stopSmoothScan, (long)address);
		
	// -------------------------------------------------------------------------
	// Copy data back to Java

	// copying data
	for (i = 0; i < numberOfScans; i++)
	{
		smoothed[i] = values[i*len2 + 2];
	}

	// release the smoothed array
	env->ReleaseFloatArrayElements(j_smoothed, (jfloat *)smoothed, 0);
	//free(values);

	// returns the smoothed values
	return j_smoothed;
}

/* =============================================================================
* Frees the Allocated space on the graphics card
*
* @param env:       JNI environment
* @param thisObj:   reference to "this" Java object
* @param address:   address of the device variables
*/
JNIEXPORT void JNICALL Java_at_tugraz_genome_lda_quantification_SavGolJNI_FreesNative(JNIEnv *env, jobject thisObj,
	jlong address)
{
	struct gh_values *ptr = (struct gh_values *)address;

	free(ptr->h_values);
	gpuErrchk(cudaFree(ptr->d_values));
	gpuErrchk(cudaFree(ptr->precalcPow));

	free(ptr);

	return;
}
