#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <highgui.h>
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
//#include <opencv2/opencv.hpp>
#include <opencv2/features2d/features2d.hpp>
#include "opencv2/nonfree/nonfree.hpp"
#include<opencv2/legacy/legacy.hpp>
#include<vector>
#define NTAPS_LUMA        8 ///< Number of taps for luma
#define NTAPS_CHROMA      4 ///< Number of taps for chroma
#define IF_INTERNAL_PREC 14 
#define IF_FILTER_PREC    6
#define IF_INTERNAL_OFFS (1<<(IF_INTERNAL_PREC-1))
#define AFIF_FILTER_PREC    8 ///< Log2 of sum of affine filter taps
#define LUMA_AFB16_FILTER_POSITIONS   64 ///< 64 phase
#define CHROMA_AFB16_FILTER_POSITIONS 64 ///< 64 phase

#define AFFINE_MIN_BLOCK_SIZE       1           ///< Minimum affine MC block size
#define AFFINE_MV_PRECISION         0.125       ///< MV precision for affine block
typedef       short          Pel;
using namespace cv;
using namespace std;

//int  bit_DepthY_ = 8;
long long MSEReferenceDirectly(Mat residue, int block_size, int CuPelX, int CuPelY);
double PSNRReferenceDirectly(Mat residue);

void xExtendPicCompBorder(Pel* piTxt, int iStride, int iWidth, int iHeight, int iMarginX, int iMarginY);

long long CalculateError(Mat Current, Pel* extendreference, int width, int height, int MarginX, int MarginY, Mat &pred, Mat H, int CuPelX, int CuPelY);

void UniFilterHorLuma(Pel *src, int srcStride, Pel *dst, int dstStride, int width, int height, int frac, bool isLast);

void filterCopy(int bitDepth, const Pel *src, int srcStride, Pel *dst, int dstStride, int width, int height, bool isFirst, bool isLast);

template<int N>
void UniFilterHor(int bitDepth, Pel *src, int srcStride, Pel *dst, int dstStride, int width, int height, bool isLast, short const *coeff);

template<int N, bool isVertical, bool isFirst, bool isLast>
void UniFilter(int bitDepth, const Pel *src, int srcStride, Pel *dst, int dstStride, int width, int height, short const *coeff);

void UniFilterVerLuma(Pel *src, int srcStride, Pel *dst, int dstStride, int width, int height, int frac, bool isFirst, bool isLast);

template<int N>
void UniFilterVer(int bitDepth, Pel *src, int srcStride, Pel *dst, int dstStride, int width, int height, bool isFirst, bool isLast, short const *coeff);

const short lumaAfB16Filter[(LUMA_AFB16_FILTER_POSITIONS)*NTAPS_LUMA] =
{
	0, 0, 0, 256, 0, 0, 0, 0,
	0, 1, -3, 256, 4, -2, 0, 0,
	0, 2, -7, 255, 8, -3, 1, 0,
	-1, 3, -10, 255, 12, -4, 1, 0,
	-1, 4, -13, 254, 16, -5, 2, -1,
	-1, 5, -16, 253, 20, -7, 2, 0,
	-1, 6, -18, 251, 25, -9, 3, -1,
	-2, 7, -21, 250, 29, -10, 4, -1,
	-2, 8, -23, 248, 34, -12, 4, -1,
	-2, 8, -25, 246, 38, -13, 5, -1,
	-2, 9, -27, 244, 43, -15, 5, -1,
	-2, 10, -30, 242, 48, -16, 6, -2,
	-2, 10, -31, 239, 52, -17, 5, 0,
	-2, 10, -32, 237, 57, -18, 6, -2,
	-2, 11, -34, 234, 63, -21, 7, -2,
	-2, 11, -35, 231, 68, -21, 6, -2,
	-3, 13, -38, 228, 74, -24, 9, -3,
	-2, 12, -38, 224, 78, -24, 7, -1,
	-3, 14, -40, 221, 84, -27, 10, -3,
	-2, 12, -39, 217, 88, -27, 8, -1,
	-3, 13, -40, 213, 94, -28, 9, -2,
	-3, 15, -43, 210, 100, -31, 11, -3,
	-3, 13, -41, 205, 104, -30, 9, -1,
	-3, 12, -41, 201, 110, -31, 9, -1,
	-3, 15, -43, 197, 116, -35, 12, -3,
	-3, 14, -43, 192, 121, -35, 12, -2,
	-2, 13, -42, 187, 126, -35, 10, -1,
	-3, 14, -43, 183, 132, -37, 12, -2,
	-2, 13, -42, 178, 137, -38, 12, -2,
	-3, 14, -42, 173, 143, -39, 12, -2,
	-3, 15, -43, 169, 148, -41, 14, -3,
	-3, 13, -41, 163, 153, -40, 13, -2,
	-3, 13, -40, 158, 158, -40, 13, -3,
	-2, 13, -40, 153, 163, -41, 13, -3,
	-3, 14, -41, 148, 169, -43, 15, -3,
	-2, 12, -39, 143, 173, -42, 14, -3,
	-2, 12, -38, 137, 178, -42, 13, -2,
	-2, 12, -37, 132, 183, -43, 14, -3,
	-1, 10, -35, 126, 187, -42, 13, -2,
	-2, 12, -35, 121, 192, -43, 14, -3,
	-3, 12, -35, 116, 197, -43, 15, -3,
	-1, 9, -31, 110, 201, -41, 12, -3,
	-1, 9, -30, 104, 205, -41, 13, -3,
	-3, 11, -31, 100, 210, -43, 15, -3,
	-2, 9, -28, 94, 213, -40, 13, -3,
	-1, 8, -27, 88, 217, -39, 12, -2,
	-3, 10, -27, 84, 221, -40, 14, -3,
	-1, 7, -24, 78, 224, -38, 12, -2,
	-3, 9, -24, 74, 228, -38, 13, -3,
	-2, 6, -21, 68, 231, -35, 11, -2,
	-2, 7, -21, 63, 234, -34, 11, -2,
	-2, 6, -18, 57, 237, -32, 10, -2,
	0, 5, -17, 52, 239, -31, 10, -2,
	-2, 6, -16, 48, 242, -30, 10, -2,
	-1, 5, -15, 43, 244, -27, 9, -2,
	-1, 5, -13, 38, 246, -25, 8, -2,
	-1, 4, -12, 34, 248, -23, 8, -2,
	-1, 4, -10, 29, 250, -21, 7, -2,
	-1, 3, -9, 25, 251, -18, 6, -1,
	0, 2, -7, 20, 253, -16, 5, -1,
	-1, 2, -5, 16, 254, -13, 4, -1,
	0, 1, -4, 12, 255, -10, 3, -1,
	0, 1, -3, 8, 255, -7, 2, 0,
	0, 0, -2, 4, 256, -3, 1, 0
};

double CalculateResidue(Mat &residue, Mat Current, int block_size, Mat pred, int CuPelX, int CuPelY);
void GetRotationandTranslation(Mat *H, int BlockNum, bool *IsHValid, Mat &R, Mat &t);
void GetRotationandTranslation(Mat H, Mat &R, Mat &t);
void DeriveTotalSSE(Mat* H, int NumBolckinWidth, int BlockNum, bool* IsHValid, double &TotalSSE);