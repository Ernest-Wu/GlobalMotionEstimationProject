#include "RefineHomography.h"
#include "CureCluster.h"       
long long MSEReferenceDirectly(Mat residue, int block_size, int CuPelX, int CuPelY)
{
/*	unsigned long long sum=0;
	for (int y = 0; y < CurrentBlock.rows; y++)
	{
		for (int x = 0; x < CurrentBlock.cols; x++)
		{
			sum += abs(CurrentBlock.data[x + y* CurrentBlock.cols] - ReferenceImage.data[(CuPelX + x) + (CuPelY + y)*ReferenceImage.cols]);// *(CurrentBlock.data[x + y* CurrentBlock.cols] - ReferenceImage.data[(CuPelX + x) + (CuPelY + y)*ReferenceImage.cols]));
		}
	}
	return double(sum)/(CurrentBlock.rows*CurrentBlock.cols);*/
	long long MSE = 0;
	for (int y = 0; y < block_size; y++)
	{
		for (int x = 0; x < block_size; x++)
		{
			MSE += residue.data[(CuPelX + x) + (CuPelY + y)* residue.cols] * residue.data[(CuPelX + x) + (CuPelY + y)* residue.cols]; // *(CurrentBlock.data[x + y* CurrentBlock.cols] - ReferenceImage.data[(CuPelX + x) + (CuPelY + y)*ReferenceImage.cols]));
		}
	}
	//return double(sum) / (block_size*block_size);
	return MSE;
}

double PSNRReferenceDirectly(Mat residue)
{
	double MSE = 0,PSNR;
	for (int y = 0; y < residue.rows; y++)
	{
		for (int x = 0; x < residue.cols; x++)
		{
			MSE += residue.data[x + y* residue.cols] * residue.data[x + y* residue.cols];
		}
	}
	long long tem = residue.cols*residue.rows;
	double t = double(255 * 255 *tem ) / MSE;
	PSNR=10 * log10(t);
	return PSNR;
}
long long CalculateError(Mat Current, Pel* extendreference, int width, int height, int MarginX, int MarginY, Mat &pred, Mat H, int CuPelX, int CuPelY)
{
	long long MSE1=0,MSE2=0;
	bool bi = false;


	Pel* dst = new Pel[pred.rows*pred.cols];
	Pel* prediction = dst;

	int refStride = width + (MarginX << 1);
	Pel* reference = extendreference + MarginX + refStride * MarginY;
	Pel *refOrg = reference + CuPelY*refStride + CuPelX;
	

	Pel* tmp = new Pel[(pred.rows + 16)*(pred.cols + 8)];
	int tmpStride = pred.rows + 16;

	int dstStride = pred.cols;
	int iBit = 6;
	int iMvShift = iBit;
	int iOffset = 8;
	//unsigned long long sum=0;
	int iMvScaleTmpHor, iMvScaleTmpVer;
	int filterSize = NTAPS_LUMA;
	int filterOrg = (filterSize >> 1) - 1;  //the original point to filter.
	
	int iHorMax = (width + iOffset - CuPelX - 1) << iMvShift;
	int iHorMin = (-pred.cols - iOffset - CuPelX + 1) << iMvShift;
	int iVerMax = (height + iOffset - CuPelY - 1) << iMvShift;
	int iVerMin = (-pred.rows - iOffset - CuPelY + 1) << iMvShift;

	double X, Y, Z;
	int x, y;
	int blockWidth = 1;
	int blockHeight = 1;
		double mv0[2], mv1[2], mv2[2], mv3[2];
	
	x = CuPelX; y = CuPelY;
    Z = 1. / (H.at<double>(2, 0)*double(x) + H.at<double>(2, 1)*double(y) + H.at<double>(2, 2));
	X = (H.at<double>(0, 0)*double(x) + H.at<double>(0, 1)*double(y) + H.at<double>(0, 2))*Z;
	Y = (H.at<double>(1, 0)*double(x) + H.at<double>(1, 1)*double(y) + H.at<double>(1, 2))*Z;
	mv0[0] = -( x - X);
	mv0[1] = -( y - Y);

	x = CuPelX + pred.cols; y = CuPelY;
	Z = 1. / (H.at<double>(2, 0)*double(x) + H.at<double>(2, 1)*double(y) + H.at<double>(2, 2));
	X = (H.at<double>(0, 0)*double(x) + H.at<double>(0, 1)*double(y) + H.at<double>(0, 2))*Z;
	Y = (H.at<double>(1, 0)*double(x) + H.at<double>(1, 1)*double(y) + H.at<double>(1, 2))*Z;
	mv1[0] = -( x - X);
	mv1[1] = -( y - Y);

	x = CuPelX; y = CuPelY + pred.rows;
	Z = 1. / (H.at<double>(2, 0)*double(x) + H.at<double>(2, 1)*double(y) + H.at<double>(2, 2));
	X = (H.at<double>(0, 0)*double(x) + H.at<double>(0, 1)*double(y) + H.at<double>(0, 2))*Z;
	Y = (H.at<double>(1, 0)*double(x) + H.at<double>(1, 1)*double(y) + H.at<double>(1, 2))*Z;
	mv2[0] = -(x - X);
	mv2[1] = -(y - Y);

	x = CuPelX + pred.cols; y = CuPelY+pred.rows;
	Z = 1. / (H.at<double>(2, 0)*double(x) + H.at<double>(2, 1)*double(y) + H.at<double>(2, 2));
	X = (H.at<double>(0, 0)*double(x) + H.at<double>(0, 1)*double(y) + H.at<double>(0, 2))*Z;
	Y = (H.at<double>(1, 0)*double(x) + H.at<double>(1, 1)*double(y) + H.at<double>(1, 2))*Z;
	mv3[0] = -(x - X);
	mv3[1] = -( y - Y);
	
	int mvWx = max(abs(mv0[0] - mv1[0]), abs(mv1[1] - mv0[1]));
	int mvWy = max(abs(mv0[0] - mv2[0]), abs(mv2[1] - mv0[1]));

	blockWidth = pred.cols;
	blockHeight = pred.rows;

	if (mvWx)
	{
		blockWidth = max((int)(pred.cols / mvWx * AFFINE_MV_PRECISION), 1);
		while (pred.cols % blockWidth)
		{
			blockWidth--;
		}
		blockWidth = max(AFFINE_MIN_BLOCK_SIZE, blockWidth);
	}
	if (mvWy)
	{
		blockHeight = max((int)(pred.rows / mvWy * AFFINE_MV_PRECISION), 1);
		while (pred.rows % blockHeight)
		{
			blockHeight--;
		}
		blockHeight = max(AFFINE_MIN_BLOCK_SIZE, blockHeight);
	}
	
	//cout << double(H.data[0]) << " " << double(H.data[1]) << " " << double(H.data[2]) << " " << endl
	//	<< double(H.data[3]) << " " << double(H.data[4]) << " " << double(H.data[5]) << " " << endl
	//	<< double(H.data[6]) << " " << double(H.data[7]) << " " << double(H.data[8]) << " " << endl;
	//cout << H.at<double>(0, 0) << " " << H.at<double>(0, 1) << " " << H.at<double>(0, 2) << " " << endl
	//	<< H.at<double>(1, 0) << " " << H.at<double>(1, 1) << " " << H.at<double>(1, 2) << " " << endl
	//	<< H.at<double>(2, 0) << " " << H.at<double>(2, 1) << " " << H.at<double>(2, 2) << " " << endl;
	for (int y = 0; y < pred.rows; y += blockHeight)
	{
		for (int x = 0; x < pred.cols; x += blockWidth)
		{
			Z = 1. / (H.at<double>(2, 0)*double(x + CuPelX) + H.at<double>(2, 1)*double(y + CuPelY) + H.at<double>(2, 2));
			X = (H.at<double>(0, 0)*double(x + CuPelX) + H.at<double>(0, 1)*double(y + CuPelY) + H.at<double>(0, 2))*Z;
			Y = (H.at<double>(1, 0)*double(x + CuPelX) + H.at<double>(1, 1)*double(y + CuPelY) + H.at<double>(1, 2))*Z;
			if (X >= width + MarginX-2*blockWidth || X<=-MarginX+2*blockWidth || Y >= height + MarginY-2*blockHeight || Y <= -MarginY+2*blockHeight)
			{
				return -1;
			}

			iMvScaleTmpHor = int((-CuPelX - x + X )*(1 << iBit) + 0.5);
			iMvScaleTmpVer = int((-CuPelY - y + Y )*(1 << iBit) + 0.5);
			//iMvScaleTmpHor = (int(-CuPelX - x + X + 0.5))*(1 << iBit);
			//iMvScaleTmpVer = (int(-CuPelY - y + Y + 0.5))*(1 << iBit);
			
			// Clip
			iMvScaleTmpHor = min(iHorMax, max(iHorMin, iMvScaleTmpHor));
			iMvScaleTmpVer = min(iVerMax, max(iVerMin, iMvScaleTmpVer));

			// get the MV in high precision
			int xFrac, yFrac, xInt, yInt, S = 1 << iMvShift;
			xInt = iMvScaleTmpHor >> iMvShift;
			yInt = iMvScaleTmpVer >> iMvShift;
			xFrac = iMvScaleTmpHor & (S - 1);
			yFrac = iMvScaleTmpVer & (S - 1);

			int refOffset = xInt + x + yInt * refStride;
			if (refOffset < (-(x + CuPelX + MarginX + (CuPelY + y + MarginY)*refStride)) || (refOffset>(width - blockWidth - NTAPS_LUMA - CuPelX - x + MarginX + (height - blockHeight - NTAPS_LUMA - CuPelY - y + MarginY)*refStride)))
			{
				return -1;
			}
			Pel *ref = refOrg + refOffset;

			if (!yFrac)
			{
				UniFilterHorLuma(ref, refStride, dst + x, dstStride, blockWidth, blockHeight, xFrac, !bi);
			}
			else if (!xFrac)
			{
				UniFilterVerLuma(ref, refStride, dst + x, dstStride, blockWidth, blockHeight, yFrac, true, !bi);
			}
			else
			{
				UniFilterHorLuma(ref - filterOrg*refStride, refStride, tmp, tmpStride, blockWidth, blockHeight + filterSize - 1, xFrac, false);
				UniFilterVerLuma(tmp + filterOrg*tmpStride, tmpStride, dst + x, dstStride, blockWidth, blockHeight, yFrac, false, !bi);
			}

		}
		dst += dstStride*blockHeight;
		refOrg += refStride*blockHeight;
	}
	for (int x = 0; x < pred.rows; x++)
	{
		for (int y = 0; y < pred.cols; y++)
		{
			pred.data[y + x*pred.cols] = ((prediction[y + x*pred.cols] < 0) ? 0 : (prediction[y + x*pred.cols] > 255) ? 255 : prediction[y + x*pred.cols]);
			MSE1 += (Current.data[(y + CuPelX) + (x + CuPelY)*Current.cols] - pred.data[y + x*pred.cols])*(Current.data[(y + CuPelX) + (x + CuPelY)*Current.cols] - pred.data[y + x*pred.cols]);// *(Current.at<unsigned char>(x, y) - pred.at<unsigned char>(x, y)));
		}
	}
      
	blockWidth = 1;
    blockHeight = 1;
	// double X, Y, Z;
	 for (int y = 0; y < pred.rows; y += blockHeight)
	{
		 for (int x = 0; x < pred.cols; x += blockWidth)
		{
			Z = 1. / (H.at<double>(2, 0)*double(x+CuPelX) + H.at<double>(2, 1)*double(y+CuPelY) + H.at<double>(2, 2));
			X = (H.at<double>(0, 0)*double(x + CuPelX) + H.at<double>(0, 1)*double(y + CuPelY) + H.at<double>(0, 2))*Z;
			Y = (H.at<double>(1, 0)*double(x + CuPelX) + H.at<double>(1, 1)*double(y + CuPelY) + H.at<double>(1, 2))*Z;
			//int x1 = int(X + 0.5);
			//int y1 = int(Y + 0.5);
			if (X>=width + MarginX || X<-MarginX || Y>=height + MarginY || Y < -MarginY)
			{
				return -1;
			}
			int x1 = int(X);
			int y1 = int(Y);
			double a = X - x1,b=Y-y1;

			prediction[x + y*pred.cols] = short(extendreference[x1 + MarginX + (y1 + MarginY)*refStride] * (1 - a)*(1 - b) + extendreference[x1 + 1 + MarginX + (y1 + MarginY)*refStride] * a*(1 - b) + extendreference[x1 + MarginX + (y1 + 1 + MarginY)*refStride] * (1 - a)*b + extendreference[x1 + 1 + MarginX + (y1 + 1 + MarginY)*refStride] * a*b);
		}
	}
	
	 for (int x = 0; x < pred.rows; x++)
	{
		 for (int y = 0; y < pred.cols; y++)
		{
			 pred.data[y + x*pred.cols] = ((prediction[y + x*pred.cols] < 0) ? 0 : (prediction[y + x*pred.cols] > 255) ? 255 : prediction[y + x*pred.cols]);
			 MSE2 += (Current.data[(y + CuPelX) + (x + CuPelY)*Current.cols] - pred.data[y + x*pred.cols])*(Current.data[(y + CuPelX) + (x + CuPelY)*Current.cols] - pred.data[y + x*pred.cols]);// *(Current.at<unsigned char>(x, y) - pred.at<unsigned char>(x, y)));
		}
	}

	 cout << MSE1 << " " << MSE2 << endl;
	 if (MSE1 < MSE2)
	 {
		 cout << "use 1/64 pixel-precision" << endl;
	 }
	 return MSE2 ;
}

/*void extendPicBorder()
{
	//if (m_bIsBorderExtended) return;

	xExtendPicCompBorder(getLumaAddr(), getStride(), getWidth(), getHeight(), m_iLumaMarginX, m_iLumaMarginY);
	//xExtendPicCompBorder(getCbAddr(), getCStride(), getWidth() >> 1, getHeight() >> 1, m_iChromaMarginX, m_iChromaMarginY);
	//xExtendPicCompBorder(getCrAddr(), getCStride(), getWidth() >> 1, getHeight() >> 1, m_iChromaMarginX, m_iChromaMarginY);

	//m_bIsBorderExtended = true;
}*/

void xExtendPicCompBorder(Pel* piTxt, int iStride, int iWidth, int iHeight, int iMarginX, int iMarginY)
{
	int   x, y;
	Pel*  pi;

	pi = piTxt;
	for (y = 0; y < iHeight; y++)
	{
		for (x = 0; x < iMarginX; x++)
		{
			pi[-iMarginX + x] = pi[0];
			pi[iWidth + x] = pi[iWidth - 1];
		}
		pi += iStride;
	}

	pi -= (iStride + iMarginX);
	for (y = 0; y < iMarginY; y++)
	{
		::memcpy(pi + (y + 1)*iStride, pi, sizeof(Pel)*(iWidth + (iMarginX << 1)));
	}

	pi -= ((iHeight - 1) * iStride);
	for (y = 0; y < iMarginY; y++)
	{
		::memcpy(pi - (y + 1)*iStride, pi, sizeof(Pel)*(iWidth + (iMarginX << 1)));
	}
}


void UniFilterHorLuma(Pel *src, int srcStride, Pel *dst, int dstStride, int width, int height, int frac, bool isLast)
{
	if (frac == 0)
	{
		filterCopy(8, src, srcStride, dst, dstStride, width, height, true, isLast);
	}
	else
	{
		UniFilterHor<NTAPS_LUMA>(8, src, srcStride, dst, dstStride, width, height, isLast, lumaAfB16Filter + frac*NTAPS_LUMA);
	}
}

void filterCopy(int bitDepth, const Pel *src, int srcStride, Pel *dst, int dstStride, int width, int height, bool isFirst, bool isLast)
{
	int row, col;

	if (isFirst == isLast)
	{
		for (row = 0; row < height; row++)
		{
			for (col = 0; col < width; col++)
			{
				dst[col] = src[col];
			}

			src += srcStride;
			dst += dstStride;
		}
	}
	else if (isFirst)
	{
		int shift = IF_INTERNAL_PREC - bitDepth;

		for (row = 0; row < height; row++)
		{
			for (col = 0; col < width; col++)
			{
				short val = src[col] << shift;
				dst[col] = val - (short)IF_INTERNAL_OFFS;
			}

			src += srcStride;
			dst += dstStride;
		}
	}
	else
	{
		int shift = IF_INTERNAL_PREC - bitDepth;
		short offset = IF_INTERNAL_OFFS;
		offset += shift ? (1 << (shift - 1)) : 0;
		short maxVal = (1 << bitDepth) - 1;
		short minVal = 0;
		for (row = 0; row < height; row++)
		{
			for (col = 0; col < width; col++)
			{
				short val = src[col];
				val = (val + offset) >> shift;
				if (val < minVal) val = minVal;
				if (val > maxVal) val = maxVal;
				dst[col] = val;
			}

			src += srcStride;
			dst += dstStride;
		}
	}
}

template<int N>
void UniFilterHor(int bitDepth, Pel *src, int srcStride, Pel *dst, int dstStride, int width, int height, bool isLast, short const *coeff)
{
	if (isLast)
	{
		UniFilter<N, false, true, true>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff);
	}
	else
	{
		UniFilter<N, false, true, false>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff);
	}
}

template<int N, bool isVertical, bool isFirst, bool isLast>
void UniFilter(int bitDepth, const Pel *src, int srcStride, Pel *dst, int dstStride, int width, int height, short const *coeff)
{
	int row, col;

	short c[8];
	c[0] = coeff[0];
	c[1] = coeff[1];
	if (N >= 4)
	{
		c[2] = coeff[2];
		c[3] = coeff[3];
	}
	if (N >= 6)
	{
		c[4] = coeff[4];
		c[5] = coeff[5];
	}
	if (N == 8)
	{
		c[6] = coeff[6];
		c[7] = coeff[7];
	}

	int cStride = (isVertical) ? srcStride : 1;
	src -= (N / 2 - 1) * cStride;

	int offset;
	short maxVal;
	int headRoom = IF_INTERNAL_PREC - bitDepth;
	int shift = AFIF_FILTER_PREC;   //affine filter precision

	if (isLast)
	{
		shift += (isFirst) ? 0 : headRoom;
		offset = 1 << (shift - 1);
		offset += (isFirst) ? 0 : (IF_INTERNAL_OFFS << (AFIF_FILTER_PREC));
		maxVal = (1 << bitDepth) - 1;
	}
	else
	{
		shift -= (isFirst) ? headRoom : 0;
		offset = (isFirst) ? -IF_INTERNAL_OFFS << shift : 0;
		maxVal = 0;
	}

	for (row = 0; row < height; row++)
	{
		for (col = 0; col < width; col++)
		{
			int sum;

			sum = src[col + 0 * cStride] * c[0];
			sum += src[col + 1 * cStride] * c[1];
			if (N >= 4)
			{
				sum += src[col + 2 * cStride] * c[2];
				sum += src[col + 3 * cStride] * c[3];
			}
			if (N >= 6)
			{
				sum += src[col + 4 * cStride] * c[4];
				sum += src[col + 5 * cStride] * c[5];
			}
			if (N == 8)
			{
				sum += src[col + 6 * cStride] * c[6];
				sum += src[col + 7 * cStride] * c[7];
			}

			// different from normal filter
			//sum = ( sum + 32 ) >> 6;

			short val = (sum + offset) >> shift;
			if (isLast)
			{
				val = (val < 0) ? 0 : val;
				val = (val > maxVal) ? maxVal : val;
			}
			dst[col] = val;
		}

		src += srcStride;
		dst += dstStride;
	}
}

void UniFilterVerLuma(Pel *src, int srcStride, Pel *dst, int dstStride, int width, int height, int frac, bool isFirst, bool isLast)
{
	if (frac == 0)
	{
		filterCopy(8, src, srcStride, dst, dstStride, width, height, isFirst, isLast);
	}
	else
	{
		UniFilterVer<NTAPS_LUMA>(8, src, srcStride, dst, dstStride, width, height, isFirst, isLast, lumaAfB16Filter + frac*NTAPS_LUMA);
	}
}

template<int N>
void UniFilterVer(int bitDepth, Pel *src, int srcStride, Pel *dst, int dstStride, int width, int height, bool isFirst, bool isLast, short const *coeff)
{
	if (isFirst && isLast)
	{
		UniFilter<N, true, true, true>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff);
	}
	else if (isFirst && !isLast)
	{
		UniFilter<N, true, true, false>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff);
	}
	else if (!isFirst && isLast)
	{
		UniFilter<N, true, false, true>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff);
	}
	else
	{
		UniFilter<N, true, false, false>(bitDepth, src, srcStride, dst, dstStride, width, height, coeff);
	}
}

double CalculateResidue(Mat &residue, Mat Current, int block_size, Mat pred, int CuPelX, int CuPelY)
{
	double MSE=0;
	for (int y = 0; y < block_size; y++)
	{
		for (int x = 0; x < block_size; x++)
		{
			residue.data[(CuPelX + x) + (CuPelY + y)*residue.cols] = abs(Current.data[(CuPelX + x) + (CuPelY + y)*  Current.cols] - pred.data[x + y* block_size]);
			MSE += residue.data[(CuPelX + x) + (CuPelY + y)*residue.cols] * residue.data[(CuPelX + x) + (CuPelY + y)*residue.cols];
		}
	}
	return MSE;
	//cout << MSE <<endl;
/*	if (CuPelX == 0 && CuPelY == 0)
	{
		for (int y = 0; y < block_size; y++)
		{
			for (int x = 0; x < block_size; x++)
			{
				cout << int(residue.data[(CuPelX + x) + (CuPelY + y)*residue.cols] )<< " ";

			}
			cout << endl;
		}
	}*/
}

void DeriveTotalSSE(Mat* H,int NumBolckinWidth, int BlockNum, bool* IsHValid, double &TotalSSE)
{
	for (int h = 0; h < 2 * BlockNum; h+=2)
	{
		if (IsHValid[h] == false&&IsHValid[h+1]==false)
		{
			if (h == 0)
			{
				if ((IsHValid[2] || IsHValid[3]) && (IsHValid[2 + 2 * NumBolckinWidth] || IsHValid[3 + 2 * NumBolckinWidth]))
				{

				}
			}
		}
	}

}
void GetRotationandTranslation(Mat *H,int BlockNum,bool *IsHValid, Mat &R, Mat &t )
{
	Mat H0 = Mat(3, 3, CV_64F, Scalar(1));
	Mat *T1 = new Mat[BlockNum*2];
	Mat *R1 = new Mat[BlockNum*2];
	double *d1 = new double[BlockNum * 2];
	Mat *n1 = new Mat[BlockNum * 2];
	int negativeNum[3] = { 0, 0, 0 };
	Mat W = Mat(3, 3, CV_64F, Scalar(1));
	Mat U = Mat(3, 3, CV_64F, Scalar(1));
	Mat VT = Mat(3, 3, CV_64F, Scalar(1));
	//cout << H[0].at<double>(0, 0) << H[0].at<double>(0, 1) << endl;
	cout <<"H"<< endl <<H[0].at<double>(0, 0) << " " << H[0].at<double>(0, 1) << " " << H[0].at<double>(0, 2) << " " << endl
		<< H[0].at<double>(1, 0) << " " << H[0].at<double>(1, 1) << " " << H[0].at<double>(1, 2) << " " << endl
		<< H[0].at<double>(2, 0) << " " << H[0].at<double>(2, 1) << " " << H[0].at<double>(2, 2) << " " << endl;
	for (int h = 0; h < 2 * BlockNum; h++)
	{
		if (IsHValid[h]==true)
		{
			H0.at<double>(0, 0) = H[h].at<double>(0, 0) / H[h].at<double>(2, 2);
			H0.at<double>(0, 1) = H[h].at<double>(0, 1) / H[h].at<double>(2, 2);
			H0.at<double>(0, 2) = H[h].at<double>(0, 2) / H[h].at<double>(2, 2);
			H0.at<double>(1, 0) = H[h].at<double>(1, 0) / H[h].at<double>(2, 2);
			H0.at<double>(1, 1) = H[h].at<double>(1, 1) / H[h].at<double>(2, 2);
			H0.at<double>(1, 2) = H[h].at<double>(1, 2) / H[h].at<double>(2, 2);
			H0.at<double>(2, 0) = H[h].at<double>(2, 0) / H[h].at<double>(2, 2);
			H0.at<double>(2, 1) = H[h].at<double>(2, 1) / H[h].at<double>(2, 2);
			H0.at<double>(2, 2) = 1.0;

			cout << "H0" << H0 << endl;

			cv::SVDecomp(H0, W, U, VT, SVD::FULL_UV);

			Mat u2 = U(Rect(0, 0, 1, 3)), u3 = U(Rect(2, 0, 1, 3));
			double sigma2 = W.at<double>(0, 0), sigma1 = W.at<double>(1, 0), sigma3 = W.at<double>(2, 0);
			double lamda = 1.0 / sigma1;
			double rou2 = lamda*sigma2, rou3 = lamda*sigma3;
			double eta1 = 1, eta2 = rou2*rou2, eta3 = rou3*rou3;
			double k = rou2 - rou3, p = rou2*rou3 - 1;
			
			Mat H0T=H0.t();
			Mat M = lamda*lamda*H0T*H0,ValuesMat,VectorsMat;
			eigen(M, ValuesMat, VectorsMat);

			Mat v2 = VectorsMat(Rect(0, 0, 1, 3)), v3 = VectorsMat(Rect(2, 0, 1, 3));

			double u = norm(v2, CV_L2), deta = norm(v3, CV_L2);

			if (p != -1)
			{
				double afa = (-k + sqrt(k*k + 4 * (p + 1))) / (2 * k*(p + 1));
				double sita = (-k - sqrt(k*k + 4 * (p + 1))) / (2 * k*(p + 1));
				Mat t0 = -(u*u2 - deta*u3) / (afa - sita);
				Mat n = -(afa*deta*u3 - sita*u*u2) / (afa - sita);

				T1[h] = t0;
				if (t0.at<double>(0, 0) < 0)
				{
					negativeNum[0]++;
				}
				if (t0.at<double>(1, 0) < 0)
				{
					negativeNum[1]++;
				}
				if (t0.at<double>(2, 0) < 0)
				{
					negativeNum[2]++;
				}
				n1[h] = n;
				Mat I = Mat(3, 3, CV_64F,Scalar(0));
				I.at<double>(0, 0) = 1;
				I.at<double>(1, 1) = 1;
				I.at<double>(2, 2) = 1;
				R1[h] = lamda*H0*((I + t0*n.t()).inv());
	
			}
		}
	}
	int maxIndex = 0;
	for (int i = 1; i < 3; i++)
	{
		if (negativeNum[i]>negativeNum[maxIndex])
		{
			maxIndex = i;
		}
	}
	
	int firstValid = 0,NumValid=0;
	for (int h = 0; h < 2 * BlockNum; h++)
	{
		if (IsHValid[h])
		{
			firstValid = h;
			break;
		}
	}
	for (int h = 0; h < 2 * BlockNum; h++)
	{
		if (IsHValid[h])
		{
			NumValid++;
			double d = -1 / T1[h].at<double>(maxIndex, 0);
			T1[h].at<double>(0, 0) = -T1[h].at<double>(0, 0)*d;
			T1[h].at<double>(1, 0) = -T1[h].at<double>(1, 0)*d;
			T1[h].at<double>(2, 0) = -T1[h].at<double>(2, 0)*d;

			d1[h] = d;


			cout << "R " << R1[h] << endl
				<< "t " << T1[h] << endl
				<< "n " << n1[h] << endl
				<< "d " << d1[h] << endl << endl;
		}

	}
	double *T2 = new double[NumValid * 2];
	double *R2 = new double[NumValid * 9];
	for (int h = 0, k1 = 0, k2 = 0 ; h < 2 * BlockNum; h++)
	{
		if (IsHValid[h])
		{
			for (int i = 0; i < 3; i++)
			{
				if (i != maxIndex)
				{
					T2[k1++] = T1[h].at<double>(i, 0);
				}
			}
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					R2[k2++] = R1[h].at<double>(i, j);
				}
			}
			for (int j = 0; j < 2; j++)
			{
				R2[k2++] = R1[h].at<double>(2, j);
			}
		}
	}

	System cure;
	cure.LoadPatterns(T2, NumValid, 2);
	cure.RunCure();
	cure.ShowClusters();
	//getch();
	System cure1;
	cure1.LoadPatterns(R2, NumValid, 8);
	cure1.RunCure();
	cure1.ShowClusters();
}

void GetRotationandTranslation(Mat H, Mat &R, Mat &t)
{
	Mat H0 = Mat(3, 3, CV_64F, Scalar(1));
	Mat n;
	double d;
	Mat W = Mat(3, 3, CV_64F, Scalar(1));
	Mat U = Mat(3, 3, CV_64F, Scalar(1));
	Mat VT = Mat(3, 3, CV_64F, Scalar(1));

	H0.at<double>(0, 0) = H.at<double>(0, 0) / H.at<double>(2, 2);
	H0.at<double>(0, 1) = H.at<double>(0, 1) / H.at<double>(2, 2);
	H0.at<double>(0, 2) = H.at<double>(0, 2) / H.at<double>(2, 2);
	H0.at<double>(1, 0) = H.at<double>(1, 0) / H.at<double>(2, 2);
	H0.at<double>(1, 1) = H.at<double>(1, 1) / H.at<double>(2, 2);
	H0.at<double>(1, 2) = H.at<double>(1, 2) / H.at<double>(2, 2);
	H0.at<double>(2, 0) = H.at<double>(2, 0) / H.at<double>(2, 2);
	H0.at<double>(2, 1) = H.at<double>(2, 1) / H.at<double>(2, 2);
	H0.at<double>(2, 2) = 1.0;

	cout << "H0" << H0 << endl;

	cv::SVDecomp(H0, W, U, VT, SVD::FULL_UV);
	Mat V = VT.t();
	Mat u2 = V(Rect(0, 0, 1, 3)), u3 = V(Rect(2, 0, 1, 3));
	double sigma2 = W.at<double>(0, 0), sigma1 = W.at<double>(1, 0), sigma3 = W.at<double>(2, 0);
	double lamda = 1.0 / sigma1;
	double rou2 = lamda*sigma2, rou3 = lamda*sigma3;
	double eta1 = 1, eta2 = rou2*rou2, eta3 = rou3*rou3;
	double k = rou2 - rou3, p = rou2*rou3 - 1;

	Mat H0T = H0.t();
	Mat M = lamda*lamda*H0T*H0, ValuesMat, VectorsMat;
	eigen(M, ValuesMat, VectorsMat);
	Mat v2 = VectorsMat(Rect(0, 0, 1, 3)), v3 = VectorsMat(Rect(2, 0, 1, 3));

	//double u = norm(v2, CV_L2), deta = norm(v3, CV_L2);

	if (p != -1)
	{
		double afa = (-k + sqrt(k*k + 4 * (p + 1))) / (2 * k*(p + 1));
		double sita = (-k - sqrt(k*k + 4 * (p + 1))) / (2 * k*(p + 1));

		double u = sqrt(afa*afa*k*k + 2 * afa*p + 1);
		double deta = sqrt(sita*sita*k*k + 2 * sita*p + 1);
		Mat t0 = -(u*u2 - deta*u3) / (afa - sita);
		n = -(afa*deta*u3 - sita*u*u2) / (afa - sita);
		cout << n << endl;
		cout << norm(t0, CV_L2) << endl;
		cout << t0.t()*t0 << endl;
		cout << k*k << endl;
		cout << norm(n, CV_L2) << endl;
		if (n.at<double>(2, 0) < 0)
		{
			n = (afa*deta*u3 - sita*u*u2) / (afa - sita);
			t0 = (u*u2 - deta*u3) / (afa - sita);
		}
		
		Mat I = Mat(3, 3, CV_64F, Scalar(0));
		I.at<double>(0, 0) = 1;
		I.at<double>(1, 1) = 1;
		I.at<double>(2, 2) = 1;
		R = lamda*H0*((I + t0*n.t()).inv());
		cout << R << endl;
		cout << R.inv() << endl;
		cout << R*R.t() << endl;
		cout << determinant(R) << endl;
		int i;
		for (i = 0; i < 3; i++)
		{
			if (t0.at<double>(i, 0)<0)
			{
				d = -1 / t0.at<double>(i, 0);
				break;
			}
		}
		if (i == 3)
		{
			d = -1 / t0.at<double>(0, 0);
		}
		
		for (int j = 0; j < 3; j++)
		{
			if (j != i)
			{
				t.at<double>(j, 0) = -t0.at<double>(j, 0)*d;
			}
			else
			{
				t.at<double>(j, 0) = 1.0;
			}
		}


		cout << "R " << R << endl
			<< "t " << t << endl
			<< "n " << n << endl
			<< "d " << d << endl << endl;
	}
}