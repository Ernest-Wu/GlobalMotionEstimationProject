
#include "RefineHomography.h"
#include "yuv.h"


int main(int argc, char **argv)
{
	int width;
	int height;
	FILE *fin = NULL;
	struct YUV_Capture cap;
	enum YUV_ReturnValue ret;
	IplImage *bgr, *image1 = NULL, *image2 = NULL;
	if (argc != 4)
	{
		fprintf(stderr, "usage: %s file.yuv width height\n", argv[0]);
		return 1;
	}
	width = atoi(argv[2]);
	height = atoi(argv[3]);
	if (width <= 0 || height <= 0)
	{
		fprintf(stderr, "error: bad frame dimensions: %d x %d\n", width, height);
		return 1;
	}
	fin = fopen(argv[1], "rb");
	if (!fin)
	{
		fprintf(stderr, "error: unable to open file: %s\n", argv[1]);
		return 1;
	}
	ret = YUV_init(fin, width, height, &cap);
	assert(ret == YUV_OK);

	bgr = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 3);
	assert(bgr);
	image1 = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);
	assert(image1);
	image2 = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 1);
	assert(image2);

	for (int poc = 0;; poc++)
	{

		ret = YUV_read(&cap);
		if (ret == YUV_EOF)
		{
			cvWaitKey(0);
			break;
		}
		else if (ret == YUV_IO_ERROR)
		{
			fprintf(stderr, "I/O error\n");
			break;
		}
		if (poc == 35)
		{
			cvCopy(cap.y, image1);
		}
		if (poc == 36)
		{
			cvCopy(cap.y, image2);
			if (image1 && image2)
			{
				cvShowImage("image1 before", image1);
				cvShowImage("image2 before", image2);
				Mat img_1(image1, 0), img_2(image2, 0);
				Mat H;

				SIFT sift1, sift2;
				vector<KeyPoint> key_points1, key_points2;
				Mat descriptors1, descriptors2, mascara;
				FlannBasedMatcher matcher;
				vector<DMatch>matches;
				vector< DMatch > good_matches;
				Mat img_matches;
				vector<Point2f> obj;
				vector<Point2f> scene;
				sift1(img_1, mascara, key_points1, descriptors1);
				sift2(img_2, mascara, key_points2, descriptors2);
				if (key_points1.size() != 0 || key_points2.size() != 0)
				{
					matcher.match(descriptors1,descriptors2,matches);
					double max_dist = 0, min_dist = 100;
					for (int i = 0; i < descriptors1.rows; i++)
					{
						double dist = matches[i].distance;
						if (dist<min_dist&&dist>0)min_dist = dist;
						if (dist>max_dist)max_dist = dist;
					}
					for (int i = 0; i < descriptors1.rows; i++)
					{
						if (matches[i].distance < 4 * min_dist)
						{
							good_matches.push_back(matches[i]);
						}
					}
					if (good_matches.size() >= 4)
					{
						drawMatches(img_1, key_points1, img_2, key_points2, good_matches, img_matches);
						imshow("Good Matches", img_matches);

						for (int i = 0; i < good_matches.size(); i++)
						{
							obj.push_back(key_points1[good_matches[i].queryIdx].pt);
							scene.push_back(key_points2[good_matches[i].trainIdx].pt);
						}

						H = findHomography(obj, scene, CV_RANSAC);
						Mat R(3, 3, CV_64F, Scalar(1));
						Mat t(3, 1, CV_64F, Scalar(1));
						GetRotationandTranslation(H, R, t);
						cvWaitKey(0);
					}
				}
			}
		}
		cvCvtColor(cap.ycrcb, bgr, CV_YCrCb2BGR);
		cvShowImage(argv[1], bgr);
		cvWaitKey(35);
	}




	return 0;
}
