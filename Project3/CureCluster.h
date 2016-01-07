/****************************************************************************
*                                                                           *
*  CURE                                                                     *
*                                                                           *
*****************************************************************************/

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string.h>
#include <math.h>
#include <conio.h>

using namespace std;

#define         SUCCESS         1
#define         FAILURE         0
#define         TRUE            1
#define         FALSE           0
//#define         NUMPATTERNS     10   // 数据点的个数
//#define         SIZEVECTOR      2     // 数据集的维数
//#define         NUMCLUSTERS     2     // 类的个数
#define         NUMPRE          3    // 代表点个数
#define         SHRINK          0.2   // 收缩率

// ***** 定义结构和类 *****
struct ClustNode {
	int          *Member;  //类中的数据项
	int          NumMembers;
	double       *Means;  //均值点
	double       *Pre[NUMPRE];//代表点
	int          PreMember[NUMPRE];//作为代表点的数据项
	int          closet;  //最近的簇
	double       MinDist;//与最近的簇之间的距离
};

class System {
private:
	int NUMPATTERNS = 0, SIZEVECTOR = 1, NUMCLUSTERS=2;
	double      *Pattern;
	ClustNode   **p;
	void        BuildClustList(); //建立初始簇链表
	ClustNode*  Merge(ClustNode *, ClustNode *); //合并两个簇
	int         MinClust();    //找到最近的两个簇
	double      dist(ClustNode*, ClustNode *);//两个簇之间的最短距离
	double      EucNorm(double *, double *);//代表点之间的距离
	void        ShrinkPre(ClustNode *);  //收缩代表点

public:
	System();
	~System();
	void LoadPatterns(double* data, int NumValid, int VectorSize);
	void RunCure();                     // 聚类的过程
	void ShowClusters();                // 显示聚类结果
};
