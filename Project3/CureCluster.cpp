#include "CureCluster.h"

System::System()
{
}
System::~System()
{
	int i;
	for (i = 0; i<NUMPATTERNS; i++){
		if (p[i] == NULL)
			continue;
		else
			delete p[i];
	}
}

void System::LoadPatterns(double* data, int NumValid, int VectorSize)
{
	NUMPATTERNS = NumValid;
	SIZEVECTOR = VectorSize;
	NUMCLUSTERS = NUMPATTERNS*2/3;
	
	
	Pattern = data;
	p = new ClustNode*[NUMPATTERNS];

}


void  System::BuildClustList(){
	int i, j;
	double MinDist, d;
	for (i = 0; i<NUMPATTERNS; i++){
		p[i] = new ClustNode;
		p[i]->NumMembers = 1;
		p[i]->Member = new int[NUMPATTERNS];
		p[i]->Means = new double[SIZEVECTOR];
		for (int l = 0; l < NUMPRE; l++)
		{
			p[i]->Pre[l] = new double[SIZEVECTOR];
		}
		p[i]->Member[0] = i;
		p[i]->PreMember[0] = i;
		for (j = 0; j<SIZEVECTOR; j++){
			p[i]->Means[j] = Pattern[i*SIZEVECTOR + j];//均值
		}
		for (j = 0; j<SIZEVECTOR; j++){
			p[i]->Pre[0][j] = Pattern[i*SIZEVECTOR + j] + (SHRINK)*(p[i]->Means[j] - Pattern[i*SIZEVECTOR + j]);//计算代表点
		}
	}
	for (i = 0; i<NUMPATTERNS; i++){
		p[i]->MinDist = 9.9e+99;
		for (j = 0; j<NUMPATTERNS; j++){
			if (i == j)
				continue;
			else{
				d = dist(p[i], p[j]);
				if (d<p[i]->MinDist){
					p[i]->closet = j;//最近的簇和距离
					p[i]->MinDist = d;
				}
			}
		}
	}
}


double System::dist(ClustNode *u, ClustNode *v){
	double dist, MinDist;
	int i, j, MaxNum1, MaxNum2;
	MinDist = 9.9e+99;
	if (u->NumMembers<NUMPRE)
		MaxNum1 = u->NumMembers;
	else
		MaxNum1 = NUMPRE;
	if (v->NumMembers<NUMPRE)
		MaxNum2 = v->NumMembers;
	else
		MaxNum2 = NUMPRE;
	for (i = 0; i<MaxNum1; i++){
		for (j = 0; j<MaxNum2; j++){
			dist = EucNorm(u->Pre[i], v->Pre[j]);
			if (dist<MinDist)
				MinDist = dist;
		}
	}
	return MinDist;
}
double System::EucNorm(double *u, double *v){
	int i;
	double dist;
	dist = 0;
	for (i = 0; i<SIZEVECTOR; i++){
		dist += (*(u + i) - *(v + i))*(*(u + i) - *(v + i));
	}
	return dist;
}
void System::ShrinkPre(ClustNode * u){
	int i, j, m;
	if (u->NumMembers<NUMPRE){
		for (i = 0; i<u->NumMembers; i++){
			m = u->PreMember[i];
			for (j = 0; j<SIZEVECTOR; j++){
				u->Pre[i][j] = Pattern[m*SIZEVECTOR + j] + (SHRINK)*(u->Means[j] - Pattern[m*SIZEVECTOR + j]);
			}
		}
	}
	else{
		for (i = 0; i<NUMPRE; i++){
			m = u->PreMember[i];
			for (j = 0; j<SIZEVECTOR; j++){
				u->Pre[i][j] = Pattern[m*SIZEVECTOR + j] + (SHRINK)*(u->Means[j] - Pattern[m*SIZEVECTOR + j]);
			}
		}
	}
}

ClustNode * System::Merge(ClustNode *u, ClustNode *v){
	int i, j, MaxPoint;
	ClustNode *w = new ClustNode;
	w->Member = new int[NUMPATTERNS];
	w->Means = new double[SIZEVECTOR];
	for (int l = 0; l < NUMPRE; l++)
	{
		w->Pre[l] = new double[SIZEVECTOR];
	}
	double MaxDist, dist;
	w->NumMembers = u->NumMembers + v->NumMembers;
	for (i = 0; i<SIZEVECTOR; i++){
		w->Means[i] = (u->Means[i] * u->NumMembers + v->Means[i] * v->NumMembers) / w->NumMembers;
	}
	for (i = 0; i<u->NumMembers; i++){
		w->Member[i] = u->Member[i];
	}
	for (i = u->NumMembers; i<w->NumMembers; i++){
		w->Member[i] = v->Member[i - u->NumMembers];
	}
	if (w->NumMembers <= NUMPRE){
		for (i = 0; i<w->NumMembers; i++){
			w->PreMember[i] = w->Member[i];
		}
	}
	else{
		for (i = 0; i<NUMPRE; i++){
			MaxDist = 0;
			for (j = 0; j<w->NumMembers; j++){
				if (i == 0){
					if (j == i)
						continue;
					dist = EucNorm(w->Means, &Pattern[w->Member[j] * SIZEVECTOR]);
				}
				else{
					if (j == i - 1)
						continue;
					dist = EucNorm(w->Pre[i - 1], &Pattern[w->Member[j] * SIZEVECTOR]);
				}
				if (dist>MaxDist){
					MaxDist = dist;
					MaxPoint = j;
				}
			}
			w->PreMember[i] = j;
		}
	}
	ShrinkPre(w);
	return w;
}
int System::MinClust(){
	int i, Index;
	double MinDist;
	MinDist = 9.9e+99;
	for (i = 0; i<NUMPATTERNS; i++){
		if (p[i] == NULL)
			continue;
		if (p[i]->MinDist<MinDist){
			MinDist = p[i]->MinDist;
			Index = i;
		}
	}
	return Index;
}
void System::RunCure(){
	int i, j, k, u, v;
	double d, MinDist;
	ClustNode *w;
	BuildClustList();
	k = NUMPATTERNS;
	for (; k>NUMCLUSTERS; k--){
		u = MinClust();
		v = p[u]->closet;
		w = Merge(p[u], p[v]);
		delete p[u];
		delete p[v];
		p[u] = w;   //将u,v移除
		p[v] = NULL;
		w->closet = 0;//0表示任意
		for (i = 0; i<NUMPATTERNS; i++){
			if (p[i] == NULL)
				continue;
			MinDist = 9.9e+99;
			for (j = 0; j<NUMPATTERNS; j++){
				if (j == i)
					continue;
				if (p[j] == NULL)
					continue;
				d = dist(p[i], p[j]);
				if (d<MinDist){
					MinDist = d;
					p[i]->MinDist = d;
					p[i]->closet = j;
				}
			}
		}
	}


}

void System::ShowClusters(){
	int i, j = 1;
	for (i = 0; i<NUMPATTERNS; i++){
		if (p[i] == NULL)
			continue;
		else
		{
			cout << "the " << j << "th cluster is:" << endl;
			j++;
			for (int k = 0; k < p[i]->NumMembers; k++)
			{
				for (int m = 0; m < SIZEVECTOR; m++)
				{
					cout << Pattern[p[i]->Member[k] * SIZEVECTOR + m] << " ";
				}
				cout << endl;
			}
			cout << endl;
			cout << endl;
		}
	}
}


