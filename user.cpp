#include <iostream>
#include "individual.h"
#include "population.h"
#include<cstdio>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<string>
#include<algorithm>
using namespace std;
#define MAX 1000000
void TestIndividual(void); //个体类测试函数
void TestPopulation(void); //种群类测试函数
void GA(void);

//////////////////////////////////////////////////////////////////////////
//
//
//////////////////////////////////////////////////////////////////////////
//BP人工神经网络
//将三位二进制数转为一位十进制数

#include<iostream>
#include<cmath>
using namespace std;

#define  innode 3  //输入结点数
#define  hidenode 10//隐含结点数
#define  outnode 1 //输出结点数
#define  trainsample 7//BP训练样本数
class BpNet
{
public:
	void train(double p[trainsample][innode ],double t[trainsample][outnode]);//Bp训练
	double p[trainsample][innode];     //输入的样本
	double t[trainsample][outnode];    //样本要输出的

	double *recognize(double *p);//Bp识别

	void writetrain(); //写训练完的权值
	void readtrain(); //读训练好的权值，这使的不用每次去训练了，只要把训练最好的权值存下来就OK

	BpNet();
	virtual ~BpNet();

public:
	void init();
	double wh[innode][hidenode];//隐含结点权值
	double wo[hidenode][outnode];//输出结点权值
	double b1[hidenode];//隐含结点阀值
	double b2[outnode];//输出结点阀值

	double rate_w; //权值学习率（输入层-隐含层)
	double rate_w1;//权值学习率 (隐含层-输出层)
	double rate_b1;//隐含层阀值学习率
	double rate_b2;//输出层阀值学习率

	double e;//误差计算
	double error;//允许的最大误差
	double result[outnode];// Bp输出
};

BpNet::BpNet()
{
	error=1.0;
	e=0.0;

	rate_w=0.9;  //权值学习率（输入层--隐含层)//学习率介于0和1之间
	rate_w1=0.9; //权值学习率 (隐含层--输出层)
	rate_b1=0.9; //隐含层阀值学习率
	rate_b2=0.9; //输出层阀值学习率
}

BpNet::~BpNet()
{

}

void winit(double w[],int n) //权值初始化
{
	for(int i=0;i<n;i++)
		w[i]=(2.0*(double)rand()/RAND_MAX)-1;//权值范围[-1,1]
}

void BpNet::init()//随机初始化权值阈值 
{
	winit((double*)wh,innode*hidenode);
	winit((double*)wo,hidenode*outnode);
	winit(b1,hidenode);
	winit(b2,outnode);
}

void BpNet::train(double p[trainsample][innode],double t[trainsample][outnode])//第一个参数是输入量，第二个参数是输入量对应的期望系统输出的输出量
{
	double pp[hidenode];//隐含结点的校正误差
	double qq[outnode];//希望输出值与实际输出值的偏差
	double yd[outnode];//希望输出值

	double x[innode]; //输入向量
	double x1[hidenode];//隐含结点状态值
	double x2[outnode];//输出结点状态值
	double o1[hidenode];//隐含层激活值
	double o2[hidenode];//输出层激活值

	for(int isamp=0;isamp<trainsample;isamp++)//循环训练一次样品
	{
		for(int i=0;i<innode;i++)
			x[i]=p[isamp][i]; //输入的样本  每个样本的三个输入量
		for(int i=0;i<outnode;i++)
			yd[i]=t[isamp][i]; //每个样本的期望输出
		//构造每个样品的输入和输出标准
		for(int j=0;j<hidenode;j++)
		{
			o1[j]=0.0;//隐藏层激活值
			for(int i=0;i<innode;i++)
				o1[j]=o1[j]+wh[i][j]*x[i];//一个样本隐含层单元输入激活值=（每个输入变量*对应到隐藏结点通道的权值之和）
			    x1[j]=1.0/(1+exp(-o1[j]-b1[j]));//此样本情况下隐含层各细胞的状态值存入数组x[j]，b1是阀值，移项后相当于权值为-1的输入
                                                //等效成的激活值就是-o1[j]-b1[j]
			//    if(o1[j]+b1[j]>0) x1[j]=1;
			//else x1[j]=0;
		}

		for(int k=0;k<outnode;k++)
		{
			o2[k]=0.0;
			for(int j=0;j<hidenode;j++)
				o2[k]=o2[k]+wo[j][k]*x1[j]; //输出层各单元输入激活值
			x2[k]=1.0/(1.0+exp(-o2[k]-b2[k])); //输出层各单元输出，状态值
			//    if(o2[k]+b2[k]>0) x2[k]=1;
			//    else x2[k]=0;
		}

		for(int k=0;k<outnode;k++)
		{
			qq[k]=(yd[k]-x2[k])*x2[k]*(1-x2[k]);
			for(int j=0;j<hidenode;j++)
			wo[j][k]+=rate_w1*qq[k]*x1[j];//权值调整量=学习率*误差对权值导数,调整量与∂E/∂Wjk成正比			                          
		}                              
      //（梯度下降法，调整量与导数的值是符号相反的，因为当导数大于0，需要减少权值才能减少误差，反之则反）
		for(int j=0;j<hidenode;j++)
		{
			pp[j]=0.0;//校正误差
			for(int k=0;k<outnode;k++)
		    pp[j]=pp[j]+qq[k]*wo[j][k];
			pp[j]=pp[j]*x1[j]*(1-x1[j]);
			for(int i=0;i<innode;i++)
				wh[i][j]+=rate_w*pp[j]*x[i];
		}

		for(int k=0;k<outnode;k++)
		{
			e+=fabs(yd[k]-x2[k])*fabs(yd[k]-x2[k]); //计算均方差
		}
		error=e/2.0;//误差为0.5均方值
		for(int k=0;k<outnode;k++)
			b2[k]=b2[k]+rate_b2*qq[k]; //下一次的隐含层和输出层之间的新阈值，和修改权值的思路一样
		for(int j=0;j<hidenode;j++)    //阈值相当于权值为-b，而它对应的固定输入量为+1
			b1[j]=b1[j]+rate_b1*pp[j]; //下一次的输入层和隐含层之间的新阈值
	}
}

double *BpNet::recognize(double *p)//作用：由输入量求网络输出量
{
	double x[innode]; //输入向量
	double x1[hidenode]; //隐含结点状态值
	double x2[outnode]; //输出结点状态值
	double o1[hidenode]; //隐含层激活值
	double o2[hidenode]; //输出层激活值

	for(int i=0;i<innode;i++)
		x[i]=p[i];
	for(int j=0;j<hidenode;j++)
	{
		o1[j]=0.0;
		for(int i=0;i<innode;i++)
			o1[j]=o1[j]+wh[i][j]*x[i]; //隐含层各单元激活值
		x1[j]=1.0/(1.0+exp(-o1[j]-b1[j])); //隐含层各单元输出
		//if(o1[j]+b1[j]>0) x1[j]=1;
		//    else x1[j]=0;
	}

	for(int k=0;k<outnode;k++)
	{
		o2[k]=0.0;
		for(int j=0;j<hidenode;j++)
			o2[k]=o2[k]+wo[j][k]*x1[j];//输出层各单元激活值
		x2[k]=1.0/(1.0+exp(-o2[k]-b2[k]));//输出层各单元输出
		//if(o2[k]+b2[k]>0) x2[k]=1;
		//else x2[k]=0;
	}

	for(int k=0;k<outnode;k++)
	{
		result[k]=x2[k];
	}
	return result;
}

void BpNet::writetrain()
{
	FILE *stream0;
	FILE *stream1;
	FILE *stream2;
	FILE *stream3;
	int i,j;
	//隐含结点权值写入
	if(( stream0 = fopen("w.txt", "w+" ))==NULL)
	{
		cout<<"创建文件失败!";
		exit(1);
	}
	for(i=0;i<innode;i++)
	{
		for(j=0;j<hidenode;j++)
		{
			fprintf(stream0, "%f\n", wh[i][j]);
		}
	}
	fclose(stream0);

	//输出结点权值写入
	if(( stream1 = fopen("w1.txt", "w+" ))==NULL)
	{
		cout<<"创建文件失败!";
		exit(1);
	}
	for(i=0;i<hidenode;i++)
	{
		for(j=0;j<outnode;j++)
		{
			fprintf(stream1, "%f\n",wo[i][j]);
		}
	}
	fclose(stream1);

	//隐含结点阀值写入
	if(( stream2 = fopen("b1.txt", "w+" ))==NULL)
	{
		cout<<"创建文件失败!";
		exit(1);
	}
	for(i=0;i<hidenode;i++)
		fprintf(stream2, "%f\n",b1[i]);
	fclose(stream2);

	//输出结点阀值写入
	if(( stream3 = fopen("b2.txt", "w+" ))==NULL)
	{
		cout<<"创建文件失败!";
		exit(1);
	}
	for(i=0;i<outnode;i++)
		fprintf(stream3, "%f\n",b2[i]);
	fclose(stream3);

}

void BpNet::readtrain()
{
	FILE *stream0;
	FILE *stream1;
	FILE *stream2;
	FILE *stream3;
	int i,j;

	//隐含结点权值读出
	if(( stream0 = fopen("w.txt", "r" ))==NULL)
	{
		cout<<"打开文件失败!";
		exit(1);
	}
	float  wx[innode][hidenode];
	for(i=0;i<innode;i++)
	{
		for(j=0;j<hidenode;j++)
		{
			fscanf(stream0, "%f", &wx[i][j]);
			wh[i][j]=wx[i][j];
		}
	}
	fclose(stream0);

	//输出结点权值读出
	if(( stream1 = fopen("w1.txt", "r" ))==NULL)
	{
		cout<<"打开文件失败!";
		exit(1);
	}
	float  wx1[hidenode][outnode];
	for(i=0;i<hidenode;i++)
	{
		for(j=0;j<outnode;j++)
		{
			fscanf(stream1, "%f", &wx1[i][j]);
			wo[i][j]=wx1[i][j];
		}
	}
	fclose(stream1);

	//隐含结点阀值读出
	if(( stream2 = fopen("b1.txt", "r" ))==NULL)
	{
		cout<<"打开文件失败!";
		exit(1);
	}
	float xb1[hidenode];
	for(i=0;i<hidenode;i++)
	{
		fscanf(stream2, "%f",&xb1[i]);
		b1[i]=xb1[i];
	}
	fclose(stream2);

	//输出结点阀值读出
	if(( stream3 = fopen("b2.txt", "r" ))==NULL)
	{
		cout<<"打开文件失败!";
		exit(1);
	}
	float xb2[outnode];
	for(i=0;i<outnode;i++)
	{
		fscanf(stream3, "%f",&xb2[i]);
		b2[i]=xb2[i];
	}
	fclose(stream3);
}


//输入样本
double X[trainsample][innode]= {
	{0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0}
};
//期望输出样本
double Y[trainsample+1][outnode]={
	{0},{0.1429},{0.2857},{0.4286},{0.5714},{0.7143},{0.8571},{1.0000}

};

int MyBP_main()
{
	BpNet bp;
	bp.init();
	int times=0;
	while(bp.error>0.0001)
	{
		bp.e=0.0;
		times++;
		bp.train(X,Y);
		cout<<"Times="<<times<<" error="<<bp.error<<endl;
	}
	cout<<"trainning complete..."<<endl;
	double m[innode]={1,1,1};
	double *r=bp.recognize(m);//r=类成员数组result的地址
	for(int i=0;i<outnode;++i)
		cout<<bp.result[i]<<" ";
	double cha[trainsample][outnode];
	double mi=100;
	double index;
	for(int i=0;i<trainsample+1;i++)
	{
		for(int j=0;j<outnode;j++)
		{
			//找差值最小的那个样本，差值是期望值与学习后的值的差
			cha[i][j]=(double)(fabs(Y[i][j]-bp.result[j]));
			if(cha[i][j]<mi)
			{
				mi=cha[i][j];
				index=i;
			}
		}
	}
	for(int i=0;i<innode;++i)
		cout<<m[i];
	cout<<" is "<<index<<endl;
	cout<<endl;
	return 0;
}
//////////////////////////////////////////////////////////////////////////
int main()
{
	do 
	{

		std::cout<<"请选择要执行的算法 0——BP人工神经网络算法    1——遗传算法"<<std::endl;
		short gen = 0;
		std::cin>>gen;
		if (gen == 0)
		{
			MyBP_main();//人工神经网络 算法
		}
		else
		{
			GA(); //遗传算法
		}

	} while (true);
}


void GA(void)
{
	//std::cout<<std::fixed;
	//std::cout.precision(2);
	//std::cout.setf(std::ios_base::showpoint);
    
    
	std::cout<<"Please enter how many generations you want the population to evolve!"<<std::endl;
	short gen = 0;
	std::cin>>gen;

	Population birds;

	for(int i = 0; i < gen; i++)
	{
		std::cout<<"the "<< i << "th generation :"<< std::endl;
		std::cout<<birds<<std::endl;//运算符已经被重载，birds对象作为重载的运算符函数的输入参数

		birds.Selection();
		std::cout<<"after selection the population is : "<< std::endl;
		std::cout<<birds<<std::endl;

		birds.Crossover();
		std::cout<<"after crossover the population is :" << std::endl;
		std::cout<<birds<<std::endl;
		

		birds.Mutation();
		std::cout<<"after mutation the population is :" << std::endl;
		std::cout<<birds<<std::endl;
		

		std::cout<<std::endl<<std::endl;

	}
	std::cout<<"After "<< gen << " evolution the best fitness Individual bird is :"<<std::endl;
	std::cout<<birds.GetbestIndividual()<<std::endl;
	std::cout<<"After "<< gen << " evolution the represetational Individual bird is :"<<std::endl;
	std::cout<<birds.GetRepresentationIndividual()<<std::endl;
	//std::cout.clear();
	//std::cout.precision(4);
	//std::cout<<0.1;
	//TestIndividual();
	//TestPopulation();
	
	/*std::cout<<"Please enter Any  Num Key!"<<std::endl;
	short numV = 0;
	std::cin>>numV;*/
}

void TestIndividual(void)
{
	Individual one;
	one.ShowIndividual();
	int fit = one.GetFitness();
	int bit = one.GetBitString();
	std::cout<<fit<<std::endl<<bit<<std::endl;
}

void TestPopulation(void)
{
	std::cout<<"Constructor:"<<std::endl;
	Population one;
	one.ShowPopulation();

	std::cout<<std::endl<<"Selection:"<<std::endl;
	one.Selection();
	one.ShowPopulation();

	std::cout<<std::endl<<"Crossover:"<<std::endl;
	one.Crossover();
	one.ShowPopulation();

	std::cout<<std::endl<<"Mutation:"<<std::endl;
	one.Mutation();
	one.ShowPopulation();
}