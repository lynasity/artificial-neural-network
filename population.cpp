#include "population.h"
#include <cstdlib>
#include <ctime>
#include <iostream>


//��Ⱥ���캯��
Population::Population()
{
	for(int i = 0; i < popuScale; i++)
	{
		popuArray[i] = Individual(i);//generate random population
		//popuArray[i].ShowIndividual();
	}

}

void Population::Selection()
{
	float selProbality[popuScale];
	float cumProbality[popuScale];
	int totalFitness = 0;
	for(int i = 0; i < popuScale; i++)
		totalFitness += popuArray[i].GetFitness();
	//std::cout<<"totalFitness:"<<totalFitness<<std::endl;
	for(int i = 0; i < popuScale; i++)
		selProbality[i] = float(popuArray[i].GetFitness())/float(totalFitness);
	cumProbality[0] = selProbality[0];
	for(int i = 1; i < popuScale; i++)
		cumProbality[i] = cumProbality[i-1] + selProbality[i];
	if(cumProbality[popuScale-1] < 1)
		cumProbality[popuScale-1] = 1;

	/*std::cout<<"selProbality:"<<std::endl;
	for(int i = 0; i < popuScale; i++)
	{
		std::cout<<i<<":";
		std::cout<<selProbality[i]<<std::endl;
	}

	std::cout<<"cumProbality:"<<std::endl;
	for(int i = 0; i < popuScale; i++)
	{
		std::cout<<i<<":";
		std::cout<<cumProbality[i]<<std::endl;
	}*/

	Individual tempIndividualArray[popuScale];//����ž���Ȼѡ������ĸ���

	srand(unsigned(time(0)));
	for(int i = 0; i < popuScale; i++)
	{
	    float temp = ((float)rand()) / (float)RAND_MAX;
		//std::cout<<"tempRand = "<<temp<<std::endl;
		int j = 0;
		for(; j < popuScale; j++)
		{
			if(temp <= cumProbality[j])
				break;
		    if(j >= popuScale)
			    std::cout<<"rand error!"<<std::endl;
		}
		tempIndividualArray[i] = popuArray[j];//j is the sign of the choosen one
	}
	for(int i = 0; i < popuScale; i++)
		popuArray[i]=tempIndividualArray[i];
}

void Population::Crossover()
{
	short crossCount = 0;
	if(short(pc * popuScale) % 2 == 0)
		crossCount = short(float(pc * popuScale) / float(2 * CROSS_MUTATE_PRO_COMP));
	else
		crossCount = short(float(pc * popuScale) / float( 2 * CROSS_MUTATE_PRO_COMP)) + 1;

	//std::cout<<"crossCount = "<<crossCount<<std::endl;

	Individual tempIndividualArray[popuScale];
	for(int i = 0; i < popuScale; i++)
		tempIndividualArray[i] = popuArray[i];

	srand((unsigned)time(0));
	for(int i = 0; i < crossCount; i++)
	{
		short temp1 = rand()%popuScale;
		short temp2 = rand()%popuScale;
		//std::cout<<"temp1 = "<<temp1<<" "<<"temp2 = "<<temp2<<std::endl;
		short crossDot = rand()%(Individual::bitCount + 1);
		short tempBitString1 = (popuArray[temp1].GetBitString() | ((1<<crossDot) - 1))
			&((popuArray[temp2].GetBitString() & ((1<<crossDot) - 1))|(((1<<Individual::bitCount)-1)<<crossDot));//my idea
		short tempBitString2 = (popuArray[temp2].GetBitString() | ((1<<crossDot) - 1))
			&((popuArray[temp1].GetBitString() & ((1<<crossDot) - 1))|(((1<<Individual::bitCount)-1)<<crossDot));
		//((popuArray[temp1].GetBitString() & ((1<<crossDot) - 1))|(((1<<Individual::bitCount)-1)>>crossDot)<<crossDot);why?
		//std::cout<<"tempBitString1 = "<<tempBitString1<<"tempBitString2 = "<<tempBitString2<<std::endl;
		tempIndividualArray[temp1] = Individual::Individual(tempBitString1,tempBitString1*tempBitString1);//ˢ��Ϊ��������������bitstrig
		tempIndividualArray[temp2] = Individual::Individual(tempBitString2,tempBitString2*tempBitString2);
	}

	for(int i = 0; i < popuScale; i++)
		popuArray[i] = tempIndividualArray[i];
}

void Population::Mutation(void)
{ 
	short mutateCount = short(float(pm) /(float) CROSS_MUTATE_PRO_COMP * popuScale * Individual::bitCount);//�ᷢ��ͻ���λ��
	//std::cout<<"mutateCount = "<< mutateCount << std::endl;

	Individual tempMutateIndividualArray[popuScale];
	for(int i = 0; i < popuScale; i++)
		tempMutateIndividualArray[i] = popuArray[i];
    
	srand((unsigned)time(0));
	for(int i = 0; i < mutateCount; i++)
	{
		short mutateIndividualIndex = rand()%popuScale;//ȷ����������ĸ���
		short mutateDot = 1 + rand()%Individual::bitCount;//ȷ�����������λ�������ҵڼ�λ��
		//std::cout<<"mutateIndividualIndex = "<< mutateIndividualIndex << " mutateDot = " << mutateDot << std::endl; 
		short tempBitString = ~(popuArray[mutateIndividualIndex].GetBitString() & (1<<(mutateDot - 1))) 
			& (popuArray[mutateIndividualIndex].GetBitString() | (1<<(mutateDot - 1)));

		tempMutateIndividualArray[mutateIndividualIndex] = Individual::Individual(tempBitString,tempBitString*tempBitString);
	}
   
	for(int i = 0; i < popuScale; i++)
	    popuArray[i] = tempMutateIndividualArray[i];
}

const Individual & Population::GetbestIndividual()const
{
	short maxflag = 0;
	for(int i = maxflag; i < (popuScale - 1); i++)
		if(popuArray[i+1] > popuArray[i])
			maxflag = i + 1;
	return popuArray[maxflag];
}
const Individual & Population::GetRepresentationIndividual()const//����ߴ����Եĸ��塣��representational
{
	short countArray[popuScale] = {0};

	for(int i = 0; i < popuScale; i++)
		for(int j = 0; j< popuScale; j++)
			if(popuArray[i] == popuArray[j])//�����==Ҳ�����ص�
				countArray[i]++;//�õ�ÿ�������ظ��ĸ��������ص���
	short maxflag = 0;//����Ĳ���ָfitness���ı�־������Ⱥ������Ӧ�����ձ����ĸ���
	for(int i = maxflag; i < (popuScale - 1); i++)
		if(countArray[i + 1] > countArray[i])
			maxflag = i + 1;

	return popuArray[maxflag];
}

std::ostream & operator<<(std::ostream & os, const Population & p)
{
	os<<"Individuals"<<std::endl;
	for(int i = 0; i < p.popuScale; i++)
		os<< i << ":"<< p.popuArray[i]<< std::endl;
	return os;
}
void Population::ShowPopulation(void)const
{
	std::cout<<"Individuals:"<<std::endl;
	for(int i = 0; i < popuScale; i++)
	{
		std::cout<<i<<":";
		popuArray[i].ShowIndividual();
	}
}
Population::~Population()
{

}