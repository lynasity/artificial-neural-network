#ifndef __POPULATION_H__
#define __POPULATION_H__

#include "individual.h"

#define CROSS_MUTATE_PRO_COMP 100 //���桢����Ļ���

class Population
{
public:
	static const short popuScale = 10; //��Ⱥ��ģ
private:
	static const short pc = 40;  //�������������ڻ���
	static const short pm = 5;   //�������������ڻ���
	Individual popuArray[popuScale];
public:
	Population();
	void Selection();  //ѡ��
	void Crossover();  //����
	void Mutation();   //����
	const Individual & GetbestIndividual()const;
	const Individual & GetRepresentationIndividual()const;
	friend std::ostream & operator<<(std::ostream &os, const Population & p);
	void ShowPopulation()const;
	~Population();
};

#endif