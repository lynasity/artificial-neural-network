#ifndef __POPULATION_H__
#define __POPULATION_H__

#include "individual.h"

#define CROSS_MUTATE_PRO_COMP 100 //交叉、变异的基数

class Population
{
public:
	static const short popuScale = 10; //种群规模
private:
	static const short pc = 40;  //交叉个体数相对于基数
	static const short pm = 5;   //变异个体数相对于基数
	Individual popuArray[popuScale];
public:
	Population();
	void Selection();  //选择
	void Crossover();  //交叉
	void Mutation();   //变异
	const Individual & GetbestIndividual()const;
	const Individual & GetRepresentationIndividual()const;
	friend std::ostream & operator<<(std::ostream &os, const Population & p);
	void ShowPopulation()const;
	~Population();
};

#endif