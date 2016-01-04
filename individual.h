#ifndef __INDIVIDUAL_H__
#define __INDIVIDUAL_H__

#include <iostream>

class Individual
{
private:
	short bitString; //个体的二进制编码
	short fitness;   //个体的适应度
public:
	static const short bitCount = 5; //个体二进制码的长度
	Individual();
	Individual(short n);
	Individual(short bitStr,int fit);
	~Individual();
	short GetFitness() const;
	short GetBitString()const;
	//void setBitString();
	//void setFitness();
	bool operator>(const Individual & indi)const;
	bool operator==(const Individual & indi)const;
	friend std::ostream & operator<< (std::ostream & os, const Individual & indi);
	void ShowIndividual(void)const;
};
#endif