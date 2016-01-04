#include "individual.h"
#include "population.h"
#include <cstdlib>
#include <ctime>
#include <iostream>


Individual::Individual()
{
	bitString = 0;
	fitness = 0;
	//std::cout<<"success!"<<std::endl;
}

Individual::Individual(short n)
{
	srand((unsigned)time(0));
	//£¿£¿£¿maybe it's used to generate the individual randomly
	bitString = n*((1<<bitCount)/Population::popuScale) + rand()%((1<<bitCount)/Population::popuScale);
	fitness = bitString * bitString;
	//std::cout<<"success!"<<std::endl;
}

Individual::Individual(short bitStr, int fit)
{
	bitString = bitStr;
	fitness = fit;
	//std::cout<<"success!"<<std::endl;
}

short Individual::GetFitness()const
{
	return fitness;
}

short Individual::GetBitString()const
{
	return bitString;
}

bool Individual::operator >(const Individual &indi) const//to see whether the fitness of the object built by the user is bigger than than this variable
{
	if(this->fitness > indi.fitness)
		return true;
	else
		return false;
}

bool Individual::operator==(const Individual & indi)const//to see whether the fitness of the object built by the user is equal to this variable
{
	if(this->fitness == indi.fitness)
		return true;
	else
		return false;
}

std::ostream & operator<<(std::ostream &os, const Individual &indi)//output the bitstring and the fitness of the object
{
	os<<"bitString = "<< indi.bitString << " " << "fitness = " << indi.fitness << std::endl;
	return os;
}

Individual::~Individual()
{
	//std::cout<<"bye!"<<std::endl;
}

void Individual::ShowIndividual(void) const
{
	std::cout<<"bitString = "<<bitString<<" ";
	std::cout<<"fitness = "<<fitness<<std::endl;
}