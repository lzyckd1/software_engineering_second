#include "Cycle.h"
#include "node.h"
Cycle::Cycle()
{

}
Cycle::Cycle(node n, double r)
{
	this->c = n;
	this->r = r;
}

node Cycle::getC()
{
	return this->c;
}

double Cycle::getR()
{
	return this->r;
}