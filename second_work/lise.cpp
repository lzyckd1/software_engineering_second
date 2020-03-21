#include "lise.h"
lise::lise()
{

}

lise::lise(node m, node n)
{
	double y1, y2, x1, x2;
	this->node1 = m;
	this->node2 = n;
	x1 = m.getX();
	y1 = m.getY();
	x2 = n.getX();
	y2 = n.getY();
	this->a = y1 - y2;
	this->b = x2 - x1;
	this->c = x1 * y2 - x2 * y1;
	this->k = 0;
	if (x1 == x2) {
		this->exitK = false;
	}
	else {
		this->exitK = true;
		this->k = (y2 - y1) / (x2 - x1);
	}
}

node lise::getNode1()
{
	return this->node1;
}

node lise::getNode2()
{
	return this->node2;
}

double lise::getA()
{
	return this->a;
}

double lise::getB()
{
	return this->b;
}

double lise::getC()
{
	return this->c;
}

double lise::getK()
{
	return this->k;
}

bool lise::getExitK()
{
	return this->exitK;
}
bool lise::judge(node n)
{
	node n1 = this->getNode1(), n2 = this->getNode2();
	double x1 = n1.getX(), x2 = n2.getX(), y1 = n1.getY(), y2 = n2.getY(), tmp;
	double x = n.getX(), y = n.getX();
	if (x1 > x2) {
		tmp = x1;
		x1 = x2;
		x2 = tmp;
	}
	if (y1 > y2) {
		tmp = y1;
		y1 = y2;
		y2 = tmp;
	}
	if (x <= x2 && x >= x1 && y <= y2 && y >= y1)
	{
		return true;
	}
	return false;
}