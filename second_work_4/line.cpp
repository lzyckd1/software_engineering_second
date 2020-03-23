#include "line.h"
line::line()
{

}

line::line(node m, node n)
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

node line::getNode1()
{
	return this->node1;
}

node line::getNode2()
{
	return this->node2;
}

double line::getA()
{
	return this->a;
}

double line::getB()
{
	return this->b;
}

double line::getC()
{
	return this->c;
}

double line::getK()
{
	return this->k;
}

bool line::getExitK()
{
	return this->exitK;
}
bool line::judge(node n)
{
	return true;
}