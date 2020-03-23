#include "rays.h"
rays::rays()
{

}
rays::rays(node m, node n)
{
	double y1, y2, x1, x2;
	this->start = m;
	this->n = n;
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
node rays::getStart()
{
	return this->start;
}
node rays::getN()
{
	return this->n;
}
double rays::getA()
{
	return this->a;
}
double rays::getB()
{
	return this->b;
}
double rays::getC()
{
	return this->c;
}
double rays::getK()
{
	return this->k;
}
bool rays::getExitK()
{
	return this->exitK;
}
bool rays::judge(node t)
{
	node n1 = this->getStart(), n2 = this->getN();
	double x1 = n1.getX(), x2 = n2.getX(), y1 = n1.getY(), y2 = n2.getY(), tmp;
	double x = t.getX(), y = t.getX();
	if (x1 <= x2) {
		if (x < x1)return false;
	}
	else {
		if (x > x1)return false;
	}
	if (y1 <= y2) {
		if (y < y1)return false;
	}
	else {
		if (y > y1)return false;
	}
	return true;
}