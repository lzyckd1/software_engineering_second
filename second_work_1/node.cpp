#include "node.h"
node::node()
{

}
node::node(double x, double y)
{
	this->x = x;
	this->y = y;
}
double node :: getX()
{
	return x;
}

double node::getY()
{
	return y;
}