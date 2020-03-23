#pragma once
#include "node.h"
class line
{
private:
	node node1, node2;
	double a, b, c, k;
	bool exitK;
public:
	line();
	line(node m, node n);
	node getNode1();
	node getNode2();
	double getA();
	double getB();
	double getC();
	double getK();
	bool getExitK();
	bool judge(node n);
};

