#pragma once
#include "node.h"
class lise
{
private:
	node node1, node2;
	double a, b, c, k;
	bool exitK;
public:
	lise();
	lise(node m, node n);
	node getNode1();
	node getNode2();
	double getA();
	double getB();
	double getC();
	double getK();
	bool getExitK();
	bool judge(node n);
	bool operator ==(const struct lise& other) const
	{
		if (node1 == other.node1 && node2 == other.node2) return true;
		else if (node1 == other.node2 && node2 == other.node1)return true;
		else return false;
	}
};

