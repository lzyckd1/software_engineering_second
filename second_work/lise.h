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
	bool operator ==(const lise& other) const
	{
		if (node1 == other.node1 && node2 == other.node2) return true;
		else if (node1 == other.node2 && node2 == other.node1)return true;
		else if (a == other.a && b == other.b && c == other.c) {
			if ((node1.x - other.node1.x) * (node2.x - other.node1.x) < 0)return true;
			if ((node1.x - other.node2.x) * (node2.x - other.node2.x) < 0)return true;
			if ((node1.y - other.node1.y) * (node2.y - other.node1.y) < 0)return true;
			if ((node1.y - other.node2.y) * (node2.y - other.node2.y) < 0)return true;
			return false;
		}
		else return false;
	}
};

