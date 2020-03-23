#pragma once
#include "node.h"
class rays
{
private:
	node start, n;
	double a, b, c, k;
	bool exitK;
public:
	rays();
	rays(node m, node n);
	node getStart();
	node getN();
	double getA();
	double getB();
	double getC();
	double getK();
	bool getExitK();
	bool judge(node n);
	bool operator ==(const struct rays& other) const
	{
		if (start == other.start && n == other.n) return true;
		else return false;
	}
};

