#pragma once
#include "node.h"
class Cycle {
private:
	node c;
	double r;
public:
	Cycle();
	Cycle(node n, double r);
	node getC();
	double getR();
	bool operator ==(const struct Cycle& other) const
	{
		if (c == other.c && r == other.r) return true;
		else return false;
	}
};