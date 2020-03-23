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
};