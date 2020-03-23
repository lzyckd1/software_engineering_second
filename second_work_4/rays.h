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
	bool operator ==(const rays& other) const
	{
		if (start == other.start) {
			if (a == other.a && b == other.b && c == other.c)
			{
				double x0 = start.x, y0 = start.y;
				double x1 = other.n.x, y1 = other.n.y;
				double x2 = n.x, y2 = n.y;
				if ((x1 > x0 && x2 > x0)||(x1 < x0&& x2 < x0)||(x1 == x0&& x2 == x0)) {
					if ((y1 > y0&& y2 > y0) || (y1 < y0 && y2 < y0) || (y1 == y0 && y2 == y0)) {
						return true;
					}
				}
			}
		}
		return false;
	}
};

