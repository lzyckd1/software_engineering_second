#pragma once
class node
{
public:
	double x, y;
	node();
	node(double x, double y);
	double getX();
	double getY();
	bool operator <(const node& other) const
	{
		if (x < other.x) return true;
		else if (x == other.x && y < other.y) return true;
		else return false;
	}
	bool operator ==(const node& other) const
	{
		if (x == other.x && y == other.y) return true;
		else return false;
	}
};

