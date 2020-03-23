#pragma once
class node
{
private:
	double x, y;
public:
	node();
	node(double x, double y);
	double getX();
	double getY();
	bool operator <(const struct node& other) const
	{
		if (x < other.x) return true;
		else if (x == other.x && y < other.y) return true;
		else return false;
	}
};

