#include "node.h"
#include "Cycle.h"
#include "line.h"
#include "lise.h"
#include "rays.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <map>
#include <string>
#include <vector>
#include <math.h>

using namespace std;
const int maxn = 1000;
map <node, int> node_set;
vector <line> lseq;
vector <lise> sseq;
vector <Cycle> cseq;
vector <rays> rseq;
int countn = 0;

void error_output(string errors)
{
	try {
		throw errors;
	}
	catch (const char* msg) {
		cerr << msg << endl;
	}
}

int string2int(string l)
{
	int ans = 0, op = 1;
	if (l[0] == '-')op = -1;
	else if (l[0] <= '9' && l[0] >= '0')ans += l[0] - '0';
	else {
		try {
			throw "wrong input!\n";
		}
		catch (const char* msg) {
			cerr << msg << endl;
		}
	}
	for (int i = 1; i < l.size(); i++)
	{
		if (!(l[i] <= '9' && l[i] >= '0'))
		{
			try {
				throw "wrong input!\n";
			}
			catch (const char* msg) {
				cerr << msg << endl;
			}
		}
		ans = ans * 10 + l[i] - '0';
	}
	ans *= op;
	return ans;
}

double string2double(string l)
{
	int ans = 0, op = 1;
	if (l[0] == '-')op = -1;
	else if (l[0] <= '9' && l[0] >= '0')ans += l[0] - '0';
	for (int i = 1; i < l.size(); i++)
	{
		if (!(l[i] <= '9' && l[i] >= '0'))
		{
			try {
				throw "wrong input!\n";
			}
			catch (const char* msg) {
				cerr << msg << endl;
			}
		}
		ans = ans * 10 + (l[i] - '0');
		//printf("%d\n", i);
	}
	ans *= op;
	return (double)ans;
}

char string2char(string l)
{
	if (l.size() != 1) {
		try {
			throw "wrong input!\n";
		}
		catch (const char* msg) {
			cerr << msg << endl;
		}
	}
	return l[0];
}

double distance_nn(node a, node b)
{
	//两点间距离的平方
	double ans, x = a.getX() - b.getX(), y = a.getY() - b.getY();
	ans = x*x+y*y;
	return ans;
}

double distance_nl(node a, line b)
{
	//点到直线距离的平方
	double ans, A = b.getA(), B = b.getB(), C = b.getC();
	ans = A*a.getX()+B*a.getY()+C;
	ans *= ans;
	ans /= A * A + B * B;
	return ans;
}

double distance_nl(node a, lise b)
{
	//点到直线距离的平方
	double ans, A = b.getA(), B = b.getB(), C = b.getC();
	ans = A * a.getX() + B * a.getY() + C;
	ans *= ans;
	ans /= A * A + B * B;
	return ans;
}

double distance_nl(node a, rays b)
{
	//点到直线距离的平方
	double ans, A = b.getA(), B = b.getB(), C = b.getC();
	ans = A * a.getX() + B * a.getY() + C;
	ans *= ans;
	ans /= A * A + B * B;
	return ans;
}

void node_set_insert(node n)
{
	if (node_set.find(n) == node_set.end()) {
		pair <node, int> p(n,1);
		node_set.insert(p);
		countn++;
		return;
	}
	else {
		node_set[n] = node_set[n] + 1;
		return;
	}
}

void node_set_delete(node n) 
{
	if (node_set.find(n) == node_set.end()) {
		
		return;
	}
	else {
		node_set[n] = node_set[n] - 1;
		if (node_set[n] == 0) {
			node_set.erase(n);
			countn--;
		}
		return;
	}
}

node findIntersectionll(line a, line b)
{
	double x = 0, y = 0, d;
	double a1 = a.getA(), b1 = a.getB(), c1 = a.getC();
	double a2 = b.getA(), b2 = b.getB(), c2 = b.getC();
	d = a1 * b2 - a2 * b1;
	if (d != 0) {
		x = (b1 * c2 - b2 * c1) / d;
		y = (a2 * c1 - a1 * c2) / d;
	}
	return node(x,y);
}

node findIntersectionlr(line a, rays b)
{
	double x = 0, y = 0, d;
	double a1 = a.getA(), b1 = a.getB(), c1 = a.getC();
	double a2 = b.getA(), b2 = b.getB(), c2 = b.getC();
	d = a1 * b2 - a2 * b1;
	if (d != 0) {
		x = (b1 * c2 - b2 * c1) / d;
		y = (a2 * c1 - a1 * c2) / d;
	}
	return node(x, y);
}

node findIntersectionrr(rays a, rays b)
{
	double x = 0, y = 0, d;
	double a1 = a.getA(), b1 = a.getB(), c1 = a.getC();
	double a2 = b.getA(), b2 = b.getB(), c2 = b.getC();
	d = a1 * b2 - a2 * b1;
	if (d != 0) {
		x = (b1 * c2 - b2 * c1) / d;
		y = (a2 * c1 - a1 * c2) / d;
	}
	return node(x, y);
}

node findIntersectionls(line a, lise b)
{
	double x = 0, y = 0, d;
	double a1 = a.getA(), b1 = a.getB(), c1 = a.getC();
	double a2 = b.getA(), b2 = b.getB(), c2 = b.getC();
	d = a1 * b2 - a2 * b1;
	if (d != 0) {
		x = (b1 * c2 - b2 * c1) / d;
		y = (a2 * c1 - a1 * c2) / d;
	}
	return node(x, y);
}

node findIntersectionrs(rays a, lise b)
{
	double x = 0, y = 0, d;
	double a1 = a.getA(), b1 = a.getB(), c1 = a.getC();
	double a2 = b.getA(), b2 = b.getB(), c2 = b.getC();
	d = a1 * b2 - a2 * b1;
	if (d != 0) {
		x = (b1 * c2 - b2 * c1) / d;
		y = (a2 * c1 - a1 * c2) / d;
	}
	return node(x, y);
}

node findIntersectionss(lise a, lise b)
{
	double x = 0, y = 0, d;
	double a1 = a.getA(), b1 = a.getB(), c1 = a.getC();
	double a2 = b.getA(), b2 = b.getB(), c2 = b.getC();
	d = a1 * b2 - a2 * b1;
	if (d != 0) {
		x = (b1 * c2 - b2 * c1) / d;
		y = (a2 * c1 - a1 * c2) / d;
	}
	return node(x, y);
}

void add_node_ll(line a, line b)
{
	node ans;
	if ((!a.getExitK()) && (!b.getExitK()))//a,b均平行于y轴
	{
		return;
	}
	else if (a.getExitK() && b.getExitK()) //a,b均存在斜率
	{
		if (a.getK() == b.getK())return;
		ans = findIntersectionll(a, b);
		node_set_insert(ans);
		return;
	}
	else {
		ans = findIntersectionll(a, b);
		node_set_insert(ans);
		return;
	}
	return;
}

void add_node_lr(line a, rays b)
{
	node ans;
	if ((!a.getExitK()) && (!b.getExitK()))//a,b均平行于y轴
	{
		return;
	}
	else if (a.getExitK() && b.getExitK()) //a,b均存在斜率
	{
		if (a.getK() == b.getK())return;
		ans = findIntersectionlr(a, b);
		if (b.judge(ans) == false)return;
		node_set_insert(ans);
		return;
	}
	else {
		ans = findIntersectionlr(a, b);
		if (b.judge(ans) == false)return;
		node_set_insert(ans);
		return;
	}
	return;
}

void add_node_rr(rays a, rays b)
{
	node ans;
	if ((!a.getExitK()) && (!b.getExitK()))//a,b均平行于y轴
	{
		return;
	}
	else if (a.getExitK() && b.getExitK()) //a,b均存在斜率
	{
		if (a.getK() == b.getK())return;
		ans = findIntersectionrr(a, b);
		if (b.judge(ans) == false || a.judge(ans) == false)return;
		node_set_insert(ans);
		return;
	}
	else {
		ans = findIntersectionrr(a, b);
		if (b.judge(ans) == false || a.judge(ans) == false)return;
		node_set_insert(ans);
		return;
	}
	return;
}

void add_node_ls(line a, lise b)
{
	node ans;
	if ((!a.getExitK()) && (!b.getExitK()))//a,b均平行于y轴
	{
		return;
	}
	else if (a.getExitK() && b.getExitK()) //a,b均存在斜率
	{
		if (a.getK() == b.getK())return;
		ans = findIntersectionls(a, b);
		if (b.judge(ans) == false)return;
		node_set_insert(ans);
		return;
	}
	else {
		ans = findIntersectionls(a, b);
		if (b.judge(ans) == false)return;
		node_set_insert(ans);
		return;
	}
	return;
}

void add_node_ss(lise a, lise b)
{
	node ans;
	if ((!a.getExitK()) && (!b.getExitK()))//a,b均平行于y轴
	{
		return;
	}
	else if (a.getExitK() && b.getExitK()) //a,b均存在斜率
	{
		if (a.getK() == b.getK())return;
		ans = findIntersectionss(a, b);
		if (b.judge(ans) == false || a.judge(ans) == false)return;
		node_set_insert(ans);
		return;
	}
	else {
		ans = findIntersectionss(a, b);
		if (b.judge(ans) == false || a.judge(ans) == false)return;
		node_set_insert(ans);
		return;
	}
	return;
}

void add_node_rs(rays a, lise b)
{
	node ans;
	if ((!a.getExitK()) && (!b.getExitK()))//a,b均平行于y轴
	{
		return;
	}
	else if (a.getExitK() && b.getExitK()) //a,b均存在斜率
	{
		if (a.getK() == b.getK())return;
		ans = findIntersectionrs(a, b);
		if (b.judge(ans) == false || a.judge(ans) == false)return;
		node_set_insert(ans);
		return;
	}
	else {
		ans = findIntersectionrs(a, b);
		if (b.judge(ans) == false || a.judge(ans) == false)return;
		node_set_insert(ans);
		return;
	}
	return;
}

void add_node_cc(Cycle a, Cycle b)
{
	double distance = distance_nn(a.getC(), b.getC());
	if (distance > (a.getR()+ b.getR())* (a.getR() + b.getR())) {
		return;
	}
	double x1 = a.getC().getX(), y1 = a.getC().getY(), r1 = a.getR();
	double x2 = b.getC().getX(), y2 = b.getC().getY(), r2 = b.getR();
	if (y1 == y2) {
		if (x1 == x2)exit(1);
		double x = ((r1 * r1 - r2 * r2) / (x2 - x1) + x1 + x2) / 2;
		double delta = sqrt(r1 * r1 - (x - x1) * (x - x1));
		node n1 = node(x, y1 + delta), n2 = node(x, y1 - delta);
		node_set_insert(n1);
		node_set_insert(n2);
		return;
	}
	else {
		double p = (x1 - x2) / (y2 - y1), t = -p * (x1 + x2) / 2 + (r1 * r1 - r2 * r2) / (y2 - y1) / 2 + (y2 - y1) / 2;
		double delta = sqrt((p * t - x1)* (p * t - x1)-(1+p*p)*(x1*x1+t*t-r1*r1));
		double ansx1 = (x1 - p * t + delta) / (1 + p * p), ansx2 = (x1 - p * t - delta) / (1 + p * p);
		node n1 = node(ansx1, (r1*r1-r2*r2)/(y2-y1)/2+p*(ansx1-(x1+x2)/2)+(y1+y2)/2),n2 = node(ansx2, (r1 * r1 - r2 * r2) / (y2 - y1) / 2 + p * (ansx2 - (x1 + x2) / 2) + (y1 + y2) / 2);
		node_set_insert(n1);
		node_set_insert(n2);
		return;
	}
	return;
}

void add_node_cl(Cycle a, line b)
{
	double distance = distance_nl(a.getC(), b);
	printf("start\n");
	if (distance > a.getR()* a.getR()) {
		return;
	}
	printf("part1\n");
	if (!b.getExitK()) {
		double x = b.getNode1().getX();
		double y = sqrt(a.getR() * a.getR() - (x - a.getC().getX()) * (x - a.getC().getX()));
		double y1 = a.getC().getY() - y, y2 = a.getC().getY() + y;
		node n1 = node(x, y1), n2 = node(x, y2);
		node_set_insert(n1);
		node_set_insert(n2);
		return;
	}
	else {
		double k = b.getK(), j = -b.getC() / b.getB();
		double x0 = a.getC().getX(), y0 = a.getC().getY(), r = a.getR();
		double t = j - y0, m = x0-k*t, n = 1 + k*k,o = x0*x0+t*t-r*r;
		if (o < 0)return;
		double x1 = (m + sqrt(m*m - n*o))/n, x2 = (m - sqrt(m * m - n * o)) / n;
		node n1 = node(x1, k * x1 + j), n2 = node(x2, k * x2 + j);
		node_set_insert(n1);
		node_set_insert(n2);
		return;
	}
	return;
}

void add_node_cs(Cycle a, lise b)
{
	double distance = distance_nl(a.getC(), b);
	if (distance > a.getR()* a.getR()) {
		return;
	}
	if (!b.getExitK()) {
		double x = b.getNode1().getX();
		double y = sqrt(a.getR() * a.getR() - (x - a.getC().getX()) * (x - a.getC().getX()));
		double y1 = a.getC().getY() - y, y2 = a.getC().getY() + y;
		node n1 = node(x, y1), n2 = node(x, y2);
		if (b.judge(n1)) {
			node_set_insert(n1);
		}
		if (b.judge(n2)) {
			node_set_insert(n2);
		}
		return;
	}
	else {
		double k = b.getK(), j = -b.getC() / b.getB();
		double x0 = a.getC().getX(), y0 = a.getC().getY(), r = a.getR();
		double t = j - y0, m = x0 - k * t, n = 1 + k * k, o = x0 * x0 + t * t - r * r;
		if (o < 0)return;
		double x1 = (m + sqrt(m * m - n * o)) / n, x2 = (m - sqrt(m * m - n * o)) / n;
		node n1 = node(x1, k * x1 + j), n2 = node(x2, k * x2 + j);
		if (b.judge(n1)) {
			node_set_insert(n1);
		}
		if (b.judge(n2)) {
			node_set_insert(n2);
		}
		return;
	}
	return;
}

void add_node_cr(Cycle a, rays b)
{
	double distance = distance_nl(a.getC(), b);
	if (distance > a.getR()* a.getR()) {
		return;
	}
	if (!b.getExitK()) {
		double x = b.getStart().getX();
		double y = sqrt(a.getR() * a.getR() - (x - a.getC().getX()) * (x - a.getC().getX()));
		double y1 = a.getC().getY() - y, y2 = a.getC().getY() + y;
		node n1 = node(x, y1), n2 = node(x, y2);
		if (b.judge(n1)) {
			node_set_insert(n1);
		}
		if (b.judge(n2)) {
			node_set_insert(n2);
		}
		return;
	}
	else {
		double k = b.getK(), j = -b.getC() / b.getB();
		double x0 = a.getC().getX(), y0 = a.getC().getY(), r = a.getR();
		double t = j - y0, m = x0 - k * t, n = 1 + k * k, o = x0 * x0 + t * t - r * r;
		if (o < 0)return;
		double x1 = (m + sqrt(m * m - n * o)) / n, x2 = (m - sqrt(m * m - n * o)) / n;
		node n1 = node(x1, k * x1 + j), n2 = node(x2, k * x2 + j);
		if (b.judge(n1)) {
			node_set_insert(n1);
		}
		if (b.judge(n2)) {
			node_set_insert(n2);
		}
		return;
	}
	return;
}

void delete_node_ll(line a, line b)
{
	node ans;
	if ((!a.getExitK()) && (!b.getExitK()))//a,b均平行于y轴
	{
		return;
	}
	else if (a.getExitK() && b.getExitK()) //a,b均存在斜率
	{
		if (a.getK() == b.getK())return;
		ans = findIntersectionll(a, b);
		node_set_delete(ans);
		return;
	}
	else {
		ans = findIntersectionll(a, b);
		node_set_delete(ans);
		return;
	}
	return;
}

void delete_node_lr(line a, rays b)
{
	node ans;
	if ((!a.getExitK()) && (!b.getExitK()))//a,b均平行于y轴
	{
		return;
	}
	else if (a.getExitK() && b.getExitK()) //a,b均存在斜率
	{
		if (a.getK() == b.getK())return;
		ans = findIntersectionlr(a, b);
		if (b.judge(ans) == false)return;
		node_set_delete(ans);
		return;
	}
	else {
		ans = findIntersectionlr(a, b);
		if (b.judge(ans) == false)return;
		node_set_delete(ans);
		return;
	}
	return;
}

void delete_node_rr(rays a, rays b)
{
	node ans;
	if ((!a.getExitK()) && (!b.getExitK()))//a,b均平行于y轴
	{
		return;
	}
	else if (a.getExitK() && b.getExitK()) //a,b均存在斜率
	{
		if (a.getK() == b.getK())return;
		ans = findIntersectionrr(a, b);
		if (b.judge(ans) == false || a.judge(ans) == false)return;
		node_set_delete(ans);
		return;
	}
	else {
		ans = findIntersectionrr(a, b);
		if (b.judge(ans) == false || a.judge(ans) == false)return;
		node_set_delete(ans);
		return;
	}
	return;
}

void delete_node_ls(line a, lise b)
{
	node ans;
	if ((!a.getExitK()) && (!b.getExitK()))//a,b均平行于y轴
	{
		return;
	}
	else if (a.getExitK() && b.getExitK()) //a,b均存在斜率
	{
		if (a.getK() == b.getK())return;
		ans = findIntersectionls(a, b);
		if (b.judge(ans) == false)return;
		node_set_delete(ans);
		return;
	}
	else {
		ans = findIntersectionls(a, b);
		if (b.judge(ans) == false)return;
		node_set_delete(ans);
		return;
	}
	return;
}

void delete_node_ss(lise a, lise b)
{
	node ans;
	if ((!a.getExitK()) && (!b.getExitK()))//a,b均平行于y轴
	{
		return;
	}
	else if (a.getExitK() && b.getExitK()) //a,b均存在斜率
	{
		if (a.getK() == b.getK())return;
		ans = findIntersectionss(a, b);
		if (b.judge(ans) == false || a.judge(ans) == false)return;
		node_set_delete(ans);
		return;
	}
	else {
		ans = findIntersectionss(a, b);
		if (b.judge(ans) == false || a.judge(ans) == false)return;
		node_set_delete(ans);
		return;
	}
	return;
}

void delete_node_rs(rays a, lise b)
{
	node ans;
	if ((!a.getExitK()) && (!b.getExitK()))//a,b均平行于y轴
	{
		return;
	}
	else if (a.getExitK() && b.getExitK()) //a,b均存在斜率
	{
		if (a.getK() == b.getK())return;
		ans = findIntersectionrs(a, b);
		if (b.judge(ans) == false || a.judge(ans) == false)return;
		node_set_delete(ans);
		return;
	}
	else {
		ans = findIntersectionrs(a, b);
		if (b.judge(ans) == false || a.judge(ans) == false)return;
		node_set_delete(ans);
		return;
	}
	return;
}

void delete_node_cc(Cycle a, Cycle b)
{
	double distance = distance_nn(a.getC(), b.getC());
	if (distance > (a.getR() + b.getR())* (a.getR() + b.getR())) {
		return;
	}
	double x1 = a.getC().getX(), y1 = a.getC().getY(), r1 = a.getR();
	double x2 = b.getC().getX(), y2 = b.getC().getY(), r2 = b.getR();
	if (y1 == y2) {
		if (x1 == x2)exit(1);
		double x = ((r1 * r1 - r2 * r2) / (x2 - x1) + x1 + x2) / 2;
		double delta = sqrt(r1 * r1 - (x - x1) * (x - x1));
		node n1 = node(x, y1 + delta), n2 = node(x, y1 - delta);
		node_set_delete(n1);
		node_set_delete(n2);
		return;
	}
	else {
		double p = (x1 - x2) / (y2 - y1), t = -p * (x1 + x2) / 2 + (r1 * r1 - r2 * r2) / (y2 - y1) / 2 + (y2 - y1) / 2;
		double delta = sqrt((p * t - x1) * (p * t - x1) - (1 + p * p) * (x1 * x1 + t * t - r1 * r1));
		double ansx1 = (x1 - p * t + delta) / (1 + p * p), ansx2 = (x1 - p * t - delta) / (1 + p * p);
		node n1 = node(ansx1, (r1 * r1 - r2 * r2) / (y2 - y1) / 2 + p * (ansx1 - (x1 + x2) / 2) + (y1 + y2) / 2), n2 = node(ansx2, (r1 * r1 - r2 * r2) / (y2 - y1) / 2 + p * (ansx2 - (x1 + x2) / 2) + (y1 + y2) / 2);
		node_set_delete(n1);
		node_set_delete(n2);
		return;
	}
	return;
}

void delete_node_cl(Cycle a, line b)
{
	double distance = distance_nl(a.getC(), b);
	printf("start\n");
	if (distance > a.getR()* a.getR()) {
		return;
	}
	printf("part1\n");
	if (!b.getExitK()) {
		double x = b.getNode1().getX();
		double y = sqrt(a.getR() * a.getR() - (x - a.getC().getX()) * (x - a.getC().getX()));
		double y1 = a.getC().getY() - y, y2 = a.getC().getY() + y;
		node n1 = node(x, y1), n2 = node(x, y2);
		node_set_delete(n1);
		node_set_delete(n2);
		return;
	}
	else {
		double k = b.getK(), j = -b.getC() / b.getB();
		double x0 = a.getC().getX(), y0 = a.getC().getY(), r = a.getR();
		double t = j - y0, m = x0 - k * t, n = 1 + k * k, o = x0 * x0 + t * t - r * r;
		if (o < 0)return;
		double x1 = (m + sqrt(m * m - n * o)) / n, x2 = (m - sqrt(m * m - n * o)) / n;
		node n1 = node(x1, k * x1 + j), n2 = node(x2, k * x2 + j);
		node_set_delete(n1);
		node_set_delete(n2);
		return;
	}
	return;
}

void delete_node_cs(Cycle a, lise b)
{
	double distance = distance_nl(a.getC(), b);
	if (distance > a.getR()* a.getR()) {
		return;
	}
	if (!b.getExitK()) {
		double x = b.getNode1().getX();
		double y = sqrt(a.getR() * a.getR() - (x - a.getC().getX()) * (x - a.getC().getX()));
		double y1 = a.getC().getY() - y, y2 = a.getC().getY() + y;
		node n1 = node(x, y1), n2 = node(x, y2);
		if (b.judge(n1)) {
			node_set_delete(n1);
		}
		if (b.judge(n2)) {
			node_set_delete(n2);
		}
		return;
	}
	else {
		double k = b.getK(), j = -b.getC() / b.getB();
		double x0 = a.getC().getX(), y0 = a.getC().getY(), r = a.getR();
		double t = j - y0, m = x0 - k * t, n = 1 + k * k, o = x0 * x0 + t * t - r * r;
		if (o < 0)return;
		double x1 = (m + sqrt(m * m - n * o)) / n, x2 = (m - sqrt(m * m - n * o)) / n;
		node n1 = node(x1, k * x1 + j), n2 = node(x2, k * x2 + j);
		if (b.judge(n1)) {
			node_set_delete(n1);
		}
		if (b.judge(n2)) {
			node_set_delete(n2);
		}
		return;
	}
	return;
}

void delete_node_cr(Cycle a, rays b)
{
	double distance = distance_nl(a.getC(), b);
	if (distance > a.getR()* a.getR()) {
		return;
	}
	if (!b.getExitK()) {
		double x = b.getStart().getX();
		double y = sqrt(a.getR() * a.getR() - (x - a.getC().getX()) * (x - a.getC().getX()));
		double y1 = a.getC().getY() - y, y2 = a.getC().getY() + y;
		node n1 = node(x, y1), n2 = node(x, y2);
		if (b.judge(n1)) {
			node_set_delete(n1);
		}
		if (b.judge(n2)) {
			node_set_delete(n2);
		}
		return;
	}
	else {
		double k = b.getK(), j = -b.getC() / b.getB();
		double x0 = a.getC().getX(), y0 = a.getC().getY(), r = a.getR();
		double t = j - y0, m = x0 - k * t, n = 1 + k * k, o = x0 * x0 + t * t - r * r;
		if (o < 0)return;
		double x1 = (m + sqrt(m * m - n * o)) / n, x2 = (m - sqrt(m * m - n * o)) / n;
		node n1 = node(x1, k * x1 + j), n2 = node(x2, k * x2 + j);
		if (b.judge(n1)) {
			node_set_delete(n1);
		}
		if (b.judge(n2)) {
			node_set_delete(n2);
		}
		return;
	}
	return;
}

void get_int(string l, double* a, double* b, double* c, double* d) {
	int ans[4] = { 0 },op = 1;
	for (int i = 1, k = 0, j, change = 0; i < l.size(); i++)
	{
		if (l[i] == '-')op = -1;
		for (j = i; l[j] <= '9' && l[j] >= '0'; j++)
		{
			ans[k] = ans[k] * 10 + l[j] - '0';
			change = 1;
		}
		if (change == 1) {
			ans[k] *= op;
			op = 1;
			k++;
			i = j;
			change = 0;
		}
	}
	*a = (double)ans[0];
	*b = (double)ans[1];
	*c = (double)ans[2];
	*d = (double)ans[3];
	return;
}

void get_int(string l, double* a, double* b, double* c) {
	int ans[3] = { 0 };
	for (int i = 1, k = 0, j, change = 0,op = 1; i < l.size(); i++)
	{
		if (l[i] == '-')op = -1;
		for (j = i; l[j] <= '9' && l[j] >= '0'; j++)
		{
			ans[k] = ans[k] * 10 + l[j] - '0';
			change = 1;
		}
		if (change == 1) {
			ans[k] *= op;
			op = 1;
			k++;
			i = j;
			change = 0;
		}

	}
	*a = (double)ans[0];
	*b = (double)ans[1];
	*c = (double)ans[2];
	return;
}

void Delete(string l)
{
	char p = l[0];
	double a, b, c, d;
	if (p == 'L') {
		int size = 0;
		get_int(l,&a,&b,&c,&d);
		node n1(a, b), n2(c, d);
		line tmp(n1, n2);
		for (vector<line>::iterator j = lseq.begin(); j != lseq.end();)
		{
			line tmp1 = *j;
			if (tmp1 == tmp) {
				j = lseq.erase(j);
				break;
			}
		}
		for (int j = 0, size = lseq.size(); j < size; j++) {
			delete_node_ll(tmp, lseq[j]);
		}
		for (int j = 0, size = sseq.size(); j < size; j++) {
			delete_node_ls(tmp, sseq[j]);
		}
		for (int j = 0, size = rseq.size(); j < size; j++) {
			delete_node_lr(tmp, rseq[j]);
		}
		for (int j = 0, size = cseq.size(); j < size; j++) {
			delete_node_cl(cseq[j], tmp);
		}
		
	}
	else if (p == 'C') {
		int size = 0;
		get_int(l, &a, &b, &c);
		node n1(a, b);
		Cycle tmp(n1, c);
		for (vector<Cycle>::iterator j = cseq.begin(); j != cseq.end();)
		{
			Cycle tmp1 = *j;
			if (tmp1 == tmp) {
				j = cseq.erase(j);
				break;
			}
		}
		for (int j = 0, size = cseq.size(); j < size; j++) {
			delete_node_cc(cseq[j], tmp);
		}
		for (int j = 0, size = lseq.size(); j < size; j++) {
			delete_node_cl(tmp, lseq[j]);
		}
		for (int j = 0, size = rseq.size(); j < size; j++) {
			delete_node_cr(tmp, rseq[j]);
		}
		for (int j = 0, size = sseq.size(); j < size; j++) {
			delete_node_cs(tmp, sseq[j]);
		}
	}
	else if (p == 'S') {
		int size = 0;
		get_int(l, &a, &b, &c, &d);
		node n1(a, b), n2(c, d);
		lise tmp(n1, n2);
		for (vector<lise>::iterator j = sseq.begin(); j != sseq.end();)
		{
			lise tmp1 = *j;
			if (tmp1 == tmp) {
				j = sseq.erase(j);
				break;
			}
		}
		for (int j = 0, size = lseq.size(); j < size; j++) {
			delete_node_ls(lseq[j], tmp);
		}
		for (int j = 0, size = sseq.size(); j < size; j++) {
			delete_node_ss(tmp, sseq[j]);
		}
		for (int j = 0, size = rseq.size(); j < size; j++) {
			delete_node_rs(rseq[j], tmp);
		}
		for (int j = 0, size = cseq.size(); j < size; j++) {
			delete_node_cs(cseq[j], tmp);
		}
		
	}
	else if (p == 'R') {
		int size = 0;
		get_int(l, &a, &b, &c, &d);
		node n1(a, b), n2(c, d);
		rays tmp(n1, n2);
		for (vector<rays>::iterator j = rseq.begin(); j != rseq.end();)
		{
			rays tmp1 = *j;
			if (tmp1 == tmp) {
				j = rseq.erase(j);
				break;
			}
		}
		for (int j = 0, size = lseq.size(); j < size; j++) {
			delete_node_lr(lseq[j], tmp);
		}
		for (int j = 0, size = sseq.size(); j < size; j++) {
			delete_node_rs(tmp, sseq[j]);
		}
		for (int j = 0, size = rseq.size(); j < size; j++) {
			delete_node_rr(tmp, rseq[j]);
		}
		for (int j = 0, size = cseq.size(); j < size; j++) {
			delete_node_cr(cseq[j], tmp);
		}
		
	}
	else {
	error_output("wrong input!\n");
	}
}

bool Add(string l) 
{
	char p = l[0];
	double a, b, c, d;
	if (p == 'L') {
		int size = 0;
		get_int(l, &a, &b, &c, &d);
		node n1(a, b), n2(c, d);
		line tmp(n1, n2);
		for (vector<line>::iterator j = lseq.begin(); j != lseq.end();)
		{
			line tmp1 = *j;
			if (tmp1 == tmp) {
				error_output("this object had already been added!\n");
				return false;
			}
		}
		for (int j = 0, size = lseq.size(); j < size; j++) {
			add_node_ll(tmp, lseq[j]);
		}
		for (int j = 0, size = sseq.size(); j < size; j++) {
			add_node_ls(tmp, sseq[j]);
		}
		for (int j = 0, size = rseq.size(); j < size; j++) {
			add_node_lr(tmp, rseq[j]);
		}
		for (int j = 0, size = cseq.size(); j < size; j++) {
			add_node_cl(cseq[j], tmp);
		}
		lseq.push_back(tmp);
	}
	else if (p == 'C') {
		int size = 0;
		get_int(l, &a, &b, &c);
		node n1(a, b);
		Cycle tmp(n1, c);
		for (vector<Cycle>::iterator j = cseq.begin(); j != cseq.end();)
		{
			Cycle tmp1 = *j;
			if (tmp1 == tmp) {
				error_output("this object had already been added!\n");
				return false;
			}
		}
		for (int j = 0, size = cseq.size(); j < size; j++) {
			add_node_cc(cseq[j], tmp);
		}
		for (int j = 0, size = lseq.size(); j < size; j++) {
			add_node_cl(tmp, lseq[j]);
		}
		for (int j = 0, size = rseq.size(); j < size; j++) {
			add_node_cr(tmp, rseq[j]);
		}
		for (int j = 0, size = sseq.size(); j < size; j++) {
			add_node_cs(tmp, sseq[j]);
		}
		cseq.push_back(tmp);
	}
	else if (p == 'S') {
		int size = 0;
		get_int(l, &a, &b, &c, &d);
		node n1(a, b), n2(c, d);
		lise tmp(n1, n2);
		for (vector<lise>::iterator j = sseq.begin(); j != sseq.end();)
		{
			lise tmp1 = *j;
			if (tmp1 == tmp) {
				error_output("this object had already been added!\n");
				return false;
			}
		}
		for (int j = 0, size = lseq.size(); j < size; j++) {
			add_node_ls(lseq[j], tmp);
		}
		for (int j = 0, size = sseq.size(); j < size; j++) {
			add_node_ss(tmp, sseq[j]);
		}
		for (int j = 0, size = rseq.size(); j < size; j++) {
			add_node_rs(rseq[j], tmp);
		}
		for (int j = 0, size = cseq.size(); j < size; j++) {
			add_node_cs(cseq[j], tmp);
		}
		sseq.push_back(tmp);
	}
	else if (p == 'R') {
		int size = 0;
		get_int(l, &a, &b, &c, &d);
		node n1(a, b), n2(c, d);
		rays tmp(n1, n2);
		for (vector<rays>::iterator j = rseq.begin(); j != rseq.end();)
		{
			rays tmp1 = *j;
			if (tmp1 == tmp) {
				error_output("this object had already been added!\n");
				return false;
			}
		}
		for (int j = 0, size = lseq.size(); j < size; j++) {
			add_node_lr(lseq[j], tmp);
		}
		for (int j = 0, size = sseq.size(); j < size; j++) {
			add_node_rs(tmp, sseq[j]);
		}
		for (int j = 0, size = rseq.size(); j < size; j++) {
			add_node_rr(tmp, rseq[j]);
		}
		for (int j = 0, size = cseq.size(); j < size; j++) {
			add_node_cr(cseq[j], tmp);
		}
		rseq.push_back(tmp);
	}
	else {
		try {
			throw "wrong input!\n";
		}
		catch (const char* msg) {
			cerr << msg << endl;
		}
	}
}

void read_file(string input_file)
{
	int n;
	double a, b, c, d;
	string s_tmp;
	char p;
	ifstream infile(input_file);
	if (!infile) {
		error_output("input file error!\n");
		
	}
	infile >> s_tmp;
	n = string2int(s_tmp);
	for (int i = 0; i < n; i++)
	{
		infile >> s_tmp;
		p = string2char(s_tmp);
		if (p == 'L') {
			int size = 0;
			infile >> s_tmp;
			a = string2double(s_tmp);
			infile >> s_tmp;
			b = string2double(s_tmp);
			infile >> s_tmp;
			c = string2double(s_tmp);
			infile >> s_tmp;
			d = string2double(s_tmp);
			node n1(a, b), n2(c, d);
			line tmp(n1, n2);
			for (int j = 0, size = lseq.size(); j < size; j++) {
				add_node_ll(tmp, lseq[j]);
			}
			for (int j = 0, size = sseq.size(); j < size; j++) {
				add_node_ls(tmp, sseq[j]);
			}
			for (int j = 0, size = rseq.size(); j < size; j++) {
				add_node_lr(tmp, rseq[j]);
			}
			for (int j = 0, size = cseq.size(); j < size; j++) {
				add_node_cl(cseq[j], tmp);
			}
			lseq.push_back(tmp);
		}
		else if (p == 'C') {
			int size = 0;
			infile >> s_tmp;
			a = string2double(s_tmp);
			infile >> s_tmp;
			b = string2double(s_tmp);
			infile >> s_tmp;
			c = string2double(s_tmp);
			node n1(a, b);
			Cycle tmp(n1, c);
			for (int j = 0, size = cseq.size(); j < size; j++) {
				add_node_cc(cseq[j], tmp);
			}
			for (int j = 0, size = lseq.size(); j < size; j++) {
				add_node_cl(tmp, lseq[j]);
			}
			for (int j = 0, size = rseq.size(); j < size; j++) {
				add_node_cr(tmp, rseq[j]);
			}
			for (int j = 0, size = sseq.size(); j < size; j++) {
				add_node_cs(tmp, sseq[j]);
			}
			cseq.push_back(tmp);
		}
		else if (p == 'S') {
			int size = 0;
			infile >> s_tmp;
			a = string2double(s_tmp);
			infile >> s_tmp;
			b = string2double(s_tmp);
			infile >> s_tmp;
			c = string2double(s_tmp);
			infile >> s_tmp;
			d = string2double(s_tmp);
			node n1(a, b), n2(c, d);
			lise tmp(n1, n2);
			for (int j = 0, size = lseq.size(); j < size; j++) {
				add_node_ls(lseq[j], tmp);
			}
			for (int j = 0, size = sseq.size(); j < size; j++) {
				add_node_ss(tmp, sseq[j]);
			}
			for (int j = 0, size = rseq.size(); j < size; j++) {
				add_node_rs(rseq[j], tmp);
			}
			for (int j = 0, size = cseq.size(); j < size; j++) {
				add_node_cs(cseq[j], tmp);
			}
			sseq.push_back(tmp);
		}
		else if (p == 'R') {
			int size = 0;
			infile >> s_tmp;
			a = string2double(s_tmp);
			infile >> s_tmp;
			b = string2double(s_tmp);
			infile >> s_tmp;
			c = string2double(s_tmp);
			infile >> s_tmp;
			d = string2double(s_tmp);
			node n1(a, b), n2(c, d);
			rays tmp(n1, n2);
			for (int j = 0, size = lseq.size(); j < size; j++) {
				add_node_lr(lseq[j], tmp);
			}
			for (int j = 0, size = sseq.size(); j < size; j++) {
				add_node_rs(tmp, sseq[j]);
			}
			for (int j = 0, size = rseq.size(); j < size; j++) {
				add_node_rr(tmp, rseq[j]);
			}
			for (int j = 0, size = cseq.size(); j < size; j++) {
				add_node_cr(cseq[j], tmp);
			}
			rseq.push_back(tmp);
		}
		else {

			try {
				throw "wrong input!\n";
			}
			catch (const char* msg) {
				cerr << msg << endl;
			}
		}
	}
	return;
}

int main(int argc, char** argv)
{
	read_file("input.txt");//argv[2]
	ofstream outfile("output.txt");   // argv[4]
	cout << countn;
	return 0;
}