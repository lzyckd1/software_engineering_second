#include "node.h"
#include "Cycle.h"
#include "line.h"
#include "lise.h"
#include "rays.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <set>
#include <vector>

using namespace std;
const int maxn = 1000;
set <node> node_set;
vector <line> lseq;
vector <lise> sseq;
vector <Cycle> cseq;
vector <rays> rseq;
int countn = 0;
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
		if (node_set.find(ans) == node_set.end()) {
			countn = countn + 1;
			node_set.insert(ans);
		}
		return;
	}
	else {
		ans = findIntersectionll(a, b);
		if (node_set.find(ans) == node_set.end()) {
			countn = countn + 1;
			node_set.insert(ans);
		}
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
		if (node_set.find(ans) == node_set.end()) {
			countn = countn + 1;
			node_set.insert(ans);
		}
		return;
	}
	else {
		ans = findIntersectionlr(a, b);
		if (b.judge(ans) == false)return;
		if (node_set.find(ans) == node_set.end()) {
			countn = countn + 1;
			node_set.insert(ans);
		}
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
		if (node_set.find(ans) == node_set.end()) {
			countn = countn + 1;
			node_set.insert(ans);
		}
		return;
	}
	else {
		ans = findIntersectionrr(a, b);
		if (b.judge(ans) == false || a.judge(ans) == false)return;
		if (node_set.find(ans) == node_set.end()) {
			countn = countn + 1;
			node_set.insert(ans);
		}
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
		if (node_set.find(ans) == node_set.end()) {
			countn = countn + 1;
			node_set.insert(ans);
		}
		return;
	}
	else {
		ans = findIntersectionls(a, b);
		if (b.judge(ans) == false)return;
		if (node_set.find(ans) == node_set.end()) {
			countn = countn + 1;
			node_set.insert(ans);
		}
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
		if (node_set.find(ans) == node_set.end()) {
			countn = countn + 1;
			node_set.insert(ans);
		}
		return;
	}
	else {
		ans = findIntersectionss(a, b);
		if (b.judge(ans) == false || a.judge(ans) == false)return;
		if (node_set.find(ans) == node_set.end()) {
			countn = countn + 1;
			node_set.insert(ans);
		}
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
		if (node_set.find(ans) == node_set.end()) {
			countn = countn + 1;
			node_set.insert(ans);
		}
		return;
	}
	else {
		ans = findIntersectionrs(a, b);
		if (b.judge(ans) == false || a.judge(ans) == false)return;
		if (node_set.find(ans) == node_set.end()) {
			countn = countn + 1;
			node_set.insert(ans);
		}
		return;
	}
	return;
}

int main(int argc, char** argv)
{
	int n;
	double a, b, c, d;
	char p;
	ifstream infile("input.txt");  // argv[2]
	ofstream outfile("output.txt");   // argv[4]
	infile >> n;
	for (int i = 0; i < n; i++)
	{
		infile >> p;
		if (p == 'L') {
			int size = 0;
			infile >> a >> b >> c >> d;
			node n1(a, b), n2(c, d);
			line tmp(n1, n2);
			for (int j = 0, size = lseq.size(); j < size; j++) {
				add_node_ll(tmp, lseq[j]);
			}
			for (int j = 0,size = sseq.size(); j < size; j++) {
				add_node_ls(tmp, sseq[j]);
			}
			for (int j = 0, size = rseq.size(); j < size; j++) {
				add_node_lr(tmp, rseq[j]);
			}
			lseq.push_back(tmp);
		}
		else if (p == 'C') {
			int size = 0;
			infile >> a >> b >> c;
			node n1(a, b);
			Cycle tmp(n1, c);
			cseq.push_back(tmp);
		}
		else if (p == 'S') {
			int size = 0;
			infile >> a >> b >> c >> d;
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
			sseq.push_back(tmp);
		}
		else if (p == 'R') {
			int size = 0;
			infile >> a >> b >> c >> d;
			node n1(a, b), n2(c, d);
			rays tmp(n1, n2);
			rseq.push_back(tmp);
			for (int j = 0, size = lseq.size(); j < size; j++) {
				add_node_lr(lseq[j], tmp);
			}
			for (int j = 0, size = sseq.size(); j < size; j++) {
				add_node_rs(tmp, sseq[j]);
			}
			for (int j = 0, size = rseq.size(); j < size; j++) {
				add_node_rr(tmp, rseq[j]);
			}
		}
		else {
			exit(1);
		}
	}
	outfile << countn;
	return 0;
}