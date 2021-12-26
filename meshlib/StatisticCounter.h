// StatisticCounter.h: interface for the StatisticCounter class.
//
//////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(STATISTICCOUNTER_H__INCLUDED)
#define STATISTICCOUNTER_H__INCLUDED

#include "DataVector.h"
#include "common.h"

class StatisticCounter  
{
public:
	enum { DATA_SIZE = 1000 };
	StatisticCounter(const string& desc = "generic");
	virtual ~StatisticCounter();
public:
	void addMark(const string& label);
	//----------------
	int getPosition() const { return position; }
	void storeCount(int ct);
	void increase() { storeCount(++counter); }
	void decrease() { storeCount(--counter); }
	void addIfDifferent(int value) { if(value != counter) storeCount(counter=value); }
	void add(int value) { storeCount(counter=value); }
protected:
	struct CountData{
		CountData() : position(0), next(0) {}
		int data[DATA_SIZE];
		int position;
		CountData* next;
	};
	struct MarkData{
		MarkData(string lab = "", int pos = 0) : position(pos), label(lab) {  }
		int position;
		string label;
	};
protected:
	int position;
	int counter;
	int stat_number;
	string description;
	CountData first;
	CountData *last;
	static int total_stat_number;
	DataVector<MarkData> *marks;
};

#endif // !defined(STATISTICCOUNTER_H__INCLUDED)
