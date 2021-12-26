// StatisticCounter.cpp: implementation of the StatisticCounter class.
//
//////////////////////////////////////////////////////////////////////

#include "StatisticCounter.h"
#include "common.h"

int StatisticCounter::total_stat_number = 0;

// init
StatisticCounter::StatisticCounter(const string& desc)
{
	description = desc;
	position = 0;
	last = &first;
	storeCount(counter = 0);
	stat_number = total_stat_number++;
	marks = new DataVector<MarkData>(100);
}

StatisticCounter::~StatisticCounter()
{
	// store to file
	ostringstream file_name;
	file_name << "stat" << stat_number << ".txt";
	ofstream file(file_name.str().c_str());

	if(file){
		file << "Statistics - " << description << endl;
		size_t i, ct = marks->countInt();
		file << "--- marks -----\n";
		for(i = 0; i < ct; i++)
			file << marks->get(i).position << '\t' << marks->get(i).label << endl;
		file << "---------------\n";
		CountData* node = &first;
		int j = 0;
		int last_value = -1;
		int prev_last_value = -2;
		for(i = 0; i < (size_t)position; i++){
			if(j >= node->position){
				node = node->next;
				j = 0;
			}
			if(node->data[j] != last_value){
				if(prev_last_value == last_value)
					file << (i-1) << '\t' << last_value << endl;
				else prev_last_value = last_value;
				file << i << '\t' << (last_value = node->data[j]) << endl;
			}else prev_last_value = last_value;
			j++;
		}
		if(last_value == prev_last_value)
			file << (position-1) << '\t' << last_value << endl;
	}
	// clean
	while(first.next){
		CountData* temp = first.next->next;
		delete first.next;
		first.next = temp;
	}
	delete marks;
}

void StatisticCounter::storeCount(int ct)
{
	++position;
	if(last->position >= DATA_SIZE){
		last->next = new CountData;
		last = last->next;
	}
	last->data[last->position++] = ct;
}

void StatisticCounter::addMark(const string& label)
{
	marks->add(MarkData(label, position));
}
