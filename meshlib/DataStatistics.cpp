/////////////////////////////////////////////////////////////////////////////
// DataStatistics.cpp
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2005-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "DataStatistics.h"
#include "common.h"

bool DataStatistics::calculate(){
	int ct = m_data.countInt();
	if(ct == 0) return (m_valid = false);
	auto it = m_data.iterator();
	assert(it.valid());
	double x = m_average = m_minimum = m_maximum = it.item();
	m_average_geometric = 1.0;
	double n_power = 1.0 / (double)ct;
	for (it.moveNext(); it.valid(); it.moveNext()) {
		double d = it.item();
		if(d < m_minimum) m_minimum = d;
		if(d > m_maximum) m_maximum = d;
		m_average += d;
		if(m_minimum > 0.0){
			x *= d;
			if(x < 0.1){
				m_average_geometric *= pow(x, n_power);
				x = 1.0;
			}
		}
	}
	if(m_minimum > 0.0) m_average_geometric *= pow(x, n_power);
	else m_average_geometric = 0.0;
	m_average /= ct;
	m_std_dev = 0.0;
	for(it = m_data.iterator(); it.valid(); it.moveNext())
		m_std_dev += sqr(m_average - it.item());
	m_std_dev = sqrt(m_std_dev/ct);
	return (m_valid = true);
}

std::vector<double> DataStatistics::getDataCountInRanges(const std::vector<double> & vranges) const
{
	int i, range_ct = (int)vranges.size();
	std::vector<double> vrec(range_ct + 1, 0.0);
	for(auto it = m_data.iterator(); it.valid(); it.moveNext()) {
		double d = it.item();
		for (i = 0; (i < range_ct) && (d >= vranges[i]); i++);
		vrec[i]++;
	}
	return vrec;
}

int DataStatistics::getDataCountInRange(double d_min, double d_max) const {
	int counter = 0;
	for (auto it = m_data.iterator(); it.valid(); it.moveNext())
		if( (it.item() >= d_min) && (it.item() < d_max)) ++counter;
	return counter;
}

/// Returns quantity of data in range (-inf, d_max)
int DataStatistics::getDataCountBottom(double d_max) const {
	int counter = 0;
	for (auto it = m_data.iterator(); it.valid(); it.moveNext())
		if( it.item() < d_max) ++counter;
	return counter;
}

/// Returns quantity of data in range [d_min,+inf)
int DataStatistics::getDataCountTop(double d_min) const {
	int counter = 0;
	for (auto it = m_data.iterator(); it.valid(); it.moveNext())
		if( it.item() >= d_min) ++counter;
	return counter;
}

void DataStatistics::logStats(const string& caption, const string& shortcut) const
{
	LOG4CPLUS_INFO(MeshLog::logger_mesh, shortcut << "-caption:\t" << caption);
	LOG4CPLUS_INFO(MeshLog::logger_mesh, shortcut << "-ave:\t" << (m_valid?m_average:0.0));
	LOG4CPLUS_INFO(MeshLog::logger_mesh, shortcut << "-stddev:\t" << (m_valid?m_std_dev:0.0));
	LOG4CPLUS_INFO(MeshLog::logger_mesh, shortcut << "-min:\t" << (m_valid?m_minimum:0.0));
	LOG4CPLUS_INFO(MeshLog::logger_mesh, shortcut << "-max:\t" << (m_valid?m_maximum:0.0));
}

void DataStatistics::logStats(const string& caption, const string& shortcut, const DataVector<double>& ranges) const
{
	logStats(caption, shortcut);
	int ct = ranges.countInt();
	if(m_valid && (ct > 0)){
		LOG4CPLUS_INFO(MeshLog::logger_mesh, 
			shortcut << "-range:\t" << "(-inf," <<
			ranges.get(0) << ")\t" << getDataCountBottom(ranges.get(0)));
		for(int i = 1; i < ct; i++){
			LOG4CPLUS_INFO(MeshLog::logger_mesh, 
				shortcut << "-range:\t" << "[" <<
				ranges.get(i-1) << "," << ranges.get(i) << ")\t" << 
				getDataCountInRange(ranges.get(i-1),ranges.get(i)));
		}
		LOG4CPLUS_INFO(MeshLog::logger_mesh, 
			shortcut << "-range:\t" << "[" <<
			ranges.get(ct-1) << ",+inf)\t" << getDataCountTop(ranges.get(ct-1)));
	}
}

