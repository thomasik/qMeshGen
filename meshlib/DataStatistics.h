/////////////////////////////////////////////////////////////////////////////
// DataStatistics.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2005-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#if !defined(DATASTATISTICS_H__INCLUDED)
#define DATASTATISTICS_H__INCLUDED

#pragma once
#include "common.h"
#include "DataList.h"
#include "DataVector.h"

#include <vector>

/**
 * Class stores data for statistical analysis
 */
class DataStatistics
{
public:
	/// Standard constructor 
	DataStatistics() : m_data(STEP), m_valid(false) {}
public:
	/// Get data count
	inline int countInt() const { return m_data.countInt(); }
	/// Performs calculations
	bool calculate();
	/// Clear data
	void clear() { m_data.clear(); m_valid = false; }
	/// Adds another data
	void add(double d) { m_data.append(d); m_valid = false; }
	/// Validity test
	inline bool valid() const { return m_valid; }
	/// Returns minimum value
	double minimum() const { assert(m_valid); return m_minimum; }
	/// Returns maximum value
	double maximum() const { assert(m_valid); return m_maximum; }
	/// Returns average
	double average() const { assert(m_valid); return m_average; }
	/// Returns geometric average
	double averageGeometric() const { assert(m_valid); return m_average_geometric; }
	/// Returns standard deviation
	double stdDev() const { assert(m_valid); return m_std_dev; } 
	/// Returns quantity of data in ranges [d_min, d_max)
	std::vector<double> getDataCountInRanges(const std::vector<double> & vranges) const;
	/// Returns quantity of data in range [d_min, d_max)
	int getDataCountInRange(double d_min, double d_max) const;
	/// Returns quantity of data in range (-inf, d_max)
	int getDataCountBottom(double d_max) const;
	/// Returns quantity of data in range [d_min,+inf)
	int getDataCountTop(double d_min) const;
	/// Stores results to log file
	void logStats(const string& caption, const string& shortcut) const;
	/// Stores results to log file
	void logStats(const string& caption, const string& shortcut, const DataVector<double>& ranges) const;
protected:
	static const unsigned int STEP = 1000;
	/// Data container
	DataCompoundList<double> m_data;
	/// Valid mark
	bool m_valid;
	/// Minimum
	double m_minimum;
	/// Maximum
	double m_maximum;
	/// Average
	double m_average;
	/// Average geometric
	double m_average_geometric;
	/// Standard deviation
	double m_std_dev;
};

#endif // !defined(DATASTATISTICS_H__INCLUDED)
