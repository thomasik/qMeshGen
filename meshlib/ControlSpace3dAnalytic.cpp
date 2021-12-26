// ControlSpace3dAnalytic.cpp: implementation of the ControlSpace3dAnalytic class.
//
//////////////////////////////////////////////////////////////////////

#include "ControlSpace3dAnalytic.h"
#include "DEquation.h"
#include "MeshLog.h"
#include "MeshData.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ControlSpace3dAnalytic::ControlSpace3dAnalytic(
	const char * default_lx_str, const char * default_ly_str, const char * default_lz_str, 
	const char * default_angle1_str, const char * default_angle2_str, const char * default_angle3_str,
		DEquationConstTable* ctable)
{
	assert(default_lx_str != nullptr);
	assert(default_ly_str != nullptr);
	assert(default_lz_str != nullptr);
	assert(default_angle1_str != nullptr);
	assert(default_angle2_str != nullptr);
	assert(default_angle3_str != nullptr);

	m_default_len[0] = new DEquation(default_lx_str, ctable);
	m_default_len[1] = new DEquation(default_ly_str, ctable);
	m_default_len[2] = new DEquation(default_lz_str, ctable);
	m_default_angle[0] = new DEquation(default_angle1_str, ctable);
	m_default_angle[1] = new DEquation(default_angle2_str, ctable);
	m_default_angle[2] = new DEquation(default_angle3_str, ctable);

	m_list = nullptr;
}

ControlSpace3dAnalytic::ControlSpace3dAnalytic(
	DEquation * lx_eq, DEquation * ly_eq, DEquation * lz_eq, 
	DEquation * angle1_eq, DEquation * angle2_eq, DEquation * angle3_eq)
{
	assert(lx_eq != nullptr);
	assert(ly_eq != nullptr);
	assert(lz_eq != nullptr);
	assert(angle1_eq != nullptr);
	assert(angle2_eq != nullptr);
	assert(angle3_eq != nullptr);

	m_default_len[0] = lx_eq;
	m_default_len[1] = ly_eq;
	m_default_len[2] = lz_eq;
	m_default_angle[0] = angle1_eq;
	m_default_angle[1] = angle2_eq;
	m_default_angle[2] = angle3_eq;

	m_list = nullptr;
}

ControlSpace3dAnalytic::~ControlSpace3dAnalytic()
{
	for(int i = 0; i < 3; i++){
		delete m_default_len[i];
		delete m_default_angle[i];
	}

	FunctionNode* node = m_list;
	while(node){
		delete node->condition[0];
		delete node->condition[1];
		for(int i = 0; i < 3; i++){
			delete node->len[i];
			delete node->angle[i];
		}
		m_list = node->next;
		delete node;
		node = m_list;
	}
}


void ControlSpace3dAnalytic::addCaseAsFirst(
		const char * test_str1, const char * test_str2, 
		const char * lx_str, const char * ly_str, const char * lz_str, 
		const char * angle1_str, const char * angle2_str, const char * angle3_str,
		DEquationConstTable* ctable)
{
	FunctionNode* item = new FunctionNode;
	item->condition[0] = new DEquation(test_str1, ctable);
	item->condition[1] = new DEquation(test_str2, ctable);
	item->len[0] = new DEquation(lx_str, ctable);
	item->len[1] = new DEquation(ly_str, ctable);
	item->len[2] = new DEquation(lz_str, ctable);
	item->angle[0] = new DEquation(angle1_str, ctable);
	item->angle[1] = new DEquation(angle2_str, ctable);
	item->angle[2] = new DEquation(angle3_str, ctable);
	item->next = m_list;
	m_list = item;
}

/// Get sizing info (matrix mode) at the given point
ControlDataMatrix3d ControlSpace3dAnalytic::getMetricAtPoint(const DPoint3d& pt) const
{
	ControlDataStretch3d data;
	bool invalid = true;
	for(FunctionNode* node = m_list; node && invalid; node = node->next){
		if(node->condition[0]->getValue(pt.x, pt.y, pt.z) >= 0 &&
			node->condition[1]->getValue(pt.x, pt.y, pt.z) >= 0)
		{
			data.lx = node->len[0]->getValue(pt.x, pt.y, pt.z);
			data.ly = node->len[1]->getValue(pt.x, pt.y, pt.z);
			data.lz = node->len[2]->getValue(pt.x, pt.y, pt.z);
			data.ax = node->angle[0]->getValue(pt.x, pt.y, pt.z);
			data.ay = node->angle[1]->getValue(pt.x, pt.y, pt.z);
			data.az = node->angle[2]->getValue(pt.x, pt.y, pt.z);
			invalid = false;
		}
	}

	if(invalid){
		data.lx = m_default_len[0]->getValue(pt.x, pt.y, pt.z);
		data.ly = m_default_len[1]->getValue(pt.x, pt.y, pt.z);
		data.lz = m_default_len[2]->getValue(pt.x, pt.y, pt.z);
		data.ax = m_default_angle[0]->getValue(pt.x, pt.y, pt.z);
		data.ay = m_default_angle[1]->getValue(pt.x, pt.y, pt.z);
		data.az = m_default_angle[2]->getValue(pt.x, pt.y, pt.z);
	}

	return DMetric3d::stretchToMatrix(data);
}
