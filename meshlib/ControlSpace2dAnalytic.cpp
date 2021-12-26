// ControlSpace2dAnalytic.cpp: implementation of the ControlSpace2dAnalytic class.
//
//////////////////////////////////////////////////////////////////////

#include "ControlSpace2dAnalytic.h"
#include "SurfaceParametric.h"
#include "DEquation.h"
#include "MeshLog.h"
#include "MeshData.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ControlSpace2dAnalytic::ControlSpace2dAnalytic(SurfaceConstPtr surface, 
										   const char * default_lx_str, 
										   const char * default_ly_str, 
										   const char * default_angle_str) : ControlSpace2d(surface)
{
	assert(default_lx_str != nullptr);
	assert(default_ly_str != nullptr);
	assert(default_angle_str != nullptr);

	m_default_lx = new DEquation(default_lx_str);
	m_default_ly = new DEquation(default_ly_str);
	m_default_angle = new DEquation(default_angle_str);

	m_list = nullptr;
}

ControlSpace2dAnalytic::ControlSpace2dAnalytic(SurfaceConstPtr surface, 
										   DEquation * lx_eq, 
										   DEquation * ly_eq, 
										   DEquation * angle_eq) : ControlSpace2d(surface)
{
	assert(lx_eq != nullptr);
	assert(ly_eq != nullptr);
	assert(angle_eq != nullptr);

	m_default_lx = lx_eq;
	m_default_ly = ly_eq;
	m_default_angle = angle_eq;

	m_list = nullptr;
}

ControlSpace2dAnalytic::~ControlSpace2dAnalytic()
{
	delete m_default_lx;
	delete m_default_ly;
	delete m_default_angle;

	FunctionNode* node = m_list;
	while(node){
		delete node->condition;
		if(node->condition2) delete node->condition2;
		delete node->lx;
		delete node->ly;
		delete node->angle;
		m_list = node->next;
		delete node;
		node = m_list;
	}
}


void ControlSpace2dAnalytic::addCaseAsFirst(const char * test_str, 
		const char * lx_str, const char * ly_str, const char * angle_str)
{
	FunctionNode* item = new FunctionNode;
	item->condition = new DEquation(test_str);
	item->condition2 = nullptr;
	item->lx = new DEquation(lx_str);
	item->ly = new DEquation(ly_str);
	item->angle = new DEquation(angle_str);
	item->next = m_list;
	m_list = item;
}

void ControlSpace2dAnalytic::addCaseAsFirst(const char * test1_str, const char * test2_str, 
		const char * lx_str, const char * ly_str, const char * angle_str)
{
	FunctionNode* item = new FunctionNode;
	item->condition = new DEquation(test1_str);
	item->condition2 = new DEquation(test2_str);
	item->lx = new DEquation(lx_str);
	item->ly = new DEquation(ly_str);
	item->angle = new DEquation(angle_str);
	item->next = m_list;
	m_list = item;
}

/// Get sizing info (matrix mode) at the given point
ControlDataMatrix2d ControlSpace2dAnalytic::getMetricAtPoint(const DPoint2d& pt) const
{
	ControlDataStretch2d data;
	bool invalid = true;
	for(FunctionNode* node = m_list; node && invalid; node = node->next){
		if(node->condition->getValue(pt.x, pt.y) < 0.0) continue;
		if(node->condition2 && node->condition2->getValue(pt.x, pt.y) < 0.0) continue;

		data.lx = node->lx->getValue(pt.x, pt.y);
		data.ly = node->ly->getValue(pt.x, pt.y);
		data.angle = node->angle->getValue(pt.x, pt.y);
		invalid = false;
	}

	if(invalid){
		data.lx = m_default_lx->getValue(pt.x, pt.y);
		data.ly = m_default_ly->getValue(pt.x, pt.y);
		data.angle = m_default_angle->getValue(pt.x, pt.y);
	}

	return DMetric2d::stretchToMatrix(data);
}
