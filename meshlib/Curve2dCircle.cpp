/////////////////////////////////////////////////////////////////////////////
// Curve2dCircle.cpp
// Circle 2D
//	[Curve2dParametric->Curve2dCircle]
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2004-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#include "Curve2dCircle.h"
#include "DPoint.h"
#include "DEquation.h"

/// t - 0-1
DPoint2d Curve2dCircle::getPoint(double t) const {
	t *= 2*PI;
	return DPoint2d(m_middle.x + m_radius*cos(t), m_middle.y + m_radius*sin(t));
}

double Curve2dCircle::getParameterInRange(const DPoint2d& pt, double ts, double t_min, double t_max) const 
{
	double t = getParameter( pt, ts );
	return std::max(t_min, std::min(t_max, t));
}

double Curve2dCircle::getCircleParameter( const DPoint2d& cmid, double cr, const DPoint2d& pt, double ts )
{
	if(cr <= 0.0) return 0.0;
	double r = cmid.distance(pt);
	double cos_dx = (pt.x - cmid.x) / r;
	if(cos_dx > 1.0) {
		cos_dx = 1.0;
	}else if(cos_dx < -1.0){
		cos_dx = -1.0;
	}
	double sin_dy = (pt.y - cmid.y) / r;
	if(sin_dy > 1.0) {
		sin_dy = 1.0;
	}else if(sin_dy < -1.0){
		sin_dy = -1.0;
	}
	double t = acos(cos_dx) / (2*PI);	// angle -> [0,PI], t -> [0-0.5]
	if(asin(sin_dy) < 0) t = 1 - t;		// angle -> [0,2PI], t -> [0-1]

	// ts - choose t nearest
	if( t > ts )
		while( (t-ts) > 0.5 ) t -= 1.0;
	else
		while( (ts-t) > 0.5 ) t += 1.0;

	return t;
}

double Curve2dCircle::getParameter( const DPoint2d& pt, double ts ) const {
	return getCircleParameter( m_middle, m_radius, pt, ts );
}

/// t - 0-1
DVector2d Curve2dCircle::getDerivative(double t) const {
	t *= 2*PI;
	return DVector2d(-m_radius*sin(t), m_radius*cos(t));
}

/// t - 0-1
DVector2d Curve2dCircle::getSecondDerivative(double t) const {
	t *= 2*PI;
	return DVector2d(-m_radius*cos(t), -m_radius*sin(t));
}

/// t - 0-1
DVector2d Curve2dCircle::getThirdDerivative(double t) const {
	t *= 2*PI;
	return DVector2d(m_radius*sin(t), -m_radius*cos(t));
}

ostream& operator<<(ostream& os, const Curve2dCircle *circle){
	os << "M=" << circle->m_middle << " R=" << circle->m_radius;
	return os;
}

istream& operator>>(istream& is, Curve2dCircle *circle){
	char z, text[200];
	double v;
	is >> ws >> z >> z;	// "M="
	is >> circle->m_middle;
	is >> ws >> z >> z >> text;	// "R="
	if(DEquation::stringToDouble(text, DEquation::v_length, v))
		circle->m_radius = v;
	else
		circle->m_radius = -1.0;

	return is;
}

/// Store XML description to stream
ostream& Curve2dCircle::storeXML(ostream& os, const string& prefix) const
{
	os << prefix << "<circle\n";
	Curve2dParametric::storeXML(os, prefix + "\t");
	os << prefix << "\t<middle> " << m_middle << " </middle>\n";
	os << prefix << "\t<radius>" << DEquation::doubleToString(m_radius, DEquation::v_length) << "</radius>\n";
	return os << prefix << "</circle>\n";
}
