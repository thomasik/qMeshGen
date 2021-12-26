/////////////////////////////////////////////////////////////////////////////
// EPSFile.cpp
// Klasa odpowiedzialna za utowrzenie pliku w formacji EPS (postscript)
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	1999/2000
//	Generacja siatek niestrukturalnych
/////////////////////////////////////////////////////////////////////////////

#include "EPSFile.h"

//////////////////////////////////////////////////////////////////////
// Konstruktor standardowy
EPSFile::EPSFile(const char* fname, double min_x, double max_x, double min_y, double max_y, int length, int width)
{
	m_file.open(fname);
	// Wspó³czynniki transformacji wspó³rzêdnych:
	double dx = max_x - min_x;
	double dy = max_y - min_y;
	if(dx > dy){
		m_ratio = width / dx;
		length = (int) (dy * m_ratio);
	}else{
		m_ratio = length / dy;
		width = (int) (dx * m_ratio);
	}
	m_dx = 5.0 - (m_ratio * min_x);
	m_dy = 5.0 - (m_ratio * min_y);
	m_line_width = std::max(width, length) / 5000.0;

	// Nag³ówek pliku
	header(length + 10, width + 10);

	m_viewport.valid = false;
}

// Konstruktor standardowy
EPSFile::EPSFile(const string& fname, double min_x, double max_x, double min_y, double max_y, int length, int width)
{
	m_file.open(fname.c_str());
	// Wspó³czynniki transformacji wspó³rzêdnych:
	double dx = max_x - min_x;
	double dy = max_y - min_y;
	if(dx > dy){
		m_ratio = width / dx;
		length = (int) (dy * m_ratio);
	}else{
		m_ratio = length / dy;
		width = (int) (dx * m_ratio);
	}
	m_dx = 5.0 - (m_ratio * min_x);
	m_dy = 5.0 - (m_ratio * min_y);
	m_line_width = std::max(width, length) / 5000.0;

	// Nag³ówek pliku
	header(length + 10, width + 10);

	m_viewport.valid = false;
}

//////////////////////////////////////////////////////////////////////
// Destruktor
EPSFile::~EPSFile()
{
	tail();
}

//////////////////////////////////////////////////////////////////////
// Nag³ówek pliku z ustaleniem obejmuj¹cego prostok¹ta
void EPSFile::header(int length, int width)
{
	m_file << "%!PS-Adobe-3.0 EPSF-3.0" << endl;
	m_file << "%%BoundingBox: 0 0 " << width << " " << length << endl;
	m_file << "%%Pages: 0" << endl;
	m_file << "%%Creator: MeshGenerator ver 3.1" << endl;
	m_file << "%%EndComments" << endl;
	m_file << "%*Font: cmr10 9.96265 9.96265 31:e000f" << endl;
	m_file << "%%EndProlog" << endl;
	m_file << "%%Page: ""mesh"" 1" << endl;
}

//////////////////////////////////////////////////////////////////////
// Koniec pliku
void EPSFile::tail()
{
	m_file << "\n%%EOF" << endl;
}

//////////////////////////////////////////////////////////////////////
// Generuje polecenia jêzyka postscript powoduj¹ce narysowanie 
//	zadanych trójk¹tów w kolorze okreœlonym przez identyfikator
//	obszaru i wspó³czynnik jakoœci,
void EPSFile::drawTriangle(const DPoint2d &a, const DPoint2d &b, const DPoint2d &c, int area_id, double quality)
{
	if(area_id < 0) return;
	if(m_viewport.valid){
		if(!m_viewport.contains(a)) return;
		if(!m_viewport.contains(b)) return;
		if(!m_viewport.contains(c)) return;
	}
	m_file << "newpath";
//	m_file.precision(2);
	m_file << (m_dx + (m_ratio * a.x)) << " ";
//	m_file.precision(2);
	m_file << (m_dy + (m_ratio * a.y)) << " moveto ", 
//	m_file.precision(2);
	m_file << (m_dx + (m_ratio * b.x)) << " ";
//	m_file.precision(2);
	m_file << (m_dy + (m_ratio * b.y)) << " lineto ";
//	m_file.precision(2);
	m_file << (m_dx + (m_ratio * c.x)) << " ";
//	m_file.precision(2);
	m_file << (m_dy + (m_ratio * c.y)) << " lineto closepath" << endl;
	double d;
	double red		= quality * modf(1.234 * area_id, &d);
	double green	= quality * modf(2.345 * area_id, &d);
	double blue		= quality * modf(3.456 * area_id, &d);
//	m_file.precision(3);
	m_file << red << " ";
//	m_file.precision(3);
	m_file << green << " ";
//	m_file.precision(3);
	m_file << blue << " setrgbcolor fill" << endl;
}

//////////////////////////////////////////////////////////////////////
// Generuje polecenia jêzyka postscript powoduj¹ce wydrukowanie
//	zadanego numeru na okreœlonej pozycji (w czarnym kolorze)
void EPSFile::drawNumber(const DPoint2d &pt, int nr)
{
	if(m_viewport.valid){
		if(!m_viewport.contains(pt)) return;
	}
	m_file << "(" << nr << ") ";
//	m_file.precision(2);
	m_file << (2 + m_dx + m_ratio * pt.x) << " ";
//	m_file.precision(2);
	m_file << (2 + m_dy + m_ratio * pt.y);
	m_file << " moveto 0.0 setgray cmr10 9.96265 fshow" << endl;
}

//////////////////////////////////////////////////////////////////////
// Generuje polecenia jêzyka postscript powoduj¹ce wydrukowanie
//	zadanej linii ³amanej
void EPSFile::drawPolyLine(DPoint2d *polyline, int ct, int type)
{
	if(m_viewport.valid){
		for(int i = 0; i < ct; i++){
			if(!m_viewport.contains(polyline[i])) return;
		}
	}
	m_file << "newpath ";
//	m_file.precision(2);
	m_file << (m_dx + (m_ratio * polyline[0].x)) << " ";
//	m_file.precision(2);
	m_file << (m_dy + (m_ratio * polyline[0].y)) << " moveto ";
	for(int i = 1; i < ct; i++){
//		m_file.precision(2);
		m_file << (m_dx + (m_ratio * polyline[i].x)) << " ";
//		m_file.precision(2);
		m_file << (m_dy + (m_ratio * polyline[i].y)) << " lineto" << endl;
	}
	double red, green;
	switch(type){
	case 0:
		red = 0.0;	green = 0.2;	break;
	case 1:
		red = 1.0;	green = 0.0;	break;
	case 2:
		red = 0.2;	green = 0.2;	break;
	default:
		red = 0.5;	green = 0.5;	break;
	}
	m_file << m_line_width << " setlinewidth ";
//	m_file.precision(3);
	m_file << red << " ";
//	m_file.precision(3);
	m_file << green << " 0.000 setrgbcolor stroke" << endl;
}

//////////////////////////////////////////////////////////////////////
// Generuje polecenia jêzyka postscript powoduj¹ce narysowanie 
//	zadanych czworok¹tów w kolorze okreœlonym przez identyfikator
//	obszaru i wspó³czynnik jakoœci,
void EPSFile::drawQuad(const DPoint2d &a, const DPoint2d &b, const DPoint2d &c, const DPoint2d &d, int area_id, double quality)
{
	if(area_id < 0) return;
	if(m_viewport.valid){
		if(!m_viewport.contains(a)) return;
		if(!m_viewport.contains(b)) return;
		if(!m_viewport.contains(c)) return;
		if(!m_viewport.contains(d)) return;
	}
	m_file << "newpath ";
//	m_file.precision(2);
	m_file << (m_dx + (m_ratio * a.x)) << " ";
//	m_file.precision(2);
	m_file << (m_dy + (m_ratio * a.y)) << " moveto ";
//	m_file.precision(2);
	m_file << (m_dx + (m_ratio * b.x)) << " ";
//	m_file.precision(2);
	m_file << (m_dy + (m_ratio * b.y)) << " lineto ";
//	m_file.precision(2);
	m_file << (m_dx + (m_ratio * c.x)) << " ";
//	m_file.precision(2);
	m_file << (m_dy + (m_ratio * c.y)) << " lineto ";
//	m_file.precision(2);
	m_file << (m_dx + (m_ratio * d.x)) << " ";
//	m_file.precision(2);
	m_file << (m_dy + (m_ratio * d.y)) << " lineto closepath" << endl;
	double df;
	double red		= quality * modf(1.234 * area_id, &df);
	double green	= quality * modf(2.345 * area_id, &df);
	double blue		= quality * modf(3.456 * area_id, &df);
//	m_file.precision(3);
	m_file << red << " ";
//	m_file.precision(3);
	m_file << green << " ";
//	m_file.precision(3);
	m_file << blue << " setrgbcolor fill" << endl;
}

//////////////////////////////////////////////////////////////////////
// Generuje polecenia jêzyka postscript powoduj¹ce wydrukowanie
//	zadanej linii 
void EPSFile::drawLine(const DPoint2d &pt1, const DPoint2d &pt2, bool marked)
{
	if(m_viewport.valid){
		if(!m_viewport.contains(pt1)) return;
		if(!m_viewport.contains(pt2)) return;
	}
	m_file << "newpath ";
//	m_file.precision(2);
	m_file << (m_dx + (m_ratio * pt1.x)) << " ";
//	m_file.precision(2);
	m_file << (m_dy + (m_ratio * pt1.y)) << " moveto ";
//	m_file.precision(2);
	m_file << (m_dx + (m_ratio * pt2.x)) << " ";
//	m_file.precision(2);
	m_file << (m_dy + (m_ratio * pt2.y)) << " lineto ";
	if(marked){
		m_file << m_line_width << " setlinewidth 0.35 0.0 0.0 setrgbcolor stroke" << endl;
	}else{
		m_file << m_line_width << " setlinewidth 0.0 setgray stroke" << endl;
	}
}

//////////////////////////////////////////////////////////////////////
// Generuje polecenia jêzyka postscript powoduj¹ce wydrukowanie
//	zadanej linii 
void EPSFile::drawLineGray(const DPoint2d &pt1, const DPoint2d &pt2, double gray)
{
	if(m_viewport.valid){
		if(!m_viewport.contains(pt1)) return;
		if(!m_viewport.contains(pt2)) return;
	}
	m_file << "newpath ";
//	m_file.precision(2);
	m_file << (m_dx + (m_ratio * pt1.x)) << " ";
//	m_file.precision(2);
	m_file << (m_dy + (m_ratio * pt1.y)) << " moveto ";
//	m_file.precision(2);
	m_file << (m_dx + (m_ratio * pt2.x)) << " ";
//	m_file.precision(2);
	m_file << (m_dy + (m_ratio * pt2.y)) << " lineto ";
	m_file << m_line_width << " setlinewidth " << gray << " setgray stroke" << endl;
}

//////////////////////////////////////////////////////////////////////
// Generuje polecenia jêzyka postscript powoduj¹ce wydrukowanie
//	zadanej linii 
void EPSFile::drawLineRGB(const DPoint2d &pt1, const DPoint2d &pt2, double r, double g, double b)
{
	if(m_viewport.valid){
		if(!m_viewport.contains(pt1)) return;
		if(!m_viewport.contains(pt2)) return;
	}
	m_file << "newpath ";
//	m_file.precision(2);
	m_file << (m_dx + (m_ratio * pt1.x)) << " ";
//	m_file.precision(2);
	m_file << (m_dy + (m_ratio * pt1.y)) << " moveto ";
//	m_file.precision(2);
	m_file << (m_dx + (m_ratio * pt2.x)) << " ";
//	m_file.precision(2);
	m_file << (m_dy + (m_ratio * pt2.y)) << " lineto ";
	m_file << m_line_width << " setlinewidth " << r << ' ' << g << ' ' << b << " setrgbcolor stroke" << endl;
}

//////////////////////////////////////////////////////////////////////
// Generuje polecenia jêzyka postscript powoduj¹ce wydrukowanie
//	zadanego punktu (kó³eczka) na okreœlonej pozycji (w czarnym kolorze)
void EPSFile::drawPoint(const DPoint2d &pt, double r, double gray_color)
{
	if(m_viewport.valid){
		if(!m_viewport.contains(pt)) return;
	}
	DPoint2d p(m_dx + m_ratio * pt.x, m_dy + m_ratio * pt.y);
	m_file << "newpath ";
//	m_file.precision(2);
	m_file << (p.x - r) << " ";
//	m_file.precision(2);
	m_file << p.y << " moveto ";
//	m_file.precision(2);
	m_file << (p.x - r) << " ";
//	m_file.precision(2);
	m_file << (p.y + 1.5*r) << " ";
//	m_file.precision(2);
	m_file << (p.x + r) << " ";
//	m_file.precision(2);
	m_file << (p.y + 1.5*r) << " ";
//	m_file.precision(2);
	m_file << (p.x + r) << " ";
//	m_file.precision(2);
	m_file << p.y << " curveto" << endl;
//	m_file.precision(2);
	m_file << (p.x + r) << " ";
//	m_file.precision(2);
	m_file << (p.y - 1.5*r) << " ";
//	m_file.precision(2);
	m_file << (p.x - r) << " ";
//	m_file.precision(2);
	m_file << (p.y - 1.5*r) << " ";
//	m_file.precision(2);
	m_file << (p.x - r) << " ";
//	m_file.precision(2);
	m_file << p.y << " curveto closepath ";
	if(gray_color < 0.0){
		const char* rgb_desc[] = {"0.9 0 0", "0 0.9 0", "0 0 0.9", "0.9 0.9 0", "0 0.9 0.9", "0.9 0 0.9"};
		int ind = (int)(-gray_color-0.5);
		assert(ind >= 0 && ind < 6);
		m_file << rgb_desc[ind] << " setrgbcolor fill" << endl;
	}else if(gray_color < 1.0) 
		m_file << gray_color << " setgray fill" << endl;
	else m_file << "0.0 setgray stroke" << endl;
}

//////////////////////////////////////////////////////////////////////
// Generuje polecenia jêzyka postscript powoduj¹ce wydrukowanie
//	zadanego punktu (kó³eczka) na okreœlonej pozycji (w czarnym kolorze)
void EPSFile::drawPointCross(const DPoint2d &pt, double r, double gray_color)
{
	if(m_viewport.valid){
		if(!m_viewport.contains(pt)) return;
	}
	DPoint2d p(m_dx + m_ratio * pt.x, m_dy + m_ratio * pt.y);
	m_file << "newpath ";
	m_file << (p.x - r) << " ";
	m_file << (p.y - r) << " moveto ";
	m_file << (p.x + r) << " ";
	m_file << (p.y + r) << " lineto" << endl;
	m_file << (p.x + r) << " ";
	m_file << (p.y - r) << " moveto ";
	m_file << (p.x - r) << " ";
	m_file << (p.y + r) << " lineto closepath " << gray_color << " setgray stroke" << endl;
}
