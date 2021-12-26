/////////////////////////////////////////////////////////////////////////////
// EPSFile.h
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	2002-
//	Generation of unstructured meshes
/////////////////////////////////////////////////////////////////////////////

#pragma once

#if !defined(EPSFILE_H__INCLUDED)
#define EPSFILE_H__INCLUDED

#include "DPoint.h"
#include "DRect.h"
#include "common.h"

/**
 * Class implements some procedures for drawing basic images in EncapsulatedPostScript format.
 * The image is stored in a file, which descriptor is given in the constructor.
 */
class EPSFile  
{
public:
	/// Standard constructor initializing the process of crating the EPS image
	EPSFile(const char* fname, double min_x, double max_x, double min_y, double max_y, int length = 150, int width = 150);
	/// Standard constructor initializing the process of crating the EPS image
	EPSFile(const string& fname, double min_x, double max_x, double min_y, double max_y, int length = 150, int width = 150);
	/// Standard denstructor finilizing the process of crating the EPS image
	~EPSFile();
public:
	/// Draws a single point (a small circle) at the given coordinates
	void drawPoint(const DPoint2d & pt, double r = 1.0, double gray_color = 0.0);
	/// Draws a single point (a small cross) at the given coordinates
	void drawPointCross(const DPoint2d & pt, double r = 1.0, double gray_color = 0.0);
	/// Draws a line segment between the two given points (with colored line if "marked")
	void drawLine(const DPoint2d & pt1, const DPoint2d & pt2, bool marked = false);
	/// Draws a line segment between the two given points (with given gray level, 0.0 == black)
	void drawLineGray(const DPoint2d & pt1, const DPoint2d & pt2, double gray = 0.0);
	/// Draws a line segment between the two given points (with given RGB color)
	void drawLineRGB(const DPoint2d & pt1, const DPoint2d & pt2, double r = 0.0, double g = 0.0, double b = 0.0);
	/// Draws a quadrangle described by the four given points (with given color hue and intensity)
	void drawQuad(const DPoint2d &a, const DPoint2d &b, const DPoint2d &c, const DPoint2d &d, int area_id, double quality);
	/// Draws a polyline described by the array of points (with color dependent on "type")
	void drawPolyLine(DPoint2d* polyline, int ct, int type);
	/// Draws a number (in text) at the given coordinates
	void drawNumber(const DPoint2d& pt, int nr);
	/// Draws a triangle described by the four given points (with given color hue and intensity)
	void drawTriangle(const DPoint2d& a, const DPoint2d& b, const DPoint2d& c, int area_id, double quality);
	/// Sets the viewport (a bounding rectangle for this image)
	void setViewport(const DRect& viewport) { m_viewport = viewport; }
	/// Clears the viewport (a bounding rectangle for this image)
	void clearViewport() { m_viewport.valid = false; }
protected:
	/// Issues finilizing postscipt commmands
	void tail();
	/// Issues initializing postscipt commmands
	void header(int length, int width);
protected:
	/// Descriptor of the text file for the EPS image
	ofstream m_file;
	/// Transformation factor between real coordinates and "document coordinates"
	double m_ratio;
	/// X-correction of coordinates to the positive quadrant of coordinate system
	double m_dx;
	/// Y-correction of coordinates to the positive quadrant of coordinate system
	double m_dy;
	/// Line width
	double m_line_width;
	/// Bounding rectangle for the drawn image
	DRect m_viewport;
};

#endif // !defined(EPSFILE_H__INCLUDED)
