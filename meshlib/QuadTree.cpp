/////////////////////////////////////////////////////////////////////////////
// QuadTree.cpp
// Klasa reprezentuj¹ca drzewo wyszukiwania geometrycznego
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	1999/2000
//	Generacja siatek niestrukturalnych
/////////////////////////////////////////////////////////////////////////////

#include "QuadTree.h"
#include "MeshData.h"
#include "DRect.h"

int QuadTree::param_qtree_threshold = 300;

//////////////////////////////////////////////////////////////////////
// Konstruktor
QuadTree::QuadTree(const DRect& rect)
{
	count = 0;
	triangle = nullptr;
	middle = rect.getMiddlePoint();
	size.x = 0.5 * rect.getDX();
	size.y = 0.5 * rect.getDY();
	parts[0] = parts[1] = parts[2] = parts[3] = nullptr;
	max_level = 0;
//	tree_count = 1;
}

//////////////////////////////////////////////////////////////////////
// Destruktor
QuadTree::~QuadTree()
{
	for(int i =0; i < 4; i++) if(parts[i]) delete parts[i];
}

//////////////////////////////////////////////////////////////////////
// Porównuje nowy trójk¹t ze wszystkimi trójk¹tami wzd³u¿ œcie¿ki
//	w drzewie wyznaczonej przez wspó³rzêdne œrodka trójk¹ta i 
//	w jeœli trzeba dokonuje uaktualnienia preferowanego trójk¹ta
//	dla danego podprostok¹ta
void QuadTree::insertTriangleLink(MeshTriangle2d *t)
{
	DataVector<DPoint2d> points(7);
	t->getRepresentativePoints(points);
	assert(points.countInt() > 0);
	for(size_t i = 0; i < points.countInt(); i++)
		insertTriangleLink(t, points[i]);
}

void QuadTree::insertTriangleLink(MeshTriangle2d *t, const DPoint2d& pt)
{
	// Scan tree
	QuadTree *tree = this;
	int level = 0;
	while(true){
		++level;
		// At each level update best triangle if required
		double d = t->distance2(tree->middle);
		if(!tree->triangle || tree->dist > d){
			tree->triangle = t;
			tree->dist = d;
		}
		if(tree->count > -1){	// If leaf
			tree->count++;
			if(tree->count > param_qtree_threshold){
				// If enough visited - make it branch
				tree->count = -1;
			}
			if(level > max_level) max_level = level;
			return;	// Done
		}
		// Select proper subtree
		int part_index = 0;
		if(pt.x > tree->middle.x) part_index++;
		if(pt.y < tree->middle.y) part_index += 2;
		if(!tree->parts[part_index]){
			double dx = tree->size.x;
			double dy = tree->size.y;
			double x1, y1;
			if(pt.x > tree->middle.x){
				x1 = tree->middle.x;
			}else{
				x1 = tree->middle.x - dx;
			}
			if(pt.y < tree->middle.y){
				y1 = tree->middle.y - dy;
			}else{
				y1 = tree->middle.y;
			}			
			tree->parts[part_index] = new QuadTree(DRect(x1, y1 + dy, x1 + dx, y1));
//			tree_count++;
		}
		tree = tree->parts[part_index];
	}
}

//////////////////////////////////////////////////////////////////////
// Usuwa zadany trójk¹t z drzewa, o ile w nim wystêpuje
void QuadTree::removeTriangleLink(const MeshTriangle2d *t)
{
	DataVector<DPoint2d> points(7);
	t->getRepresentativePoints(points);
	assert(points.countInt() > 0);
	for(size_t i = 0; i < points.countInt(); i++)
		removeTriangleLink(t, points[i]);
}

//////////////////////////////////////////////////////////////////////
// Usuwa zadany trójk¹t z drzewa, o ile w nim wystêpuje
void QuadTree::removeTriangleLink(const MeshTriangle2d *t, const DPoint2d& pt)
{
	QuadTree *tree = this;
	do{
		// At each level check if best triangle should be removed
		if(tree->triangle == t) tree->triangle = nullptr;
		if(tree->count > -1)	// If leaf
			tree->count--;
		// Select proper subtree
		int part_index = 0;
		if(pt.x > tree->middle.x) part_index++;
		if(pt.y < tree->middle.y) part_index += 2;
		tree = tree->parts[part_index];
	}while(tree);

}

//////////////////////////////////////////////////////////////////////
// Zwraca ten z trójk¹tów znajduj¹cych siê na œcie¿ce w drzewie 
//	wyznaczonej przez punkt point, który to trójk¹t jest najbli¿ej
//	zadanego punktu
MeshTriangle2d* QuadTree::getNearestTriangle(const DPoint2d& pt, MeshTriangle2d *last_triangle) const
{
	MeshTriangle2d* best_triangle = nullptr;
	double min_distance = mesh_data.relative_infinity;
	if(last_triangle){
		min_distance = last_triangle->distance2(pt);
		best_triangle = last_triangle;
	}
		
	const QuadTree* tree = this;
	while(tree){
		if(tree->triangle){ // If triangle is defined
			double distance = tree->triangle->distance2(pt);
			if(distance < min_distance){	// ... and is nearer to the given point
				best_triangle = tree->triangle;	// make it currently-best
				min_distance = distance;
			}
		}
		// Select proper subtree
		int part_index = 0;
		if(pt.x > tree->middle.x) part_index++;
		if(pt.y < tree->middle.y) part_index += 2;
		tree = tree->parts[part_index];
	}

	return best_triangle;
}
