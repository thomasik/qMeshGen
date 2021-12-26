/////////////////////////////////////////////////////////////////////////////
// OctTree.cpp
// Klasa reprezentuj¹ca drzewo wyszukiwania geometrycznego
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	1999/2000
//	Generacja siatek niestrukturalnych
/////////////////////////////////////////////////////////////////////////////

#include "OctTree.h"
#include "MeshData.h"
#include "MeshPoint3d.h"

int OctTree::param_octtree_threshold = 1000;

//////////////////////////////////////////////////////////////////////
// Konstruktor
OctTree::OctTree(const DBox& box)
{
	count = 0;
	tetrahedron = nullptr;
	middle = box.getMiddlePoint();
	size.x = 0.5 * box.getDX();
	size.y = 0.5 * box.getDY();
	size.z = 0.5 * box.getDZ();
	for(int i = 0; i < 8; i++) parts[i] = nullptr;
	max_level = 0;
//	tree_count = 1;
}

//////////////////////////////////////////////////////////////////////
// Destruktor
OctTree::~OctTree()
{
	for(int i =0; i < 8; i++) if(parts[i]) delete parts[i];
//	if(max_level > 0){	// Czyli korzeñ
//		LOG("\n*** OctTree -> max_level = %d\n", max_level);
//		LOG("                 tree_count = %d\n", tree_count);
//	}
}

void OctTree::insertTetrahedronLink(MeshTetrahedron *tetra)
{
	DataVector<DPoint3d> points(15);
	tetra->getRepresentativePoints(points);
	assert(points.countInt() > 0);
	for(size_t i = 0; i < points.countInt(); i++)
		insertTetrahedronLink(tetra, points[i]);
}

//////////////////////////////////////////////////////////////////////
// Porównuje nowy trójk¹t ze wszystkimi trójk¹tami wzd³u¿ œcie¿ki
//	w drzewie wyznaczonej przez wspó³rzêdne œrodka trójk¹ta i 
//	w jeœli trzeba dokonuje uaktualnienia preferowanego trójk¹ta
//	dla danego podprostok¹ta
void OctTree::insertTetrahedronLink(MeshTetrahedron *tetra, const DPoint3d& pt)
{
	// Przegl¹danie drzewa
	OctTree *tree = this;
	int level = 0;
	while(true){
		level++;
		// Na ka¿dym poziomie uaktualnienie najbli¿szego czworoœcianu w miarê potrzeby
		double d = tetra->distance2(pt);
		if(!tree->tetrahedron || tree->dist > d){
			tree->tetrahedron = tetra;
			tree->dist = d;
		}
		if(tree->count > -1){	// Jeœli liœæ
			tree->count++;
			if(tree->count > param_octtree_threshold){
				// Jeœli by³ odpowiednio czêsto odwiedzany - przekszta³æ w ga³¹Ÿ
				tree->count = -1;
			}
			if(level > max_level) max_level = level;
			return;	// Koniec
		}
		// Wybierz odpowiednie poddrzewo
		int part_index = 0;
		if(pt.x > tree->middle.x) part_index += 1;
		if(pt.y > tree->middle.y) part_index += 2;
		if(pt.z > tree->middle.z) part_index += 4;
		if(!tree->parts[part_index]){
			double dx = tree->size.x;
			double dy = tree->size.y;
			double dz = tree->size.z;
			double x1, y1, z1;
			if(pt.x > tree->middle.x)
				x1 = tree->middle.x;
			else x1 = tree->middle.x - dx;
			if(pt.y > tree->middle.y)
				y1 = tree->middle.y;
			else y1 = tree->middle.y - dy;
			if(pt.z > tree->middle.z)
				z1 = tree->middle.z;
			else z1 = tree->middle.z - dz;
			tree->parts[part_index] = new OctTree(DBox(x1, x1 + dx, y1, y1 + dy, z1, z1 + dz));
//			tree_count++;
		}
		tree = tree->parts[part_index];
	}
}

void OctTree::removeTetrahedronLink(const MeshTetrahedron *tetra)
{
	DataVector<DPoint3d> points(15);
	tetra->getRepresentativePoints(points);
	assert(points.countInt() > 0);
	for(size_t i = 0; i < points.countInt(); i++)
		removeTetrahedronLink(tetra, points[i]);
}

//////////////////////////////////////////////////////////////////////
// Usuwa zadany trójk¹t z drzewa, o ile w nim wystêpuje
void OctTree::removeTetrahedronLink(const MeshTetrahedron *tetra, const DPoint3d& pt)
{
	OctTree *tree = this;
	do{
		// Na ka¿dym poziomie uaktualnienie najbli¿szego trójk¹ta
		if(tree->tetrahedron == tetra) tree->tetrahedron = nullptr;
		if(tree->count > -1){	// Jeœli liœæ
			tree->count--;
		}
		// Wybierz odpowiednie poddrzewo
		int part_index = 0;
		if(pt.x > tree->middle.x) part_index += 1;
		if(pt.y > tree->middle.y) part_index += 2;
		if(pt.z > tree->middle.z) part_index += 4;
		tree = tree->parts[part_index];
	}while(tree);

}

//////////////////////////////////////////////////////////////////////
// Zwraca ten z trójk¹tów znajduj¹cych siê na œcie¿ce w drzewie 
//	wyznaczonej przez punkt point, który to trójk¹t jest najbli¿ej
//	zadanego punktu
MeshTetrahedron* OctTree::getNearestTetrahedron(const DPoint3d& pt, MeshTetrahedron *last_tetrahedron) const
{
	MeshTetrahedron* best_tetrahedron = nullptr;
	double min_distance = mesh_data.relative_infinity;
	if(last_tetrahedron){
		min_distance = last_tetrahedron->distance2(pt);
		best_tetrahedron = last_tetrahedron;
	}
	int counter = 0;
	const OctTree* tree = this;
	while(tree){
		++counter;
		if(tree->tetrahedron){ // Jeœli trójk¹t jest zdefiniowany
			double distance = tree->tetrahedron->distance2(pt);
			if(distance < min_distance){	// ... i jest bli¿ej danego punktu
				best_tetrahedron = tree->tetrahedron;	// niech stanie siê aktualnie najlepszym
				min_distance = distance;
			}
		}
		// Wybierz odpowiednie poddrzewo
		int part_index = 0;
		if(pt.x > tree->middle.x) part_index += 1;
		if(pt.y > tree->middle.y) part_index += 2;
		if(pt.z > tree->middle.z) part_index += 4;
		tree = tree->parts[part_index];
	}
//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "TT: " << counter);
	return best_tetrahedron;
}

///////////////////////////////////////////////////////////////

#define OCT_TREE_MESH_POINTS_THRESHOLD 10

OctTreeMeshPoints::OctTreeMeshPoints(const DBox& box, OctTreeMeshPoints* father) 
	: m_points(OCT_TREE_MESH_POINTS_THRESHOLD), m_father(father)
{
	m_middle = box.getMiddlePoint();
	m_size.x = 0.5 * box.getDX();
	m_size.y = 0.5 * box.getDY();
	m_size.z = 0.5 * box.getDZ();
	m_parts[0] = m_parts[1] = nullptr;
}

//////////////////////////////////////////////////////////////////////
// Destruktor
OctTreeMeshPoints::~OctTreeMeshPoints()
{
	if(m_parts[0]) delete m_parts[0];
	if(m_parts[1]) delete m_parts[1];
}

void OctTreeMeshPoints::insertMeshPointLink(MeshPoint3d *point)
{
	// Scan the tree
	OctTreeMeshPoints *tree = findLeaf(point->getCoordinates());

	while(tree->m_points.countInt() >= OCT_TREE_MESH_POINTS_THRESHOLD){
		// create sub-trees
		if(tree->m_size.x >= tree->m_size.y && tree->m_size.x >= tree->m_size.z){
			tree->m_split = 0;
			tree->m_parts[0] = new OctTreeMeshPoints(DBox(
				tree->m_middle.x - tree->m_size.x, tree->m_middle.x, 
				tree->m_middle.y - tree->m_size.y, tree->m_middle.y + tree->m_size.y, 
				tree->m_middle.z - tree->m_size.z, tree->m_middle.z + tree->m_size.z), tree);
			tree->m_parts[1] = new OctTreeMeshPoints(DBox(
				tree->m_middle.x, tree->m_middle.x + tree->m_size.x, 
				tree->m_middle.y - tree->m_size.y, tree->m_middle.y + tree->m_size.y, 
				tree->m_middle.z - tree->m_size.z, tree->m_middle.z + tree->m_size.z), tree);
		}else if(tree->m_size.y >= tree->m_size.x && tree->m_size.y >= tree->m_size.z){
			tree->m_split = 1;
			tree->m_parts[0] = new OctTreeMeshPoints(DBox(
				tree->m_middle.x - tree->m_size.x, tree->m_middle.x + tree->m_size.x, 
				tree->m_middle.y - tree->m_size.y, tree->m_middle.y, 
				tree->m_middle.z - tree->m_size.z, tree->m_middle.z + tree->m_size.z), tree);
			tree->m_parts[1] = new OctTreeMeshPoints(DBox(
				tree->m_middle.x - tree->m_size.x, tree->m_middle.x + tree->m_size.x, 
				tree->m_middle.y, tree->m_middle.y + tree->m_size.y, 
				tree->m_middle.z - tree->m_size.z, tree->m_middle.z + tree->m_size.z), tree);
		}else{
			tree->m_split = 2;
			tree->m_parts[0] = new OctTreeMeshPoints(DBox(
				tree->m_middle.x - tree->m_size.x, tree->m_middle.x + tree->m_size.x, 
				tree->m_middle.y - tree->m_size.y, tree->m_middle.y + tree->m_size.y, 
				tree->m_middle.z - tree->m_size.z, tree->m_middle.z ), tree);
			tree->m_parts[1] = new OctTreeMeshPoints(DBox(
				tree->m_middle.x - tree->m_size.x, tree->m_middle.x + tree->m_size.x, 
				tree->m_middle.y - tree->m_size.y, tree->m_middle.y + tree->m_size.y, 
				tree->m_middle.z, tree->m_middle.z + tree->m_size.z), tree);
		}
		// split points
		for(size_t i = 0; i < tree->m_points.countInt(); i++){
			MeshPoint3d * mp = tree->m_points[i];
			switch(tree->m_split){
				case 0: tree->m_parts[ (mp->getCoordinates().x < tree->m_middle.x) ? 0 : 1 ]->m_points.add(mp); break;
				case 1: tree->m_parts[ (mp->getCoordinates().y < tree->m_middle.y) ? 0 : 1 ]->m_points.add(mp); break;
				case 2: tree->m_parts[ (mp->getCoordinates().z < tree->m_middle.z) ? 0 : 1 ]->m_points.add(mp); break;
				default: assert(false);
			}
		}
		tree->m_points.clear();
		// step down
		tree = tree->findLeaf(point->getCoordinates());
	}
	// insert
	tree->m_points.add(point);
}

void OctTreeMeshPoints::removeMeshPointLink(MeshPoint3d *point)
{
	// Scan the tree
	OctTreeMeshPoints *tree = findLeaf(point->getCoordinates());
	// check and remove
	tree->m_points.remove(point);
}

OctTreeMeshPoints* OctTreeMeshPoints::findLeaf(const DPoint3d& pt)
{
	OctTreeMeshPoints *tree = this;
//	int level = 0;
	while(tree->m_parts[0]){ // while not leaf
//		++level;
		switch(tree->m_split){
			case 0: tree = tree->m_parts[ pt.x < tree->m_middle.x ? 0 : 1 ]; break;
			case 1: tree = tree->m_parts[ pt.y < tree->m_middle.y ? 0 : 1 ]; break;
			case 2: tree = tree->m_parts[ pt.z < tree->m_middle.z ? 0 : 1 ]; break;
			default:
				assert(false);
		}
		assert(tree);
	}
	return tree;
}

const OctTreeMeshPoints* OctTreeMeshPoints::findLeaf(const DPoint3d& pt) const
{
	const OctTreeMeshPoints *tree = this;
	while(tree->m_parts[0]){ // while not leaf
		switch(tree->m_split){
			case 0: tree = tree->m_parts[ pt.x < tree->m_middle.x ? 0 : 1 ]; break;
			case 1: tree = tree->m_parts[ pt.y < tree->m_middle.y ? 0 : 1 ]; break;
			case 2: tree = tree->m_parts[ pt.z < tree->m_middle.z ? 0 : 1 ]; break;
			default:
				assert(false);
		}
		assert(tree);
	}
	return tree;
}

bool OctTreeMeshPoints::anyMeshPointInProximity(Metric3dContext& mc, const DPoint3d& pt, double max_metric_distance) const
{
	// find leaf
	const OctTreeMeshPoints *tree = findLeaf(pt);
	const DMPoint3d dpt = mc.transformRStoMS(pt);
	// check in-points
	const double max_metric_distance2 = sqr(max_metric_distance);
	if(tree->anyMeshPointInProximityForLeaf(mc, dpt, max_metric_distance2)) return true;
	// check uptree
	double h_max = max_metric_distance * mc.maxLength();
	while(tree->m_father){
		OctTreeMeshPoints* other_tree = (tree->m_father->m_parts[0] == tree) 
			? tree->m_father->m_parts[1] : tree->m_father->m_parts[0];
		if(other_tree->anyMeshPointInProximityForNearLeaves(mc, pt, dpt, h_max, max_metric_distance2)) return true;
		tree = tree->m_father;
	}
	return false;
}

bool OctTreeMeshPoints::anyMeshPointInProximityForLeaf(Metric3dContext& mc, const DMPoint3d& dpt, double max_metric_distance2) const
{
	for(size_t i = 0; i < m_points.countInt(); i++)
		if(dpt.distance2(m_points.get(i)->getMetricCoordinates(mc)) < max_metric_distance2) 
			return true;
	return false;
}

bool OctTreeMeshPoints::anyMeshPointInProximityForNearLeaves(Metric3dContext& mc, 
		const DPoint3d& pt, const DMPoint3d& dpt, double h_max, double max_metric_distance2) const
{
	double safe_dist2 = sqr(h_max + m_size.length());
	if(pt.distance2(m_middle) > safe_dist2) return false;
	if(m_parts[0]){ // if not leaf
		if(m_parts[0]->anyMeshPointInProximityForNearLeaves(mc, pt, dpt, h_max, max_metric_distance2)) return true;
		else return m_parts[1]->anyMeshPointInProximityForNearLeaves(mc, pt, dpt, h_max, max_metric_distance2);
	}else{
		return anyMeshPointInProximityForLeaf(mc, dpt, max_metric_distance2);
	}
}
