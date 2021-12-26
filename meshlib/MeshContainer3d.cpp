/////////////////////////////////////////////////////////////////////////////
// MeshContainer3d.cpp
// Klasa odpowiedzialna za przechowywanie punktów i elementów siatki
/////////////////////////////////////////////////////////////////////////////
//	Tomasz Jurczyk,	1999/2000
//	Generacja siatek niestrukturalnych
/////////////////////////////////////////////////////////////////////////////


#define USE_OPENMP_HERE

#include "MeshContainer3d.h"
#include "MeshPoint3d.h"
#include "DataContainer.h"
#include "Curve2dParametric.h"
#include "Curve3dParametric.h"
#include "MeshEdge3d.h"
#include "MeshDomainEdge3d.h"
#include "MeshTriangle3d.h"
#include "MeshDomainVolume.h"
#include "MeshDomainSurface.h"
#include "OctTree.h"
#include "MeshViewSet.h"
#include "ControlSpace2dIdentity.h"
#include "ControlSpace3dMatrixUniform.h"
#include "ControlSpace3dOctree.h"
#include "DMetric3d.h"
#include "MeshBoundaryCondition.h"
#include "MeshContainer2d.h"
#include "ControlSpace2d.h"
#include "MeshGenerator2d.h"
#include "DataHashTable.h"
#include "DTriangle.h"
#include "DTetrahedron.h"

/////////////////////////////////////////////////////////////////////////////
// Standardowy konstruktor (part_size) okreœla mno¿nik rozmiaru tablic, które
//	zwiêkszane s¹ w razie potrzeby krokowo o t¹ w³aœnie wartoœæ
MeshContainer3d::MeshContainer3d(int part_size, 
		MeshContainer3d* total_volume_mesh, 
		MeshContainer3dSurface* total_surface_mesh) :
	m_part_size(part_size), m_last_tetrahedron(nullptr),
	m_constraining(MeshData::CONSTRAIN_NONE), m_discretization_state(0)
{
	m_points = new DataContainer<MeshPoint3d>(part_size);
	m_blocks = new DataContainer<MeshBlock>(part_size);

	if(total_volume_mesh) total_mdv.setMesh(total_volume_mesh);
	if(total_surface_mesh) total_mdv.setSurfaceMesh(total_surface_mesh);
}

/////////////////////////////////////////////////////////////////////////////
// Destruktor
MeshContainer3d::~MeshContainer3d()
{
	deleteAll();
	delete m_blocks;
	delete m_points;
}

/////////////////////////////////////////////////////////////////////////////
// Usuwa dane dotycz¹ce przechowywanej geometrii obszaru
void MeshContainer3d::deleteAll()
{
	m_blocks->deleteAll();
	m_points->deleteAll();
}

//////////////////////////////////////////////////////////////////////
// Zwraca wspó³rzêdne prostok¹ta obejmuj¹cego wszystkie punkty 
DBox MeshContainer3d::getBoundingBox() const
{
	DBox box;

	int count = m_points->countInt();
	for(int i = 0; i < count; i++)
		box.addPoint(m_points->getDataAt(i)->getCoordinates());

	count = m_blocks->countInt();
	for(int i = 0; i < count; i++)
		box.addBox(m_blocks->getDataAt(i)->getBoundingBox());

	box.addBox(total_mdv.getBoundingBox());

	return box;
}

void MeshContainer3d::clearDiscretization1d()
{
	// 1d
	int i, count = m_points->countInt();
	for(i = 0; i < count; i++){
		MeshPoint3d* point = m_points->getDataAt(i);
		int rank = point->getRank();
		for(int j = 0; j < rank; j++){
			MeshDomainEdge3d* edge = (MeshDomainEdge3d*)point->getEdge(j);
			if(edge->getType() == EDGE_DOMAIN_3D)
				edge->clearDiscretization();
		}
	}
}

void MeshContainer3d::clearDiscretization2d()
{
	int i, j, count = m_blocks->countInt();
	for(i = 0; i < count; i++){
		MeshDomainVolume* volume = (MeshDomainVolume*)m_blocks->getDataAt(i);
		if(volume->getType() == BLOCK_DOMAIN){			
			int face_count = volume->getFaceCount();
			for(j = 0; j < face_count; j++){
				MeshDomainSurface* surface = (MeshDomainSurface*)volume->getFace(j);
				surface->clearDiscretization();
			}
		}
	}
}

void MeshContainer3d::clearCS2d()
{
	int i, j, count = m_blocks->countInt();
	for(i = 0; i < count; i++){
		MeshDomainVolume* volume = (MeshDomainVolume*)m_blocks->getDataAt(i);
		if(volume->getType() == BLOCK_DOMAIN){			
			int face_count = volume->getFaceCount();
			for(j = 0; j < face_count; j++){
				MeshDomainSurface* surface = (MeshDomainSurface*)volume->getFace(j);
				surface->clearControlSpace();
			}
		}
	}
}

void MeshContainer3d::clearDiscretization3d()
{
	int i, count = m_blocks->countInt();
	for(i = 0; i < count; i++){
		MeshDomainVolume* volume = (MeshDomainVolume*)m_blocks->getDataAt(i);
		if(volume->getType() == BLOCK_DOMAIN){
			volume->clearDiscretization();
		}
	}
}

void MeshContainer3d::clearCS3d()
{
	int i, count = m_blocks->countInt();
	for(i = 0; i < count; i++){
		MeshDomainVolume* volume = (MeshDomainVolume*)m_blocks->getDataAt(i);
		if(volume->getType() == BLOCK_DOMAIN){
			volume->clearControlSpace();
		}
	}
}

void MeshContainer3d::clearDiscretization()
{
	clearDiscretization3d();
	clearDiscretization2d();
	clearDiscretization1d();
}

void MeshContainer3d::clearCS()
{
	clearCS3d();
	clearCS2d();
}

MeshContainer2d* MeshContainer3d::getFirst2dBoundary()
{
	int block_count = getBlocksCount();
	for(int i = 0; i < block_count; i++){
		MeshDomainVolume* volume = (MeshDomainVolume*)getBlockAt(i);
		assert(volume && (volume->getType() == BLOCK_DOMAIN));
		int face_count = volume->getFaceCount();
		for(int j = 0; j < face_count; j++){
			MeshDomainSurface* domain_surface = (MeshDomainSurface*)volume->getFace(j);
			assert(domain_surface && (domain_surface->getType() == FACE_DOMAIN));
			MeshContainer2d* mesh = domain_surface->getBoundary();
			if(mesh) return mesh;
		}
	}
	return nullptr;
}

bool MeshContainer3d::isValid() const
{
	/// TODO -> remove, placed here only temporary
	return true;

	int bct = getBlocksCount();
	for(int i = 0; i < bct; i++){
		MeshBlock* block = getBlockAt(i);
		if(!block){
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "invalid> empty block?");
			return false;
		}
		if(m_constraining == MeshData::CONSTRAIN_DONE && block->getAreaID() < 0){
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "invalid> constrained block with negative area id");
			return false;
		}
		if(block->isInverted()){
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "invalid> inverted block, vol = " << block->getVolumeNoMetric());
			return false;
		}
		int fct = block->getFaceCount();
		if(fct < 4 || fct > 6){
			LOG4CPLUS_INFO(MeshLog::logger_mesh, "invalid> block with " << fct << " faces?");
			return false;
		}
		for(int j = 0; j < fct; j++){
			MeshFace* face = block->getFace(j);
			if(!face){
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "invalid> empty face?");
				return false;
			}
			MeshTetrahedron* tetrahedron = (MeshTetrahedron*)face->getBlock(0);
			if(tetrahedron && !tetrahedron->properOrientation(face->getPoint(0), face->getPoint(1), face->getPoint(2))){
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "invalid> face with wrong oriented 0-block");
				return false;
			}
			tetrahedron = (MeshTetrahedron*)face->getBlock(1);
			if(tetrahedron && !tetrahedron->properOrientation(face->getPoint(2), face->getPoint(1), face->getPoint(0))){
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "invalid> face with wrong oriented 1-block");
				return false;
			}
			int fect = face->getEdgeCount();
			if(fect < 3 || fect > 4){
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "invalid> face with " << fect << " edges?");
				return false;
			}
			if(face->getBlockIndex(block) < 0){
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "invalid> missing face->block link");
				return false;
			}
			if(!face->isBorder() && (face->getBlock(0) == nullptr || face->getBlock(1) == nullptr)){
				LOG4CPLUS_INFO(MeshLog::logger_mesh, "invalid> non-border face not adjacent to 2 blocks");
				return false;
			}
			//--------
			if(!face->isBorder()){
				MeshTetrahedron* t0 = (MeshTetrahedron*)face->getBlock(0);
				MeshTetrahedron* t1 = (MeshTetrahedron*)face->getBlock(1);
				MeshPoint3d* mpt0 = t0->getOppositePoint(face);
				MeshPoint3d* mpt1 = t1->getOppositePoint(face);
				double vol0 = DTriangle3d::orient3d(
					face->getPoint(0)->getCoordinates(),
					face->getPoint(1)->getCoordinates(), 
					face->getPoint(2)->getCoordinates(),
					mpt0->getCoordinates());
				double vol1 = DTriangle3d::orient3d(
					face->getPoint(0)->getCoordinates(),
					face->getPoint(1)->getCoordinates(), 
					face->getPoint(2)->getCoordinates(),
					mpt1->getCoordinates());
				if(vol0 * vol1 > 0.0){
					LOG4CPLUS_INFO(MeshLog::logger_mesh, "invalid> non-border face has two blocks on same side physically");
					LOG4CPLUS_INFO(MeshLog::logger_mesh, "t0 inverted = " << t0->isInverted());
					LOG4CPLUS_INFO(MeshLog::logger_mesh, "volume0 = " << vol0);
					LOG4CPLUS_INFO(MeshLog::logger_mesh, "t1 inverted = " << t1->isInverted());
					LOG4CPLUS_INFO(MeshLog::logger_mesh, "volume1 = " << vol1);
					return false;
				}
			}
			//--------
			if(face->isBorder()){
				int ect = face->getEdgeCount();
				for(int k = 0; k < ect; k++)
					if(!face->getEdge(k)->isBorder()){
						LOG4CPLUS_INFO(MeshLog::logger_mesh, "invalid> border face adjacent to non-border edge");
						return false;
					}
			}
		}
	}

	// check stranded nodes
	int pct = getPointsCount();
	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = getPointAt(i);
		int point_rank = point->getRank();
		if(point_rank == 0){
			LOG4CPLUS_INFO(MeshLog::logger_mesh, 
				"invalid> stranded point! - index=" << point->getIndex() << 
				(point->isBorder()?" boundary":" inner"));	
			return false;
		}
		// ... stranded edges
		for(int j = 0; j < point_rank; j++){
			MeshEdge3d* edge = point->getEdge(j);
			int edge_rank = edge->getFaceCount();
			if(edge_rank == 0){
				LOG4CPLUS_INFO(MeshLog::logger_mesh, 
					"invalid> stranded edge! - " << (edge->isBorder()?"boundary":"inner"));	
				return false;
			}
			// ... stranded faces
			for(int k = 0; k < edge_rank; k++){
				MeshFace* face = edge->getFaceAt(k);
				if(!face->isBounded()){
					LOG4CPLUS_INFO(MeshLog::logger_mesh, 
						"invalid> stranded face! - " << (face->isBorder()?"boundary":"inner"));	
					return false;
				}
				
			}			
		}
	}
	// check stranded faces

	return true;
}

void MeshContainer3d::clearSearchTree()
{ 
	m_oct_tree.reset();
}

void MeshContainer3d::setSearchTree(std::shared_ptr<OctTree> tree)
{
	m_oct_tree = tree;
}

int MeshContainer3d::getMaxSearchTreeLevel() const
{
	return m_oct_tree ? m_oct_tree->getMaxLevel() : 0;
}



int MeshContainer3d::addMeshTetrahedron(MeshTetrahedron *tetrahedron)
{
	m_last_tetrahedron = tetrahedron;
	if(m_oct_tree) m_oct_tree->insertTetrahedronLink(tetrahedron);
	return m_blocks->addDataItem(tetrahedron);
}

MeshTetrahedron* MeshContainer3d::removeMeshTetrahedron(MeshTetrahedron* tetrahedron)
{ 
	if(m_oct_tree) m_oct_tree->removeTetrahedronLink(tetrahedron);
	if(m_last_tetrahedron == tetrahedron) m_last_tetrahedron = nullptr;
	return (MeshTetrahedron*)m_blocks->removeDataItem(tetrahedron->getIndex());
}

MeshTetrahedron* MeshContainer3d::getNearTetrahedron(const DPoint3d& pt)
{ 
	if(m_oct_tree)
		return m_oct_tree->getNearestTetrahedron(pt, m_last_tetrahedron);
	else return m_last_tetrahedron; 
}

MeshViewSet* MeshContainer3d::getDebugViewSet(const MeshBlock* el1, 
				const MeshBlock* el2, double radius) const
{
	if(MeshViewSet::param_show_visualization == 0) return nullptr;
	if(!el1) return nullptr;

	double max_len = 0.0;
	int edge_count = el1->getEdgeCount();
	for(int i = 0; i < edge_count; i++)
		max_len = std::max(max_len, el1->getEdge(i)->getLengthNoMetric());
	if(el2){
		edge_count = el2->getEdgeCount();
		for(int i = 0; i < edge_count; i++)
			max_len = std::max(max_len, el2->getEdge(i)->getLengthNoMetric());
	}
	const DPoint3d middle = el1->getMiddlePoint();
	double r2 = sqr(radius*max_len);

	int bct = getBlocksCount();
	int pct = getPointsCount();

	ofstream of("debug_view.log");

	MeshViewSet *view_set = new MeshViewSet(pct, 3*bct, 0, bct);
	view_set->setMesh(this);

	// blocks
	for(int i = 0; i < bct; i++){
		MeshBlock* el = getBlockAt(i);
		if(middle.distance2(el->getMiddlePoint()) > r2) continue;
		if(el == el1)
			view_set->addBlockWithEdges(el, 1);
		else if(el == el2)
			view_set->addBlockWithEdges(el, 2);
		else
			view_set->addBlockWithEdges(el, 0);
		// log
		int epct = el->getPointCount();
		of << "Block index=" << el->getIndex() << " pts=" << epct << " [";
		for(int j = 0; j < epct; j++) of << el->getPoint(j)->getIndex() << " ";
		of << "]" << endl;
		of << "\tinverted=" << (el->isInverted()?"true":"false") << endl;
		of << "\tarea-0123=" << DTetrahedron::volume(
			el->getPoint(0)->getCoordinates(),
			el->getPoint(1)->getCoordinates(),
			el->getPoint(2)->getCoordinates(),
			el->getPoint(3)->getCoordinates()) << endl;
	}

	// points
	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = getPointAt(i);
		if(middle.distance2(point->getCoordinates()) > r2) continue;
		view_set->addPoint(point);
	}

	return view_set;
}

MeshViewSet* MeshContainer3d::getDebugViewSet(const MeshPoint3d* pt1, 
				const MeshPoint3d* pt2, double radius) const
{
	if(MeshViewSet::param_show_visualization == 0) return nullptr;
	if(!pt1) return nullptr;

	const DPoint3d middle = pt1->getCoordinates();
	double r2 = 0.0;
	if(pt2) r2 = middle.distance2(pt2->getCoordinates());
	else{
		int rank = pt1->getRank();
		for(int i = 0; i < rank; i++)
			r2 = std::max(r2, sqr(pt1->getEdge(i)->getLengthNoMetric()));
	}
	r2 *= sqr(radius);

	int bct = getBlocksCount();
	int pct = getPointsCount();

	ofstream of("debug_view.log");

	MeshViewSet *view_set = new MeshViewSet(pct, 3*bct, 0, bct);
	view_set->setMesh(this);

	// blocks
	for(int i = 0; i < bct; i++){
		MeshBlock* el = getBlockAt(i);
		if(middle.distance2(el->getMiddlePoint()) > r2) continue;
		view_set->addBlockWithEdges(el, 0);
		// log
		int epct = el->getPointCount();
		of << "Block index=" << el->getIndex() << " pts=" << epct << " [";
		for(int j = 0; j < epct; j++) of << el->getPoint(j)->getIndex() << " ";
		of << "]" << endl;
		of << "\tinverted=" << (el->isInverted()?"true":"false") << endl;
		of << "\tarea-0123=" << DTetrahedron::volume(
			el->getPoint(0)->getCoordinates(),
			el->getPoint(1)->getCoordinates(),
			el->getPoint(2)->getCoordinates(),
			el->getPoint(3)->getCoordinates()) << endl;
	}

	// others
	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = getPointAt(i);
		if(middle.distance2(point->getCoordinates()) > r2) continue;
		if(point == pt1)
			view_set->addPoint(point, 3);
		else if(point == pt2)
			view_set->addPoint(point, 4);
		else
			view_set->addPoint(point, 0);
	}

	return view_set;
}

MeshViewSet* MeshContainer3d::getDebugViewSet(const MeshContainer3d* boundary) const
{
	if(MeshViewSet::param_show_visualization == 0) return nullptr;

	//	ofstream of("debug_view.log");

	int pct = boundary->getPointsCount();
	MeshViewSet *view_set = new MeshViewSet(pct, 6*pct, 6*pct, 0);
	view_set->setMesh(this);

	int bct = boundary->getBlocksCount();
	for(int i = 0; i < bct; i++){
		const MeshBlock* block = boundary->getBlockAt(i);
		int fct = block->getFaceCount();
		for(int j = 0; j < fct; j++){
			const MeshFace* face = block->getFace(j);
			const MeshPoint3d* bpt0 = face->getPoint(0);
			const MeshPoint3d* bpt1 = face->getPoint(1);
			const MeshPoint3d* bpt2 = face->getPoint(2);
			const MeshPoint3d* mpt0 = getPointAt(bpt0->getIndex() + 8);
			const MeshPoint3d* mpt1 = getPointAt(bpt1->getIndex() + 8);
			const MeshPoint3d* mpt2 = getPointAt(bpt2->getIndex() + 8);
			const MeshFace* mface = mpt0->getFaceToPoints(mpt1, mpt2);
			view_set->addFace(face, mface ? 1 : 0, 0.95);
		}
		int fpct = block->getPointCount();
		for(int j = 0; j < fpct; j++){
			const MeshPoint3d* point = block->getPoint(j);
			int rank = point->getRank();
			for(int k = 0; k < rank; k++){
				const MeshEdge3d* edge = point->getEdge(k);
				if(edge->getPointIndex(point) > 0) continue;
				const MeshPoint3d* bpt0 = edge->getMeshPoint(0);
				const MeshPoint3d* bpt1 = edge->getMeshPoint(1);
				const MeshPoint3d* mpt0 = getPointAt(bpt0->getIndex() + 8);
				const MeshPoint3d* mpt1 = getPointAt(bpt1->getIndex() + 8);
				const MeshEdge3d* medge = mpt0->getEdgeToPoint(mpt1);
				if(!medge) view_set->addEdge(edge);
			}
		}
	}
	int mpct = getPointsCount();
	for(int j = pct+8; j < mpct; j++){
		view_set->addPoint(getPointAt(j), 100);
//	for(int j = 0; j < mpct; j++){
//		view_set->addPoint(getPointAt(j));
	}

	return view_set;
}

MeshViewSet* MeshContainer3d::getViewSetWithVisibleBlocksOnly(
	const DataVector<MeshViewSet::ClipPlane> * clip_planes,
	MeshViewSet* view_set, int part_id) const
{
	if (MeshViewSet::param_show_visualization == 0) return nullptr;

	int bct = getBlocksCount();
	if (bct == 0) return view_set;

	if (view_set) {
		view_set->prepareFreePlace(0, 3 * bct, 0, bct);
	}
	else {
		view_set = new MeshViewSet(0, 3 * bct, 0, bct);
	}

	view_set->setMesh(this);

	DataVector<int> hidden_blocks(bct, 0);

	// 1. check clip
	if (clip_planes != nullptr && clip_planes->notEmpty()) {
		for (int i = 0; i < bct; i++) {
			for (int j = 0; j < clip_planes->countInt(); j++) {
				if (clip_planes->get(j).clipped(getBlockAt(i)->getMiddlePoint())) {
					hidden_blocks[i] = 1;
					break;
				}
			}
		}
	}

	// 2. check sight-line
	for (int i = 0; i < bct; i++) {
		if (hidden_blocks[i] > 0) continue;
		MeshBlock* block = getBlockAt(i);
		hidden_blocks[i] = 2;
		int fct = block->getFaceCount();
		for (int j = 0; j < fct; j++) {
			// if boundary or adjacent to clip-plane
			MeshBlock* fblock = block->getNeighbour(j);
			if (!fblock || hidden_blocks[fblock->getIndex()] == 1) {
				hidden_blocks[i] = 0; // make visible
				break;
			}
		}
	}

	// blocks
	for (int i = 0; i < bct; i++) {
		if (hidden_blocks[i] > 0) continue; // hidden
		MeshBlock* block = getBlockAt(i);
		int index = 0;
		if (part_id == -2 || part_id == block->getAreaID()) {
			view_set->addBlock(block);
		}
	}

	return view_set;
}

MeshViewSet* MeshContainer3d::getViewSetWithVisibleFacesOnly(
	const DataVector<MeshViewSet::ClipPlane> * clip_planes,
	MeshViewSet* view_set, int id, int show_part_id) const
{
	if (MeshViewSet::param_show_visualization == 0) return nullptr;

	int bct = getBlocksCount();
	if (bct == 0) return view_set;

	if (view_set) {
		view_set->prepareFreePlace(0, 3 * bct, 3*bct, 0);
	}
	else {
		view_set = new MeshViewSet(0, 3 * bct, 3*bct, 0);
	}

	view_set->setMesh(this);

	DataVector<bool> hidden_blocks(bct, false);

	// 1. check clip
	if (clip_planes != nullptr && clip_planes->notEmpty()) {
		for (int i = 0; i < bct; i++) {
			for (int j = 0; j < clip_planes->countInt(); j++) {
				if (clip_planes->get(j).clipped(getBlockAt(i)->getMiddlePoint())) {
					hidden_blocks[i] = true;
					break;
				}
			}
		}
	}

	// visible faces
	for (int i = 0; i < bct; i++) {
		if (hidden_blocks[i]) continue; // hidden
		MeshBlock* block = getBlockAt(i);
		int fct = block->getFaceCount();
		for (int j = 0; j < fct; j++) {
			MeshBlock* fblock = block->getNeighbour(j);
			if (!fblock || hidden_blocks[fblock->getIndex()]) {
				if (show_part_id == -2 || show_part_id == block->getAreaID()) {
					MeshFace* face = block->getFace(j);
					view_set->addFaceWithEdges(face, id, 1.0, 
						//face->outsideOrientation(block)
						face->getBlockIndex(block)==1
					);
				}
			}
		}
	}

	return view_set;
}

MeshViewSet* MeshContainer3d::getViewSet(MeshViewSet* view_set, int part_id) const
{
	if(MeshViewSet::param_show_visualization == 0) return nullptr;

	int bct = getBlocksCount();
	int pct = getPointsCount();
	if(bct + pct == 0) return view_set;

	if(view_set){
		view_set->prepareFreePlace(pct, 3*bct, 0, bct);
	}else{
		view_set = new MeshViewSet(pct, 3*bct, 0, bct);
	}

	view_set->setMesh(this);

	DataVector<bool> active_points(pct, false);

	// blocks
	for(int i = 0; i < bct; i++){
		MeshBlock* block = getBlockAt(i);
		int index = 0;
		if(part_id == -2 || part_id == block->getAreaID()){
			view_set->addBlock(block);
			int bpct = block->getPointCount();
			for(int j = 0; j < bpct; j++){
				active_points[block->getPoint(j)->getIndex()] = true;
			}
		}
	}

	// others
	for(int i = 0; i < pct; i++){
		if(!active_points[i]) continue;
		MeshPoint3d* point = getPointAt(i);
		view_set->addPoint(point);
		int rank = point->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge3d* edge = point->getEdge(j);
			if(edge->getPointIndex(point) == 0){
				int edge_rank = edge->getFaceCount();
				DataVector<MeshBlock*> adjacent_blocks(edge_rank*2);
				for(int k = 0; k < edge_rank; k++){
					MeshFace* face = edge->getFaceAt(k);
					MeshBlock* block = face->getBlock(0);
					if(block && (part_id == -2 || part_id == block->getAreaID())) 
						adjacent_blocks.addIfNew(block);
					block = face->getBlock(1);
					if(block && (part_id == -2 || part_id == block->getAreaID())) 
						adjacent_blocks.addIfNew(block);
					if(face->getEdgeIndex(edge) == 0 && !face->isBounded() && part_id == -2){
						view_set->addFace(face);
					}
				}
				if(!adjacent_blocks.empty()){
					view_set->addEdge(edge);
				}
			}
		}
	}

	return view_set;
}

MeshViewSet* MeshContainer3d::getBoundaryEdgesViewSet(MeshViewSet* view_set) const
{
	if(MeshViewSet::param_show_visualization == 0) return nullptr;

	int pct = getPointsCount();
	if(pct == 0) return view_set;

	if(view_set){
		view_set->prepareFreePlace(pct, 15*pct, 3*pct, 0);
	}else{
		view_set = new MeshViewSet(pct, 15*pct, 3*pct, 0);
	}

	view_set->setMesh(this);

	Metric3dContext mc(m_control);

	DataStatistics data;

	// edges (+ points)
	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = getPointAt(i);
		if(!point->isBorder()) continue;
		view_set->addPoint(point);
		int rank = point->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge3d* edge = point->getEdge(j);
			if(!edge->isBorder()) continue;
			if(edge->getPointIndex(point) == 0){
				double len = edge->getLengthQuality(mc, false);
				data.add(len);
				const DPoint3d pt0 = edge->getPoint(0.0);
				const DPoint3d pt1 = edge->getPoint(1.0);
				if(len <= 0.7) view_set->addEdge(pt0, pt1, 1);
				else if(len >= 1.5) view_set->addEdge(pt0, pt1, 2);
				else view_set->addEdge(pt0, pt1, 0);
			}
		}
	}

	if(data.calculate()){
		LOG4CPLUS_DEBUG(MeshLog::logger_console, 
			"Average length of boundary edges: " << data.average());
		LOG4CPLUS_DEBUG(MeshLog::logger_console, "Std dev ========================: " << data.stdDev());
		LOG4CPLUS_DEBUG(MeshLog::logger_console, "Minimum ========================: " << data.minimum());
		LOG4CPLUS_DEBUG(MeshLog::logger_console, "Maximum ========================: " << data.maximum());
	}

	return view_set;
}

void MeshContainer3d::countQuality(Metric3dContext& mc, int criterion){
	int bct = getBlocksCount();
	for(int i = 0; i < bct; i++)
		getBlockAt(i)->countQuality(mc, criterion);
}

void MeshContainer3d::countMetricDifferenceQuality(){
	int bct = getBlocksCount();
	Metric3dContext mc(m_control);
	for(int i = 0; i < bct; i++)
		getBlockAt(i)->countMetricDiffQuality(mc);
}


void MeshContainer3d::setControlSpace(CS3dPtr space)
{
	m_control = space;
}

int MeshContainer3d::getMetricEdgeRange(Metric3dContext& mc, int ranges_count, double ranges[], int stats[])
{
	for(int i = 0; i <= ranges_count; i++) stats[i] = 0;
	CS3dPtr space = getControlSpace();
	if(!space) return 0;

	int pct = getPointsCount();
	int count = 0;
	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = getPointAt(i);
		int rank = point->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge3d* edge = point->getEdge(j);
			if(edge->getPointIndex(point) == 0){
				++count;
				mc.countMetricAtPoint(edge->getPoint(0.5));
				double len = edge->getMeshPoint(0)->getMetricCoordinates(mc).distance(
					edge->getMeshPoint(1)->getMetricCoordinates(mc));
				int k;
				for(k = 0; (k < ranges_count) && (len < ranges[k]); k++);
				++stats[k];
			}
		}
	}
	return count;
}

int MeshContainer3d::getTetrahedraQualityRange(Metric3dContext& mc, int ranges_count, double ranges[], int stats[], 
	MeshData::StatData stat_data[], int count_data[])
{
	// 0 - inner blocks
	// 1 - blocks incident to boundary via vertex only
	// 2 - blocks incident to boundary via edge only
	// 3 - blocks incident to boundary via face
	// 4 - all blocks
	const int BLOCKS_TYPES = 4;
	for(int i = 0; i <= BLOCKS_TYPES ; i++){
		stat_data[i].minimum = 1.0;
		stat_data[i].maximum = 0.0;
		stat_data[i].average = 1.0;
		if(i < BLOCKS_TYPES) count_data[i] = 0;
	}

	for(int i = 0; i <= ranges_count; i++) stats[i] = 0;
	int bct = getBlocksCount();
	if(bct <= 0) return 0;

	for(int i = 0; i < bct; i++){
		MeshBlock* block = getBlockAt(i);
		if(block->hasBoundaryFace())
			++count_data[3];
		else if(block->hasBoundaryEdge())
			++count_data[2];
		else if(block->hasBoundaryVertex())
			++count_data[1];
		else
			++count_data[0];
	}
	double n_power = 1.0 / (double)bct;
	double x = 1.0;
	int count = 0;
	double nt_power[BLOCKS_TYPES];
	double xt[BLOCKS_TYPES];
//	int countt[BLOCKS_TYPES];
	for(int i = 0; i < BLOCKS_TYPES; i++){
		nt_power[i] = 1.0 / (double)count_data[i];
		xt[i] = 1.0;
//		countt[i] = 0;
	}

	for(int i = 0; i < bct; i++){
		MeshBlock* block = getBlockAt(i);
		if(block->getType() != BLOCK_TETRA) continue;
		++count;
		double quality = block->countQuality(mc, MeshData::QUALITY3D_MEAN_RATIO);

		if(quality < stat_data[BLOCKS_TYPES].minimum) stat_data[BLOCKS_TYPES].minimum = quality;
		if(quality > stat_data[BLOCKS_TYPES].maximum) stat_data[BLOCKS_TYPES].maximum = quality;
		if(stat_data[BLOCKS_TYPES].minimum > 0.0){
			x *= quality;
			if(x < 0.1){
				stat_data[BLOCKS_TYPES].average *= pow(x, n_power);
				x = 1.0;
			}
		}
		int block_type = 0;
		if(block->hasBoundaryFace())
			block_type = 3;
		else if(block->hasBoundaryEdge())
			block_type = 2;
		else if(block->hasBoundaryVertex())
			block_type = 1;
		if(quality < stat_data[block_type].minimum) stat_data[block_type].minimum = quality;
		if(quality > stat_data[block_type].maximum) stat_data[block_type].maximum = quality;
		if(stat_data[block_type].minimum > 0.0){
			xt[block_type] *= quality;
			if(xt[block_type] < 0.1){
				stat_data[block_type].average *= pow(xt[block_type], nt_power[block_type]);
				xt[block_type] = 1.0;
			}
		}
		int k;
		for(k = 0; (k < ranges_count) && (quality < ranges[k]); k++);
		++stats[k];
	}

	if(stat_data[BLOCKS_TYPES].minimum > 0.0)
		stat_data[BLOCKS_TYPES].average *= pow(x, n_power);
	else
		stat_data[BLOCKS_TYPES].average = 0.0;

	for(int i = 0; i < BLOCKS_TYPES; i++){
		if(stat_data[i].minimum > 0.0)
			stat_data[i].average *= pow(xt[i], nt_power[i]);
		else
			stat_data[i].average = 0.0;
	}

	return count;
}

MeshData::StatData MeshContainer3d::getTetrahedraQualityStats(Metric3dContext& mc) const
{
	MeshData::StatData stat_data;

	int bct = getBlocksCount();
	if(bct <= 0) return stat_data;

	stat_data.minimum = 1.0;
	stat_data.maximum = 0.0;
	stat_data.average = 1.0;

	double n_power = 1.0 / (double)bct;
	double x = 1.0;
	int count = 0;

	for(int i = 0; i < bct; i++){
		MeshBlock* block = getBlockAt(i);
		if(block->getType() != BLOCK_TETRA) continue;
		++count;
		double quality = block->countQuality(mc, MeshData::QUALITY3D_MEAN_RATIO);

		if(quality < stat_data.minimum) stat_data.minimum = quality;
		if(quality > stat_data.maximum) stat_data.maximum = quality;
		if(stat_data.minimum > 0.0){
			x *= quality;
			if(x < 0.1){
				stat_data.average *= pow(x, n_power);
				x = 1.0;
			}
		}
	}

	if(stat_data.minimum > 0.0)
		stat_data.average *= pow(x, n_power);
	else
		stat_data.average = 0.0;

	return stat_data;
}

int MeshContainer3d::getBlocksCountForFaces(int face_count) const
{
	int count = 0;
	int bcount = m_blocks->countInt();
	for(int i = 0; i < bcount; i++)
		if(getBlockAt(i)->getFaceCount() == face_count) ++count;
	return count;
}

/// Stores 3d mesh into tetmesh file
bool MeshContainer3d::storeTetmesh(const string& fname, int id) const
{
	int pct = getPointsCount();
	int bct = getBlocksCount();
	if (pct < 1 || bct < 1) {
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Empty mesh - nothing to write!");
		return false;
	}
	//	int tct = getBlocksCountForFaces(4);

	string tet_fname = fname;
	if (id >= 0) tet_fname += "-" + to_string(id);
	tet_fname += ".tetmesh";

	ofstream ftet(tet_fname.c_str());

	ftet << "Vertices " << pct << endl << fixed;
	ftet.precision(6);
	for (int i = 0; i < pct; i++) {
		MeshPoint3d* mpt = getPointAt(i);
		DPoint3d pt = mpt->getCoordinates();
		ftet << pt.x << "\t" << pt.y << "\t" << pt.z << "\t0" << endl;
	}

	const int rct[] = { 1, 0, 2, 3 }; // for switching orientation of tetrahedra (sequence of vertices)

	ftet << "Tetrahedra " << bct << endl;
	for (int i = 0; i < bct; i++) {
		MeshBlock* block = getBlockAt(i);
		int points_ct = block->getPointCount();
		assert(points_ct == 4);

		for (int j = 0; j < points_ct; j++) {
			ftet << "\t" << block->getPoint(rct[j])->getIndex();
		}
		ftet << endl;
	}

	return true;
}

/// Stores 3d mesh into tetmesh file
bool MeshContainer3d::storeFacesOFF(const string& fname, int id) const
{
	int pct = getPointsCount();
	int fct = 0;
	for (auto it = getFirstFace(); it.isValid(); it.nextFace()) { fct++; }
	if (pct < 1 || fct < 1) {
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Empty mesh - nothing to write!");
		return false;
	}

	string off_fname = fname;
	if (id >= 0) off_fname += "-" + to_string(id);
	off_fname += ".off";

	ofstream foff(off_fname.c_str());

	foff << "OFF\n" << pct << " " << fct << " 0" << endl << fixed;
	foff.precision(6);
	for (int i = 0; i < pct; i++) {
		MeshPoint3d* mpt = getPointAt(i);
		DPoint3d pt = mpt->getCoordinates();
		pt.storeSimple(foff);
	}

	for (auto it = getFirstFace(); it.isValid(); it.nextFace()) {
		auto face = it.getFace();
		int fpct = face->getPointCount();
		foff << fpct;
		for (int i = 0; i < fpct; i++)
			foff << ' ' << face->getPoint(i)->getIndex();
		foff << endl;
	}

	return foff.good();
}

/// Stores 3d mesh into text file
bool MeshContainer3d::storeTxt(const string& fname, int id) const
{
	int pct = getPointsCount();
	int bct = getBlocksCount();
	if(pct < 1 || bct < 1){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Empty mesh - nothing to write!");
		return false;
	}
//	int tct = getBlocksCountForFaces(4);

	ostringstream fname_points, fname_blocks;
	fname_points << fname;
	if(id > 0) fname_points << "-" << id << "-p.txt";
	fname_blocks << fname;
	if(id > 0) fname_blocks << "-" << id << "-b.txt";

	ofstream fp(fname_points.str().c_str());
	ofstream fb(fname_blocks.str().c_str());

	fp << pct << endl << fixed;
	fp.precision(6);
	for(int i = 0; i < pct; i++){
		MeshPoint3d* mpt = getPointAt(i);
		DPoint3d pt = mpt->getCoordinates();
		fp << i << "\t" << pt.x << "\t" << pt.y << "\t" << pt.z;
		if(mpt->availableTag(TagExtended::TAG_BOUNDARY_POINT))
			mpt = (MeshPoint3d*) mpt->getPtrTag(TagExtended::TAG_BOUNDARY_POINT);
		MeshBoundaryCondition* mbc = (MeshBoundaryCondition*) mpt->getPtrTag(TagExtended::TAG_BOUNDARY_COND);
		if(mbc)
			fp << "\t" << mbc->getCondition();
		fp << endl;
	}

	const int rct[] = {1, 0, 2, 3}; // for switching orientation of tetrahedra (sequence of vertices)

	fb << bct << endl;
	for(int i = 0; i < bct; i++){
		MeshBlock* block = getBlockAt(i);
		int points_ct = block->getPointCount();
		assert(points_ct == 4);

		fb << i << "\t" << points_ct;
		for(int j = 0; j < points_ct; j++){
			fb << "\t" << block->getPoint(rct[j])->getIndex();
		}
		fb << endl;
	}

	return true;
}

bool MeshContainer3d::storePJM(const string& fname, int grid_type, int id) const
{
	int pct = getPointsCount();
	int bct = getBlocksCount();
	if(pct < 1 || bct < 1){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Empty mesh - nothing to write!");
		return false;
	}
	int tct = getBlocksCountForFaces(4);
/*	if(tct > 0 && qct > 0){
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Mixed mesh (tri and quad) - don't know how to write!");
		return false;
	}
*/
//	MeshBlock* first_block = getBlockAt(0);
	int edge_inner_nodes_count = 0; // first_block->getEdge(0)->getInnerPointsCount(); // TODO for inner nodes

	if(grid_type == 25){
		LOG4CPLUS_INFO(MeshLog::logger_console, "Storing grid GR_3D_T4");
		if(tct != bct){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Mesh contains blocks other than tetrahedra! Leaving."); return false; }
		if(edge_inner_nodes_count > 0){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Non-zero number of inner nodes for edges! Leaving."); return false; }
	}else if(grid_type == 26){
		LOG4CPLUS_INFO(MeshLog::logger_console, "Storing grid GR_3D_T10");
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Not implemented yet - no format for grd data given...");
		return false;
//		if(tct != bct){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Mesh contains blocks other than tetrahedra! Leaving."); return false; }
//		if(edge_inner_nodes_count != 1){ LOG4CPLUS_ERROR(MeshLog::logger_console,   "Number of inner nodes for edges different from 1! Leaving."); return false; }
	}else{
		LOG4CPLUS_ERROR(MeshLog::logger_console,   "Unknown grid type! Leaving."); 
		return false;
	}

	ostringstream fname_desc, fname_data;
	fname_data << fname;
	if(id > 0) fname_data << '-' << id;
	fname_desc << fname;
	if(id > 0) fname_desc << '-' << id;
	fname_desc << ".desc";

	ofstream fdata(fname_data.str().c_str());
	ofstream fdesc(fname_desc.str().c_str());

	fdata << "GRID_TYPE " << grid_type << endl;
	DBox box = getBoundingBox();
	fdata << "DIMENSIONS\t" << box.getDX() << '\t' << box.getDY() << '\t' << box.getDZ() << endl;
//	fdata << "GRID_LEVELS\t0\t2" << endl << endl;

	fdata << "POINTS " << pct << endl << endl;
	for(int i = 0; i < pct; i++){
		DPoint3d pt = getPointAt(i)->getCoordinates();
		fdata << (abs(pt.x) < mesh_data.relative_small_number ? 0.0 : pt.x) << "\t" 
			  << (abs(pt.y) < mesh_data.relative_small_number ? 0.0 : pt.y) << "\t" 
			  << (abs(pt.z) < mesh_data.relative_small_number ? 0.0 : pt.z) <<endl;
	}

	const int rct[] = {1, 0, 2, 3}; // for switching orientation of tetrahedra (sequence of vertices)

	fdata << endl << endl << "ELEMENTS \t" << bct << '\t' << pct << endl << endl;
	for(int i = 0; i < bct; i++){
		MeshBlock* block = getBlockAt(i);
		int points_ct = block->getPointCount();
		assert(points_ct == 4);

		for(int j = 0; j < points_ct; j++){
			fdata << block->getPoint(rct[j])->getIndex() << '\t';
/*
			MeshEdge2d* edge = element->getEdge(j);
			int ipct = edge->getInnerPointsCount();
			if(ipct != edge_inner_nodes_count)
				LOG4CPLUS_WARN(MeshLog::logger_console, "Different number of inner nodes for some edges!!!");
			for(int k = 0; k < ipct; k++)
				fdata << edge->getInnerMeshPoint(k)->getIndex() << '\t';
*/
		}
		int faces_ct = block->getFaceCount();
		assert(faces_ct == 4); // for now its only possible number
//		int pjm_hash_table[4] = {3, 2, 1, 0}; // different sequence of faces for PJM format
		int pjm_hash_table[4] = {0, 1, 2, 3}; // different sequence of faces for PJM format // <- 8.11.2008
		for(int j = 0; j < faces_ct; j++){
			MeshFace* face = block->getFace(pjm_hash_table[rct[j]]);
			MeshBlock* other_block = face->getOtherBlock(block);
			if(other_block){
				fdata << '\t' << other_block->getIndex();
			}else{
				MeshBoundaryCondition* condition = 
					(MeshBoundaryCondition*)face->getPtrTag(TagExtended::TAG_BOUNDARY_COND);
				if(condition) fdata << "\t" << condition->getCondition();
				else fdata << "\t?";
			}
		}
		fdata << "\t" << block->getAreaID() << endl;
	}

	fdata << endl << "COLOR 1" << endl << bct << '\t';
	for(int i = 0; i < bct; i++) fdata << i << ' ';
	fdata << endl;

	fdesc << pct << " nodes" << endl;
	fdesc << tct << " tetrahedra" << endl;
	fdesc << (bct-tct) << " other blocks" << endl << endl;
	ControlSpace3dIdentity csi;
	Metric3dContext mc(&csi);
	MeshData::StatData quality = getTetrahedraQualityStats(mc);
	fdesc << "Quality (mean ratio with metric - ave): " << quality.average << endl;
	fdesc << "Quality (mean ratio with metric - min): " << quality.minimum << endl;
//	MeshData::StatData angles = getAngleQuality(surface);
//	fdesc << "Inner angle (min)  : " << angles.minimum*180.0/PI << "o" << endl;

	return true;
}

/*
//////////////////////////////////////////////////////////////////////
// Wstawia okreœlon¹ liczbê punktów poœrednich do wszystkich krawêdzi
//	elementów wystêpuj¹cych w siatce
void MeshContainer3d::addEdgeInnerPoints(int count, bool use_shapes)
{
	// W trakcie dodawania punktów wewnêtrznych ta wartoœæ bêdzie
	//	zwiêkszana, ale przegl¹danie nowych punktów by³oby bezu¿yteczne
	int pct = m_points->countInt();	
									
	int id = pct;
	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = m_points->getDataAt(i);
		int ect = point->getRank();
		if(ect < 2) continue;
		for(int k = 0; k < ect; k++){	// For all incident edges
			MeshEdge3d* edge3d = point->getEdge(k);
			if(edge3d->getPointIndex(point) == 0){	// Each edge only once
				MeshPoint3d **points, **old_points;
				if(count > 0){
					points = new MeshPoint3d*[count];	
				}else{
					points = nullptr;
				}
				// Tworzenie tablicy nowych punktów
				for(int n = 0; n < count; n++){
					DPoint3d pt = edge3d->getAspectPoint(((n+1) / (double)(count+1)));
					m_points->addDataItem(points[n] = new MeshPoint3d(pt));
					points[n]->addEdgeLink(edge3d);
					points[n]->setAsInnerNode();
				}
				// Zamiana tablic
				int old_ct = edge3d->addInnerPoints(count, points, &old_points);
				// Usuwanie starych punktów (o ile istniej¹)
				for(int n = 0; n < old_ct; n++)
					delete removeMeshPoint(old_points[n]->getIndex());
				if(old_ct > 0) delete[] old_points;
				if(old_ct > count) pct -= (old_ct - count);
			}
		}
	}
//	renumeratePoints();
}
*/

/// Returns number of faces with the given number of edges
int MeshContainer3d::getFaceCount(int ect) const
{
	int bct = m_blocks->countInt();
	int count = 0;
	for(int i = 0; i < bct; i++){
		MeshBlock* block = getBlockAt(i);
		int fct = block->getFaceCount();
		for(int j = 0; j < fct; j++){
			MeshFace* face = block->getFace(j);
			if(face->getEdgeCount() == ect && 
				face->getBlockIndex(block) == 0)
				++count;
		}
	}
	return count;
}

DPoint3d MeshContainer3d::getInertialCenter() const
{
	DPoint3d pt;
	int pct = m_points->countInt();
	if(!pct) return pt;
	for(int i = 0; i < pct; i++)
	    pt.add(getPointAt(i)->getCoordinates());
	pt /= (double)pct;
	return pt;
}

int MeshContainer3d::logMetricQuality2d()
{
	int mesh_count = 0;
	for(IteratorMesh2d it = getFirstValidMesh2d(); it.isValid(); it.nextValidMesh()){
		MeshContainer2d* mesh = it.getMesh();
		++mesh_count;
		Metric2dContext mc(mesh->getControlSpace());
		DataStatistics stats;

		LOG4CPLUS_INFO(MeshLog::logger_mesh, "===> Patch " << mesh_count);

		if(MeshGenerator2d::statMetricQuality(mc, mesh, stats, MeshData::QUALITY_SPACE))
			stats.logStats("MetricQuality - Triangle Area", "MQ-TS");
		if(MeshGenerator2d::statMetricQuality(mc, mesh, stats, MeshData::QUALITY_CIRCLE_SPACE))
			stats.logStats("MetricQuality - Circumcircle Area", "MQ-CS");
		if(mesh->statMetricEdgeLength(mc, stats, false))
			stats.logStats("MetricQuality - Edge Length", "MQ-EL");
	}

	return mesh_count;
}

/// Returns statistical information about minimum inner angles (dihedral)
bool MeshContainer3d::statMinDihedralAngles(DataStatistics& stats) const
{
	int bct = getBlocksCount();
	for(int i = 0; i < bct; i++){
		const MeshTetrahedron* tetra = (MeshTetrahedron*)getBlockAt(i);
		if(tetra->getType() != BLOCK_TETRA) continue; // tetrahedra only
			//  - 3D shape quality in metric
		double qt = asin(tetra->getMinDihedralAngleSin());
		stats.add(qt);
	}
	return bct > 0;
}

bool MeshContainer3d::statMeanRatio(Metric3dContext& mc, DataStatistics& stats, bool ext_metric) const
{
	int ect = getBlocksCount();
	for(int i = 0; i < ect; i++){
		const MeshTetrahedron* tetra = (MeshTetrahedron*)getBlockAt(i);
		if(tetra->getType() != BLOCK_TETRA) continue; // tetrahedra only
			//  - 3D shape quality in metric
		double qt = tetra->getMeanRatio(mc, ext_metric);
		stats.add(qt);
	}
	return true;
}

#if defined(USE_OPENMP_HERE) && defined(_OPENMP)

bool MeshContainer3d::statMetricEdgeLength(Metric3dContext& mc, DataStatistics& stats, bool ext_metric) const
{
	DataVector<MeshEdge3d*> edge_list(getPointsCount());
	for (IteratorEdge3d it = getFirstEdge3d(); it.isValid(); it.nextEdge())
		edge_list.add(it.getEdge());
	int ect = edge_list.countInt();

	Metric3dContext mc_local(mc);

	#pragma omp parallel for firstprivate(mc_local) shared(stats,ext_metric)
	for (int i = 0; i < ect; i++) {
		double q = edge_list[i]->getLengthQuality(mc_local, ext_metric);
		#pragma omp critical(estats_update)
		{
			stats.add(q);
		}
	}

	return true;
}

#else

bool MeshContainer3d::statMetricEdgeLength(Metric3dContext& mc, DataStatistics& stats, bool ext_metric) const
{
	for(IteratorEdge3d it = getFirstEdge3d(); it.isValid(); it.nextEdge())
		stats.add(it.getEdge()->getLengthQuality(mc, ext_metric));
	return true;
}

#endif // USE_OPENMP_HERE && _OPENMP

bool MeshContainer3d::statMetricEdgeLengthNoOMP(Metric3dContext& mc, DataStatistics& stats, bool ext_metric) const
{
	for (IteratorEdge3d it = getFirstEdge3d(); it.isValid(); it.nextEdge())
		stats.add(it.getEdge()->getLengthQuality(mc, ext_metric));
	return true;
}

/// Returns statistical information about number of blocks adjacent to each edge
bool MeshContainer3d::statEdgeBlockAdjacency(DataStatistics& stats, bool inner_only) const
{
	for(IteratorEdge3d it = getFirstEdge3d(); it.isValid(); it.nextEdge()){
		MeshEdge3d* edge = it.getEdge();
		if(inner_only && edge->isBorder()) continue;
		int fct = edge->getFaceCount();
		DataHashTable<MeshBlock*> unique_blocks(2*fct, nullptr);
		for(int i = 0; i < fct; i++){
			MeshFace* face = edge->getFaceAt(i);
			if(face->getBlock(0)) unique_blocks.insert(face->getBlock(0));
			if(face->getBlock(1)) unique_blocks.insert(face->getBlock(1));
		}
		stats.add(unique_blocks.countInt());
	}
	return stats.countInt() > 0;
}

/// Returns statistical information about number of blocks adjacent to each vertex
bool MeshContainer3d::statVertexBlockAdjacency(DataStatistics& stats, bool inner_only) const
{
	int pct = getPointsCount();
	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = getPointAt(i);
		if(inner_only && point->isBorder()) continue;
		DataVector<MeshBlock*> blocks(10*point->getRank());
		point->adjacentBlocks(blocks);
		stats.add((int)blocks.countInt());
	}
	return stats.countInt() > 0;
}

void MeshContainer3d::removeAllTags(TagExtended::TagType tag_type)
{
	int pct = getPointsCount();
	#pragma omp parallel for
	for(int i = 0; i < pct ; i++){ // clear all tags for points, edges and faces
		MeshPoint3d* point = getPointAt(i);
		point->removeTag(tag_type);
		int rank = point->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge3d* edge = point->getEdge(j);
			if(edge->getPointIndex(point) == 0){
				edge->removeTag(tag_type);
				int fcount = edge->getFaceCount();
				for(int k = 0; k < fcount; k++)
					edge->getFaceAt(k)->removeTag(tag_type);
			}
		}
	}
}

void MeshContainer3d::removeAllTags()
{
	int pct = getPointsCount();
	#pragma omp parallel for
	for(int i = 0; i < pct ; i++){ // clear all tags for points, edges, faces and blocks
		MeshPoint3d* point = getPointAt(i);
		point->removeAllTags();
		int rank = point->getRank();
		for(int j = 0; j < rank; j++){
			MeshEdge3d* edge = point->getEdge(j);
			if(edge->getPointIndex(point) == 0){
				edge->removeAllTags();
				int fcount = edge->getFaceCount();
				for(int k = 0; k < fcount; k++)
					edge->getFaceAt(k)->removeAllTags();
			}
		}
	}
	int bct = getBlocksCount();
	#pragma omp parallel for
	for(int i = 0; i < bct; i++)
		getBlockAt(i)->removeAllTags();
}

bool MeshContainer3d::storeMatlabFile(const string& fname, const DPoint3d& clip, double quality_clip) const
{
	ofstream os(fname.c_str());

	os << "clear all; figure(1); clf; hold on; axis equal; axis off; view(3);" << endl;
	os << "c(1,:) = [0.9,0.9,0.9];" << endl;
	os << "c(2,:) = [0.7,0.7,0.7];" << endl;
	os << "colormap(c);" << endl;
	os << "load 'mesh_vertices.txt' -ascii" << endl;
	os << "load 'mesh_faces.txt' -ascii" << endl;
	os << "load 'mesh_ci.txt' -ascii" << endl;
	os << "patch('Vertices',mesh_vertices,'faces',mesh_faces,'FaceVertexCData',mesh_ci(:,1),'FaceColor','flat');" << endl;

	// store vertices
	ofstream ofsv("mesh_vertices.txt");
	int pct = getPointsCount();
	for(int j = 0; j < pct; j++){
		const DPoint3d& pt = getPointAt(j)->getCoordinates();
		ofsv << pt.x << ' ' << pt.y << ' ' << pt.z << endl;
	}

	// check blocks
	int bct = getBlocksCount();
	DataVector<bool> hidden(bct, true);
	for(int i = 0; i < bct; i++){
		MeshBlock* block = getBlockAt(i);
		if(quality_clip < 1.05 && block->getQuality() > quality_clip) continue;
		int bpct = block->getPointCount();
		for(int j = 0; j < bpct; j++){
			const DPoint3d& bpt = block->getPoint(j)->getCoordinates();
			if(bpt.x > clip.x || bpt.y > clip.y || bpt.z > clip.z){
				hidden[i] = true; break;
			}
		}
	}

	ofstream ofsf("mesh_faces.txt");
	ofstream ofsc("mesh_ci.txt");

	// store exterior faces of visible (not-clipped) blocks
	for(int i = 0; i < bct; i++){
		if(hidden[i]) continue;
		MeshBlock* block = getBlockAt(i);
		// store faces
		int bfct = block->getFaceCount();
		for(int j = 0; j < bfct; j++){
			const MeshFace* face = block->getFace(j);
			const MeshBlock* bl0 = face->getBlock(0);
			const MeshBlock* bl1 = face->getBlock(1);
			if(bl0 && bl1 && hidden[bl0->getIndex()] && hidden[bl1->getIndex()]) 
				continue; // inner face
			int fpct = face->getEdgeCount();
			for(int k = 0; k < fpct; k++)
				ofsf << (1+ face->getPoint(k)->getIndex()) << ' ';
			ofsf << endl;
			ofsc << (face->isBorder() ? 2 : 1) << endl;
		}
	}

	return true;
}

bool MeshContainer3d::statMetricDifference(Metric3dContext& mc, DataStatistics& stats) const
{
	int bct = getBlocksCount();

	Metric3dContext mc_local(mc);

	#pragma omp parallel for firstprivate(mc_local) shared(stats)
	for(int i = 0; i < bct; i++){
		const MeshTetrahedron* tetra = (MeshTetrahedron*)getBlockAt(i);
		if(tetra->getType() != BLOCK_TETRA) continue; // tetrahedra only
			//  - 2D metric from simplex (three points)
		double diff = tetra->countMetricDiff(mc_local);

		#pragma omp critical(mstats_update)
		{
			stats.add(diff);
		}
	}

	return true;
}

bool MeshContainer3d::statMetricDifferenceNoOMP(Metric3dContext& mc, DataStatistics& stats) const
{
	int bct = getBlocksCount();

	for (int i = 0; i < bct; i++) {
		const MeshTetrahedron* tetra = (MeshTetrahedron*)getBlockAt(i);
		if (tetra->getType() != BLOCK_TETRA) continue; // tetrahedra only
													   //  - 2D metric from simplex (three points)
		double diff = tetra->countMetricDiff(mc);
		stats.add(diff);
	}

	return true;
}

void MeshContainer3d::logQuality(bool with_openmp) const
{
	if(with_openmp)
		START_CLOCK("MC3d::logQualityOMP");
	else
		START_CLOCK("MC3d::logQuality");

	DataStatistics mean_ratio_stats;
	DataStatistics edge_stats;
	DataStatistics nonconf_stats;

	unsigned int np = getPointsCount();
	unsigned int nt = getBlocksCount();

//	LOG4CPLUS_INFO(MeshLog::logger_mesh, "logQuality: np=" << np << ", nt=" << nt);

	assert(m_control);
	Metric3dContext mc(m_control);
	//statMeanRatio(mc, mean_ratio_stats, false);

	if (with_openmp) {
		statMetricEdgeLength(mc, edge_stats, false);
		statMetricDifference(mc, nonconf_stats);
	}
	else {
		statMetricEdgeLengthNoOMP(mc, edge_stats, false);
		statMetricDifferenceNoOMP(mc, nonconf_stats);
	}

	auto logger_cstat = Logger::getInstance(LOG4CPLUS_TEXT("mesh.stat.cstat"));

	LOG4CPLUS_INFO(logger_cstat, "NP\t" << np);
	LOG4CPLUS_INFO(logger_cstat, "NT\t" << nt);
	//if(mean_ratio_stats.calculate()){
	//	LOG4CPLUS_INFO(MeshLog::logger_console, "Metric mean ratio (ave)", mean_ratio_stats.average());
	//	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "MetricMeanRatio-ave\t" << mean_ratio_stats.average());
	//	LOG4CPLUS_INFO(MeshLog::logger_console, "Metric mean ratio (min)", mean_ratio_stats.minimum());
	//	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "MetricMeanRatio-min\t" << mean_ratio_stats.minimum());
	//	LOG4CPLUS_INFO(MeshLog::logger_console, "------------------(dev)", mean_ratio_stats.stdDev());
	//	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, "MetricMeanRatio-dev\t" << mean_ratio_stats.stdDev());
	//}
	if(edge_stats.calculate()){
		LOG4CPLUS_INFO(logger_cstat, "MetricEdgeLength-ave\t" << edge_stats.average());
		LOG4CPLUS_INFO(logger_cstat, "MetricEdgeLength-dev\t" << edge_stats.stdDev());
		int low_ct = edge_stats.getDataCountBottom(0.8);
		int hi_ct = edge_stats.getDataCountTop(1.25);
		int all_ct = edge_stats.countInt();
		int good_ct = all_ct - low_ct - hi_ct;
		double e_good_ratio = ((double)good_ct / all_ct);
		LOG4CPLUS_INFO(logger_cstat, "MetricEdgeLength-good-ratio-[0.8, 1.25]\t" << e_good_ratio);
	}
	if(nonconf_stats.calculate()){
		LOG4CPLUS_INFO(logger_cstat, "MetricDiff-ave\t" << nonconf_stats.average());
		LOG4CPLUS_INFO(logger_cstat, "MetricDiff-dev\t" << nonconf_stats.stdDev());
		int good_ct = nonconf_stats.getDataCountBottom(2.0);
		int all_ct = nonconf_stats.countInt();
		double m_good_ratio = ((double)good_ct / all_ct);
		LOG4CPLUS_INFO(logger_cstat, "MetricDiff-good-ratio-[1.0, 2.0]\t" << m_good_ratio);
		}

	if (with_openmp)
		STOP_CLOCK("MC3d::logQualityOMP");
	else
		STOP_CLOCK("MC3d::logQuality");
}

/// Set freepoints list
void MeshContainer3d::setFreePoints(
	std::shared_ptr<DataVector<std::shared_ptr<MeshPoint3d>>> fp_list)
{
	m_freepoints = fp_list;
}

/// Set boundary-condition list
void MeshContainer3d::setBoundaryConditions(
	std::shared_ptr<DataVector<std::shared_ptr<MeshBoundaryCondition>>> mbc_list)
{
	m_mbc = mbc_list;
}

/// update boundary flags for outer faces/edges/nodes
void MeshContainer3d::markOuterBoundary()
{
	int bct = getBlocksCount();
	// ... first clear all border info
	for(int i = 0; i < bct; i++){
		MeshBlock* block = getBlockAt(i);
		int fct = block->getFaceCount();
		for(int j = 0; j < fct; j++){
			block->getFace(j)->setBorderWithEdgesAndPoints(TagBorder::NONE);
		}
	}
	// ... and then set border accordingly
	for(int i = 0; i < bct; i++){
		MeshBlock* block = getBlockAt(i);
		int bid = block->getAreaID();
		int fct = block->getFaceCount();
		for(int j = 0; j < fct; j++){
			MeshBlock* other_block = block->getNeighbour(j);
			if(other_block == nullptr || other_block->getAreaID() != bid){
				block->getFace(j)->setBorderWithEdgesAndPoints(TagBorder::OUTER);
			}
		}
	}
}

/// mark boundary edges for boundary faces with sharp angle
void MeshContainer3d::markSharpBoundaryEdges(double fmax)
{
	int pct = getPointsCount();
	if(pct < 1) return;
	DataVector<int> rpoints(pct, 0);

	// Identify ridge edges
	for(IteratorEdge3d it = getFirstEdge3d(); it.isValid(); it.nextEdge()){
		MeshEdge3d* edge = it.getEdge();
		if(!edge->isBorder()) continue;
		int fct = edge->getFaceCount();
		DataVector<MeshFace*> bfaces;
		for(int i = 0; i < fct; i++){
			MeshFace* face = edge->getFaceAt(i);
			if(face->isBorder()) bfaces.add(face);
		}
		if(bfaces.countInt() > 2){
			edge->setBorderFlags(TagBorder::RIDGE);
			rpoints[edge->getMeshPoint(0)->getIndex()]++;
			rpoints[edge->getMeshPoint(1)->getIndex()]++;
			continue;
		}else if(bfaces.countInt() < 2) continue;
		// if exactly == 2, check angles
		DVector3d dn0 = bfaces[0]->getNormalVector();
		DVector3d dn1 = bfaces[1]->getNormalVector();
		// ... orientation
		if(bfaces[0]->getBlockIndex(nullptr) != bfaces[1]->getBlockIndex(nullptr))
			dn0.reverse();
		// ... check scalar product
		if(dn0.scalarProduct(dn1) < fmax){
			edge->setBorderFlags(TagBorder::RIDGE);
			rpoints[edge->getMeshPoint(0)->getIndex()]++;
			rpoints[edge->getMeshPoint(1)->getIndex()]++;
		}
	}

	// Identify ridge/corner points
	for(int i = 0; i < pct; i++){
		if(rpoints[i] == 1 || rpoints[i] > 2) 
			getPointAt(i)->setBorderFlags(TagBorder::RIDGE | TagBorder::CORNER);
		else if(rpoints[i] == 2){
			getPointAt(i)->setBorderFlags(TagBorder::RIDGE);
		}
	}
}

/// Gather volume meshes from all domain volumes in this model and combine them in the total_mdv
bool MeshContainer3d::gatherVolumeMeshes()
{
	total_mdv.clearDiscretization();
	int bct = getBlocksCount();
	if(bct == 0) return false;
	for(int i = 0; i < bct; i++){
		MeshDomainVolume* mdv = (MeshDomainVolume*) getBlockAt(i);
		if(mdv->getType() != BLOCK_DOMAIN) continue;
		total_mdv.combineVolumeMeshFrom(mdv);
	}
	return total_mdv.getMesh() != nullptr;
}

/// Removes all local surface definitions and counters
void MeshContainer3d::clearLocalShapes()
{
	if(m_local_surfaces.countInt() > 0 || m_local_curves.countInt() > 0){
		m_local_surfaces.clear();
		m_local_curves.clear();
		// clear info in points' tags
		int pct = getPointsCount();
		for(int i = 0; i < pct; i++){
			getPointAt(i)->clearLocalShapes();
		}
		// and in faces' tags
		for(IteratorFace it = getFirstFace(); it.isValid(); it.nextFace()){
			it.getFace()->clearSurfaceData();
		}
		// and in edges' tags
		for(IteratorEdge3d it = getFirstEdge3d(); it.isValid(); it.nextEdge())
			it.getEdge()->clearLocalCurve();
	}
}

/// Tries to move points to fit theirs ascribed local surfaces
int MeshContainer3d::movePointsToLocalShape(Metric3dContext& mc, 
			TagExtended::TagType tag_type, int tag_value, int forbid_tag_value)
{
	if(m_local_surfaces.empty() && m_local_curves.empty()) return 0;

	int moved_count = 0;
	int pct = getPointsCount();
	for(int i = 0; i < pct; i++){
		MeshPoint3d* point = getPointAt(i);
		if(!point->isBorder()) continue;
		if(tag_type != TagExtended::TAG_NONE && !point->hasAnyIntFlags(tag_type, tag_value)) 
			continue;
		if(forbid_tag_value > 0 && point->hasAnyIntFlags(tag_type, forbid_tag_value)) 
			continue;
		if(point->moveToLocalShape(mc)) moved_count++;
	}
	return moved_count;
}

/// Add new local surface
void MeshContainer3d::addLocalSurface(const SurfacePtr & surface)
{
	m_local_surfaces.add(surface);
}

/// Add new local curve
void MeshContainer3d::addLocalCurve(const Curve3dPtr & curve)
{
	m_local_curves.add(curve);
}

DataVector<MeshEdge3d*> MeshContainer3d::getEdges3d(TagExtended::TagType tag_type, int tag_value)
{
	DataVector<MeshEdge3d*> edges(getBlocksCount());
	for (IteratorEdge3d it = getFirstEdge3d(tag_type, tag_value); it.isValid(); it.nextEdge())
		edges.add(it.getEdge());

	return edges;
}

DataVector<const MeshEdge3d*> MeshContainer3d::getEdges3d(TagExtended::TagType tag_type, int tag_value) const
{
	DataVector<const MeshEdge3d*> edges(getBlocksCount());
	for (IteratorEdge3d it = getFirstEdge3d(tag_type, tag_value); it.isValid(); it.nextEdge())
		edges.add(it.getEdge());

	return edges;
}

DataVector<MeshDomainSurface*>  MeshContainer3d::getValidBoundaries2d()
{
	DataVector<MeshDomainSurface*> dsurfaces(4 * getBlocksCount());

	for (IteratorBoundary2d it = getFirstValidBoundary2d(); it.isValid(); it.nextValidBoundary())
		dsurfaces.add(it.getDomainSurface());

	return dsurfaces;
}

DataVector<MeshContainer2d*> MeshContainer3d::getValidMeshes2d()
{
	DataVector<MeshContainer2d*> meshes(getBlocksCount());

	for (IteratorMesh2d it = getFirstValidMesh2d(); it.isValid(); it.nextValidMesh())
		meshes.add(it.getMesh());

	return meshes;
}
