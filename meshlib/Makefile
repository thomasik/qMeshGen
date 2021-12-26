CXX       = g++
DIR       = linux
INST_DIR  = ../linux
CXXFLAGS  = -DUNIX_MODE -DNDEBUG -std=c++14 -pedantic -I../libboard -I../log4cplus-linux/include -Wall -W -Wconversion -Wshadow -Wcast-qual -Wwrite-strings -Wno-unknown-pragmas -Wno-char-subscripts
LDFLAGS   = -O2

MESHLIB   = meshlib.a
SRCS      = ControlSpace2dAdaptive.cpp ControlSpace2dAnalytic.cpp\
	 ControlSpace2d.cpp ControlSpace2dKdTreeL.cpp\
	 ControlSpace2dKdTreeV.cpp ControlSpace2dMatrixUniform.cpp\
	 ControlSpace2dMesh.cpp ControlSpace2dQuadTree.cpp\
	 ControlSpace3dAdaptive.cpp ControlSpace3dAnalytic.cpp\
	 ControlSpace3d.cpp ControlSpace3dInternal.cpp\
	 ControlSpace3dKdTree.cpp ControlSpace3dKdTreeL.cpp\
	 ControlSpace3dKdTreeLi.cpp ControlSpace3dKdTreeV.cpp\
	 ControlSpace3dMatrixUniform.cpp ControlSpace3dOctree.cpp\
	 ControlSpace3dOctreeLB.cpp ControlSpace3dOctreeL.cpp\
	 ControlSpace3dOctreeV.cpp Curve2dAnalytic.cpp\
	 Curve2dBSpline.cpp Curve2dCircle.cpp\
	 Curve2dLinearQuadric.cpp Curve2dParametric.cpp\
	 Curve2dSegment.cpp Curve3dBSpline.cpp\
	 Curve3dLinearQuadric.cpp Curve3dParametric.cpp\
	 Curve3dSegment.cpp Curve3dSurfaceParametric.cpp\
	 DataKdSplit2dContField.cpp DataKdSplit2d.cpp\
	 DataKdSplit2dSingleSource.cpp DataKdSplit3dContField.cpp\
	 DataKdSplit3d.cpp DataKdSplit3dSingleSource.cpp\
	 DataStatistics.cpp DEquation.cpp\
	 DHesjan.cpp DLeastSquaresFitting.cpp\
	 DLinearQuadric.cpp DLine.cpp\
	 DMatrix.cpp DMetric2d.cpp\
	 DMetric3d.cpp DOrientedBox.cpp\
	 DPlanarQuadric.cpp DPlane.cpp\
	 DPoint.cpp DQuad.cpp\
	 DQuadric.cpp DRect.cpp\
	 DSegment.cpp DSphere.cpp\
	 DSurfaceFitting.cpp DTetrahedron.cpp\
	 DTriangle.cpp DVector.cpp\
	 EPSFile.cpp FrontEdge.cpp\
	 FrontFace.cpp FrontLines.cpp\
	 GeometricPredicates.cpp IteratorEdge2d.cpp\
	 IteratorEdge3d.cpp IteratorFace.cpp\
	 IteratorMesh2d.cpp IteratorMesh3d.cpp\
	 MarchingCubes.cpp MeshArea.cpp\
	 MeshBlock.cpp MeshBRep.cpp\
	 MeshBRepRSurf.cpp MeshBRepXML.cpp\
	 MeshContainer2d.cpp MeshContainer3d.cpp\
	 MeshContainer3dSurface.cpp MeshData.cpp\
	 MeshDecompositionUTC.cpp MeshDomainEdge3d.cpp\
	 MeshDomainSurface.cpp MeshDomainVolume.cpp\
	 MeshEdge2d.cpp MeshEdge2dCurve.cpp\
	 MeshEdge3d.cpp MeshElement.cpp\
	 MeshFace.cpp MeshGenerator1d.cpp\
	 MeshGenerator2d.cpp MeshGenerator2dQMorph.cpp\
	 MeshGenerator2dQuad.cpp MeshGenerator3dAdapt.cpp\
	 MeshGenerator3d.cpp MeshGenerator3dDelaunayBoundary.cpp\
	 MeshGenerator3dDirectBoundary.cpp MeshGenerator3dQuality.cpp\
	 MeshGenerator3dSurface.cpp MeshGrain.cpp\
	 MeshLog.cpp MeshModel.cpp\
	 MeshPoint2d.cpp MeshPoint3d.cpp\
	 MeshPoly3d.cpp MeshQuad2d.cpp\
	 MeshQuad3d.cpp MeshSpecialRoutinesDAT.cpp\
	 MeshSpecialRoutinesPaszynski.cpp MeshSpecialRoutinesUTC.cpp\
	 MeshSplit3d.cpp MeshTetrahedron.cpp\
	 MeshTriangle2d.cpp MeshTriangle3d.cpp\
	 MeshViewSet.cpp Metric2dContext.cpp\
	 Metric3dContext.cpp ModelMultiSphere.cpp\
	 nrutil.cpp OctTree.cpp\
	 Quad12.cpp QuadTree.cpp\
	 StatisticCounter.cpp SurfaceAnalytic.cpp\
	 SurfaceBSplineCylindrical.cpp SurfaceBSplinePlanar.cpp\
	 SurfaceCorrected.cpp SurfaceCylinder.cpp\
	 SurfaceDomain.cpp SurfaceDomainFaces.cpp\
	 SurfaceDomainFacesPlanar.cpp SurfaceDomainHull.cpp\
	 SurfaceMonoAnalytic.cpp SurfaceMulti.cpp\
	 SurfaceParametric.cpp SurfacePlanarQuadric.cpp\
	 SurfacePlane.cpp SurfaceQuadricQTree.cpp\
	 SurfaceSphere.cpp tinystr.cpp\
	 tinyxml.cpp tinyxmlerror.cpp\
	 tinyxmlparser.cpp

OBJS      =  $(DIR)/ControlSpace2dAdaptive.o $(DIR)/ControlSpace2dAnalytic.o\
	 $(DIR)/ControlSpace2d.o $(DIR)/ControlSpace2dKdTreeL.o\
	 $(DIR)/ControlSpace2dKdTreeV.o $(DIR)/ControlSpace2dMatrixUniform.o\
	 $(DIR)/ControlSpace2dMesh.o $(DIR)/ControlSpace2dQuadTree.o\
	 $(DIR)/ControlSpace3dAdaptive.o $(DIR)/ControlSpace3dAnalytic.o\
	 $(DIR)/ControlSpace3d.o $(DIR)/ControlSpace3dInternal.o\
	 $(DIR)/ControlSpace3dKdTree.o $(DIR)/ControlSpace3dKdTreeL.o\
	 $(DIR)/ControlSpace3dKdTreeLi.o $(DIR)/ControlSpace3dKdTreeV.o\
	 $(DIR)/ControlSpace3dMatrixUniform.o $(DIR)/ControlSpace3dOctree.o\
	 $(DIR)/ControlSpace3dOctreeLB.o $(DIR)/ControlSpace3dOctreeL.o\
	 $(DIR)/ControlSpace3dOctreeV.o $(DIR)/Curve2dAnalytic.o\
	 $(DIR)/Curve2dBSpline.o $(DIR)/Curve2dCircle.o\
	 $(DIR)/Curve2dLinearQuadric.o $(DIR)/Curve2dParametric.o\
	 $(DIR)/Curve2dSegment.o $(DIR)/Curve3dBSpline.o\
	 $(DIR)/Curve3dLinearQuadric.o $(DIR)/Curve3dParametric.o\
	 $(DIR)/Curve3dSegment.o $(DIR)/Curve3dSurfaceParametric.o\
	 $(DIR)/DataKdSplit2dContField.o $(DIR)/DataKdSplit2d.o\
	 $(DIR)/DataKdSplit2dSingleSource.o $(DIR)/DataKdSplit3dContField.o\
	 $(DIR)/DataKdSplit3d.o $(DIR)/DataKdSplit3dSingleSource.o\
	 $(DIR)/DataStatistics.o $(DIR)/DEquation.o\
	 $(DIR)/DHesjan.o $(DIR)/DLeastSquaresFitting.o\
	 $(DIR)/DLinearQuadric.o $(DIR)/DLine.o\
	 $(DIR)/DMatrix.o $(DIR)/DMetric2d.o\
	 $(DIR)/DMetric3d.o $(DIR)/DOrientedBox.o\
	 $(DIR)/DPlanarQuadric.o $(DIR)/DPlane.o\
	 $(DIR)/DPoint.o $(DIR)/DQuad.o\
	 $(DIR)/DQuadric.o $(DIR)/DRect.o\
	 $(DIR)/DSegment.o $(DIR)/DSphere.o\
	 $(DIR)/DSurfaceFitting.o $(DIR)/DTetrahedron.o\
	 $(DIR)/DTriangle.o $(DIR)/DVector.o\
	 $(DIR)/EPSFile.o $(DIR)/FrontEdge.o\
	 $(DIR)/FrontFace.o $(DIR)/FrontLines.o\
	 $(DIR)/GeometricPredicates.o $(DIR)/IteratorEdge2d.o\
	 $(DIR)/IteratorEdge3d.o $(DIR)/IteratorFace.o\
	 $(DIR)/IteratorMesh2d.o $(DIR)/IteratorMesh3d.o\
	 $(DIR)/MarchingCubes.o $(DIR)/MeshArea.o\
	 $(DIR)/MeshBlock.o $(DIR)/MeshBRep.o\
	 $(DIR)/MeshBRepRSurf.o $(DIR)/MeshBRepXML.o\
	 $(DIR)/MeshContainer2d.o $(DIR)/MeshContainer3d.o\
	 $(DIR)/MeshContainer3dSurface.o $(DIR)/MeshData.o\
	 $(DIR)/MeshDecompositionUTC.o $(DIR)/MeshDomainEdge3d.o\
	 $(DIR)/MeshDomainSurface.o $(DIR)/MeshDomainVolume.o\
	 $(DIR)/MeshEdge2d.o $(DIR)/MeshEdge2dCurve.o\
	 $(DIR)/MeshEdge3d.o $(DIR)/MeshElement.o\
	 $(DIR)/MeshFace.o $(DIR)/MeshGenerator1d.o\
	 $(DIR)/MeshGenerator2d.o $(DIR)/MeshGenerator2dQMorph.o\
	 $(DIR)/MeshGenerator2dQuad.o $(DIR)/MeshGenerator3dAdapt.o\
	 $(DIR)/MeshGenerator3d.o $(DIR)/MeshGenerator3dDelaunayBoundary.o\
	 $(DIR)/MeshGenerator3dDirectBoundary.o $(DIR)/MeshGenerator3dQuality.o\
	 $(DIR)/MeshGenerator3dSurface.o $(DIR)/MeshGrain.o\
	 $(DIR)/MeshLog.o $(DIR)/MeshModel.o\
	 $(DIR)/MeshPoint2d.o $(DIR)/MeshPoint3d.o\
	 $(DIR)/MeshPoly3d.o $(DIR)/MeshQuad2d.o\
	 $(DIR)/MeshQuad3d.o $(DIR)/MeshSpecialRoutinesDAT.o\
	 $(DIR)/MeshSpecialRoutinesPaszynski.o $(DIR)/MeshSpecialRoutinesUTC.o\
	 $(DIR)/MeshSplit3d.o $(DIR)/MeshTetrahedron.o\
	 $(DIR)/MeshTriangle2d.o $(DIR)/MeshTriangle3d.o\
	 $(DIR)/MeshViewSet.o $(DIR)/Metric2dContext.o\
	 $(DIR)/Metric3dContext.o $(DIR)/ModelMultiSphere.o\
	 $(DIR)/nrutil.o $(DIR)/OctTree.o\
	 $(DIR)/Quad12.o $(DIR)/QuadTree.o\
	 $(DIR)/StatisticCounter.o $(DIR)/SurfaceAnalytic.o\
	 $(DIR)/SurfaceBSplineCylindrical.o $(DIR)/SurfaceBSplinePlanar.o\
	 $(DIR)/SurfaceCorrected.o $(DIR)/SurfaceCylinder.o\
	 $(DIR)/SurfaceDomain.o $(DIR)/SurfaceDomainFaces.o\
	 $(DIR)/SurfaceDomainFacesPlanar.o $(DIR)/SurfaceDomainHull.o\
	 $(DIR)/SurfaceMonoAnalytic.o $(DIR)/SurfaceMulti.o\
	 $(DIR)/SurfaceParametric.o $(DIR)/SurfacePlanarQuadric.o\
	 $(DIR)/SurfacePlane.o $(DIR)/SurfaceQuadricQTree.o\
	 $(DIR)/SurfaceSphere.o $(DIR)/tinystr.o\
	 $(DIR)/tinyxml.o $(DIR)/tinyxmlerror.o\
	 $(DIR)/tinyxmlparser.o

default: $(INST_DIR)/$(MESHLIB)

$(INST_DIR)/$(MESHLIB): $(OBJS)
	ar rcs $(INST_DIR)/$(MESHLIB) $(OBJS) 

clean:
	rm -f $(OBJS) $(INST_DIR)/$(MESHLIB)

$(DIR)/ControlSpace2dAdaptive.o: ControlSpace2dAdaptive.cpp  ControlSpace2dAdaptive.h SurfaceParametric.h Curve2dParametric.h Curve2dSegment.h DEquation.h MeshData.h MeshContainer2d.h MeshPoint2d.h MeshEdge2d.h MeshEdge2dCurve.h EPSFile.h ControlSpace3d.h ControlSpace3dAdaptive.h Metric2dContext.h MeshViewSet.h
	${CXX} $(CXXFLAGS) -c ControlSpace2dAdaptive.cpp -o $(DIR)/ControlSpace2dAdaptive.o

$(DIR)/ControlSpace2dAnalytic.o: ControlSpace2dAnalytic.cpp  ControlSpace2dAnalytic.h SurfaceParametric.h DEquation.h MeshLog.h MeshData.h
	${CXX} $(CXXFLAGS) -c ControlSpace2dAnalytic.cpp -o $(DIR)/ControlSpace2dAnalytic.o

$(DIR)/ControlSpace2d.o: ControlSpace2d.cpp  ControlSpace2d.h SurfaceParametric.h DEquation.h ControlSpace2dAdaptive.h
	${CXX} $(CXXFLAGS) -c ControlSpace2d.cpp -o $(DIR)/ControlSpace2d.o

$(DIR)/ControlSpace2dKdTreeL.o: ControlSpace2dKdTreeL.cpp  ControlSpace2dKdTreeL.h DataKdTree.h
	${CXX} $(CXXFLAGS) -c ControlSpace2dKdTreeL.cpp -o $(DIR)/ControlSpace2dKdTreeL.o

$(DIR)/ControlSpace2dKdTreeV.o: ControlSpace2dKdTreeV.cpp  ControlSpace2dKdTreeV.h DataKdTree.h
	${CXX} $(CXXFLAGS) -c ControlSpace2dKdTreeV.cpp -o $(DIR)/ControlSpace2dKdTreeV.o

$(DIR)/ControlSpace2dMatrixUniform.o: ControlSpace2dMatrixUniform.cpp  ControlSpace2dMatrixUniform.h common.h MeshPoint2d.h MeshEdge2d.h SurfaceParametric.h MeshContainer2d.h Curve2dSegment.h EPSFile.h ControlSpace3dAdaptive.h
	${CXX} $(CXXFLAGS) -c ControlSpace2dMatrixUniform.cpp -o $(DIR)/ControlSpace2dMatrixUniform.o

$(DIR)/ControlSpace2dMesh.o: ControlSpace2dMesh.cpp  ControlSpace2dMesh.h ControlSpace2dIdentity.h MeshContainer2d.h MeshTriangle2d.h MeshGenerator2d.h MeshData.h SurfaceParametric.h EPSFile.h Curve2dSegment.h DTriangle.h ControlSpace3dAdaptive.h
	${CXX} $(CXXFLAGS) -c ControlSpace2dMesh.cpp -o $(DIR)/ControlSpace2dMesh.o

$(DIR)/ControlSpace2dQuadTree.o: ControlSpace2dQuadTree.cpp  ControlSpace2dQuadTree.h MeshPoint2d.h MeshEdge2d.h SurfaceParametric.h MeshContainer2d.h Curve2dSegment.h EPSFile.h MeshViewSet.h ReparameterizationData.h SurfacePlane.h DTriangle.h DSegment.h ControlSpace3d.h ControlSpace3dAdaptive.h MeshGenerator2d.h
	${CXX} $(CXXFLAGS) -c ControlSpace2dQuadTree.cpp -o $(DIR)/ControlSpace2dQuadTree.o

$(DIR)/ControlSpace3dAdaptive.o: ControlSpace3dAdaptive.cpp  ControlSpace3dAdaptive.h ControlSpace2dAdaptive.h MeshContainer3d.h MeshContainer3dSurface.h MeshPoint3d.h MeshEdge3d.h MeshFace.h SurfacePlane.h Metric3dContext.h DataHashTable.h Curve3dParametric.h DPlane.h
	${CXX} $(CXXFLAGS) -c ControlSpace3dAdaptive.cpp -o $(DIR)/ControlSpace3dAdaptive.o

$(DIR)/ControlSpace3dAnalytic.o: ControlSpace3dAnalytic.cpp  ControlSpace3dAnalytic.h DEquation.h MeshLog.h MeshData.h
	${CXX} $(CXXFLAGS) -c ControlSpace3dAnalytic.cpp -o $(DIR)/ControlSpace3dAnalytic.o

$(DIR)/ControlSpace3d.o: ControlSpace3d.cpp  ControlSpace3d.h MeshData.h
	${CXX} $(CXXFLAGS) -c ControlSpace3d.cpp -o $(DIR)/ControlSpace3d.o

$(DIR)/ControlSpace3dInternal.o: ControlSpace3dInternal.cpp  ControlSpace3dInternal.h DTriangle.h
	${CXX} $(CXXFLAGS) -c ControlSpace3dInternal.cpp -o $(DIR)/ControlSpace3dInternal.o

$(DIR)/ControlSpace3dKdTree.o: ControlSpace3dKdTree.cpp  ControlSpace3dKdTree.h DataHashTable.h DLine.h DMetric3d.h DMatrix.h DataKdTree.h DataOctree.h MeshGenerator2d.h MeshContainer3d.h ControlSpace3dAnalytic.h MeshGenerator3d.h MeshGenerator3dQuality.h MeshBRepXML.h MeshViewSet.h ControlSpace3dKdTreeL.h ControlSpace3dKdTreeV.h ControlSpace3dKdTreeLi.h ControlSpace3dOctreeL.h ControlSpace3dOctreeV.h ControlSpace3dOctreeLB.h ControlSpace3dOctree.h DSphere.h SurfaceAnalytic.h SurfacePlane.h DataStatistics.h
	${CXX} $(CXXFLAGS) -c ControlSpace3dKdTree.cpp -o $(DIR)/ControlSpace3dKdTree.o

$(DIR)/ControlSpace3dKdTreeL.o: ControlSpace3dKdTreeL.cpp  ControlSpace3dKdTreeL.h DataKdTree.h
	${CXX} $(CXXFLAGS) -c ControlSpace3dKdTreeL.cpp -o $(DIR)/ControlSpace3dKdTreeL.o

$(DIR)/ControlSpace3dKdTreeLi.o: ControlSpace3dKdTreeLi.cpp  ControlSpace3dKdTreeLi.h
	${CXX} $(CXXFLAGS) -c ControlSpace3dKdTreeLi.cpp -o $(DIR)/ControlSpace3dKdTreeLi.o

$(DIR)/ControlSpace3dKdTreeV.o: ControlSpace3dKdTreeV.cpp  ControlSpace3dKdTreeV.h DataKdTree.h
	${CXX} $(CXXFLAGS) -c ControlSpace3dKdTreeV.cpp -o $(DIR)/ControlSpace3dKdTreeV.o

$(DIR)/ControlSpace3dMatrixUniform.o: ControlSpace3dMatrixUniform.cpp  ControlSpace3dMatrixUniform.h ControlSpace2dAdaptive.h DataVector.h common.h
	${CXX} $(CXXFLAGS) -c ControlSpace3dMatrixUniform.cpp -o $(DIR)/ControlSpace3dMatrixUniform.o

$(DIR)/ControlSpace3dOctree.o: ControlSpace3dOctree.cpp  ControlSpace3dOctree.h MeshPoint3d.h MeshEdge3d.h MeshFace.h MeshContainer3d.h ControlSpace2dAdaptive.h common.h DTriangle.h DSegment.h DTetrahedron.h MeshViewSet.h MeshContainer3dSurface.h
	${CXX} $(CXXFLAGS) -c ControlSpace3dOctree.cpp -o $(DIR)/ControlSpace3dOctree.o

$(DIR)/ControlSpace3dOctreeLB.o: ControlSpace3dOctreeLB.cpp  ControlSpace3dOctreeLB.h
	${CXX} $(CXXFLAGS) -c ControlSpace3dOctreeLB.cpp -o $(DIR)/ControlSpace3dOctreeLB.o

$(DIR)/ControlSpace3dOctreeL.o: ControlSpace3dOctreeL.cpp  ControlSpace3dOctreeL.h
	${CXX} $(CXXFLAGS) -c ControlSpace3dOctreeL.cpp -o $(DIR)/ControlSpace3dOctreeL.o

$(DIR)/ControlSpace3dOctreeV.o: ControlSpace3dOctreeV.cpp  ControlSpace3dOctreeV.h DataKdTree.h
	${CXX} $(CXXFLAGS) -c ControlSpace3dOctreeV.cpp -o $(DIR)/ControlSpace3dOctreeV.o

$(DIR)/Curve2dAnalytic.o: Curve2dAnalytic.cpp  Curve2dAnalytic.h DEquation.h MeshLog.h MeshData.h
	${CXX} $(CXXFLAGS) -c Curve2dAnalytic.cpp -o $(DIR)/Curve2dAnalytic.o

$(DIR)/Curve2dBSpline.o: Curve2dBSpline.cpp  Curve2dSegment.h Curve2dBSpline.h MeshData.h DMetric2d.h SurfaceParametric.h MeshViewSet.h
	${CXX} $(CXXFLAGS) -c Curve2dBSpline.cpp -o $(DIR)/Curve2dBSpline.o

$(DIR)/Curve2dCircle.o: Curve2dCircle.cpp  Curve2dCircle.h DPoint.h DEquation.h
	${CXX} $(CXXFLAGS) -c Curve2dCircle.cpp -o $(DIR)/Curve2dCircle.o

$(DIR)/Curve2dLinearQuadric.o: Curve2dLinearQuadric.cpp  Curve2dLinearQuadric.h DPoint.h DEquation.h
	${CXX} $(CXXFLAGS) -c Curve2dLinearQuadric.cpp -o $(DIR)/Curve2dLinearQuadric.o

$(DIR)/Curve2dParametric.o: Curve2dParametric.cpp  Curve2dParametric.h MeshData.h ControlSpace2d.h DRect.h SurfaceParametric.h DEquation.h Metric2dContext.h
	${CXX} $(CXXFLAGS) -c Curve2dParametric.cpp -o $(DIR)/Curve2dParametric.o

$(DIR)/Curve2dSegment.o: Curve2dSegment.cpp  Curve2dSegment.h MeshData.h ControlSpace2d.h DRect.h SurfaceParametric.h DEquation.h DMatrix.h
	${CXX} $(CXXFLAGS) -c Curve2dSegment.cpp -o $(DIR)/Curve2dSegment.o

$(DIR)/Curve3dBSpline.o: Curve3dBSpline.cpp  Curve3dSegment.h Curve3dBSpline.h MeshData.h DMetric3d.h Metric3dContext.h DLeastSquaresFitting.h DLinearQuadric.h Curve2dBSpline.h MeshViewSet.h
	${CXX} $(CXXFLAGS) -c Curve3dBSpline.cpp -o $(DIR)/Curve3dBSpline.o

$(DIR)/Curve3dLinearQuadric.o: Curve3dLinearQuadric.cpp  Curve3dLinearQuadric.h DPoint.h DEquation.h
	${CXX} $(CXXFLAGS) -c Curve3dLinearQuadric.cpp -o $(DIR)/Curve3dLinearQuadric.o

$(DIR)/Curve3dParametric.o: Curve3dParametric.cpp  Curve3dParametric.h MeshData.h DRect.h DEquation.h Metric3dContext.h DLeastSquaresFitting.h DLine.h DSphere.h DLinearQuadric.h SurfaceParametric.h Curve2dSegment.h Curve2dCircle.h Curve2dLinearQuadric.h Curve3dSurfaceParametric.h
	${CXX} $(CXXFLAGS) -c Curve3dParametric.cpp -o $(DIR)/Curve3dParametric.o

$(DIR)/Curve3dSegment.o: Curve3dSegment.cpp  Curve3dSegment.h MeshData.h DRect.h Metric3dContext.h
	${CXX} $(CXXFLAGS) -c Curve3dSegment.cpp -o $(DIR)/Curve3dSegment.o

$(DIR)/Curve3dSurfaceParametric.o: Curve3dSurfaceParametric.cpp  Curve3dSurfaceParametric.h MeshData.h DRect.h Metric3dContext.h Curve2dParametric.h SurfaceParametric.h DEquation.h
	${CXX} $(CXXFLAGS) -c Curve3dSurfaceParametric.cpp -o $(DIR)/Curve3dSurfaceParametric.o

$(DIR)/DataKdSplit2dContField.o: DataKdSplit2dContField.cpp  DataKdSplit2dContField.h
	${CXX} $(CXXFLAGS) -c DataKdSplit2dContField.cpp -o $(DIR)/DataKdSplit2dContField.o

$(DIR)/DataKdSplit2d.o: DataKdSplit2d.cpp  DataKdSplit2d.h
	${CXX} $(CXXFLAGS) -c DataKdSplit2d.cpp -o $(DIR)/DataKdSplit2d.o

$(DIR)/DataKdSplit2dSingleSource.o: DataKdSplit2dSingleSource.cpp  DataKdSplit2dSingleSource.h
	${CXX} $(CXXFLAGS) -c DataKdSplit2dSingleSource.cpp -o $(DIR)/DataKdSplit2dSingleSource.o

$(DIR)/DataKdSplit3dContField.o: DataKdSplit3dContField.cpp  DataKdSplit3dContField.h
	${CXX} $(CXXFLAGS) -c DataKdSplit3dContField.cpp -o $(DIR)/DataKdSplit3dContField.o

$(DIR)/DataKdSplit3d.o: DataKdSplit3d.cpp  DataKdSplit3d.h
	${CXX} $(CXXFLAGS) -c DataKdSplit3d.cpp -o $(DIR)/DataKdSplit3d.o

$(DIR)/DataKdSplit3dSingleSource.o: DataKdSplit3dSingleSource.cpp  DataKdSplit3dSingleSource.h ControlSpace3dKdTree.h
	${CXX} $(CXXFLAGS) -c DataKdSplit3dSingleSource.cpp -o $(DIR)/DataKdSplit3dSingleSource.o

$(DIR)/DataStatistics.o: DataStatistics.cpp  DataStatistics.h common.h
	${CXX} $(CXXFLAGS) -c DataStatistics.cpp -o $(DIR)/DataStatistics.o

$(DIR)/DEquation.o: DEquation.cpp  common.h DEquation.h
	${CXX} $(CXXFLAGS) -c DEquation.cpp -o $(DIR)/DEquation.o

$(DIR)/DHesjan.o: DHesjan.cpp  DHesjan.h common.h
	${CXX} $(CXXFLAGS) -c DHesjan.cpp -o $(DIR)/DHesjan.o

$(DIR)/DLeastSquaresFitting.o: DLeastSquaresFitting.cpp  common.h DLeastSquaresFitting.h DRect.h DMatrix.h DataMatrix.h DVectorN.h DMatrixN.h DQuadric.h DLine.h DPlane.h DPlanarQuadric.h DLinearQuadric.h SurfacePlane.h SurfaceCylinder.h DSphere.h MeshViewSet.h
	${CXX} $(CXXFLAGS) -c DLeastSquaresFitting.cpp -o $(DIR)/DLeastSquaresFitting.o

$(DIR)/DLinearQuadric.o: DLinearQuadric.cpp  common.h DLinearQuadric.h DVectorN.h DPoint.h
	${CXX} $(CXXFLAGS) -c DLinearQuadric.cpp -o $(DIR)/DLinearQuadric.o

$(DIR)/DLine.o: DLine.cpp  common.h DLine.h DMatrix.h
	${CXX} $(CXXFLAGS) -c DLine.cpp -o $(DIR)/DLine.o

$(DIR)/DMatrix.o: DMatrix.cpp  DMatrix.h
	${CXX} $(CXXFLAGS) -c DMatrix.cpp -o $(DIR)/DMatrix.o

$(DIR)/DMetric2d.o: DMetric2d.cpp  DMetric2d.h DMetric3d.h DEquation.h MeshData.h SurfaceParametric.h SurfacePlane.h ControlSpace2dAdaptive.h ControlSpace3dAdaptive.h Curve2dParametric.h
	${CXX} $(CXXFLAGS) -c DMetric2d.cpp -o $(DIR)/DMetric2d.o

$(DIR)/DMetric3d.o: DMetric3d.cpp  DMetric3d.h MeshData.h EPSFile.h DEquation.h ControlSpace2dAdaptive.h ControlSpace3dAdaptive.h
	${CXX} $(CXXFLAGS) -c DMetric3d.cpp -o $(DIR)/DMetric3d.o

$(DIR)/DOrientedBox.o: DOrientedBox.cpp  DOrientedBox.h DEquation.h MeshViewSet.h Metric3dContext.h
	${CXX} $(CXXFLAGS) -c DOrientedBox.cpp -o $(DIR)/DOrientedBox.o

$(DIR)/DPlanarQuadric.o: DPlanarQuadric.cpp  common.h DPlanarQuadric.h DVectorN.h DPoint.h MeshViewSet.h DEquation.h
	${CXX} $(CXXFLAGS) -c DPlanarQuadric.cpp -o $(DIR)/DPlanarQuadric.o

$(DIR)/DPlane.o: DPlane.cpp  common.h DPlane.h DVector.h DPoint.h DMatrix.h
	${CXX} $(CXXFLAGS) -c DPlane.cpp -o $(DIR)/DPlane.o

$(DIR)/DPoint.o: DPoint.cpp  common.h DPoint.h DEquation.h
	${CXX} $(CXXFLAGS) -c DPoint.cpp -o $(DIR)/DPoint.o

$(DIR)/DQuad.o: DQuad.cpp  DQuad.h common.h DTriangle.h
	${CXX} $(CXXFLAGS) -c DQuad.cpp -o $(DIR)/DQuad.o

$(DIR)/DQuadric.o: DQuadric.cpp  common.h DQuadric.h DVectorN.h DPoint.h
	${CXX} $(CXXFLAGS) -c DQuadric.cpp -o $(DIR)/DQuadric.o

$(DIR)/DRect.o: DRect.cpp  DRect.h DEquation.h MeshViewSet.h Metric3dContext.h
	${CXX} $(CXXFLAGS) -c DRect.cpp -o $(DIR)/DRect.o

$(DIR)/DSegment.o: DSegment.cpp  common.h DSegment.h GeometricPredicates.h
	${CXX} $(CXXFLAGS) -c DSegment.cpp -o $(DIR)/DSegment.o

$(DIR)/DSphere.o: DSphere.cpp  DSphere.h DMatrix.h SurfaceCurvature.h DEquation.h MeshViewSet.h
	${CXX} $(CXXFLAGS) -c DSphere.cpp -o $(DIR)/DSphere.o

$(DIR)/DSurfaceFitting.o: DSurfaceFitting.cpp  common.h DSurfaceFitting.h Metric3dContext.h MeshPoint3d.h SurfaceParametric.h MeshViewSet.h
	${CXX} $(CXXFLAGS) -c DSurfaceFitting.cpp -o $(DIR)/DSurfaceFitting.o

$(DIR)/DTetrahedron.o: DTetrahedron.cpp  DTetrahedron.h common.h DMatrix.h
	${CXX} $(CXXFLAGS) -c DTetrahedron.cpp -o $(DIR)/DTetrahedron.o

$(DIR)/DTriangle.o: DTriangle.cpp  DTriangle.h common.h
	${CXX} $(CXXFLAGS) -c DTriangle.cpp -o $(DIR)/DTriangle.o

$(DIR)/DVector.o: DVector.cpp  common.h DVector.h DPoint.h DEquation.h
	${CXX} $(CXXFLAGS) -c DVector.cpp -o $(DIR)/DVector.o

$(DIR)/EPSFile.o: EPSFile.cpp  EPSFile.h
	${CXX} $(CXXFLAGS) -c EPSFile.cpp -o $(DIR)/EPSFile.o

$(DIR)/FrontEdge.o: FrontEdge.cpp  common.h FrontEdge.h MeshPoint2d.h MeshEdge2d.h MeshElement.h ControlSpace2d.h
	${CXX} $(CXXFLAGS) -c FrontEdge.cpp -o $(DIR)/FrontEdge.o

$(DIR)/FrontFace.o: FrontFace.cpp  FrontFace.h MeshPoint3d.h MeshEdge3d.h MeshFace.h MeshBlock.h Metric3dContext.h DTriangle.h MeshViewSet.h
	${CXX} $(CXXFLAGS) -c FrontFace.cpp -o $(DIR)/FrontFace.o

$(DIR)/FrontLines.o: FrontLines.cpp  FrontLines.h FrontEdge.h MeshEdge2d.h MeshPoint2d.h
	${CXX} $(CXXFLAGS) -c FrontLines.cpp -o $(DIR)/FrontLines.o

$(DIR)/GeometricPredicates.o: GeometricPredicates.cpp  GeometricPredicates.h
	${CXX} $(CXXFLAGS) -c GeometricPredicates.cpp -o $(DIR)/GeometricPredicates.o

$(DIR)/IteratorEdge2d.o: IteratorEdge2d.cpp  IteratorEdge2d.h MeshContainer2d.h MeshEdge2d.h
	${CXX} $(CXXFLAGS) -c IteratorEdge2d.cpp -o $(DIR)/IteratorEdge2d.o

$(DIR)/IteratorEdge3d.o: IteratorEdge3d.cpp  IteratorEdge3d.h MeshContainer3d.h MeshContainer3dSurface.h MeshEdge3d.h
	${CXX} $(CXXFLAGS) -c IteratorEdge3d.cpp -o $(DIR)/IteratorEdge3d.o

$(DIR)/IteratorFace.o: IteratorFace.cpp  IteratorFace.h MeshContainer3d.h MeshBlock.h MeshFace.h MeshDomainVolume.h
	${CXX} $(CXXFLAGS) -c IteratorFace.cpp -o $(DIR)/IteratorFace.o

$(DIR)/IteratorMesh2d.o: IteratorMesh2d.cpp  IteratorMesh2d.h MeshContainer3d.h MeshDomainVolume.h MeshDomainSurface.h common.h
	${CXX} $(CXXFLAGS) -c IteratorMesh2d.cpp -o $(DIR)/IteratorMesh2d.o

$(DIR)/IteratorMesh3d.o: IteratorMesh3d.cpp  IteratorMesh3d.h MeshContainer3d.h MeshDomainVolume.h common.h
	${CXX} $(CXXFLAGS) -c IteratorMesh3d.cpp -o $(DIR)/IteratorMesh3d.o

$(DIR)/MarchingCubes.o: MarchingCubes.cpp  MarchingCubes.h common.h MeshContainer3dSurface.h MeshPoint3d.h MeshEdge3d.h MeshTriangle3d.h MeshViewSet.h ControlSpace2dIdentity.h MeshDomainVolume.h MeshGenerator3dSurface.h Metric3dContext.h
	${CXX} $(CXXFLAGS) -c MarchingCubes.cpp -o $(DIR)/MarchingCubes.o

$(DIR)/MeshArea.o: MeshArea.cpp  DPoint.h MeshPoint2d.h MeshPoint3d.h MeshEdge2d.h MeshArea.h Curve2dParametric.h SurfaceParametric.h MeshData.h MeshViewSet.h
	${CXX} $(CXXFLAGS) -c MeshArea.cpp -o $(DIR)/MeshArea.o

$(DIR)/MeshBlock.o: MeshBlock.cpp  MeshBlock.h MeshPoint3d.h DRect.h MeshFace.h MeshEdge3d.h MeshData.h MeshTriangle3d.h MeshQuad3d.h MeshElement.h MeshTetrahedron.h DMetric3d.h
	${CXX} $(CXXFLAGS) -c MeshBlock.cpp -o $(DIR)/MeshBlock.o

$(DIR)/MeshBRep.o: MeshBRep.cpp  MeshBRep.h MeshData.h MeshContainer2d.h MeshContainer3d.h MeshPoint2d.h MeshEdge2d.h MeshEdge2dCurve.h MeshPoint3d.h MeshTriangle2d.h MeshQuad2d.h MeshDomainSurface.h MeshDomainVolume.h MeshArea.h MeshTetrahedron.h Curve2dAnalytic.h Curve2dCircle.h Curve2dBSpline.h SurfaceParametric.h SurfaceTranslated.h SurfacePlane.h SurfaceAnalytic.h SurfaceBSplinePlanar.h SurfaceMulti.h MeshEdge3d.h MeshDomainEdge3d.h ControlSpace2dMatrixUniform.h ControlSpace2dMesh.h ControlSpace2dQuadTree.h ControlSpace2dAnalytic.h ControlSpace2dIdentity.h ControlSpace3dAnalytic.h ControlSpace3dMatrixUniform.h ControlSpace3dOctree.h ControlSpace3dInternal.h DataHashTable.h MeshBoundaryCondition.h MeshViewSet.h DPlane.h DLeastSquaresFitting.h MeshContainer3dSurface.h MeshGenerator3dSurface.h DataVector.h
	${CXX} $(CXXFLAGS) -c MeshBRep.cpp -o $(DIR)/MeshBRep.o

$(DIR)/MeshBRepRSurf.o: MeshBRepRSurf.cpp  MeshBRepRSurf.h MeshData.h DPoint.h DataHashTable.h MeshBoundaryCondition.h MeshViewSet.h
	${CXX} $(CXXFLAGS) -c MeshBRepRSurf.cpp -o $(DIR)/MeshBRepRSurf.o

$(DIR)/MeshBRepXML.o: MeshBRepXML.cpp  MeshBRepXML.h MeshData.h DPoint.h tinyxml.h Curve2dCircle.h Curve2dAnalytic.h Curve2dBSpline.h Curve3dParametric.h Curve3dBSpline.h Curve3dSegment.h Curve3dSurfaceParametric.h SurfacePlane.h SurfaceAnalytic.h SurfaceBSplinePlanar.h SurfaceCorrected.h SurfaceCylinder.h SurfaceDomainHull.h DPlanarQuadric.h SurfacePlanarQuadric.h MeshBlock.h MeshDomainVolume.h MeshContainer3d.h MeshTriangle3d.h MeshQuad3d.h MeshPoly3d.h MeshViewSet.h MeshGenerator3dSurface.h MeshContainer3dSurface.h Metric3dContext.h
	${CXX} $(CXXFLAGS) -c MeshBRepXML.cpp -o $(DIR)/MeshBRepXML.o

$(DIR)/MeshContainer2d.o: MeshContainer2d.cpp  MeshContainer2d.h MeshEdge2d.h MeshEdge3d.h DataContainer.h Curve2dParametric.h MeshArea.h SurfaceParametric.h ControlSpace2dIdentity.h ControlSpace2dMatrixUniform.h ControlSpace2dMesh.h ControlSpace2dQuadTree.h ControlSpace3d.h QuadTree.h SurfacePlane.h EPSFile.h MeshViewSet.h MeshBoundaryCondition.h MeshGenerator2d.h IteratorEdge2d.h DataHashTable.h DTriangle.h DQuad.h
	${CXX} $(CXXFLAGS) -c MeshContainer2d.cpp -o $(DIR)/MeshContainer2d.o

$(DIR)/MeshContainer3d.o: MeshContainer3d.cpp  MeshContainer3d.h MeshPoint3d.h DataContainer.h Curve2dParametric.h Curve3dParametric.h MeshEdge3d.h MeshDomainEdge3d.h MeshTriangle3d.h MeshDomainVolume.h MeshDomainSurface.h OctTree.h MeshViewSet.h ControlSpace2dIdentity.h ControlSpace3dMatrixUniform.h ControlSpace3dOctree.h DMetric3d.h MeshBoundaryCondition.h MeshContainer2d.h ControlSpace2d.h MeshGenerator2d.h DataHashTable.h DTriangle.h DTetrahedron.h
	${CXX} $(CXXFLAGS) -c MeshContainer3d.cpp -o $(DIR)/MeshContainer3d.o

$(DIR)/MeshContainer3dSurface.o: MeshContainer3dSurface.cpp  MeshContainer3dSurface.h MeshViewSet.h MeshEdge3d.h MeshFace.h MeshDomainVolume.h Metric3dContext.h DLeastSquaresFitting.h MeshGenerator3dSurface.h DataList.h DataStatistics.h SurfaceParametric.h SurfacePlane.h SurfacePlanarQuadric.h SurfaceCylinder.h SurfaceBSplinePlanar.h SurfaceBSplineCylindrical.h Curve3dParametric.h Curve3dBSpline.h Curve2dSegment.h Curve3dSurfaceParametric.h DPlane.h DLine.h DTriangle.h DPlanarQuadric.h SurfaceDomainHull.h SurfaceDomainFaces.h SurfaceDomainFacesPlanar.h ControlSpace3dAdaptive.h ControlSpace2dIdentity.h DOrientedBox.h DLinearQuadric.h Curve3dLinearQuadric.h SurfaceQuadricQTree.h
	${CXX} $(CXXFLAGS) -c MeshContainer3dSurface.cpp -o $(DIR)/MeshContainer3dSurface.o

$(DIR)/MeshData.o: MeshData.cpp  MeshData.h DPoint.h ControlSpace2d.h ControlSpace3d.h ControlSpace2dMatrixUniform.h ControlSpace3dMatrixUniform.h ControlSpace2dQuadTree.h ControlSpace3dOctree.h ControlSpace2dMesh.h MeshLog.h QuadTree.h OctTree.h MeshGenerator1d.h MeshGenerator2d.h MeshGenerator2dQuad.h MeshGenerator2dQMorph.h MeshGenerator3d.h MeshGenerator3dDelaunayBoundary.h MeshGenerator3dQuality.h MeshGenerator3dSurface.h DHesjan.h Metric2dContext.h Metric3dContext.h FrontEdge.h MeshViewSet.h MeshGrain.h MeshPoint3d.h MeshContainer3dSurface.h
	${CXX} $(CXXFLAGS) -c MeshData.cpp -o $(DIR)/MeshData.o

$(DIR)/MeshDecompositionUTC.o: MeshDecompositionUTC.cpp  MeshDecompositionUTC.h MeshContainer2d.h MeshContainer3d.h MeshPoint3d.h MeshPoint2d.h MeshArea.h MeshFace.h MeshBlock.h MeshDomainVolume.h DataVector.h MeshViewSet.h MeshTriangle3d.h MeshQuad3d.h MeshEdge2d.h MeshEdge3d.h SurfacePlane.h ControlSpace2dMatrixUniform.h ControlSpace3dOctree.h MeshGenerator2d.h DataHashTable.h
	${CXX} $(CXXFLAGS) -c MeshDecompositionUTC.cpp -o $(DIR)/MeshDecompositionUTC.o

$(DIR)/MeshDomainEdge3d.o: MeshDomainEdge3d.cpp  MeshDomainEdge3d.h ControlSpace3dAdaptive.h DMetric3d.h
	${CXX} $(CXXFLAGS) -c MeshDomainEdge3d.cpp -o $(DIR)/MeshDomainEdge3d.o

$(DIR)/MeshDomainSurface.o: MeshDomainSurface.cpp  MeshDomainSurface.h MeshDomainVolume.h MeshPoint3d.h MeshEdge3d.h MeshDomainEdge3d.h SurfaceParametric.h MeshContainer2d.h SurfacePlane.h MeshArea.h MeshPoint2d.h MeshEdge2d.h MeshQuad2d.h QuadTree.h MeshGenerator2d.h MeshGenerator2dQuad.h MeshGenerator2dQMorph.h ControlSpace2d.h ControlSpace2dMatrixUniform.h EPSFile.h DEquation.h DTriangle.h DQuad.h MeshViewSet.h MeshTriangle3d.h MeshQuad3d.h IteratorFace.h MeshContainer3d.h MeshContainer3dSurface.h MeshGenerator3dSurface.h
	${CXX} $(CXXFLAGS) -c MeshDomainSurface.cpp -o $(DIR)/MeshDomainSurface.o

$(DIR)/MeshDomainVolume.o: MeshDomainVolume.cpp  MeshDomainVolume.h MeshDomainSurface.h MeshEdge2d.h MeshFace.h MeshPoint2d.h MeshPoint3d.h MeshEdge3d.h MeshContainer2d.h MeshContainer3d.h MeshContainer3dSurface.h SurfaceParametric.h MeshTriangle3d.h MeshQuad3d.h MeshGenerator2d.h MeshGenerator3dSurface.h MeshGenerator3d.h MeshGenerator3dAdapt.h MeshGenerator3dDelaunayBoundary.h MeshGenerator3dDirectBoundary.h MeshGenerator3dQuality.h MeshBoundaryCondition.h ControlSpace3d.h ControlSpace2d.h ControlSpace2dAdaptive.h ControlSpace3dMatrixUniform.h ControlSpace3dOctree.h ControlSpace3dKdTreeL.h ControlSpace3dKdTreeV.h OctTree.h DataHashTable.h IteratorEdge3d.h MeshViewSet.h
	${CXX} $(CXXFLAGS) -c MeshDomainVolume.cpp -o $(DIR)/MeshDomainVolume.o

$(DIR)/MeshEdge2d.o: MeshEdge2d.cpp  DPoint.h MeshEdge2d.h MeshPoint3d.h MeshPoint2d.h MeshElement.h MeshTriangle2d.h SurfaceParametric.h Curve2dSegment.h ControlSpace2d.h MeshViewSet.h MeshEdge3d.h
	${CXX} $(CXXFLAGS) -c MeshEdge2d.cpp -o $(DIR)/MeshEdge2d.o

$(DIR)/MeshEdge2dCurve.o: MeshEdge2dCurve.cpp  MeshEdge2dCurve.h Curve2dParametric.h SurfaceParametric.h ControlSpace2dAdaptive.h MeshPoint2d.h MeshViewSet.h
	${CXX} $(CXXFLAGS) -c MeshEdge2dCurve.cpp -o $(DIR)/MeshEdge2dCurve.o

$(DIR)/MeshEdge3d.o: MeshEdge3d.cpp  MeshEdge3d.h MeshPoint3d.h MeshEdge2d.h MeshFace.h MeshDomainSurface.h MeshContainer2d.h MeshArea.h Curve2dParametric.h Curve3dParametric.h SurfaceParametric.h MeshBoundaryCondition.h MeshViewSet.h MeshTetrahedron.h Metric3dContext.h
	${CXX} $(CXXFLAGS) -c MeshEdge3d.cpp -o $(DIR)/MeshEdge3d.o

$(DIR)/MeshElement.o: MeshElement.cpp  MeshElement.h MeshPoint2d.h MeshEdge2d.h SurfaceParametric.h MeshViewSet.h MeshPoint3d.h SurfacePlane.h DTriangle.h
	${CXX} $(CXXFLAGS) -c MeshElement.cpp -o $(DIR)/MeshElement.o

$(DIR)/MeshFace.o: MeshFace.cpp  MeshData.h MeshFace.h MeshPoint3d.h MeshViewSet.h MeshTetrahedron.h MeshBoundaryCondition.h MeshEdge3d.h Metric3dContext.h DTriangle.h MeshGenerator3dSurface.h SurfaceParametric.h MeshTriangle3d.h MeshGenerator2d.h MeshContainer2d.h MeshTriangle2d.h DPlane.h GeometricPredicates.h
	${CXX} $(CXXFLAGS) -c MeshFace.cpp -o $(DIR)/MeshFace.o

$(DIR)/MeshGenerator1d.o: MeshGenerator1d.cpp  MeshGenerator1d.h MeshContainer2d.h MeshContainer3d.h MeshPoint2d.h MeshPoint3d.h MeshEdge2d.h MeshDomainEdge3d.h SurfaceParametric.h MeshDomainSurface.h MeshDomainVolume.h MeshArea.h ControlSpace2d.h ControlSpace2dAdaptive.h ControlSpace3d.h Metric3dContext.h
	${CXX} $(CXXFLAGS) -c MeshGenerator1d.cpp -o $(DIR)/MeshGenerator1d.o

$(DIR)/MeshGenerator2d.o: MeshGenerator2d.cpp  MeshGenerator1d.h MeshGenerator2d.h MeshGenerator3d.h MeshContainer2d.h MeshContainer3d.h MeshDomainVolume.h MeshDomainSurface.h MeshDomainEdge3d.h MeshTriangle2d.h QuadTree.h ControlSpace2d.h ControlSpace2dIdentity.h ControlSpace2dMatrixUniform.h ControlSpace2dMesh.h ControlSpace2dProjected.h ControlSpace2dQuadTree.h ControlSpace2dKdTreeL.h ControlSpace2dKdTreeV.h SurfaceParametric.h MeshViewSet.h DEquation.h MeshEdge3d.h DataVector.h DataList.h IteratorEdge2d.h IteratorEdge3d.h DTriangle.h OctTree.h ControlSpace3dOctree.h ControlSpace3dKdTreeL.h ControlSpace3dKdTreeV.h DLeastSquaresFitting.h DPlane.h SurfacePlane.h MeshArea.h Metric2dContext.h
	${CXX} $(CXXFLAGS) -c MeshGenerator2d.cpp -o $(DIR)/MeshGenerator2d.o

$(DIR)/MeshGenerator2dQMorph.o: MeshGenerator2dQMorph.cpp  MeshGenerator2d.h MeshGenerator2dQuad.h MeshGenerator2dQMorph.h MeshContainer2d.h FrontEdge.h FrontLines.h MeshEdge2d.h MeshPoint2d.h MeshTriangle2d.h MeshQuad2d.h SurfaceParametric.h ControlSpace2d.h DEquation.h DataVector.h DTriangle.h
	${CXX} $(CXXFLAGS) -c MeshGenerator2dQMorph.cpp -o $(DIR)/MeshGenerator2dQMorph.o

$(DIR)/MeshGenerator2dQuad.o: MeshGenerator2dQuad.cpp  MeshGenerator2dQuad.h MeshContainer2d.h MeshContainer3d.h MeshGenerator2d.h MeshDomainVolume.h MeshDomainSurface.h FrontLines.h MeshTriangle2d.h MeshQuad2d.h SurfaceParametric.h MeshGenerator2dQMorph.h DTriangle.h DQuad.h
	${CXX} $(CXXFLAGS) -c MeshGenerator2dQuad.cpp -o $(DIR)/MeshGenerator2dQuad.o

$(DIR)/MeshGenerator3dAdapt.o: MeshGenerator3dAdapt.cpp  MeshGenerator3dAdapt.h MeshContainer3d.h MeshContainer3dSurface.h MeshTetrahedron.h MeshPoint3d.h MeshEdge3d.h MeshFace.h MeshTriangle3d.h MeshBlock.h DataVector.h DataList.h DataHeapVector.h MeshData.h MeshLog.h ControlSpace3dAdaptive.h ControlSpace3dOctree.h MeshGenerator3d.h MeshGenerator3dSurface.h MeshGenerator3dQuality.h MeshGenerator3dDelaunayBoundary.h MeshSplit3d.h MeshViewSet.h SurfacePlane.h SurfacePlanarQuadric.h SurfaceBSplinePlanar.h IteratorEdge3d.h IteratorFace.h ControlSpace2dAdaptive.h GeometricPredicates.h DTetrahedron.h DTriangle.h DPlane.h DPlanarQuadric.h DLeastSquaresFitting.h
	${CXX} $(CXXFLAGS) -c MeshGenerator3dAdapt.cpp -o $(DIR)/MeshGenerator3dAdapt.o

$(DIR)/MeshGenerator3d.o: MeshGenerator3d.cpp  MeshGenerator3d.h MeshContainer3d.h MeshDomainVolume.h MeshDomainSurface.h MeshContainer2d.h OctTree.h MeshTetrahedron.h MeshPoint2d.h MeshPoint3d.h MeshEdge3d.h MeshFace.h MeshTriangle3d.h MeshBlock.h DataVector.h DataList.h DataHeapVector.h MeshData.h MeshLog.h MeshBoundaryCondition.h ControlSpace3d.h ControlSpace3dOctree.h ControlSpace2dIdentity.h MeshGenerator2d.h MeshGenerator3dQuality.h MeshGenerator3dDelaunayBoundary.h MeshViewSet.h SurfacePlane.h IteratorEdge3d.h IteratorFace.h ControlSpace2dAdaptive.h GeometricPredicates.h DTetrahedron.h DTriangle.h MeshDomainEdge3d.h
	${CXX} $(CXXFLAGS) -c MeshGenerator3d.cpp -o $(DIR)/MeshGenerator3d.o

$(DIR)/MeshGenerator3dDelaunayBoundary.o: MeshGenerator3dDelaunayBoundary.cpp  MeshGenerator3dDelaunayBoundary.h MeshGenerator3d.h MeshContainer3d.h MeshContainer3dSurface.h MeshDomainVolume.h OctTree.h MeshTetrahedron.h MeshPoint2d.h MeshPoint3d.h MeshEdge3d.h MeshFace.h MeshTriangle3d.h MeshBlock.h DataVector.h DataList.h DataHeapVector.h MeshData.h MeshLog.h MeshBoundaryCondition.h ControlSpace3d.h ControlSpace3dOctree.h MeshGenerator2d.h MeshGenerator3dQuality.h MeshViewSet.h SurfacePlane.h IteratorEdge3d.h IteratorFace.h ControlSpace2dAdaptive.h GeometricPredicates.h DTetrahedron.h DTriangle.h DSegment.h
	${CXX} $(CXXFLAGS) -c MeshGenerator3dDelaunayBoundary.cpp -o $(DIR)/MeshGenerator3dDelaunayBoundary.o

$(DIR)/MeshGenerator3dDirectBoundary.o: MeshGenerator3dDirectBoundary.cpp  MeshGenerator3dDirectBoundary.h DPoint.h DTetrahedron.h MeshContainer3d.h MeshEdge3d.h MeshFace.h MeshTriangle3d.h MeshTetrahedron.h Metric3dContext.h ControlSpace3dAdaptive.h MeshDomainVolume.h MeshGenerator3d.h MeshViewSet.h DataHashTable.h FrontFace.h DataContainer.h MeshContainer3dSurface.h
	${CXX} $(CXXFLAGS) -c MeshGenerator3dDirectBoundary.cpp -o $(DIR)/MeshGenerator3dDirectBoundary.o

$(DIR)/MeshGenerator3dQuality.o: MeshGenerator3dQuality.cpp  MeshGenerator3dQuality.h MeshContainer3d.h MeshDomainVolume.h OctTree.h MeshTetrahedron.h MeshPoint3d.h MeshEdge3d.h MeshFace.h MeshBlock.h DataVector.h MeshData.h MeshLog.h MeshBoundaryCondition.h ControlSpace3d.h MeshGenerator3d.h MeshGenerator3dAdapt.h MeshGenerator3dDelaunayBoundary.h DataStatistics.h DPoint.h MeshViewSet.h DataStatistics.h GeometricPredicates.h
	${CXX} $(CXXFLAGS) -c MeshGenerator3dQuality.cpp -o $(DIR)/MeshGenerator3dQuality.o

$(DIR)/MeshGenerator3dSurface.o: MeshGenerator3dSurface.cpp  MeshGenerator3dSurface.h MeshContainer3d.h MeshContainer3dSurface.h MeshFace.h MeshPoint3d.h MeshEdge3d.h DataVector.h DataHashTable.h MeshDomainVolume.h DataMatrix.h DLeastSquaresFitting.h SurfaceParametric.h SurfacePlane.h SurfacePlanarQuadric.h Curve3dParametric.h DTriangle.h Metric3dContext.h MeshTriangle3d.h MeshQuad3d.h MeshPoly3d.h ControlSpace3dOctree.h ControlSpace2dIdentity.h DPlanarQuadric.h DLine.h DPlane.h DTriangle.h MeshGenerator2d.h MeshContainer2d.h MeshTriangle2d.h SurfaceAnalytic.h DataHeapVector.h DataPtrMatrix.h MeshViewSet.h
	${CXX} $(CXXFLAGS) -c MeshGenerator3dSurface.cpp -o $(DIR)/MeshGenerator3dSurface.o

$(DIR)/MeshGrain.o: MeshGrain.cpp  MeshGrain.h MeshViewSet.h DMetric3d.h DataMatrix.h SurfacePlane.h SurfaceCylinder.h SurfaceBSplinePlanar.h SurfaceBSplineCylindrical.h SurfaceMulti.h SurfaceCorrected.h DataHashTable.h DLeastSquaresFitting.h DLine.h Curve2dSegment.h Curve2dBSpline.h DSegment.h DataList.h MeshPoint3d.h MeshGenerator3d.h MeshGenerator3dDelaunayBoundary.h MeshContainer3d.h ControlSpace2dIdentity.h ControlSpace3d.h Metric3dContext.h DPlane.h
	${CXX} $(CXXFLAGS) -c MeshGrain.cpp -o $(DIR)/MeshGrain.o

$(DIR)/MeshLog.o: MeshLog.cpp  MeshLog.h
	${CXX} $(CXXFLAGS) -c MeshLog.cpp -o $(DIR)/MeshLog.o

$(DIR)/MeshModel.o: MeshModel.cpp  MeshModel.h common.h IteratorMesh2d.h IteratorMesh3d.h MeshContainer2d.h MeshContainer3d.h MeshContainer3dSurface.h MeshDomainVolume.h MeshDomainSurface.h DEquation.h MeshGenerator1d.h MeshGenerator2d.h MeshGenerator2dQuad.h MeshGenerator3d.h MeshGenerator3dQuality.h MeshSpecialRoutinesDAT.h MeshSpecialRoutinesPaszynski.h MeshBRepXML.h ControlSpace2dAdaptive.h ControlSpace2dMesh.h ControlSpace3dAdaptive.h SurfaceParametric.h MeshEdge3d.h SurfaceAnalytic.h SurfacePlane.h ControlSpace2dQuadTree.h EPSFile.h MeshDecompositionUTC.h DataStatistics.h MeshGrain.h MeshViewSet.h Metric2dContext.h MeshElement.h MeshTriangle2d.h MeshSpecialRoutinesUTC.h MeshSplit3d.h ControlSpace2dIdentity.h DLeastSquaresFitting.h MeshGenerator3dAdapt.h ModelMultiSphere.h MarchingCubes.h ControlSpace3dKdTree.h DSegment.h DTriangle.h
	${CXX} $(CXXFLAGS) -c MeshModel.cpp -o $(DIR)/MeshModel.o

$(DIR)/MeshPoint2d.o: MeshPoint2d.cpp  MeshPoint2d.h MeshEdge2d.h MeshElement.h MeshPoint3d.h DHesjan.h SurfaceParametric.h MeshViewSet.h ControlSpace2d.h Metric2dContext.h
	${CXX} $(CXXFLAGS) -c MeshPoint2d.cpp -o $(DIR)/MeshPoint2d.o

$(DIR)/MeshPoint3d.o: MeshPoint3d.cpp  MeshPoint2d.h MeshPoint3d.h MeshEdge3d.h MeshFace.h MeshViewSet.h DataHashTable.h Metric3dContext.h SurfaceParametric.h Curve3dParametric.h ControlSpace3dAdaptive.h DataStatistics.h
	${CXX} $(CXXFLAGS) -c MeshPoint3d.cpp -o $(DIR)/MeshPoint3d.o

$(DIR)/MeshPoly3d.o: MeshPoly3d.cpp  MeshPoly3d.h MeshEdge3d.h MeshTriangle3d.h DTriangle.h DQuad.h DPlane.h DSegment.h DLeastSquaresFitting.h SurfaceParametric.h SurfacePlane.h ControlSpace2d.h ControlSpace2dIdentity.h MeshPoint2d.h MeshContainer2d.h MeshGenerator2d.h MeshArea.h MeshEdge2d.h MeshTriangle2d.h MeshQuad3d.h MeshViewSet.h Metric3dContext.h
	${CXX} $(CXXFLAGS) -c MeshPoly3d.cpp -o $(DIR)/MeshPoly3d.o

$(DIR)/MeshQuad2d.o: MeshQuad2d.cpp  MeshQuad2d.h ControlSpace2d.h MeshPoint2d.h MeshEdge2d.h Metric2dContext.h DTriangle.h DQuad.h
	${CXX} $(CXXFLAGS) -c MeshQuad2d.cpp -o $(DIR)/MeshQuad2d.o

$(DIR)/MeshQuad3d.o: MeshQuad3d.cpp  MeshQuad3d.h MeshEdge3d.h MeshTriangle3d.h DTriangle.h MeshViewSet.h DQuad.h
	${CXX} $(CXXFLAGS) -c MeshQuad3d.cpp -o $(DIR)/MeshQuad3d.o

$(DIR)/MeshSpecialRoutinesDAT.o: MeshSpecialRoutinesDAT.cpp  MeshData.h MeshSpecialRoutinesDAT.h MeshContainer2d.h MeshContainer3d.h MeshContainer3dSurface.h MeshPoint2d.h MeshPoint3d.h MeshEdge2d.h MeshEdge2dCurve.h MeshEdge3d.h MeshDomainEdge3d.h MeshArea.h MeshFace.h MeshDomainSurface.h MeshDomainVolume.h MeshQuad2d.h Curve2dAnalytic.h SurfacePlane.h ControlSpace2dAnalytic.h ControlSpace2dMatrixUniform.h ControlSpace2dMesh.h ControlSpace2dIdentity.h MeshTriangle3d.h Quad12.h MeshGenerator1d.h MeshGenerator2d.h MeshGenerator2dQuad.h ControlSpace3dOctree.h DHesjan.h DataHashTable.h DSegment.h DTriangle.h DLine.h
	${CXX} $(CXXFLAGS) -c MeshSpecialRoutinesDAT.cpp -o $(DIR)/MeshSpecialRoutinesDAT.o

$(DIR)/MeshSpecialRoutinesPaszynski.o: MeshSpecialRoutinesPaszynski.cpp  MeshData.h MeshSpecialRoutinesPaszynski.h MeshContainer3d.h MeshPoint3d.h MeshEdge3d.h MeshTriangle3d.h MeshTetrahedron.h MeshDomainSurface.h MeshDomainVolume.h
	${CXX} $(CXXFLAGS) -c MeshSpecialRoutinesPaszynski.cpp -o $(DIR)/MeshSpecialRoutinesPaszynski.o

$(DIR)/MeshSpecialRoutinesUTC.o: MeshSpecialRoutinesUTC.cpp  MeshData.h MeshSpecialRoutinesUTC.h MeshContainer2d.h MeshContainer3d.h MeshBlock.h MeshTetrahedron.h MeshTriangle3d.h MeshEdge3d.h MeshViewSet.h ControlSpace3dOctree.h MeshGenerator3d.h MeshGenerator3dAdapt.h MeshGenerator3dQuality.h MeshSplit3d.h MeshContainer3dSurface.h MeshGenerator3dSurface.h MeshDomainVolume.h DLeastSquaresFitting.h SurfacePlane.h GeometricPredicates.h DPlane.h
	${CXX} $(CXXFLAGS) -c MeshSpecialRoutinesUTC.cpp -o $(DIR)/MeshSpecialRoutinesUTC.o

$(DIR)/MeshSplit3d.o: MeshSplit3d.cpp  MeshData.h MeshSplit3d.h MeshContainer3d.h MeshEdge3d.h MeshFace.h MeshTriangle3d.h MeshViewSet.h
	${CXX} $(CXXFLAGS) -c MeshSplit3d.cpp -o $(DIR)/MeshSplit3d.o

$(DIR)/MeshTetrahedron.o: MeshTetrahedron.cpp  MeshTetrahedron.h MeshTriangle3d.h MeshData.h ControlSpace3d.h MeshViewSet.h Metric3dContext.h MeshGenerator3d.h MeshEdge3d.h DTriangle.h DTetrahedron.h GeometricPredicates.h
	${CXX} $(CXXFLAGS) -c MeshTetrahedron.cpp -o $(DIR)/MeshTetrahedron.o

$(DIR)/MeshTriangle2d.o: MeshTriangle2d.cpp  MeshTriangle2d.h common.h MeshEdge2d.h MeshData.h MeshGenerator2d.h ControlSpace2d.h DMetric2d.h MeshViewSet.h DTriangle.h
	${CXX} $(CXXFLAGS) -c MeshTriangle2d.cpp -o $(DIR)/MeshTriangle2d.o

$(DIR)/MeshTriangle3d.o: MeshTriangle3d.cpp  MeshTriangle3d.h MeshEdge3d.h DTriangle.h Metric3dContext.h SurfaceParametric.h MeshViewSet.h
	${CXX} $(CXXFLAGS) -c MeshTriangle3d.cpp -o $(DIR)/MeshTriangle3d.o

$(DIR)/MeshViewSet.o: MeshViewSet.cpp  MeshViewSet.h MeshPoint2d.h MeshEdge2d.h MeshPoint3d.h MeshEdge3d.h MeshElement.h MeshFace.h MeshBlock.h SurfaceParametric.h MeshContainer2d.h MeshContainer3d.h MeshTetrahedron.h DLeastSquaresFitting.h SurfacePlane.h DTriangle.h DPlane.h
	${CXX} $(CXXFLAGS) -c MeshViewSet.cpp -o $(DIR)/MeshViewSet.o

$(DIR)/Metric2dContext.o: Metric2dContext.cpp  Metric2dContext.h ControlSpace2d.h ControlSpace2dAdaptive.h SurfaceParametric.h MeshPoint2d.h MeshTriangle2d.h ReparameterizationData.h
	${CXX} $(CXXFLAGS) -c Metric2dContext.cpp -o $(DIR)/Metric2dContext.o

$(DIR)/Metric3dContext.o: Metric3dContext.cpp  Metric3dContext.h ControlSpace3d.h MeshPoint3d.h MeshTetrahedron.h ControlSpace3dAdaptive.h
	${CXX} $(CXXFLAGS) -c Metric3dContext.cpp -o $(DIR)/Metric3dContext.o

$(DIR)/ModelMultiSphere.o: ModelMultiSphere.cpp  ModelMultiSphere.h ControlSpace2dAdaptive.h MeshPoint3d.h
	${CXX} $(CXXFLAGS) -c ModelMultiSphere.cpp -o $(DIR)/ModelMultiSphere.o

$(DIR)/nrutil.o: nrutil.cpp  nrutil.h common.h
	${CXX} $(CXXFLAGS) -c nrutil.cpp -o $(DIR)/nrutil.o

$(DIR)/OctTree.o: OctTree.cpp  OctTree.h MeshData.h MeshPoint3d.h
	${CXX} $(CXXFLAGS) -c OctTree.cpp -o $(DIR)/OctTree.o

$(DIR)/Quad12.o: Quad12.cpp  Quad12.h MeshQuad2d.h MeshEdge2d.h MeshPoint2d.h nrutil.h
	${CXX} $(CXXFLAGS) -c Quad12.cpp -o $(DIR)/Quad12.o

$(DIR)/QuadTree.o: QuadTree.cpp  QuadTree.h MeshData.h DRect.h
	${CXX} $(CXXFLAGS) -c QuadTree.cpp -o $(DIR)/QuadTree.o

$(DIR)/StatisticCounter.o: StatisticCounter.cpp  StatisticCounter.h common.h
	${CXX} $(CXXFLAGS) -c StatisticCounter.cpp -o $(DIR)/StatisticCounter.o

$(DIR)/SurfaceAnalytic.o: SurfaceAnalytic.cpp  SurfaceAnalytic.h MeshLog.h
	${CXX} $(CXXFLAGS) -c SurfaceAnalytic.cpp -o $(DIR)/SurfaceAnalytic.o

$(DIR)/SurfaceBSplineCylindrical.o: SurfaceBSplineCylindrical.cpp  SurfaceBSplineCylindrical.h DataVector.h DataMatrix.h DLeastSquaresFitting.h DPoint.h DVector.h MeshViewSet.h DMatrix.h DQuadric.h Curve2dBSpline.h DEquation.h SurfaceMulti.h DPlane.h MeshGrain.h
	${CXX} $(CXXFLAGS) -c SurfaceBSplineCylindrical.cpp -o $(DIR)/SurfaceBSplineCylindrical.o

$(DIR)/SurfaceBSplinePlanar.o: SurfaceBSplinePlanar.cpp  SurfaceBSplinePlanar.h DataVector.h DataMatrix.h DLeastSquaresFitting.h DPoint.h DVector.h MeshViewSet.h DMatrix.h DQuadric.h Curve2dBSpline.h DEquation.h SurfaceMulti.h DPlane.h MeshGrain.h
	${CXX} $(CXXFLAGS) -c SurfaceBSplinePlanar.cpp -o $(DIR)/SurfaceBSplinePlanar.o

$(DIR)/SurfaceCorrected.o: SurfaceCorrected.cpp  SurfaceCorrected.h DEquation.h
	${CXX} $(CXXFLAGS) -c SurfaceCorrected.cpp -o $(DIR)/SurfaceCorrected.o

$(DIR)/SurfaceCylinder.o: SurfaceCylinder.cpp  SurfaceCylinder.h DLeastSquaresFitting.h MeshFace.h MeshPoint3d.h
	${CXX} $(CXXFLAGS) -c SurfaceCylinder.cpp -o $(DIR)/SurfaceCylinder.o

$(DIR)/SurfaceDomain.o: SurfaceDomain.cpp  SurfaceDomain.h MeshViewSet.h
	${CXX} $(CXXFLAGS) -c SurfaceDomain.cpp -o $(DIR)/SurfaceDomain.o

$(DIR)/SurfaceDomainFaces.o: SurfaceDomainFaces.cpp  SurfaceDomainFaces.h GeometricPredicates.h MeshViewSet.h SurfaceParametric.h DTriangle.h DLine.h Metric3dContext.h DEquation.h
	${CXX} $(CXXFLAGS) -c SurfaceDomainFaces.cpp -o $(DIR)/SurfaceDomainFaces.o

$(DIR)/SurfaceDomainFacesPlanar.o: SurfaceDomainFacesPlanar.cpp  SurfaceDomainFacesPlanar.h GeometricPredicates.h MeshViewSet.h SurfaceParametric.h DTriangle.h DSegment.h DLine.h Metric3dContext.h DEquation.h
	${CXX} $(CXXFLAGS) -c SurfaceDomainFacesPlanar.cpp -o $(DIR)/SurfaceDomainFacesPlanar.o

$(DIR)/SurfaceDomainHull.o: SurfaceDomainHull.cpp  SurfaceDomainHull.h GeometricPredicates.h MeshViewSet.h SurfaceParametric.h DTriangle.h DLine.h Metric3dContext.h DEquation.h
	${CXX} $(CXXFLAGS) -c SurfaceDomainHull.cpp -o $(DIR)/SurfaceDomainHull.o

$(DIR)/SurfaceMonoAnalytic.o: SurfaceMonoAnalytic.cpp  SurfaceMonoAnalytic.h MeshLog.h
	${CXX} $(CXXFLAGS) -c SurfaceMonoAnalytic.cpp -o $(DIR)/SurfaceMonoAnalytic.o

$(DIR)/SurfaceMulti.o: SurfaceMulti.cpp  SurfaceMulti.h SurfacePlane.h
	${CXX} $(CXXFLAGS) -c SurfaceMulti.cpp -o $(DIR)/SurfaceMulti.o

$(DIR)/SurfaceParametric.o: SurfaceParametric.cpp  SurfaceParametric.h common.h DEquation.h MeshData.h Curve2dParametric.h DMetric3d.h Curve2dSegment.h DMatrix.h Metric3dContext.h ControlSpace2dAdaptive.h ControlSpace3dAdaptive.h MeshViewSet.h DLeastSquaresFitting.h DPlane.h DPlanarQuadric.h SurfacePlane.h SurfacePlanarQuadric.h SurfaceBSplinePlanar.h SurfaceDomain.h SurfaceCylinder.h SurfaceSphere.h MeshGenerator3dSurface.h MeshPoint3d.h MeshFace.h DataStatistics.h
	${CXX} $(CXXFLAGS) -c SurfaceParametric.cpp -o $(DIR)/SurfaceParametric.o

$(DIR)/SurfacePlanarQuadric.o: SurfacePlanarQuadric.cpp  SurfacePlanarQuadric.h DEquation.h Curve2dParametric.h DMatrix.h DPlane.h
	${CXX} $(CXXFLAGS) -c SurfacePlanarQuadric.cpp -o $(DIR)/SurfacePlanarQuadric.o

$(DIR)/SurfacePlane.o: SurfacePlane.cpp  SurfacePlane.h DEquation.h Curve2dParametric.h DMatrix.h DPlane.h DataVector.h MeshFace.h MeshPoint3d.h
	${CXX} $(CXXFLAGS) -c SurfacePlane.cpp -o $(DIR)/SurfacePlane.o

$(DIR)/SurfaceQuadricQTree.o: SurfaceQuadricQTree.cpp  SurfaceQuadricQTree.h DPlanarQuadric.h MeshFace.h DLeastSquaresFitting.h MeshContainer3dSurface.h GeometricPredicates.h Metric3dContext.h MeshViewSet.h Metric2dContext.h ControlSpace2dIdentity.h ControlSpace2dProjected.h
	${CXX} $(CXXFLAGS) -c SurfaceQuadricQTree.cpp -o $(DIR)/SurfaceQuadricQTree.o

$(DIR)/SurfaceSphere.o: SurfaceSphere.cpp  SurfaceSphere.h DLeastSquaresFitting.h DataVector.h MeshPoint3d.h
	${CXX} $(CXXFLAGS) -c SurfaceSphere.cpp -o $(DIR)/SurfaceSphere.o

$(DIR)/tinystr.o: tinystr.cpp  tinystr.h
	${CXX} $(CXXFLAGS) -c tinystr.cpp -o $(DIR)/tinystr.o

$(DIR)/tinyxml.o: tinyxml.cpp  tinyxml.h
	${CXX} $(CXXFLAGS) -c tinyxml.cpp -o $(DIR)/tinyxml.o

$(DIR)/tinyxmlerror.o: tinyxmlerror.cpp  tinyxml.h
	${CXX} $(CXXFLAGS) -c tinyxmlerror.cpp -o $(DIR)/tinyxmlerror.o

$(DIR)/tinyxmlparser.o: tinyxmlparser.cpp  tinyxml.h
	${CXX} $(CXXFLAGS) -c tinyxmlparser.cpp -o $(DIR)/tinyxmlparser.o

