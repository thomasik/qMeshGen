#include <QtGui>
#include <QtOpenGL>

#undef UNICODE
#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
using namespace log4cplus;
#define UNICODE

#include <math.h>

#include "ModelGlWidget.h"

#include "MeshViewSet.h"
#include "MeshContainer3d.h"
#include "MeshDomainVolume.h"
#include "MeshTetrahedron.h"
#include "MeshEdge3d.h"

ModelGLWidget::ModelGLWidget(QWidget *parent)
    : QGLWidget(parent), m_view_set(nullptr), m_view_mesh(nullptr), 
		m_initialized(false), 
		clip_quality(1.1), view_mode(VIEW_BLOCKS | VIEW_FACES | VIEW_EDGES | VIEW_NODES), 
		nav_mode(0),
//		xRot(0.0f), yRot(0.0f), xTrans(0.0f), yTrans(0.0f), 
		last_zoom(1.0f), total_zoom(1.0f)
{
	setFormat(QGLFormat(QGL::DoubleBuffer | QGL::DepthBuffer));

    QAction* viewPrefAct = new QAction(tr("&Preferences"), this);
    viewPrefAct->setShortcut(tr("Ctrl+P"));
    viewPrefAct->setStatusTip(tr("Mesh view preferences"));
    connect(viewPrefAct, SIGNAL(triggered()), this, SLOT(showPrefDlg()));
	addAction(viewPrefAct);

    QAction* resetAct = new QAction(tr("&Reset"), this);
    resetAct->setShortcut(tr("Ctrl+R"));
    resetAct->setStatusTip(tr("Reset positioning"));
	connect(resetAct, SIGNAL(triggered()), this, SLOT(resetPositioning()));
	addAction(resetAct);

    QAction* clearAct = new QAction(tr("&Clear"), this);
    clearAct->setShortcut(tr("Ctrl+L"));
    clearAct->setStatusTip(tr("Clear view"));
	connect(clearAct, SIGNAL(triggered()), this, SLOT(clearView()));
	addAction(clearAct);

    QAction* storeMFAct = new QAction(tr("Store &Matlab"), this);
    storeMFAct->setShortcut(tr("Ctrl+M"));
    storeMFAct->setStatusTip(tr("Store matlab file for graphics export"));
	connect(storeMFAct, SIGNAL(triggered()), this, SLOT(storeMatlabFile()));
	addAction(storeMFAct);

    QAction* storeEPSAct = new QAction(tr("Store &EPS"), this);
    storeEPSAct->setShortcut(tr("Ctrl+E"));
    storeEPSAct->setStatusTip(tr("Store eps file for graphics export"));
	connect(storeEPSAct, SIGNAL(triggered()), this, SLOT(storeEPSFile()));
	addAction(storeEPSAct);

	setContextMenuPolicy(Qt::ActionsContextMenu);

	viewPrefDlg = new ViewPrefDialog(this, this);

}

ModelGLWidget::~ModelGLWidget()
{
    if(m_view_set) 
		delete m_view_set;
}

QSize ModelGLWidget::minimumSizeHint() const
{
    return QSize(100, 100);
}

QSize ModelGLWidget::sizeHint() const
{
    return QSize(500, 500);
}

	//"uniform lowp vec4 color_white  = vec4(1.0, 1.0, 1.0, 1.0);"
	//"uniform lowp vec4 color_black  = vec4(0.0, 0.0, 0.0, 1.0);"
	//"uniform lowp vec4 color_dpurple = vec4(0.6, 0.0, 0.6, 1.0);"
	//"uniform lowp vec4 color_dblue   = vec4(0.0, 0.0, 0.6, 1.0);"
	//"uniform lowp vec4 color_dcyan   = vec4(0.0, 0.6, 0.6, 1.0);"
	//"uniform lowp vec4 color_dgreen  = vec4(0.0, 0.6, 0.0, 1.0);"
	//"uniform lowp vec4 color_dyellow = vec4(0.6, 0.6, 0.0, 1.0);"
	//"uniform lowp vec4 color_dred    = vec4(0.6, 0.0, 0.0, 1.0);"
	//"uniform lowp vec4 color_purple = vec4(0.9, 0.0, 0.9, 1.0);"
	//"uniform lowp vec4 color_blue   = vec4(0.0, 0.0, 0.9, 1.0);"
	//"uniform lowp vec4 color_cyan   = vec4(0.0, 0.9, 0.9, 1.0);"
	//"uniform lowp vec4 color_green  = vec4(0.0, 0.9, 0.0, 1.0);"
	//"uniform lowp vec4 color_yellow = vec4(0.9, 0.9, 0.0, 1.0);"
	//"uniform lowp vec4 color_red    = vec4(0.9, 0.0, 0.0, 1.0);"

void ModelGLWidget::initializeGL()
{
	if(m_initialized) return;

//	qglClearColor(Qt::black);

	glDepthFunc( GL_LEQUAL );
	glEnable( GL_DEPTH_TEST );

	glPointSize( 3.0 );
	glLineWidth( 1.0 );

	m_bf_render.init(this);
	m_ep_render.init(this);
	m_label_render.init(this);

	m_initialized = true;
}

void ModelGLWidget::paintGL()
{
	if(view_mode & VIEW_WHITE){
		qglClearColor(Qt::white);
	}else{
		qglClearColor(Qt::black);
	}

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	draw();
}

void ModelGLWidget::resizeGL(int w, int h)
{
	//LOG4CPLUS_INFO(MeshLog::logger_mesh, "ModelGLWidget::resizeGL(" << w << " x " << h  << ")");
	glViewport (0, 0, w, h);

	if(m_bounding_box.valid){
		double dl = 0.4 * m_bounding_box.getDiameter();
		//double dlx(0.55 * m_bounding_box.getDX());
		//double dly(0.55 * m_bounding_box.getDY());
		double dlx = dl;
		double dly = dl;
		//double dlz( 100 * m_bounding_box.getDZ());
		double dlz = 100 * dl;
		if (w <= h)
			dly *= h/(GLfloat)w;
		else
			dlx *= w/(GLfloat)h;

		if(dlz < dlx) dlz = dlx;
		if(dlz < dly) dlz = dly;

		m_proj_matrix.setToIdentity();
		m_proj_matrix.ortho(-dlx, dlx, -dly, dly, -dlz, dlz);

		axis_ox = QVector3D( 2*dlx / w, 0.0, 0.0);
		axis_oy = QVector3D( 0.0, 2*dly / h, 0.0);

		//m_proj_matrix.frustum + perspective ?
	}
}

void ModelGLWidget::dumpMatrixInfo()
{
	QMatrix4x4 tm;
	tm.rotate(total_rot * last_rot);
	QMatrix3x3 nm = tm.normalMatrix();
	tm.scale(total_zoom * last_zoom);
	tm.translate(total_trans + last_trans);

	QMatrix4x4 modelview = m_proj_matrix * tm;

	//	m_view_set->storeEPSFile(trans_matrix, mvmatrix, view_mode, 
	// m_view_set->storeEPSFile(modelview.data(), tm.data(), view_mode,

//#undef UNICODE
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "*** QtMeshGen *** modelview:");
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "		float mv_m[] = {");
	for (int i = 0; i < 4; i++) {
		ostringstream line;
		for (int j = 0; j < 4; j++)
			line << modelview(i, j) << ((j < 3 || i < 3) ? ", " : " };");
		LOG4CPLUS_INFO(MeshLog::logger_mesh, line.str());
	}
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "*** QtMeshGen *** tm:");
	LOG4CPLUS_INFO(MeshLog::logger_mesh, "		float tm_m[] = {");
	for (int i = 0; i < 4; i++) {
		ostringstream line;
		for (int j = 0; j < 4; j++)
			line  << tm(i, j) << ((j < 3 || i < 3) ? ", " : " };");
		LOG4CPLUS_INFO(MeshLog::logger_mesh, line.str());
	}
//#define UNICODE

	if (viewPrefDlg->clipPlanes.notEmpty()) {
		for (int i = 1; i < viewPrefDlg->clipPlanes.countInt(); i++) {
			auto cp = viewPrefDlg->clipPlanes[i];
			auto nv = cp->getNormalVector();
			auto line = tr("MeshViewSet::ClipPlane cp(FVector3d(%1,%2,%3), %4);").arg(nv.x).arg(nv.y).arg(nv.z).arg(cp->getParamD());
			LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, line.toStdString());
		}
	}
}

void ModelGLWidget::mousePressEvent(QMouseEvent *event)
{
	lastPos = event->pos();
	if(event->buttons() & Qt::LeftButton){
		if(event->modifiers() & Qt::ShiftModifier)
			if(viewPrefDlg->selected_clip_plane > 0)
				nav_mode = 4; // translate clip plane
			else
				nav_mode = 1; // translate view set
		else if(event->modifiers() & Qt::ControlModifier)
			nav_mode = 2; // scale
		else{
			if(viewPrefDlg->selected_clip_plane > 0)
				nav_mode = 5; // rotate clip plane
			else
				nav_mode = 3; // rotate view set
		}
	}
	else if (event->buttons() & Qt::RightButton) {
		dumpMatrixInfo();
	}
}

void ModelGLWidget::mouseReleaseEvent(QMouseEvent *  /*event*/)
{
	if(nav_mode == 1){ // translate view set
		total_trans += last_trans;
		last_trans = QVector3D();
		updateGL();
	}else if(nav_mode == 2){ // scale view set
		total_zoom *= last_zoom;
		last_zoom = 1.0f;
		updateGL();
	}else if(nav_mode == 3){ // rotate view set
		total_rot *= last_rot;
		last_rot = QQuaternion();
		axis_rot = last_axis_rot * axis_rot;
		last_axis_rot = QQuaternion();
		updateGL();
	}else if(nav_mode == 4){ // translate clip plane
		prepareGLScene();
		updateGL();
	}else if(nav_mode == 5){ // rotate clip plane
		prepareGLScene();
		updateGL();
	}

	nav_mode = 0;
}

void ModelGLWidget::mouseMoveEvent(QMouseEvent *event)
{
	if(nav_mode == 1){ // translate view set
		int dx = event->x() - lastPos.x();
		int dy = lastPos.y() - event->y(); // OpenGL renders with (0,0) on bottom, mouse reports with (0,0) on top

		QVector3D axis_ox_rot = axis_rot.rotatedVector(axis_ox);
		QVector3D axis_oy_rot = axis_rot.rotatedVector(axis_oy);

		last_trans = ( axis_ox_rot * dx + axis_oy_rot * dy ) / total_zoom;

		updateGL();
	}else if(nav_mode == 2){ // scale view set
		QPoint center(width()/2, height()/2);
		float start_dist = sqrtf(sqr(center.x() - lastPos.x()) + sqr(center.y() - lastPos.y()));
		float curr_dist  = sqrtf(sqr(center.x() - event->x())  + sqr(center.y() - event->y()));
		if(start_dist > 0.0){
			last_zoom = curr_dist / start_dist;
			updateGL();
		}
	}else if(nav_mode == 3){ // rotate view set
		GLfloat dx = GLfloat(event->x() - lastPos.x()) / width();
		GLfloat dy = GLfloat(event->y() - lastPos.y()) / height();

		QVector3D axis_ox_rot = axis_rot.rotatedVector(axis_ox);
		QVector3D axis_oy_rot = axis_rot.rotatedVector(axis_oy);

		last_rot = QQuaternion::fromAxisAndAngle(axis_ox_rot, 180.0 * dy / total_zoom) 
				 * QQuaternion::fromAxisAndAngle(axis_oy_rot, 180.0 * dx / total_zoom);

		last_axis_rot = QQuaternion::fromAxisAndAngle(axis_oy_rot, -180.0 * dx / total_zoom)
						* QQuaternion::fromAxisAndAngle(axis_ox_rot, -180.0 * dy / total_zoom) ;


		updateGL();
	}else if(nav_mode == 4){ // translate clip plane
		int dx = event->x() - lastPos.x();
		int dy = lastPos.y() - event->y(); // OpenGL renders with (0,0) on bottom, mouse reports with (0,0) on top

		QVector3D axis_ox_rot = axis_rot.rotatedVector(axis_ox);
		QVector3D axis_oy_rot = axis_rot.rotatedVector(axis_oy);

		QVector3D trans = ( axis_ox_rot * dx + axis_oy_rot * dy ) / total_zoom;

		viewPrefDlg->clipPlanes[viewPrefDlg->selected_clip_plane]->translate(DVector3d( trans.x(), trans.y(), trans.z() ));
		viewPrefDlg->clipPlanes[viewPrefDlg->selected_clip_plane]->setTableRow(
			viewPrefDlg->clipTable, viewPrefDlg->selected_clip_plane);

		lastPos = event->pos();

		updateGL();
	}else if(nav_mode == 5){ // rotate clip plane
		GLfloat dx = GLfloat(event->x() - lastPos.x()) / width();
		GLfloat dy = GLfloat(event->y() - lastPos.y()) / height();

		QVector3D axis_ox_rot = axis_rot.rotatedVector(axis_ox);
		QVector3D axis_oy_rot = axis_rot.rotatedVector(axis_oy);

		QQuaternion rot = QQuaternion::fromAxisAndAngle(axis_ox_rot, 180.0 * dy / total_zoom) 
				 * QQuaternion::fromAxisAndAngle(axis_oy_rot, 180.0 * dx / total_zoom);


		DVector3d vn = viewPrefDlg->clipPlanes[viewPrefDlg->selected_clip_plane]->getNormalVector();

		QVector3D rot_vn = rot.rotatedVector( QVector3D( vn.x, vn.y, vn.z ) );

		viewPrefDlg->clipPlanes[viewPrefDlg->selected_clip_plane]->setNormalVector( DVector3d( rot_vn.x(), rot_vn.y(), rot_vn.z() ) );
		viewPrefDlg->clipPlanes[viewPrefDlg->selected_clip_plane]->setTableRow(
			viewPrefDlg->clipTable, viewPrefDlg->selected_clip_plane);

		lastPos = event->pos();

		updateGL();
	}
}

void ModelGLWidget::wheelEvent(QWheelEvent *event)
{
	static const float ZOOM_RATIO = 1.1f;

	if(nav_mode == 0){
		GLfloat z = (event->delta() > 0.0) ? ZOOM_RATIO : 1/ZOOM_RATIO;
		total_zoom *= z;

		updateGL();
	}
}

bool ModelGLWidget::BlockFaceRender::init(QObject * parent)
{
	if(program != nullptr) return true; // alredy initialized?

	QOpenGLShader *vshader = new QOpenGLShader(QOpenGLShader::Vertex, parent);
	bool vok = vshader->compileSourceCode(
		"attribute highp vec4 vertex;"
		"attribute mediump vec3 normal;"
		"attribute lowp float colorIndex;"
		"uniform mediump mat4 mvMatrix;"
		"uniform mediump mat3 normalMatrix;"
		"varying mediump vec4 color;"
		"void main(void)"
		"{"
		"    vec3 trNormal = normalMatrix * normal;"
		"	 vec3 toLight  = normalize(vec3(0.0, -0.5, -1.0));"
		"    float angle = max(dot(trNormal, toLight), 0.0);"
		"	 if( trNormal.z <= 0.0 ) color = vec4(0.3, 0.3, 0.3, 1.0);"
		"    else if( colorIndex > 250.0 ) color = vec4(0.3, 0.3, 0.3, 1.0);"
		"    else if( colorIndex > 239.0 ) color = mix( vec4(1.0, 1.0, 1.0, 1.0), vec4(0.0, 0.0, 0.0, 1.0), (250.0-colorIndex) * 0.1);"
		"    else if( colorIndex > 220.0 ) color = mix( vec4(0.9, 0.0, 0.9, 1.0), vec4(0.0, 0.0, 0.9, 1.0), (240.0-colorIndex) * 0.05);"
		"    else if( colorIndex > 200.0 ) color = mix( vec4(0.0, 0.0, 0.9, 1.0), vec4(0.0, 0.9, 0.9, 1.0), (220.0-colorIndex) * 0.05);"
		"    else if( colorIndex > 180.0 ) color = mix( vec4(0.0, 0.9, 0.9, 1.0), vec4(0.0, 0.9, 0.0, 1.0), (200.0-colorIndex) * 0.05);"
		"    else if( colorIndex > 160.0 ) color = mix( vec4(0.0, 0.9, 0.0, 1.0), vec4(0.9, 0.9, 0.0, 1.0), (180.0-colorIndex) * 0.05);"
		"    else if( colorIndex > 120.0 ) color = mix( vec4(0.9, 0.9, 0.0, 1.0), vec4(0.9, 0.0, 0.0, 1.0), (160.0-colorIndex) * 0.025);"
		"    else if( colorIndex > 100.0 ) color = mix( vec4(0.6, 0.0, 0.6, 1.0), vec4(0.0, 0.0, 0.6, 1.0), (120.0-colorIndex) * 0.05);"
		"    else if( colorIndex >  80.0 ) color = mix( vec4(0.0, 0.0, 0.6, 1.0), vec4(0.0, 0.6, 0.6, 1.0), (100.0-colorIndex) * 0.05);"
		"    else if( colorIndex >  60.0 ) color = mix( vec4(0.0, 0.6, 0.6, 1.0), vec4(0.0, 0.6, 0.0, 1.0), ( 80.0-colorIndex) * 0.05);"
		"    else if( colorIndex >  40.0 ) color = mix( vec4(0.0, 0.6, 0.0, 1.0), vec4(0.6, 0.6, 0.0, 1.0), ( 60.0-colorIndex) * 0.05);"
		"    else color = mix( vec4(0.6, 0.6, 0.0, 1.0), vec4(0.6, 0.0, 0.0, 1.0), (40.0-colorIndex) * 0.025);"
		"    color.xyz = color.xyz * 0.2 + color.xyz * 0.8 * angle;"
		"    color = sqrt(clamp(color, 0.0, 1.0));"
		"    gl_Position = mvMatrix * vertex;"
		"}");
	if(!vok) LOG4CPLUS_WARN(MeshLog::logger_console, 
		"block/face vshader compile: " << vshader->log().toStdString());

	QOpenGLShader *fshader = new QOpenGLShader(QOpenGLShader::Fragment, parent);
	bool fok = fshader->compileSourceCode(
		"varying mediump vec4 color;"
		"void main(void)"
		"{"
		"    gl_FragColor = color;"
		"}");
	if(!fok) LOG4CPLUS_WARN(MeshLog::logger_console, 
		"block/face fshader compile: " << fshader->log().toStdString());

	if(!vok || !fok) return false;

	program = new QOpenGLShaderProgram(parent);
	program->addShader(vshader);
	program->addShader(fshader);
	bool pok = program->link();
	if(!pok){
		LOG4CPLUS_WARN(MeshLog::logger_console, 
			"block/face gl program link: " << program->log().toStdString());
		return false;
	}

	vertexAttr = program->attributeLocation("vertex");
	normalAttr = program->attributeLocation("normal");
	colorAttr  = program->attributeLocation("colorIndex");
	matrixUniform = program->uniformLocation("mvMatrix");
	normalMatrixUniform = program->uniformLocation("normalMatrix");
	
	return true;
}

void ModelGLWidget::BlockFaceRender::renderBlocks(const QMatrix4x4& modelview, const QMatrix3x3& nm) const
{
	if(program == nullptr) return;
	if(block_vertices.isEmpty()) return;

	glEnable(GL_CULL_FACE);

	program->bind();
	program->setUniformValue(matrixUniform, modelview);
	program->setUniformValue(normalMatrixUniform, nm);

	program->enableAttributeArray(normalAttr);
	program->enableAttributeArray(vertexAttr);
	program->enableAttributeArray(colorAttr);

	program->setAttributeArray(vertexAttr, block_vertices.constData());
	program->setAttributeArray(normalAttr, block_normals.constData());
	program->setAttributeArray(colorAttr, GL_FLOAT, block_colors.constData(), 1);

	glDrawArrays(GL_TRIANGLES, 0, block_vertices.size());

	program->disableAttributeArray(normalAttr);
	program->disableAttributeArray(vertexAttr);
	program->disableAttributeArray(colorAttr);

	program->release();

	glDisable(GL_CULL_FACE);
}

void ModelGLWidget::BlockFaceRender::renderFaces(const QMatrix4x4& modelview, const QMatrix3x3& nm) const
{
	if(program == nullptr) return;
	if(face_vertices.isEmpty()) return;

	program->bind();
	program->setUniformValue(matrixUniform, modelview);
	program->setUniformValue(normalMatrixUniform, nm);

	program->enableAttributeArray(normalAttr);
	program->enableAttributeArray(vertexAttr);
	program->enableAttributeArray(colorAttr);

	program->setAttributeArray(vertexAttr, face_vertices.constData());
	program->setAttributeArray(normalAttr, face_normals.constData());
	program->setAttributeArray(colorAttr, GL_FLOAT, face_colors.constData(), 1);

	glDrawArrays(GL_TRIANGLES, 0, face_vertices.size());

	program->disableAttributeArray(normalAttr);
	program->disableAttributeArray(vertexAttr);
	program->disableAttributeArray(colorAttr);

	program->release();
}

bool ModelGLWidget::EdgePointRender::init(QObject * parent)
{
	if( program != nullptr ) return true;

	QOpenGLShader *vshader = new QOpenGLShader(QOpenGLShader::Vertex, parent);
	bool vok = vshader->compileSourceCode(
		"attribute highp vec4 vertex;"
		"attribute lowp float colorIndex;"
		"uniform mediump mat4 mvMatrix;"
		"varying mediump vec4 color;"
		"void main(void)"
		"{"
		"    if( colorIndex > 250.0 ) color = vec4(0.3, 0.3, 0.3, 1.0);"
		"    else if( colorIndex > 239.0 ) color = mix( vec4(1.0, 1.0, 1.0, 1.0), vec4(0.0, 0.0, 0.0, 1.0), (250.0-colorIndex) * 0.1);"
		"    else if( colorIndex > 220.0 ) color = mix( vec4(0.9, 0.0, 0.9, 1.0), vec4(0.0, 0.0, 0.9, 1.0), (240.0-colorIndex) * 0.05);"
		"    else if( colorIndex > 200.0 ) color = mix( vec4(0.0, 0.0, 0.9, 1.0), vec4(0.0, 0.9, 0.9, 1.0), (220.0-colorIndex) * 0.05);"
		"    else if( colorIndex > 180.0 ) color = mix( vec4(0.0, 0.9, 0.9, 1.0), vec4(0.0, 0.9, 0.0, 1.0), (200.0-colorIndex) * 0.05);"
		"    else if( colorIndex > 160.0 ) color = mix( vec4(0.0, 0.9, 0.0, 1.0), vec4(0.9, 0.9, 0.0, 1.0), (180.0-colorIndex) * 0.05);"
		"    else if( colorIndex > 120.0 ) color = mix( vec4(0.9, 0.9, 0.0, 1.0), vec4(0.9, 0.0, 0.0, 1.0), (160.0-colorIndex) * 0.025);"
		"    else if( colorIndex > 100.0 ) color = mix( vec4(0.6, 0.0, 0.6, 1.0), vec4(0.0, 0.0, 0.6, 1.0), (120.0-colorIndex) * 0.05);"
		"    else if( colorIndex >  80.0 ) color = mix( vec4(0.0, 0.0, 0.6, 1.0), vec4(0.0, 0.6, 0.6, 1.0), (100.0-colorIndex) * 0.05);"
		"    else if( colorIndex >  60.0 ) color = mix( vec4(0.0, 0.6, 0.6, 1.0), vec4(0.0, 0.6, 0.0, 1.0), ( 80.0-colorIndex) * 0.05);"
		"    else if( colorIndex >  40.0 ) color = mix( vec4(0.0, 0.6, 0.0, 1.0), vec4(0.6, 0.6, 0.0, 1.0), ( 60.0-colorIndex) * 0.05);"
		"    else color = mix( vec4(0.6, 0.6, 0.0, 1.0), vec4(0.6, 0.0, 0.0, 1.0), (40.0-colorIndex) * 0.025);"
		"    gl_Position = mvMatrix * vertex;"
		"}");
	if(!vok) LOG4CPLUS_WARN(MeshLog::logger_console, 
		"edge vshader compile: " << vshader->log().toStdString());

	QOpenGLShader *fshader = new QOpenGLShader(QOpenGLShader::Fragment, parent);
	bool fok = fshader->compileSourceCode(
		"varying mediump vec4 color;"
		"void main(void)"
		"{"
		"    gl_FragColor = color;"
		"}");
	if(!fok) LOG4CPLUS_WARN(MeshLog::logger_console, 
		"edge fshader compile: " << fshader->log().toStdString());

	if(!vok || !fok) return false;

	program = new QOpenGLShaderProgram(parent);
	program->addShader(vshader);
	program->addShader(fshader);
	bool pok = program->link();
	if(!pok){
		LOG4CPLUS_WARN(MeshLog::logger_console, 
			"edge gl program link: " << program->log().toStdString());
		return false;
	}

	vertexAttr = program->attributeLocation("vertex");
	colorAttr  = program->attributeLocation("colorIndex");
	matrixUniform = program->uniformLocation("mvMatrix");
	
	return true;
}

void ModelGLWidget::EdgePointRender::renderEdges(const QMatrix4x4& modelview) const
{
	if(program == nullptr) return;
	if(edge_vertices.isEmpty()) return;

	program->bind();
	program->setUniformValue(matrixUniform, modelview);

	program->enableAttributeArray(vertexAttr);
	program->enableAttributeArray(colorAttr);

	program->setAttributeArray(vertexAttr, edge_vertices.constData());
	program->setAttributeArray(colorAttr, GL_FLOAT, edge_colors.constData(), 1);
	glDrawArrays(GL_LINES, 0, edge_vertices.size());

	program->disableAttributeArray(vertexAttr);
	program->disableAttributeArray(colorAttr);

	program->release();
}

QVector3D ModelGLWidget::EdgePointRender::drawNumber(
	int number, const QVector3D & pos, GLfloat color, 
	const QVector3D& svx, const QVector3D& svy, 
	QVector<QVector3D> & nvertices, QVector<GLfloat> & ncolors) const
{
	QVector3D npos = pos;

	if( number < 0 ) {
		// draw '-'
		QVector3D mpos = pos + svy * (N_CY/2);
		nvertices << mpos << (mpos + svx * N_CX);
		ncolors << color << color;
		npos += svx * (N_CX + N_CS);
		number = -number;
	}

	if( number > 9 ) {
		npos = drawNumber(number / 10, npos, color, svx, svy, nvertices, ncolors);
		number %= 10;
	}

	QVector3D vx = svx * N_CX;
	QVector3D vy = svy * N_CY;
	QVector3D vmy = svy * (N_CY/2);
	// draw number ...
	if(number == 0){
		QVector3D nposx = npos + vx;
		QVector3D nposy = npos + vy;
		QVector3D nposxy = nposx + vy;
		nvertices << npos << nposx << nposy << nposxy;
		ncolors << color << color << color << color;
		nvertices << npos << nposy << nposx << nposxy;
		ncolors << color << color << color << color;
	}else if(number == 1){
		QVector3D nposy = npos + vy;
		nvertices << npos << nposy;
		ncolors << color << color;
		return npos + (svx * (N_CS+1));
	}else if(number == 2){
		QVector3D nposx = npos + vx;
		QVector3D nposy = npos + vy;
		QVector3D nposxy = nposx + vy;
		QVector3D nposmy = npos + vmy;
		QVector3D nposxmy = nposx + vmy;
		nvertices << nposy << nposxy << nposxy << nposxmy;
		ncolors << color << color << color << color;
		nvertices << nposxmy << nposmy << nposmy << npos << npos << nposx;
		ncolors << color << color << color << color << color << color;
	}else if(number == 3){
		QVector3D nposx = npos + vx;
		QVector3D nposy = npos + vy;
		QVector3D nposxy = nposx + vy;
		QVector3D nposmy = npos + vmy;
		QVector3D nposxmy = nposx + vmy;
		nvertices << nposy << nposxy << nposxy << nposx;
		ncolors << color << color << color << color;
		nvertices << nposx << npos << nposmy << nposxmy;
		ncolors << color << color << color << color;
	}else if(number == 4){
		QVector3D nposx = npos + vx;
		QVector3D nposy = npos + vy;
		QVector3D nposmy = npos + vmy;
		QVector3D nposxmy = nposx + vmy;
		nvertices << nposy << nposmy << nposmy << nposxmy << nposxmy << nposx;
		ncolors << color << color << color << color << color << color;
	}else if(number == 5){
		QVector3D nposx = npos + vx;
		QVector3D nposy = npos + vy;
		QVector3D nposxy = nposx + vy;
		QVector3D nposmy = npos + vmy;
		QVector3D nposxmy = nposx + vmy;
		nvertices << nposxy << nposy << nposy << nposmy;
		ncolors << color << color << color << color;
		nvertices << nposmy << nposxmy << nposxmy << nposx << nposx << npos;
		ncolors << color << color << color << color << color << color;
	}else if(number == 6){
		QVector3D nposx = npos + vx;
		QVector3D nposy = npos + vy;
		QVector3D nposxy = nposx + vy;
		QVector3D nposmy = npos + vmy;
		QVector3D nposxmy = nposx + vmy;
		nvertices << nposxy << nposy << nposy << npos;
		ncolors << color << color << color << color;
		nvertices << npos << nposx << nposx << nposxmy << nposxmy << nposmy;
		ncolors << color << color << color << color << color << color;
	}else if(number == 7){
		QVector3D nposx = npos + vx;
		QVector3D nposy = npos + vy;
		QVector3D nposxy = nposx + vy;
		nvertices << nposy << nposxy << nposxy << nposx;
		ncolors << color << color << color << color;
	}else if(number == 8){
		QVector3D nposx = npos + vx;
		QVector3D nposy = npos + vy;
		QVector3D nposxy = nposx + vy;
		QVector3D nposmy = npos + vmy;
		QVector3D nposxmy = nposx + vmy;
		nvertices << nposxy << nposy << nposy << npos;
		ncolors << color << color << color << color;
		nvertices << npos << nposx << nposx << nposxy << nposxmy << nposmy;
		ncolors << color << color << color << color << color << color;
	}else if(number == 9){
		QVector3D nposx = npos + vx;
		QVector3D nposy = npos + vy;
		QVector3D nposxy = nposx + vy;
		QVector3D nposmy = npos + vmy;
		QVector3D nposxmy = nposx + vmy;
		nvertices << npos << nposx << nposx << nposxy;
		ncolors << color << color << color << color;
		nvertices << nposxy << nposy << nposy << nposmy << nposmy << nposxmy;
		ncolors << color << color << color << color << color << color;
	}
	return npos + (svx * (N_CX+N_CS));
}

void ModelGLWidget::EdgePointRender::renderPoints(const QMatrix4x4& modelview, const QVector3D& svx, const QVector3D& svy) const
{
	if(program == nullptr) return;
	if(point_vertices.isEmpty()) return;

	program->bind();
	program->setUniformValue(matrixUniform, modelview);

	program->enableAttributeArray(vertexAttr);
	program->enableAttributeArray(colorAttr);

	program->setAttributeArray(vertexAttr, point_vertices.constData());
	program->setAttributeArray(colorAttr, GL_FLOAT, point_colors.constData(), 1);
	glDrawArrays(GL_POINTS, 0, point_vertices.size());

	if(any_number){
		QVector<QVector3D> number_vertices;
		QVector<GLfloat>   number_colors;
		for(int i = 0; i < point_vertices.size(); i++){
			if( point_ids[i] == EMPTY_ID) continue;
			drawNumber(point_ids[i], point_vertices[i] + svx*2 + svy, point_colors[i], svx, svy, number_vertices, number_colors);
		}
		program->setAttributeArray(vertexAttr, number_vertices.constData());
		program->setAttributeArray(colorAttr, GL_FLOAT, number_colors.constData(), 1);
		glDrawArrays(GL_LINES, 0, number_vertices.size());
	}

	program->disableAttributeArray(vertexAttr);
	program->disableAttributeArray(colorAttr);

	program->release();
}

void ModelGLWidget::EdgePointRender::renderClipPlane(const QMatrix4x4& modelview, const ClipPlane* plane, double diameter) const
{
	if(program == nullptr) return;

	QVector<QVector3D> clip_vertices;
	QVector<GLfloat> clip_colors;
	plane->drawGL(clip_vertices, clip_colors, diameter);

	if(clip_vertices.isEmpty()) return;

	program->bind();
	program->setUniformValue(matrixUniform, modelview);

	program->enableAttributeArray(vertexAttr);
	program->enableAttributeArray(colorAttr);

	program->setAttributeArray(vertexAttr, clip_vertices.constData());
	program->setAttributeArray(colorAttr, GL_FLOAT, clip_colors.constData(), 1);
	glDrawArrays(GL_LINES, 0, clip_vertices.size());

	program->disableAttributeArray(vertexAttr);
	program->disableAttributeArray(colorAttr);

	program->release();
}

bool ModelGLWidget::LabelRender::init(QObject * parent)
{
	if( program != nullptr ) return true;

	if( FT_Init_FreeType( &ft )) {
		LOG4CPLUS_WARN(MeshLog::logger_console, "error initing freetype library");
		return false;
	}
	if( FT_New_Face( ft, "FreeSans.ttf", 0, &face) ) {
		LOG4CPLUS_WARN(MeshLog::logger_console, "error opening font FreeSans.ttf");
		return false;
	}

	FT_Set_Pixel_Sizes( face, 0, 14 );

	QOpenGLShader *vshader = new QOpenGLShader(QOpenGLShader::Vertex, parent);
	bool vok = vshader->compileSourceCode(
		"attribute highp vec4 vertex;"
		"attribute mediump vec2 texCoord;"
		"attribute lowp float colorIndex;"
		"uniform mediump mat4 mvMatrix;"
		"varying mediump vec4 color;"
		"varying mediump vec2 texc;"
		"void main(void)"
		"{"
		"    if( colorIndex > 250.0 ) color = vec4(0.3, 0.3, 0.3, 1.0);"
		"    else if( colorIndex > 239.0 ) color = mix( vec4(1.0, 1.0, 1.0, 1.0), vec4(0.0, 0.0, 0.0, 1.0), (250.0-colorIndex) * 0.1);"
		"    else if( colorIndex > 220.0 ) color = mix( vec4(0.9, 0.0, 0.9, 1.0), vec4(0.0, 0.0, 0.9, 1.0), (240.0-colorIndex) * 0.05);"
		"    else if( colorIndex > 200.0 ) color = mix( vec4(0.0, 0.0, 0.9, 1.0), vec4(0.0, 0.9, 0.9, 1.0), (220.0-colorIndex) * 0.05);"
		"    else if( colorIndex > 180.0 ) color = mix( vec4(0.0, 0.9, 0.9, 1.0), vec4(0.0, 0.9, 0.0, 1.0), (200.0-colorIndex) * 0.05);"
		"    else if( colorIndex > 160.0 ) color = mix( vec4(0.0, 0.9, 0.0, 1.0), vec4(0.9, 0.9, 0.0, 1.0), (180.0-colorIndex) * 0.05);"
		"    else if( colorIndex > 120.0 ) color = mix( vec4(0.9, 0.9, 0.0, 1.0), vec4(0.9, 0.0, 0.0, 1.0), (160.0-colorIndex) * 0.025);"
		"    else if( colorIndex > 100.0 ) color = mix( vec4(0.6, 0.0, 0.6, 1.0), vec4(0.0, 0.0, 0.6, 1.0), (120.0-colorIndex) * 0.05);"
		"    else if( colorIndex >  80.0 ) color = mix( vec4(0.0, 0.0, 0.6, 1.0), vec4(0.0, 0.6, 0.6, 1.0), (100.0-colorIndex) * 0.05);"
		"    else if( colorIndex >  60.0 ) color = mix( vec4(0.0, 0.6, 0.6, 1.0), vec4(0.0, 0.6, 0.0, 1.0), ( 80.0-colorIndex) * 0.05);"
		"    else if( colorIndex >  40.0 ) color = mix( vec4(0.0, 0.6, 0.0, 1.0), vec4(0.6, 0.6, 0.0, 1.0), ( 60.0-colorIndex) * 0.05);"
		"    else color = mix( vec4(0.6, 0.6, 0.0, 1.0), vec4(0.6, 0.0, 0.0, 1.0), (40.0-colorIndex) * 0.025);"
		"	 texc = texCoord;"
		"    gl_Position = mvMatrix * vertex;"
		"}");
	if(!vok) LOG4CPLUS_WARN(MeshLog::logger_console, 
		"label vshader compile: " << vshader->log().toStdString());

	QOpenGLShader *fshader = new QOpenGLShader(QOpenGLShader::Fragment, parent);
	bool fok = fshader->compileSourceCode(
		"varying mediump vec4 color;"
		"varying mediump vec2 texc;"
		"uniform sampler2D texture;"
		"void main(void)"
		"{"
		"    gl_FragColor = vec4(1, 1, 1, texture2D(texture, texc).a) * color;"
		"}");
	if(!fok) LOG4CPLUS_WARN(MeshLog::logger_console, 
		"label fshader compile: " << fshader->log().toStdString());

	if(!vok || !fok) return false;

	program = new QOpenGLShaderProgram(parent);
	program->addShader(vshader);
	program->addShader(fshader);
	bool pok = program->link();
	if(!pok){
		LOG4CPLUS_WARN(MeshLog::logger_console, 
			"label gl program link: " << program->log().toStdString());
		return false;
	}

	vertexAttr = program->attributeLocation("vertex");
	colorAttr  = program->attributeLocation("colorIndex");
	texcAttr   = program->attributeLocation("texCoord");
	matrixUniform = program->uniformLocation("mvMatrix");
	texUniform    = program->uniformLocation("texture");

	program->bind();
	program->setUniformValue(texUniform, 0);

	//glActiveTexture(GL_TEXTURE0);
	glGenTextures(1, &tex);
	glBindTexture( GL_TEXTURE_2D, tex);
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
	glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
	glPixelStorei( GL_UNPACK_ALIGNMENT, 1);
	
	return true;
}
void ModelGLWidget::LabelRender::renderText(
		const char* text, float color, 
		const QVector3D& pos, const QVector3D& svx, const QVector3D& svy) const
{
	QVector3D text_pos = pos;

	for(const char* p = text; *p; p++) {
		if(FT_Load_Char(face, *p, FT_LOAD_RENDER)) continue;

		glTexImage2D( GL_TEXTURE_2D, 0, GL_ALPHA, face->glyph->bitmap.width, face->glyph->bitmap.rows,
						0, GL_ALPHA, GL_UNSIGNED_BYTE, face->glyph->bitmap.buffer);
 
		QVector3D gpos = text_pos + svx * face->glyph->bitmap_left + svy * face->glyph->bitmap_top;

		QVector<QVector3D> gvertices;
		QVector<GLfloat> gcolors;
		QVector<GLfloat> gtex;

		gvertices << gpos;
		gcolors << color;
		gtex << 0 << 0;
		gvertices << (gpos + svx * face->glyph->bitmap.width);
		gcolors << color;
		gtex << 1 << 0;
		gvertices << (gpos - svy * face->glyph->bitmap.rows);
		gcolors << color;
		gtex << 0 << 1;
		gvertices << (gpos + svx * face->glyph->bitmap.width - svy * face->glyph->bitmap.rows);
		gcolors << color;
		gtex << 1 << 1;

		program->setAttributeArray(vertexAttr, gvertices.constData());
		program->setAttributeArray(colorAttr, GL_FLOAT, gcolors.constData(), 1);
		program->setAttributeArray(texcAttr, GL_FLOAT, gtex.constData(), 2);
 
	    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
 
		text_pos += svx * (face->glyph->advance.x >> 6);
		text_pos += svy * (face->glyph->advance.y >> 6);
	}
}

void ModelGLWidget::LabelRender::render(const QMatrix4x4& modelview, const QVector3D& svx, const QVector3D& svy) const
{
	if(program == nullptr) return;
	if(vertices.isEmpty()) return;

	glEnable( GL_BLEND );
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

	program->bind();
	program->setUniformValue(matrixUniform, modelview);

	program->enableAttributeArray(vertexAttr);
	program->enableAttributeArray(colorAttr);
	program->enableAttributeArray(texcAttr);

	for(int i = 0; i < vertices.size(); i++){
		renderText("+", colors[i], vertices[i] - svx * 4 - svy * 4, svx, svy);
		//renderText("*", colors[i], vertices[i] - svx * 2 - svy * 9, svx, svy);
		const char* text = labels[i].c_str();
		if(*text != '\0')
			renderText(text, colors[i], vertices[i] + svx - svy, svx, svy);
	}

	program->disableAttributeArray(vertexAttr);
	program->disableAttributeArray(colorAttr);
	program->disableAttributeArray(texcAttr);

	program->release();

	glDisable( GL_BLEND );

}

void ModelGLWidget::draw()
{
	if( !m_initialized )
		initializeGL();

	if( ! m_bounding_box.valid ) return;

	QMatrix4x4 tm;
	tm.rotate( total_rot * last_rot );
	QMatrix3x3 nm = tm.normalMatrix();
	tm.scale( total_zoom * last_zoom );
	tm.translate( total_trans + last_trans );

	QMatrix4x4 modelview = m_proj_matrix * tm;

	//	m_view_set->storeEPSFile(trans_matrix, mvmatrix, view_mode, 
	// m_view_set->storeEPSFile(modelview.data(), tm.data(), view_mode,
	//LOG4CPLUS_INFO(MeshLog::logger_mesh, "*** QtMeshGen *** modelview:");
	//for (int i = 0; i < 4; i++) {
	//	for (int j = 0; j < 4; j++)
	//		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, modelview(i, j) << "\t";
	//	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, endl;
	//}
	//LOG4CPLUS_INFO(MeshLog::logger_mesh, "*** QtMeshGen *** tm:");
	//for (int i = 0; i < 4; i++) {
	//	for (int j = 0; j < 4; j++)
	//		LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, tm(i, j) << "\t";
	//	LOG4CPLUS_INFO(MeshLog::logger_mesh_stat, endl;
	//}
	//if (viewPrefDlg->clipPlanes.notEmpty()) {
	//	for (int i = 0; i < viewPrefDlg->clipPlanes.countInt(); i++) {
	//		auto cp = viewPrefDlg->clipPlanes[i];
	//		LOG4CPLUS_INFO(MeshLog::logger_mesh, "clip#" << i << "\t" << cp->getNormalVector() << "\t" << cp->getParamD());
	//	}
	//}

	float sf = 1.0f / total_zoom;
	QVector3D axis_ox_rot = axis_rot.rotatedVector(axis_ox) * sf;
	QVector3D axis_oy_rot = axis_rot.rotatedVector(axis_oy) * sf;

	// blocks+faces
	m_bf_render.renderBlocks(modelview, nm);
	m_bf_render.renderFaces(modelview, nm);
	// edges+points
	m_ep_render.renderEdges(modelview);
	m_ep_render.renderPoints(modelview, axis_ox_rot, axis_oy_rot);
	if(viewPrefDlg->selected_clip_plane > 0)
		m_ep_render.renderClipPlane(modelview, 
			viewPrefDlg->clipPlanes[viewPrefDlg->selected_clip_plane].get(), 
			m_bounding_box.getDiameter());

	// numbered points + labels
	m_label_render.render(modelview, axis_ox_rot, axis_oy_rot);

	glFlush();
}

void ModelGLWidget::setViewSet(MeshViewSet* set, bool with_reset)
{
	if(m_view_set) delete m_view_set;
	m_view_set = set;
	m_view_mesh = nullptr;

	makeCurrent();
	prepareGLScene();

	resizeGL(width(), height());

	if(with_reset) resetPositioning();

	updateGL();
}

void ModelGLWidget::setViewMesh(MeshContainer3d* mesh)
{
	if(m_view_set) delete m_view_set;
	m_view_set = nullptr;
	m_view_mesh = mesh;

	makeCurrent();
	prepareGLScene();

	resizeGL(width(), height());
	resetPositioning();
	updateGL();
}

void ModelGLWidget::BlockFaceRender::addTriangleFace(const QVector3D& pt1, const QVector3D& pt2, const QVector3D& pt3, unsigned char ci)
{
	face_vertices << pt1 << pt2 << pt3;
	QVector3D n = QVector3D::normal( pt2-pt1, pt3-pt1 );
	for(int i = 0; i < 3; i++){
		face_normals << n;
		face_colors << ci;
	}
}

void ModelGLWidget::BlockFaceRender::addBlockTriangleFace(const QVector3D& pt1, const QVector3D& pt2, const QVector3D& pt3, unsigned char ci)
{
	block_vertices << pt1 << pt2 << pt3;
	QVector3D n = QVector3D::normal( pt2-pt1, pt3-pt1 );
	for(int i = 0; i < 3; i++){
		block_normals << n;
		block_colors << ci;
	}
}

void ModelGLWidget::BlockFaceRender::addFaceVertex(const QVector3D& pt, const QVector3D& n, unsigned char ci)
{
	face_vertices << pt;
	face_normals << n;
	face_colors << ci;
}

void ModelGLWidget::BlockFaceRender::addFaceVertex(const FPoint3d& pt, const FVector3d& n, unsigned char ci)
{
	addFaceVertex( QVector3D(pt.x, pt.y, pt.z), QVector3D(n.x, n.y, n.z), ci);
}

void ModelGLWidget::BlockFaceRender::addBlockVertex(const QVector3D& pt, const QVector3D& n, unsigned char ci)
{
	block_vertices << pt;
	block_normals << n;
	block_colors << ci;
}

void ModelGLWidget::BlockFaceRender::addBlockVertex(const FPoint3d& pt, const FVector3d& n, unsigned char ci)
{
	addBlockVertex( QVector3D(pt.x, pt.y, pt.z), QVector3D(n.x, n.y, n.z), ci);
}

void ModelGLWidget::EdgePointRender::addEdgeVertex(const QVector3D& pt, unsigned char ci)
{
	edge_vertices << pt;
	edge_colors << ci;
}

void ModelGLWidget::EdgePointRender::addEdgeVertex(const FPoint3d& pt, unsigned char ci)
{
	addEdgeVertex( QVector3D(pt.x, pt.y, pt.z), ci);
}

void ModelGLWidget::EdgePointRender::addPointVertex(const QVector3D& pt, unsigned char ci, int id)
{
	point_vertices << pt;
	point_colors << ci;
	point_ids << id;
	any_number |= (id != EMPTY_ID);
}

void ModelGLWidget::EdgePointRender::addPointVertex(const FPoint3d& pt, unsigned char ci, int id)
{
	addPointVertex( QVector3D(pt.x, pt.y, pt.z), ci, id);
}

void ModelGLWidget::LabelRender::addLabel(const QVector3D& pt, const string& s, unsigned char ci)
{
	vertices << pt;
	labels << s;
	colors << ci;
}

void ModelGLWidget::LabelRender::addLabel(const FPoint3d& pt, const string& s, unsigned char ci)
{
	addLabel( QVector3D(pt.x, pt.y, pt.z), s, ci);
}

void ModelGLWidget::prepareGLScene()
{
	m_bf_render.clear();
	m_ep_render.clear();
	m_label_render.clear();

	m_bounding_box.clear();

	if( m_view_set || m_view_mesh ) {
		markHidden();
		drawGL();
	}else{
		// pyramid, just for test...
		QVector3D pt1( 0.0, 0.0, 0.0 );
		QVector3D pt2( 0.2, 0.0, 0.0 );
		QVector3D pt3( 0.0, 0.4, 0.0 );
		QVector3D pt4( 0.0, 0.0, 0.8 );

		m_bounding_box.addPoint( DPoint3d( 0.0, 0.0, 0.0) );
		m_bounding_box.addPoint( DPoint3d( 0.2, 0.0, 0.0) );
		m_bounding_box.addPoint( DPoint3d( 0.0, 0.4, 0.0) );
		m_bounding_box.addPoint( DPoint3d( 0.0, 0.0, 0.8) );

		m_bf_render.addTriangleFace(pt1, pt2, pt4, ICOLOR_DYELLOW);
		m_bf_render.addTriangleFace(pt2, pt3, pt4, ICOLOR_DRED);
		m_bf_render.addTriangleFace(pt3, pt1, pt4, ICOLOR_DGREEN);
		m_bf_render.addTriangleFace(pt1, pt3, pt2, ICOLOR_DBLUE);

		// ------
		m_label_render.addLabel(pt4, "ala ma 3 kotki?", ICOLOR_WHITE);
		m_ep_render.addPointVertex(pt3, ICOLOR_WHITE, 123450);
		m_ep_render.addPointVertex(pt2, ICOLOR_WHITE, -678910);
	}

	if(m_bounding_box.valid){
		m_center = m_bounding_box.getMiddlePoint();
		total_trans = QVector3D(-m_center.x, -m_center.y, -m_center.z);
	}
}

void ModelGLWidget::drawGL()
{
	if(!m_view_set && !m_view_mesh) return;

	const unsigned char col_faces[] = { 
		ICOLOR_DBLUE, ICOLOR_DGREEN, ICOLOR_DRED, ICOLOR_DPURPLE, ICOLOR_DYELLOW, ICOLOR_DCYAN };
	const unsigned char col_points[] = { 
		ICOLOR_WHITE-1, ICOLOR_YELLOW, ICOLOR_CYAN, ICOLOR_RED, ICOLOR_GREEN, ICOLOR_BLUE };
	const unsigned char col_edges[] = { 
		ICOLOR_WHITE-2, ICOLOR_YELLOW, ICOLOR_CYAN, ICOLOR_PURPLE, ICOLOR_RED };
	const unsigned char col_edges_x[] = { 
		ICOLOR_DYELLOW, ICOLOR_YELLOW, ICOLOR_CYAN, ICOLOR_PURPLE, ICOLOR_RED };

	DataVector<MeshContainer3d*> available_meshes;
	if(!m_view_set){
		int dvct = m_view_mesh->getBlocksCount();
		if(dvct > 0){
			MeshBlock* block0 = m_view_mesh->getBlockAt(0);
			if(block0->getType() != BLOCK_DOMAIN) dvct = 1; // single 3d mesh
		}

		MeshContainer3d* mesh3d = m_view_mesh->getTotalMDV()->getMesh();

		if(mesh3d){
			available_meshes.add(mesh3d);
		}else{
			for(int k = 0; k < dvct; k++){
				MeshBlock* block = m_view_mesh->getBlockAt(k);
				if(block->getType() == BLOCK_DOMAIN){
					mesh3d = ((MeshDomainVolume*)block)->getMesh();
					if(mesh3d) 
						available_meshes.add(mesh3d);
				}else
					available_meshes.add(m_view_mesh);
			}
		}
	}

	// blocks
	if(view_mode & VIEW_BLOCKS){
		if(m_view_set){
			int count = m_view_set->m_blocks.count();

			MeshViewSet::PolygonFill mode = m_view_set->getPolygonFillMode();
			for(int i = 0; i < count; i++){
				auto data = m_view_set->m_blocks[i];
				if(data->hidden) continue;

				for(int j = 0; j < data->pts.count(); j++)
					m_bounding_box.addPoint( data->pts[j] );

				assert( data->indices.notEmpty() );
				int ci = ICOLOR_DGRAY;
				if(mode == MeshViewSet::FILL_QUALITY){
					ci = (data->quality < 0.0) ? ICOLOR_DRED : (int) ( data->quality * ICOLOR_DPURPLE );
				}else{
					assert(mode == MeshViewSet::FILL_AREA);
					if(data->part < 0) ci = ICOLOR_DORANGE;
					else if(data->area_id >= 0) {
						if(data->area_id < 6) ci = col_faces[data->area_id];
						else ci = genIColor(data->area_id);
					}
				}

				for(int j = 0; j < data->normals.count(); j++){
					for(int k = 0; k < 3; k++){
						m_bf_render.addBlockVertex(data->pts[data->indices[3*j+k]], data->normals[j], ci);
					}
				}
			}
		}else{ // m_view_mesh
			MeshViewSet::PolygonFill mode = MeshViewSet::FILL_AREA;

			for(int k = 0; k < available_meshes.count(); k++){
				MeshContainer3d* mesh3d = available_meshes[k];

				int bct = mesh3d->getBlocksCount();

				for(int i = 0; i < bct; i++){
					MeshTetrahedron* block = (MeshTetrahedron*)mesh3d->getBlockAt(i);
					if(block->getType() != BLOCK_TETRA || block->availableTag(TagExtended::TAG_HIDDEN)) continue;

					for(int j = 0; j < block->getPointCount(); j++)
						m_bounding_box.addPoint( block->getPoint(j)->getCoordinates() );

					int ci = ICOLOR_DGRAY;
					if(mode == MeshViewSet::FILL_QUALITY){
						double q = block->getQuality();
						ci = (q < 0.0) ? ICOLOR_DRED : (int) ( q * ICOLOR_DPURPLE );
					}else{
						assert(mode == MeshViewSet::FILL_AREA);
						int area_id = block->getAreaID();
						if(area_id < 0) ci = ICOLOR_DORANGE;
						else if(area_id < 6) ci = col_faces[area_id];
						else ci = genIColor(area_id);
					}

					int fct = block->getFaceCount();
					const FPoint3d middle = block->getMiddlePoint();
					for(int j = 0; j < fct; j++){
						MeshFace* face = block->getFace(j);
						int fpct = face->getPointCount();
						assert( fpct == 3 );
						DVector3d dn = face->getNormalVector();
						auto item = face->getViewData(1.0, face->getBlockIndex(block) == 0);
						if( item ){
							DataVector<FPoint3d> item_pts(item->pts.count());
							for(int k = 0; k < item->pts.count(); k++)
								item_pts.add( FPoint3d( middle, item->pts[k], MeshViewSet::param_shrink) );
							for(int k = 0; k < item->indices.count(); k++)
								m_bf_render.addBlockVertex( item_pts[ item->indices[k] ], item->normal, ci );
						}
					}
				}
			}
		}
	}

	if(view_mode & VIEW_FACES){
		// faces
		if(m_view_set){
			int count = m_view_set->m_faces.count();
			MeshViewSet::PolygonFill mode = m_view_set->getPolygonFillMode();
			for(int i = 0; i < count; i++){
				auto data = m_view_set->m_faces[i];
				if(data->hidden) continue;

				for(int j = 0; j < data->pts.count(); j++)
					m_bounding_box.addPoint( data->pts[j] );

				// draw
				assert( data->indices.notEmpty() );
				if(mode == MeshViewSet::FILL_NODES){
					for(int j = 0; j < data->indices.count(); j+=3){
						for(int k = 0; k < 3; k++){
							float q = data->wpts[data->indices[j+k]];
							int ci = (q < 0.0) ? ICOLOR_DRED : (int) ( q * ICOLOR_DBLUE );
							m_bf_render.addFaceVertex(data->pts[data->indices[j+k]], data->normal, ci);
						}
					}
				}else{
					int ci = ICOLOR_DGRAY;
					if(mode == MeshViewSet::FILL_QUALITY){
						ci = (data->quality < 0.0) ? ICOLOR_DRED : (int) ( data->quality * ICOLOR_DBLUE );
					}else{
						assert(mode == MeshViewSet::FILL_AREA);
						if(data->part > 0) ci = ICOLOR_DORANGE;
						else if(data->area_id >= 0) {
							if(data->area_id < 6) ci = col_faces[data->area_id];
							else ci = genIColor(data->area_id);
						}
					}

					for(int j = 0; j < data->indices.count(); j+=3){
						for(int k = 0; k < 3; k++){
							m_bf_render.addFaceVertex(data->pts[data->indices[j+k]], data->normal, ci);
						}
					}
				}
			}
		}
	}

	if(view_mode & VIEW_EDGES){
		// edges
		if(m_view_set){
			int count = m_view_set->m_edges.countInt();
			for(int i = 0; i < count; i++){
				auto data = m_view_set->m_edges[i];
				if(data->hidden) continue;

				m_bounding_box.addPoint( data->pt1 );
				m_bounding_box.addPoint( data->pt2 );

				int ci = ICOLOR_DGRAY;
				if(data->part > 0) ci = col_edges[data->part % 5];
				else if(data->part < 0) ci = (view_mode & VIEW_WHITE) ? ICOLOR_GRAY : ICOLOR_DGRAY;
				else{
					assert(data->border < 4);
					ci = (view_mode & VIEW_WHITE) ? col_edges_x[data->border] : col_edges[data->border];
				}

				m_ep_render.addEdgeVertex( data->pt1, ci);
				m_ep_render.addEdgeVertex( data->pt2, ci);
			}
		}else{
			for(int k = 0; k < available_meshes.count(); k++){
				MeshContainer3d* mesh3d = available_meshes[k];
				
				for(IteratorEdge3d it = mesh3d->getFirstEdge3d(); it.isValid(); it.nextEdge()){
					MeshEdge3d* edge = it.getEdge();
					if(edge->availableTag(TagExtended::TAG_HIDDEN)) continue;

					const DPoint3d& pt1 = edge->getMeshPoint(0)->getCoordinates();
					const DPoint3d& pt2 = edge->getMeshPoint(1)->getCoordinates();

					m_bounding_box.addPoint( pt1 );
					m_bounding_box.addPoint( pt2 );

					int mat_id = 0;
					if(edge->isBorder(TagBorder::FIXED)) mat_id = 3;
					else if(edge->isBorder(TagBorder::RIDGE)) mat_id = 2;
					else if(edge->isBorder()) mat_id = 1;

					int ci = (view_mode & VIEW_WHITE) ? col_edges_x[mat_id] : col_edges[mat_id];

					m_ep_render.addEdgeVertex( pt1, ci);
					m_ep_render.addEdgeVertex( pt2, ci);
				}
			}
		}
	}

	// vertices (without numbers)
	if(view_mode & VIEW_NODES){
		if(m_view_set){
			int count = m_view_set->m_points.count();

			for(int i = 0; i < count; i++){
				auto data = m_view_set->m_points[i];
				if(data->hidden) continue;

				m_bounding_box.addPoint( data->pt );

				int ci = ICOLOR_DGRAY;
				if(data->part > 0) ci = col_points[data->part % 5];
				else if(data->part < 0) ci = (view_mode & VIEW_WHITE) ? ICOLOR_GRAY : ICOLOR_DGRAY;
				else{
					assert(data->border < 4);
					ci = (view_mode & VIEW_WHITE) ? col_edges_x[data->border] : col_edges[data->border];
				}

				if(data->numbered && (view_mode & VIEW_LABELS))
					m_ep_render.addPointVertex( data->pt, ci, data->id);	//m_label_render.addLabel( data->pt, to_string(data->id), ci);
				else
					m_ep_render.addPointVertex( data->pt, ci);
			}
		}else{
			for(int k = 0; k < available_meshes.count(); k++){
				MeshContainer3d* mesh3d = available_meshes[k];

				int count = mesh3d->getPointsCount();

				for(int i = 0; i < count; i++){
					MeshPoint3d* point = mesh3d->getPointAt(i);
					if(point->availableTag(TagExtended::TAG_HIDDEN)) continue;

					m_bounding_box.addPoint( point->getCoordinates() );

					int mat_id = 0;
					if(point->isBorder(TagBorder::CORNER | TagBorder::FIXED)) mat_id = 3;
					else if(point->isBorder(TagBorder::RIDGE)) mat_id = 2;
					else if(point->isBorder()) mat_id = 1;

					int ci = (view_mode & VIEW_WHITE) ? col_edges_x[mat_id] : col_edges[mat_id];

					if(view_mode & VIEW_LABELS)
						m_ep_render.addPointVertex( point->getCoordinates(), ci, point->getIndex());
					else
						m_ep_render.addPointVertex( point->getCoordinates(), ci);

					//m_label_render.addLabel( point->getCoordinates(), to_string(point->getIndex()), ci);
				}
			}
		}
	}

	if(view_mode & VIEW_LABELS){
		if(m_view_set && m_view_set->m_labels.count() > 0){
			for(int i = 0; i < m_view_set->m_labels.count(); i++){
				const MeshViewSet::LabelInfo &li = m_view_set->m_labels[i];
				m_bounding_box.addPoint( li.pt );
				m_label_render.addLabel( li.pt, li.label, ICOLOR_RED);
			}
		}
	}
}

void ModelGLWidget::markHidden()
{
	if(!m_view_set && !m_view_mesh) return;

	DataVector<MeshContainer3d*> available_meshes;
	if(!m_view_set){
		int dvct = m_view_mesh->getBlocksCount();
		if(dvct > 0){
			MeshBlock* block0 = m_view_mesh->getBlockAt(0);
			if(block0->getType() != BLOCK_DOMAIN) dvct = 1; // single 3d mesh
		}

		MeshContainer3d* mesh3d = m_view_mesh->getTotalMDV()->getMesh();
		if(mesh3d){
			available_meshes.add(mesh3d);
		}else{
			for(int k = 0; k < dvct; k++){
				MeshBlock* block = m_view_mesh->getBlockAt(k);
				if(block->getType() == BLOCK_DOMAIN){
					mesh3d = ((MeshDomainVolume*)block)->getMesh();
					if(mesh3d) 
						available_meshes.add(mesh3d);
				}else
					available_meshes.add(m_view_mesh);
			}
		}
	}

	// blocks

	if(view_mode & VIEW_BLOCKS){
		if(m_view_set){
			int count = m_view_set->m_blocks.count();
			// 1. check clip
			for(int i = 0; i < count; i++){
				auto data = m_view_set->m_blocks[i];
				data->hidden = 0;
				if(clip_quality < 1.05 && data->quality > clip_quality){
					data->hidden = 1;
					continue;
				}
				for(int k = 1; k < viewPrefDlg->clipPlanes.count() && !data->hidden; k++){
					for(int j = 0; j < data->pts.count(); j++){
						if(viewPrefDlg->clipPlanes[k]->clipped(data->pts[j])){
							data->hidden = 1; 
							break;
						}
					}
				}
			}
			// 2. check sight-line
			for(int i = 0; i < count; i++){
				auto data = m_view_set->m_blocks[i];
				if(data->hidden) continue;
				if(data->adjacent.count() == data->normals.count()){ // if has all neighbours
					data->hidden = 2;
					for(int j = 0; j < data->adjacent.count(); j++){
						// if boundary or adjacent to clip-plane
						if(data->adjacent[j]->hidden == 1){
							data->hidden = 0; // make visible
							break;
						}
					}
				}
			}
		}else{

			for(int k = 0; k < available_meshes.count(); k++){
				MeshContainer3d* mesh3d = available_meshes[k];

				int bct = mesh3d->getBlocksCount();
				for(int i = 0; i < bct; i++){
					MeshTetrahedron* block = (MeshTetrahedron*)mesh3d->getBlockAt(i);
					// 1. check clip
					if(clip_quality < 1.05 && block->getQuality() > clip_quality){
						block->setIntTag(TagExtended::TAG_HIDDEN, 1);
					}else{
						bool hidden = false;
						for(int j = 1; j < viewPrefDlg->clipPlanes.count() && !hidden; j++){
							if(viewPrefDlg->clipPlanes[j]->clipped(block->getMiddlePoint())){
								hidden = true;
								block->setIntTag(TagExtended::TAG_HIDDEN, 1);
							}
						}
						if(!hidden) block->removeTag(TagExtended::TAG_HIDDEN);
					}
				}
			}

			for(int k = 0; k < available_meshes.count(); k++){
				MeshContainer3d* mesh3d = available_meshes[k];

				int bct = mesh3d->getBlocksCount();
				for(int i = 0; i < bct; i++){
					MeshTetrahedron* block = (MeshTetrahedron*)mesh3d->getBlockAt(i);
					// 2. check sight-line
					if(block->availableTag(TagExtended::TAG_HIDDEN)) continue;
					block->setIntTag(TagExtended::TAG_HIDDEN, 2);
					int fct = block->getFaceCount();
					for(int j = 0; j < fct; j++){
					// if boundary or adjacent to clip-plane
						MeshBlock* fblock = block->getNeighbour(j);
						if(!fblock || fblock->checkIntTag(TagExtended::TAG_HIDDEN, 1)){
							block->removeTag(TagExtended::TAG_HIDDEN); // make visible
							break;
						}
					}
				}
			}
		}
	}else{ // all hidden
		if(m_view_set){
			int count = m_view_set->m_blocks.count();
			for(int i = 0; i < count; i++)
				m_view_set->m_blocks.get(i)->hidden = 1;
		}
	}

	if(view_mode & VIEW_FACES){
		// faces
		if(m_view_set){
			int count = m_view_set->m_faces.count();
			for(int i = 0; i < count; i++){
				auto data = m_view_set->m_faces[i];
				// check clip
				if(clip_quality < 1.05 && data->quality > clip_quality){
					data->hidden = 1;
					continue;
				}
				data->hidden = 0;
				for(int k = 1; k < viewPrefDlg->clipPlanes.count() && !data->hidden; k++){
					for(int j = 0; j < data->pts.count(); j++){
						if(viewPrefDlg->clipPlanes[k]->clipped(data->pts[j])){
							data->hidden = 1; 
							break;
						}
					}
				}
			}
		}
	}else{ // all hidden
		if(m_view_set){
			int count = m_view_set->m_faces.count();
			for(int i = 0; i < count; i++)
				m_view_set->m_faces.get(i)->hidden = 1;
		}
	}

	if(view_mode & VIEW_EDGES){
		// edges
		if(m_view_set){
			int count = m_view_set->m_edges.countInt();
			for(int i = 0; i < count; i++){
				auto data = m_view_set->m_edges[i];
				data->hidden = 0;
				if(data->adjacent_blocks.count() == 0){	// check clip
					for(int k = 1; k < viewPrefDlg->clipPlanes.count(); k++){
						if(viewPrefDlg->clipPlanes[k]->clipped(data->pt1)){
							data->hidden = 1;
							break;
						}else if(viewPrefDlg->clipPlanes[k]->clipped(data->pt2)){
							data->hidden = 1;
							break;
						}
					}
				}else{	// check if any adjacent block is visible
					data->hidden = 1;
					for(int j = 0; j < data->adjacent_blocks.count(); j++){
						if(!data->adjacent_blocks[j]->hidden){
							data->hidden = 0;
							break;
						}
					}
				}
			}
		}else{

			for(int k = 0; k < available_meshes.count(); k++){
				MeshContainer3d* mesh3d = available_meshes[k];

				for(IteratorEdge3d it = mesh3d->getFirstEdge3d(); it.isValid(); it.nextEdge()){
					it.getEdge()->setIntTag(TagExtended::TAG_HIDDEN, 1);
				}
				int bct = mesh3d->getBlocksCount();
				for(int i = 0; i < bct; i++){
					MeshTetrahedron* block = (MeshTetrahedron*)mesh3d->getBlockAt(i);
					if(!block->availableTag(TagExtended::TAG_HIDDEN)){
						int ect = block->getEdgeCount();
						for(int j = 0; j < ect; j++)
							block->getEdge(j)->removeTag(TagExtended::TAG_HIDDEN);
					}
				}
			}
		}
	}else{ // all hidden
		if(m_view_set){
			int count = m_view_set->m_edges.countInt();
			for(int i = 0; i < count; i++)
				m_view_set->m_edges.get(i)->hidden = 1;
		}
	}

	// vertices
	if(view_mode & VIEW_NODES){
		if(m_view_set){
			int count = m_view_set->m_points.count();
			for(int i = 0; i < count; i++){
				auto data = m_view_set->m_points[i];

				data->hidden = 0;
				for(int k = 1; k < viewPrefDlg->clipPlanes.count() && !data->hidden; k++){
					if(viewPrefDlg->clipPlanes[k]->clipped(data->pt)){
						data->hidden = 1;
						break;
					}
				}
			}
		}else{

			for(int k = 0; k < available_meshes.count(); k++){
				MeshContainer3d* mesh3d = available_meshes[k];
				for(int i = 0; i < mesh3d->getPointsCount(); i++){
					MeshPoint3d* point = mesh3d->getPointAt(i);
					bool hidden = false;
					for(int j = 1; j < viewPrefDlg->clipPlanes.count(); j++){
						if(viewPrefDlg->clipPlanes[j]->clipped(point->getCoordinates())){
							hidden = true;
							break;
						}
					}
					if(hidden)
						point->setIntTag(TagExtended::TAG_HIDDEN, 1);
					else
						point->removeTag(TagExtended::TAG_HIDDEN);
				}
			}
		}
	}else{
		if(m_view_set){
			int count = m_view_set->m_points.count();
			for(int i = 0; i < count; i++)
				m_view_set->m_points.get(i)->hidden = 1;
		}
	}
}

int ModelGLWidget::genIColor(int area_id)
{
	return (1103515245 * area_id + 12345) % 120;
}

void ModelGLWidget::resetPositioning()
{
	last_rot = total_rot = QQuaternion();
	last_axis_rot = axis_rot = QQuaternion();
	last_trans = QVector3D();
	total_trans = QVector3D(-m_center.x, -m_center.y, -m_center.z);
	last_zoom = total_zoom = 1.0f;
	nav_mode = 0;

	updateGL();
}

void ModelGLWidget::setViewMode(int m)
{
	view_mode = m;
	redraw(true);
}

void ModelGLWidget::redraw(bool recreate)
{
	if(recreate){
		makeCurrent();
		prepareGLScene();
		resizeGL(width(), height());
	}
	updateGL();
}

void ModelGLWidget::showPrefDlg()
{
	viewPrefDlg->setViewMode(view_mode);
	viewPrefDlg->show();
}

void ModelGLWidget::clearView()
{
	setViewSet(nullptr);
}

void ModelGLWidget::storeMatlabFile()
{
	if(!m_view_set) return;

	static QString dir = ".";
	QString fullFileName = QFileDialog::getSaveFileName(this,
						tr("Choose a file name"), dir,
						tr("Matlab File (*.m)"));
	if(fullFileName.isEmpty()) return;
	dir = fullFileName;
	m_view_set->storeMatlabFile(
		QFileInfo(fullFileName).absolutePath().toStdString(), 
		QFileInfo(fullFileName).fileName().toStdString());
}

void ModelGLWidget::storeEPSFile()
{
	if(!m_view_set) return;

	static QString dir = ".";
	QString fullFileName = QFileDialog::getSaveFileName(this,
						tr("Choose a file name"), dir,
						tr("EPS File (*.eps)"));
	if(fullFileName.isEmpty()) return;
	dir = fullFileName;

	QMatrix4x4 tm;
	tm.rotate( total_rot * last_rot );
	tm.scale( total_zoom * last_zoom );
	tm.translate( total_trans + last_trans );

	QMatrix4x4 modelview = m_proj_matrix * tm;

	//double mvmatrix[16];
	//double projmatrix[16];
	//double trans_matrix[16];
	//makeCurrent();
	//glPushMatrix();
	//	glMultMatrixf(tMatrix);
	//	glGetDoublev (GL_MODELVIEW_MATRIX, mvmatrix);
	//	glGetDoublev (GL_PROJECTION_MATRIX, projmatrix);
	//glPopMatrix();
	//MultiplyMatrices4by4OpenGL(trans_matrix, projmatrix, mvmatrix);

//	m_view_set->storeEPSFile(trans_matrix, mvmatrix, view_mode, 
	m_view_set->storeEPSFile(modelview.data(), tm.data(), view_mode, 
//	m_view_set->storeEPSFile(trans_matrix, view_mode, 
		QFileInfo(fullFileName).absolutePath().append(
			tr("/")).append(QFileInfo(fullFileName).fileName()).toStdString()
	);
}

ViewPrefDialog::ViewPrefDialog(QWidget *parent, ModelGLWidget* glw)
    : QDialog(parent), glWidget(glw), selected_clip_plane(0)
{
	chkViewBlocks = new QCheckBox(tr("&blocks"));
	chkViewFaces  = new QCheckBox(tr("&faces"));
	chkViewEdges  = new QCheckBox(tr("&edges"));
	chkViewNodes  = new QCheckBox(tr("&nodes"));
	chkViewLabels = new QCheckBox(tr("&labels"));

	chkViewWhite  = new QCheckBox(tr("&white"));
	chkStoreBackFaces  = new QCheckBox(tr("ba&ck faces"));

	setViewMode(glw->getViewMode());

	QVBoxLayout *vbox = new QVBoxLayout;
	vbox->addWidget(chkViewBlocks);
	vbox->addWidget(chkViewFaces);
	vbox->addWidget(chkViewEdges);
	vbox->addWidget(chkViewNodes);
	vbox->addWidget(chkViewLabels);
	vbox->addStretch(1);
	vbox->addWidget(chkViewWhite);
	vbox->addWidget(chkStoreBackFaces);
	vbox->addStretch(1);

	QGroupBox* view_box = new QGroupBox(tr("View entities"));
    view_box->setLayout(vbox);

	vbox = new QVBoxLayout;
//	QListWidget *lw = new QListWidget;

	clipTable = new QTableWidget(this);
	clipTable->setColumnCount(5);
	clipTable->setHorizontalHeaderLabels(
		QStringList() << tr("Plane") << tr("A") << tr("B") << tr("C") << tr("D"));

	clipTable->setRowCount(1);
	clipPlanes.add(
		std::make_shared<ClipPlane>(DPoint3d(0.0, 0.0, 0.0), DVector3d(0.0, 0.0, 0.0), tr("none")));
	clipPlanes[0]->setTableRow(clipTable, 0);

	vbox->addWidget(clipTable);

	for(int i = 0; i < 5; i++)
			clipTable->resizeColumnToContents(i);
	connect(clipTable, SIGNAL(itemSelectionChanged()), this, SLOT(otherClipPlaneSelected()));

	QHBoxLayout *clipBtnLayout = new QHBoxLayout;
	QPushButton *addClipPlaneButton = new QPushButton(tr("Add OXY"));
	connect(addClipPlaneButton, SIGNAL(clicked()), this, SLOT(addClipPlaneOXY()));
	QPushButton *delClipPlaneButton = new QPushButton(tr("Delete"));
	connect(delClipPlaneButton, SIGNAL(clicked()), this, SLOT(delSelectedClipPlane()));
	QPushButton *applyClipPlaneButton = new QPushButton(tr("Apply"));
	connect(applyClipPlaneButton, SIGNAL(clicked()), this, SLOT(applyEditedClipPlane()));
	clipBtnLayout->addWidget(addClipPlaneButton);
	clipBtnLayout->addWidget(delClipPlaneButton);
	clipBtnLayout->addWidget(applyClipPlaneButton);
//	clipBtnLayout->addStretch();

	vbox->addLayout(clipBtnLayout);

	QGroupBox* clip_planes_box = new QGroupBox(tr("Clip planes"));
	clip_planes_box->setLayout(vbox);

	redrawButton = new QPushButton(tr("&Redraw"));
	closeButton  = new QPushButton(tr("Close"));

	connect(redrawButton, SIGNAL(clicked()), this, SLOT(redrawGL()));
	connect(closeButton, SIGNAL(clicked()), this, SLOT(close()));

	QHBoxLayout *bottomLayout = new QHBoxLayout;
	bottomLayout->addStretch();
	bottomLayout->addWidget(redrawButton);
	bottomLayout->addWidget(closeButton);
	bottomLayout->addStretch();

	QHBoxLayout *chLayout = new QHBoxLayout;
	chLayout->addWidget(view_box);
	chLayout->addWidget(clip_planes_box);

	QVBoxLayout *cvLayout = new QVBoxLayout;
	cvLayout->addLayout(chLayout);
	cvLayout->addLayout(bottomLayout);

	// apply layout
	setLayout(cvLayout);

	setWindowTitle(tr("Mesh view preferences"));
}

void ViewPrefDialog::addClipPlaneOXY()
{
	int row_ct = clipTable->rowCount();
	clipTable->setRowCount(row_ct+1);

	auto plane = std::make_shared<ClipPlane>(
		glWidget->getCenter(), DVector3d(0.0, 0.0, 1.0), 
		tr("clip plane [%1]").arg(ClipPlane::counter++));
	clipPlanes.add(plane);
	plane->setTableRow(clipTable, row_ct);

	for(int i = 0; i < 5; i++)
			clipTable->resizeColumnToContents(i);

	clipTable->selectRow(row_ct);
	selected_clip_plane = row_ct;

	glWidget->redraw(true);
}

void ViewPrefDialog::delSelectedClipPlane()
{
	int row = clipTable->currentRow();
	if(row > 0){
		clipPlanes.removeOrderedAt(row);
		clipTable->removeRow(row);
		if(row == selected_clip_plane){ 
			selected_clip_plane = 0;
		}
		glWidget->redraw(true);
	}
}

void ViewPrefDialog::applyEditedClipPlane()
{
	int row = clipTable->currentRow();
	if(row < 1) return;

	QTableWidgetItem *item = clipTable->item(row, 0);
	clipPlanes[row]->setName(item->text());

	bool ok;
	DVector3d vn = clipPlanes[row]->getNormalVector();
	// vx
	item = clipTable->item(row, 1);
	double v = item->text().toDouble(&ok);
	if(ok) vn.x = v;
	else item->setText(tr("%1").arg(vn.x));
	// vy
	item = clipTable->item(row, 2);
	v = item->text().toDouble(&ok);
	if(ok) vn.y = v;
	else item->setText(tr("%1").arg(vn.y));
	// vz
	item = clipTable->item(row, 3);
	v = item->text().toDouble(&ok);
	if(ok) vn.z = v;
	else item->setText(tr("%1").arg(vn.z));
	// D
	item = clipTable->item(row, 4);
	v = item->text().toDouble(&ok);
	if(!ok) item->setText(tr("%1").arg(v = clipPlanes[row]->getParamD()));

	clipPlanes[row]->setData(vn, v);

	selected_clip_plane = row;
	glWidget->redraw(true);
}

void ViewPrefDialog::otherClipPlaneSelected()
{
	int row = clipTable->currentRow();
	if(row != selected_clip_plane){
		selected_clip_plane = row;
		glWidget->redraw();
	}
}

void ViewPrefDialog::setViewMode(int mode)
{
	last_mode = mode;
	chkViewBlocks->setChecked(mode & ModelGLWidget::VIEW_BLOCKS);
	chkViewFaces->setChecked (mode & ModelGLWidget::VIEW_FACES);
	chkViewEdges->setChecked (mode & ModelGLWidget::VIEW_EDGES);
	chkViewNodes->setChecked (mode & ModelGLWidget::VIEW_NODES);
	chkViewLabels->setChecked (mode & ModelGLWidget::VIEW_LABELS);
	chkViewWhite->setChecked (mode & ModelGLWidget::VIEW_WHITE);
	chkStoreBackFaces->setChecked (mode & ModelGLWidget::STORE_BACK_FACES);
}

void ViewPrefDialog::redrawGL()
{
	last_mode = 0;
	if(chkViewBlocks->isChecked()) last_mode |= ModelGLWidget::VIEW_BLOCKS;
	if(chkViewFaces->isChecked()) last_mode |= ModelGLWidget::VIEW_FACES;
	if(chkViewEdges->isChecked()) last_mode |= ModelGLWidget::VIEW_EDGES;
	if(chkViewNodes->isChecked()) last_mode |= ModelGLWidget::VIEW_NODES;
	if(chkViewLabels->isChecked()) last_mode |= ModelGLWidget::VIEW_LABELS;
	if(chkViewWhite->isChecked()) last_mode |= ModelGLWidget::VIEW_WHITE;
	if(chkStoreBackFaces->isChecked()) last_mode |= ModelGLWidget::STORE_BACK_FACES;

	glWidget->setViewMode(last_mode);
}

void ViewPrefDialog::closeEvent(QCloseEvent *event)
{
	int mode = 0;
	if(chkViewBlocks->isChecked()) mode |= ModelGLWidget::VIEW_BLOCKS;
	if(chkViewFaces->isChecked())  mode |= ModelGLWidget::VIEW_FACES;
	if(chkViewEdges->isChecked())  mode |= ModelGLWidget::VIEW_EDGES;
	if(chkViewNodes->isChecked())  mode |= ModelGLWidget::VIEW_NODES;
	if(chkViewLabels->isChecked())  mode |= ModelGLWidget::VIEW_LABELS;
	if(chkViewWhite->isChecked())  mode |= ModelGLWidget::VIEW_WHITE;
	if(chkStoreBackFaces->isChecked())  mode |= ModelGLWidget::STORE_BACK_FACES;

	if(mode != last_mode){
		last_mode = mode;
		selected_clip_plane = 0;
		glWidget->setViewMode(mode);
	}else if (selected_clip_plane > 0){
		selected_clip_plane = 0;
		glWidget->redraw();
	}

	event->accept();
}

ClipPlane::ClipPlane(const DPoint3d& pt, const DVector3d& vt, const QString& name)
	: m_pt(pt), m_vt(vt), m_d(-(pt.x * vt.x + pt.y * vt.y + pt.z * vt.z)), m_name(name)
{ }

void ClipPlane::setTableRow(QTableWidget* table, int row) const
{
	QTableWidgetItem *item = new QTableWidgetItem(m_name);
	table->setItem(row, 0, item);
	item->setFlags(Qt::ItemIsEnabled);

	if(m_name == QString("none")){
		for(int i = 1; i < 5; i++){
			item = new QTableWidgetItem(QString("-"));
			item->setFlags(Qt::NoItemFlags);
			table->setItem(row, i, item);
		}
	}else{
		item = new QTableWidgetItem(QString("%1").arg(m_vt.x));
		item->setFlags(Qt::ItemIsEditable | Qt::ItemIsEnabled);
		table->setItem(row, 1, item);

		item = new QTableWidgetItem(QString("%1").arg(m_vt.y));
		item->setFlags(Qt::ItemIsEditable | Qt::ItemIsEnabled);
		table->setItem(row, 2, item);

		item = new QTableWidgetItem(QString("%1").arg(m_vt.z));
		item->setFlags(Qt::ItemIsEditable | Qt::ItemIsEnabled);
		table->setItem(row, 3, item);

		item = new QTableWidgetItem(QString("%1").arg(m_d));
		item->setFlags(Qt::ItemIsEditable | Qt::ItemIsEnabled);
		table->setItem(row, 4, item);
	}
}

bool ClipPlane::clipped(const DPoint3d& pt) const
{
	return (pt.x * m_vt.x + pt.y * m_vt.y + pt.z * m_vt.z + m_d) > 0.0;
}

bool ClipPlane::clipped(const FPoint3d& pt) const
{
	return (pt.x * m_vt.x + pt.y * m_vt.y + pt.z * m_vt.z + m_d) > 0.0;
}

void ClipPlane::translate(const DVector3d& v)
{
	m_pt += v;
	m_d = -(m_pt.x * m_vt.x + m_pt.y * m_vt.y + m_pt.z * m_vt.z);
}

void ClipPlane::setNormalVector(const DVector3d& v)
{
	m_vt = v;
	m_d = -(m_pt.x * m_vt.x + m_pt.y * m_vt.y + m_pt.z * m_vt.z);
}

void ClipPlane::setData(const DVector3d& v, double d)
{
	m_vt = v;
	m_d = d;
	// calculate m_pt;
	double dist = m_pt.x * m_vt.x + m_pt.y * m_vt.y + m_pt.z * m_vt.z + m_d;
	m_pt += m_vt * (dist / m_vt.length());
}

void ClipPlane::drawGL(QVector<QVector3D> & clip_vertices, QVector<GLfloat> & clip_colors, float diameter) const
{
	const unsigned char ci = 140; // ICOLOR_ORANGE ?
	static const int GRID_COUNT = 10;

	double vt_x = fabs(m_vt.x);
	double vt_y = fabs(m_vt.y);
	double vt_z = fabs(m_vt.z);

	QVector3D vt0(0.0, 0.0, 1.0);
	if(vt_x < vt_y && vt_x < vt_z){
		vt0 = QVector3D(1.0, 0.0, 0.0);
	}else if(vt_y < vt_z){
		vt0 = QVector3D(0.0, 1.0, 0.0);
	}

	QVector3D mvt(m_vt.x, m_vt.y, m_vt.z);
	QVector3D vt1 = QVector3D::crossProduct(mvt, vt0);
	vt0 = QVector3D::crossProduct(mvt, vt1);

	double vlen = diameter / GRID_COUNT;
	vt0 *= vlen / vt0.length();
	vt1 *= vlen / vt1.length();

	QVector3D pt0( m_pt.x, m_pt.y, m_pt.z );

	QVector3D pt00 = pt0 - vt0*(GRID_COUNT/2) - vt1*(GRID_COUNT/2);
	QVector3D pt10 = pt0 + vt0*(GRID_COUNT/2) - vt1*(GRID_COUNT/2);
	QVector3D pt01 = pt0 - vt0*(GRID_COUNT/2) + vt1*(GRID_COUNT/2);
	for(int i = 0; i <= GRID_COUNT; i++){
		QVector3D pt1 = pt00 + vt0 * i;
		QVector3D pt2 = pt01 + vt0 * i;

		clip_vertices << pt1 << pt2;
		clip_colors << ci << ci;

		pt1 = pt00 + vt1 * i;
		pt2 = pt10 + vt1 * i;

		clip_vertices << pt1 << pt2;
		clip_colors << ci << ci;
	}

	clip_vertices << pt0 << (pt0 + mvt.normalized() * vlen);
	clip_colors << ci << ci;
}

int ClipPlane::counter = 1;
