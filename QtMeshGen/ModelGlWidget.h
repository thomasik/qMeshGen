
#ifndef MODELGLWIDGET_H
#define MODELGLWIDGET_H

#include <QGLWidget>
#include <QDialog>
#include <QOpenGLShaderProgram>

#include <ft2build.h>
#include FT_FREETYPE_H

QT_BEGIN_NAMESPACE
class QAction;
class QPushButton;
class QCheckBox;
class QTableWidget;
QT_END_NAMESPACE

#include <memory>
#include "DRect.h"
#include "DataVector.h"

class MeshViewSet;
class ModelGLWidget;
class MeshContainer3d;

class ClipPlane
{
public:
	ClipPlane(const DPoint3d& pt, const DVector3d& vt, const QString& name);
	void setTableRow(QTableWidget* table, int row) const;
	bool clipped(const DPoint3d& pt) const;
	bool clipped(const FPoint3d& pt) const;
	void drawGL(QVector<QVector3D> & clip_vertices, QVector<GLfloat> & clip_colors, float diameter) const;
	void translate(const DVector3d& v);
	DVector3d getNormalVector() const { return m_vt; }
	DPoint3d getCenterPoint() const { return m_pt; }
	void setNormalVector(const DVector3d& v);
	void setName(const QString& name) { m_name = name; }
	void setData(const DVector3d& v, double d);
	double getParamD() const { return m_d; }
public:
	static int counter;
private:
	DPoint3d m_pt;
	DVector3d m_vt;
	double m_d;
	QString m_name;
};

class ViewPrefDialog : public QDialog
{
    Q_OBJECT

public:
    ViewPrefDialog(QWidget *parent, ModelGLWidget* glw);
	void setViewMode(int mode);

	int selected_clip_plane;
	DataVector<std::shared_ptr<ClipPlane>> clipPlanes;
	QTableWidget *clipTable;

protected:
	void closeEvent(QCloseEvent *event);

private slots:
	void redrawGL();
	void addClipPlaneOXY();
	void delSelectedClipPlane();
	void otherClipPlaneSelected();
	void applyEditedClipPlane();

private:
	ModelGLWidget* glWidget;
	int last_mode;

	QPushButton *closeButton;
	QPushButton *redrawButton;

	QCheckBox* chkViewBlocks;
	QCheckBox* chkViewFaces;
	QCheckBox* chkViewEdges;
	QCheckBox* chkViewNodes;
	QCheckBox* chkViewLabels;
	QCheckBox* chkViewWhite;
	QCheckBox* chkStoreBackFaces;
};

class ModelGLWidget : public QGLWidget
{
    Q_OBJECT

public:
    ModelGLWidget(QWidget *parent = 0);
    ~ModelGLWidget();

	enum {VIEW_BLOCKS = 1, VIEW_FACES = 2, VIEW_EDGES = 4, VIEW_NODES = 8, VIEW_WHITE = 16, STORE_BACK_FACES = 32, VIEW_LABELS = 64 };
	enum { ICOLOR_WHITE = 250, ICOLOR_GRAY = 245, ICOLOR_DGRAY = 243, ICOLOR_BLACK = 240, 
		ICOLOR_PURPLE = 239, ICOLOR_DPURPLE = 120, ICOLOR_BLUE =220, ICOLOR_DBLUE = 100, 
		ICOLOR_CYAN = 200, ICOLOR_DCYAN = 80, ICOLOR_GREEN = 180, ICOLOR_DGREEN = 60, 
		ICOLOR_YELLOW = 160, ICOLOR_DYELLOW = 40, ICOLOR_ORANGE = 140, ICOLOR_DORANGE = 20, 
		ICOLOR_RED = 121, ICOLOR_DRED = 0 };

    QSize minimumSizeHint() const;
    QSize sizeHint() const;

public:
	void setViewSet(MeshViewSet* set, bool with_reset = true);
	void setViewMesh(MeshContainer3d* mesh);
	void prepareGLScene();
	void drawGL();
	int getViewMode() const { return view_mode; }
	void setViewMode(int m);
	void redraw(bool recreate = false);
	DPoint3d getCenter() const { return m_center; }
	void markHidden();
private slots:
    void clearView();
    void resetPositioning();
    void showPrefDlg();
	void storeMatlabFile();
	void storeEPSFile();

protected:
    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);

	void dumpMatrixInfo();

	void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
	void wheelEvent(QWheelEvent *event);

private:
	void draw();

	static int genIColor(int area_id);

private:
	DBox m_bounding_box;
	DPoint3d m_center;
	MeshViewSet* m_view_set;
	MeshContainer3d* m_view_mesh;

	bool m_initialized;

	class BlockFaceRender {
	public:
		BlockFaceRender() : program(nullptr) {}
		bool init(QObject * parent = nullptr);
		void clear() { 
			face_vertices.clear(); face_normals.clear(); face_colors.clear(); 
			block_vertices.clear(); block_normals.clear(); block_colors.clear();
		}
		void renderBlocks(const QMatrix4x4& modelview, const QMatrix3x3& nm) const;
		void renderFaces(const QMatrix4x4& modelview, const QMatrix3x3& nm) const;
	public:
		void addFaceVertex(const QVector3D& pt, const QVector3D& n, unsigned char ci);
		void addFaceVertex(const FPoint3d& pt, const FVector3d& n, unsigned char ci);
		void addTriangleFace(const QVector3D& pt1, const QVector3D& pt2, const QVector3D& pt3, unsigned char ci);
		void addBlockVertex(const QVector3D& pt, const QVector3D& n, unsigned char ci);
		void addBlockVertex(const FPoint3d& pt, const FVector3d& n, unsigned char ci);
		void addBlockTriangleFace(const QVector3D& pt1, const QVector3D& pt2, const QVector3D& pt3, unsigned char ci);
	private:
		QOpenGLShaderProgram *program;
	    int colorAttr;
		int vertexAttr;
		int normalAttr;
		int matrixUniform;
		int normalMatrixUniform;
		QVector<QVector3D> face_vertices;
		QVector<QVector3D> face_normals;
		QVector<GLfloat>   face_colors;
		QVector<QVector3D> block_vertices;
		QVector<QVector3D> block_normals;
		QVector<GLfloat>   block_colors;
	} m_bf_render;

	class EdgePointRender {
		static const int EMPTY_ID = -100;
		static const int N_CX = 5;
		static const int N_CY = 10;
		static const int N_CS = 2;
	public:
		EdgePointRender() : program(nullptr), any_number(false) {}
		bool init(QObject * parent = nullptr);
		void clear() { 
			edge_vertices.clear(); edge_colors.clear(); 
			point_vertices.clear(); point_colors.clear(); point_ids.clear();
			any_number = false;
		}
		void renderEdges(const QMatrix4x4& modelview) const;
		void renderPoints(const QMatrix4x4& modelview, const QVector3D& svx, const QVector3D& svy) const;
		void renderClipPlane(const QMatrix4x4& modelview, const ClipPlane* plane, double diameter) const;
		QVector3D drawNumber(int number, const QVector3D & pos, GLfloat color, 
			const QVector3D& svx, const QVector3D& svy, 
			QVector<QVector3D> & nvertices, QVector<GLfloat> & ncolors) const;
	public:
		void addEdgeVertex(const QVector3D& pt, unsigned char ci);
		void addEdgeVertex(const FPoint3d& pt, unsigned char ci);
		void addPointVertex(const QVector3D& pt, unsigned char ci, int id = EMPTY_ID);
		void addPointVertex(const FPoint3d& pt, unsigned char ci, int id = EMPTY_ID);
	private:
		QOpenGLShaderProgram *program;
		bool any_number;
	    int colorAttr;
		int vertexAttr;
		int matrixUniform;
		QVector<QVector3D> edge_vertices;
		QVector<GLfloat>   edge_colors;
		QVector<QVector3D> point_vertices;
		QVector<GLfloat>   point_colors;
		QVector<int>	   point_ids;
	} m_ep_render;

	class LabelRender {
	public:
		LabelRender() : program(nullptr) {}
		bool init(QObject * parent = nullptr);
		void clear() { vertices.clear(); colors.clear(); labels.clear(); }
		void renderText(const char* text, float color, 
			const QVector3D& pos, const QVector3D& svx, const QVector3D& svy) const;
		void render(const QMatrix4x4& modelview, const QVector3D& svx, const QVector3D& svy) const;
	public:
		void addLabel(const QVector3D& pt, const string& s, unsigned char ci);
		void addLabel(const FPoint3d& pt, const string& s, unsigned char ci);
	private:
		FT_Library ft;
		FT_Face face;
	private:
		QOpenGLShaderProgram *program;
	    int colorAttr;
		int vertexAttr;
		int texcAttr;
		int matrixUniform;
		int texUniform;
		GLuint tex;
		QVector<QVector3D> vertices;
		QVector<GLfloat>   colors;
		QVector<string>	   labels;
	} m_label_render;

	QMatrix4x4 m_proj_matrix;

 	double clip_quality;
	int view_mode;
	ViewPrefDialog *viewPrefDlg;

	int nav_mode;
	QPoint lastPos;
	QVector3D axis_ox, axis_oy;
	GLfloat last_zoom, total_zoom;
	QVector3D last_trans, total_trans;
	QQuaternion last_rot, total_rot, axis_rot, last_axis_rot;
};

#endif
