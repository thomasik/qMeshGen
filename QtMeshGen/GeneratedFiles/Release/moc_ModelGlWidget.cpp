/****************************************************************************
** Meta object code from reading C++ file 'ModelGlWidget.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.9.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../ModelGlWidget.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'ModelGlWidget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.9.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_ViewPrefDialog_t {
    QByteArrayData data[7];
    char stringdata0[106];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_ViewPrefDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_ViewPrefDialog_t qt_meta_stringdata_ViewPrefDialog = {
    {
QT_MOC_LITERAL(0, 0, 14), // "ViewPrefDialog"
QT_MOC_LITERAL(1, 15, 8), // "redrawGL"
QT_MOC_LITERAL(2, 24, 0), // ""
QT_MOC_LITERAL(3, 25, 15), // "addClipPlaneOXY"
QT_MOC_LITERAL(4, 41, 20), // "delSelectedClipPlane"
QT_MOC_LITERAL(5, 62, 22), // "otherClipPlaneSelected"
QT_MOC_LITERAL(6, 85, 20) // "applyEditedClipPlane"

    },
    "ViewPrefDialog\0redrawGL\0\0addClipPlaneOXY\0"
    "delSelectedClipPlane\0otherClipPlaneSelected\0"
    "applyEditedClipPlane"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_ViewPrefDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       5,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   39,    2, 0x08 /* Private */,
       3,    0,   40,    2, 0x08 /* Private */,
       4,    0,   41,    2, 0x08 /* Private */,
       5,    0,   42,    2, 0x08 /* Private */,
       6,    0,   43,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void ViewPrefDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        ViewPrefDialog *_t = static_cast<ViewPrefDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->redrawGL(); break;
        case 1: _t->addClipPlaneOXY(); break;
        case 2: _t->delSelectedClipPlane(); break;
        case 3: _t->otherClipPlaneSelected(); break;
        case 4: _t->applyEditedClipPlane(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject ViewPrefDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_ViewPrefDialog.data,
      qt_meta_data_ViewPrefDialog,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *ViewPrefDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *ViewPrefDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_ViewPrefDialog.stringdata0))
        return static_cast<void*>(const_cast< ViewPrefDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int ViewPrefDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 5)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 5;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 5)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 5;
    }
    return _id;
}
struct qt_meta_stringdata_ModelGLWidget_t {
    QByteArrayData data[7];
    char stringdata0[83];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_ModelGLWidget_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_ModelGLWidget_t qt_meta_stringdata_ModelGLWidget = {
    {
QT_MOC_LITERAL(0, 0, 13), // "ModelGLWidget"
QT_MOC_LITERAL(1, 14, 9), // "clearView"
QT_MOC_LITERAL(2, 24, 0), // ""
QT_MOC_LITERAL(3, 25, 16), // "resetPositioning"
QT_MOC_LITERAL(4, 42, 11), // "showPrefDlg"
QT_MOC_LITERAL(5, 54, 15), // "storeMatlabFile"
QT_MOC_LITERAL(6, 70, 12) // "storeEPSFile"

    },
    "ModelGLWidget\0clearView\0\0resetPositioning\0"
    "showPrefDlg\0storeMatlabFile\0storeEPSFile"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_ModelGLWidget[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       5,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   39,    2, 0x08 /* Private */,
       3,    0,   40,    2, 0x08 /* Private */,
       4,    0,   41,    2, 0x08 /* Private */,
       5,    0,   42,    2, 0x08 /* Private */,
       6,    0,   43,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void ModelGLWidget::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        ModelGLWidget *_t = static_cast<ModelGLWidget *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->clearView(); break;
        case 1: _t->resetPositioning(); break;
        case 2: _t->showPrefDlg(); break;
        case 3: _t->storeMatlabFile(); break;
        case 4: _t->storeEPSFile(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject ModelGLWidget::staticMetaObject = {
    { &QGLWidget::staticMetaObject, qt_meta_stringdata_ModelGLWidget.data,
      qt_meta_data_ModelGLWidget,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *ModelGLWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *ModelGLWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_ModelGLWidget.stringdata0))
        return static_cast<void*>(const_cast< ModelGLWidget*>(this));
    return QGLWidget::qt_metacast(_clname);
}

int ModelGLWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QGLWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 5)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 5;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 5)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 5;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
