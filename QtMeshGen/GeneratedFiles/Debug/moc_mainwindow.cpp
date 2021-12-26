/****************************************************************************
** Meta object code from reading C++ file 'mainwindow.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.9.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../mainwindow.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.9.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_MeshKernelThread_t {
    QByteArrayData data[1];
    char stringdata0[17];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_MeshKernelThread_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_MeshKernelThread_t qt_meta_stringdata_MeshKernelThread = {
    {
QT_MOC_LITERAL(0, 0, 16) // "MeshKernelThread"

    },
    "MeshKernelThread"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_MeshKernelThread[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       0,    0, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

       0        // eod
};

void MeshKernelThread::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    Q_UNUSED(_o);
    Q_UNUSED(_id);
    Q_UNUSED(_c);
    Q_UNUSED(_a);
}

const QMetaObject MeshKernelThread::staticMetaObject = {
    { &QThread::staticMetaObject, qt_meta_stringdata_MeshKernelThread.data,
      qt_meta_data_MeshKernelThread,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *MeshKernelThread::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *MeshKernelThread::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_MeshKernelThread.stringdata0))
        return static_cast<void*>(const_cast< MeshKernelThread*>(this));
    return QThread::qt_metacast(_clname);
}

int MeshKernelThread::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QThread::qt_metacall(_c, _id, _a);
    return _id;
}
struct qt_meta_stringdata_MeshViewDialog_t {
    QByteArrayData data[5];
    char stringdata0[58];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_MeshViewDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_MeshViewDialog_t qt_meta_stringdata_MeshViewDialog = {
    {
QT_MOC_LITERAL(0, 0, 14), // "MeshViewDialog"
QT_MOC_LITERAL(1, 15, 15), // "continuePressed"
QT_MOC_LITERAL(2, 31, 0), // ""
QT_MOC_LITERAL(3, 32, 14), // "skipAllPressed"
QT_MOC_LITERAL(4, 47, 10) // "logPressed"

    },
    "MeshViewDialog\0continuePressed\0\0"
    "skipAllPressed\0logPressed"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_MeshViewDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   29,    2, 0x08 /* Private */,
       3,    0,   30,    2, 0x08 /* Private */,
       4,    0,   31,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void MeshViewDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        MeshViewDialog *_t = static_cast<MeshViewDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->continuePressed(); break;
        case 1: _t->skipAllPressed(); break;
        case 2: _t->logPressed(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject MeshViewDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_MeshViewDialog.data,
      qt_meta_data_MeshViewDialog,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *MeshViewDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *MeshViewDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_MeshViewDialog.stringdata0))
        return static_cast<void*>(const_cast< MeshViewDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int MeshViewDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 3)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 3;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 3)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 3;
    }
    return _id;
}
struct qt_meta_stringdata_MeshPropertiesDialog_t {
    QByteArrayData data[3];
    char stringdata0[28];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_MeshPropertiesDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_MeshPropertiesDialog_t qt_meta_stringdata_MeshPropertiesDialog = {
    {
QT_MOC_LITERAL(0, 0, 20), // "MeshPropertiesDialog"
QT_MOC_LITERAL(1, 21, 5), // "apply"
QT_MOC_LITERAL(2, 27, 0) // ""

    },
    "MeshPropertiesDialog\0apply\0"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_MeshPropertiesDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   19,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,

       0        // eod
};

void MeshPropertiesDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        MeshPropertiesDialog *_t = static_cast<MeshPropertiesDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->apply(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject MeshPropertiesDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_MeshPropertiesDialog.data,
      qt_meta_data_MeshPropertiesDialog,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *MeshPropertiesDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *MeshPropertiesDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_MeshPropertiesDialog.stringdata0))
        return static_cast<void*>(const_cast< MeshPropertiesDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int MeshPropertiesDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 1)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 1;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 1)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 1;
    }
    return _id;
}
struct qt_meta_stringdata_ConsoleWidget_t {
    QByteArrayData data[11];
    char stringdata0[128];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_ConsoleWidget_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_ConsoleWidget_t qt_meta_stringdata_ConsoleWidget = {
    {
QT_MOC_LITERAL(0, 0, 13), // "ConsoleWidget"
QT_MOC_LITERAL(1, 14, 14), // "addMeshCommand"
QT_MOC_LITERAL(2, 29, 0), // ""
QT_MOC_LITERAL(3, 30, 3), // "cmd"
QT_MOC_LITERAL(4, 34, 17), // "enableApplyButton"
QT_MOC_LITERAL(5, 52, 4), // "text"
QT_MOC_LITERAL(6, 57, 12), // "applyClicked"
QT_MOC_LITERAL(7, 70, 14), // "cmdListClicked"
QT_MOC_LITERAL(8, 85, 16), // "QListWidgetItem*"
QT_MOC_LITERAL(9, 102, 4), // "item"
QT_MOC_LITERAL(10, 107, 20) // "cmdListDoubleClicked"

    },
    "ConsoleWidget\0addMeshCommand\0\0cmd\0"
    "enableApplyButton\0text\0applyClicked\0"
    "cmdListClicked\0QListWidgetItem*\0item\0"
    "cmdListDoubleClicked"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_ConsoleWidget[] = {

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
       1,    1,   39,    2, 0x0a /* Public */,
       4,    1,   42,    2, 0x08 /* Private */,
       6,    0,   45,    2, 0x08 /* Private */,
       7,    1,   46,    2, 0x08 /* Private */,
      10,    1,   49,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void, QMetaType::QString,    3,
    QMetaType::Void, QMetaType::QString,    5,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 8,    9,
    QMetaType::Void, 0x80000000 | 8,    9,

       0        // eod
};

void ConsoleWidget::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        ConsoleWidget *_t = static_cast<ConsoleWidget *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->addMeshCommand((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 1: _t->enableApplyButton((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 2: _t->applyClicked(); break;
        case 3: _t->cmdListClicked((*reinterpret_cast< QListWidgetItem*(*)>(_a[1]))); break;
        case 4: _t->cmdListDoubleClicked((*reinterpret_cast< QListWidgetItem*(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObject ConsoleWidget::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_ConsoleWidget.data,
      qt_meta_data_ConsoleWidget,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *ConsoleWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *ConsoleWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_ConsoleWidget.stringdata0))
        return static_cast<void*>(const_cast< ConsoleWidget*>(this));
    return QWidget::qt_metacast(_clname);
}

int ConsoleWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
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
struct qt_meta_stringdata_PropertiesWidget_t {
    QByteArrayData data[5];
    char stringdata0[48];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_PropertiesWidget_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_PropertiesWidget_t qt_meta_stringdata_PropertiesWidget = {
    {
QT_MOC_LITERAL(0, 0, 16), // "PropertiesWidget"
QT_MOC_LITERAL(1, 17, 18), // "updateMeshProperty"
QT_MOC_LITERAL(2, 36, 0), // ""
QT_MOC_LITERAL(3, 37, 3), // "row"
QT_MOC_LITERAL(4, 41, 6) // "column"

    },
    "PropertiesWidget\0updateMeshProperty\0"
    "\0row\0column"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_PropertiesWidget[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    2,   19,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void, QMetaType::Int, QMetaType::Int,    3,    4,

       0        // eod
};

void PropertiesWidget::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        PropertiesWidget *_t = static_cast<PropertiesWidget *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->updateMeshProperty((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        default: ;
        }
    }
}

const QMetaObject PropertiesWidget::staticMetaObject = {
    { &QDockWidget::staticMetaObject, qt_meta_stringdata_PropertiesWidget.data,
      qt_meta_data_PropertiesWidget,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *PropertiesWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *PropertiesWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_PropertiesWidget.stringdata0))
        return static_cast<void*>(const_cast< PropertiesWidget*>(this));
    return QDockWidget::qt_metacast(_clname);
}

int PropertiesWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDockWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 1)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 1;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 1)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 1;
    }
    return _id;
}
struct qt_meta_stringdata_MainWindow_t {
    QByteArrayData data[17];
    char stringdata0[224];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_MainWindow_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_MainWindow_t qt_meta_stringdata_MainWindow = {
    {
QT_MOC_LITERAL(0, 0, 10), // "MainWindow"
QT_MOC_LITERAL(1, 11, 4), // "open"
QT_MOC_LITERAL(2, 16, 0), // ""
QT_MOC_LITERAL(3, 17, 14), // "openRecentFile"
QT_MOC_LITERAL(4, 32, 7), // "saveGRD"
QT_MOC_LITERAL(5, 40, 7), // "saveTXT"
QT_MOC_LITERAL(6, 48, 5), // "about"
QT_MOC_LITERAL(7, 54, 15), // "updateModelInfo"
QT_MOC_LITERAL(8, 70, 16), // "showMeshItemView"
QT_MOC_LITERAL(9, 87, 16), // "showMeshItemStat"
QT_MOC_LITERAL(10, 104, 13), // "clearMeshItem"
QT_MOC_LITERAL(11, 118, 24), // "statTreeSelectionChanged"
QT_MOC_LITERAL(12, 143, 17), // "triangulateAction"
QT_MOC_LITERAL(13, 161, 10), // "quadAction"
QT_MOC_LITERAL(14, 172, 16), // "tetrahedraAction"
QT_MOC_LITERAL(15, 189, 16), // "parseGrainAction"
QT_MOC_LITERAL(16, 206, 17) // "meshingPropAction"

    },
    "MainWindow\0open\0\0openRecentFile\0saveGRD\0"
    "saveTXT\0about\0updateModelInfo\0"
    "showMeshItemView\0showMeshItemStat\0"
    "clearMeshItem\0statTreeSelectionChanged\0"
    "triangulateAction\0quadAction\0"
    "tetrahedraAction\0parseGrainAction\0"
    "meshingPropAction"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_MainWindow[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      15,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   89,    2, 0x08 /* Private */,
       3,    0,   90,    2, 0x08 /* Private */,
       4,    0,   91,    2, 0x08 /* Private */,
       5,    0,   92,    2, 0x08 /* Private */,
       6,    0,   93,    2, 0x08 /* Private */,
       7,    0,   94,    2, 0x08 /* Private */,
       8,    0,   95,    2, 0x08 /* Private */,
       9,    0,   96,    2, 0x08 /* Private */,
      10,    0,   97,    2, 0x08 /* Private */,
      11,    0,   98,    2, 0x08 /* Private */,
      12,    0,   99,    2, 0x08 /* Private */,
      13,    0,  100,    2, 0x08 /* Private */,
      14,    0,  101,    2, 0x08 /* Private */,
      15,    0,  102,    2, 0x08 /* Private */,
      16,    0,  103,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void MainWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        MainWindow *_t = static_cast<MainWindow *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->open(); break;
        case 1: _t->openRecentFile(); break;
        case 2: _t->saveGRD(); break;
        case 3: _t->saveTXT(); break;
        case 4: _t->about(); break;
        case 5: _t->updateModelInfo(); break;
        case 6: _t->showMeshItemView(); break;
        case 7: _t->showMeshItemStat(); break;
        case 8: _t->clearMeshItem(); break;
        case 9: _t->statTreeSelectionChanged(); break;
        case 10: _t->triangulateAction(); break;
        case 11: _t->quadAction(); break;
        case 12: _t->tetrahedraAction(); break;
        case 13: _t->parseGrainAction(); break;
        case 14: _t->meshingPropAction(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject MainWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_MainWindow.data,
      qt_meta_data_MainWindow,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *MainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *MainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_MainWindow.stringdata0))
        return static_cast<void*>(const_cast< MainWindow*>(this));
    if (!strcmp(_clname, "MeshViewExt"))
        return static_cast< MeshViewExt*>(const_cast< MainWindow*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int MainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 15)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 15;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 15)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 15;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
