TARGET = trimesh_project
DEPENDPATH += . 
INCLUDEPATH += . ../..
CONFIG += console stl
TEMPLATE = app
SOURCES += trimesh_project.cpp ../../wrap/ply/plylib.cpp
QMAKE_LFLAGS += -O2
# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle
