TARGET = mesh_closest
DEPENDPATH += . 
INCLUDEPATH += . ../..
CONFIG += console stl
TEMPLATE = app
SOURCES += main.cpp ../../wrap/ply/plylib.cpp
QMAKE_LFLAGS += -O2
# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle
