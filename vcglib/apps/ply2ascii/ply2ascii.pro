TARGET = ply2ascii
DEPENDPATH += . 
INCLUDEPATH += . ../..
CONFIG += console stl 
TEMPLATE = app
SOURCES += ply2ascii.cpp ../../wrap/ply/plylib.cpp
QMAKE_LFLAGS += -O2
# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle
