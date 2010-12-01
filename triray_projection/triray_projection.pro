TARGET = triray_project
DEPENDPATH += . 
INCLUDEPATH += . ../..
CONFIG += console stl 
TEMPLATE = app
SOURCES += main.cpp ../../wrap/ply/plylib.cpp

# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle
