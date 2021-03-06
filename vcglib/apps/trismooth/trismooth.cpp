// Author: Stefan Schlager
// Date: 15 September 2010
#
#include <string.h>
#include <vector>
using namespace std;
#include <stdio.h>
#include <cstddef>

// VCG headers for triangular mesh processing
#include<vcg/simplex/edge/base.h>
#include<vcg/simplex/vertex/base.h>
#include<vcg/simplex/face/base.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/edges.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/intersection.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/space/index/spatial_hashing.h>
#include <vcg/complex/algorithms/closest.h>
#include <vcg/complex/algorithms/smooth.h>

// VCG File Format Importer/Exporter
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <wrap/io_trimesh/export_ply.h>
#include <vcg/complex/algorithms/update/color.h>

using namespace vcg;

class MyFace;
class MyEdge;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>		::AsVertexType,
Use<MyEdge>			::AsEdgeType,
Use<MyFace>			::AsFaceType>{};
class MyEdge : public Edge<MyUsedTypes>{};
class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::BitFlags, vertex::Normal3f, vertex::Mark,vertex::Color4b, vertex::Qualityf>{};
class MyFace    : public Face  <MyUsedTypes, face::VertexRef,face::BitFlags,face::Mark, face::Normal3f> {};
class MyMesh : public tri::TriMesh< vector<MyVertex>, vector<MyFace > >{};
typedef MyMesh::ScalarType ScalarType;

// Uncomment only one of the two following lines to test different data structures
//typedef vcg::GridStaticPtr<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid;
//typedef vcg::SpatialHashTable<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid;

int main(int argc,char ** argv){
  char filename[256];
  int iteration = 10;
  float lambda = 0.5;
  bool lcheck=false;
  bool icheck=false, laplace=false, lapang=false, lapHC=false;
  bool noout=false;
  if (argc< 3){
    printf("\n");
    printf("    Smooth a mesh using a taubin smooth\n");
    printf("    Usage: trismooth <input.mesh> <output.ply>\n");
    printf("       <input>             any common mesh file (any common mesh file).\n");
    printf("       [-it <integer>]     smoothing iterations (default is 10).\n");
    printf("       [-l <float>]        lambda value (default is 0.5).\n");
    printf("       [--laplace]         use laplacian smoothing.\n");
    printf("       [--lapHC]           use laplacian HC smoothing.\n");
    printf("       <output>            smoothed mesh with taubin smooth applied (PLY Format).\n");
    return 0;
  }
  for (int i = 1; i < argc; i++) {
    
    if (i + 1 != argc) {// Check that we haven't finished parsing already
      if (strcmp("-l", argv[i]) == 0) {
	lambda = atof(argv[i + 1]);
	if (i==argc-2)
	  noout=true;
      }
      if (strcmp("-it", argv[i]) == 0) {
	iteration=atoi(argv[i + 1]);
	if (i== (argc-2))
	  noout=true;
      }
    }
  }
  for (int i = 1; i < argc; i++) {
    if (i != argc) {
      if (strcmp("--laplace", argv[i]) == 0) {
	laplace=true;
	if (i== (argc-1))
	  noout=true;
      }		  
    }
    if (strcmp("--lapang", argv[i]) == 0) {
      lapang=true;
      if (i== (argc-1))
	noout=true;
    }
    if (strcmp("--lapHC", argv[i]) == 0) {
      lapHC=true;
      if (i== (argc-1))
	noout=true;
    }
  }
  if (noout==true) {
    printf("Error: please specify output file\n");
    return 0;
  }
  strcpy(filename, argv[argc-1]);
  MyMesh mesh;
  
  //--------------------------------------------------------------------------------------//
  //
  //                                   OPENING THE FILES  
  //
  //--------------------------------------------------------------------------------------//

  int err2 = tri::io::Importer<MyMesh>::Open(mesh,argv[1]);
  if(err2) {
    printf("Error in reading %s: '%s'\n",argv[1],tri::io::Importer<MyMesh>::ErrorMsg(err2));
    exit(-1);  
  }
  
  tri::UpdateBounding<MyMesh>::Box(mesh);
  //tri::UpdateNormals<MyMesh>::PerFaceNormalized(mesh);
  //tri::UpdateNormals<MyMesh>::PerVertexAngleWeighted(mesh);
  //tri::UpdateNormals<MyMesh>::NormalizeVertex(mesh);
  if (laplace == true) {
    tri::Smooth<MyMesh>::VertexCoordLaplacian(mesh,iteration,false);
  } else if (lapang == true) {
    // tri::Smooth<MyMesh>::VertexCoordLaplacianCurvatureFlow(mesh,iteration, lambda);
    tri::Smooth<MyMesh>::VertexCoordLaplacianAngleWeighted(mesh,iteration,lambda);
  } else if (lapHC == true) {
    // tri::Smooth<MyMesh>::VertexCoordLaplacianCurvatureFlow(mesh,iteration, lambda);
    tri::Smooth<MyMesh>::VertexCoordLaplacianHC(mesh,iteration);
  } else {
    tri::Smooth<MyMesh>::VertexCoordTaubin(mesh,iteration,lambda,-0.53);
  }  
  //tri::UpdateNormals<MyMesh>::PerFaceNormalized(mesh);
  tri::UpdateNormals<MyMesh>::PerVertexAngleWeighted(mesh);
  tri::UpdateNormals<MyMesh>::NormalizeVertex(mesh);
  tri::io::ExporterPLY<MyMesh>::Save(mesh,filename,tri::io::Mask::IOM_VERTNORMAL, false); // in ASCII
  return 0;
}
