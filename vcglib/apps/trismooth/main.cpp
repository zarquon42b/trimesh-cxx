// Author: Stefan Schlager
// Date: 15 September 2010
#
#include <string.h>
#include <vector>
using namespace std;
#include <stdio.h>

// VCG headers for triangular mesh processing
#include<vcg/simplex/edge/base.h>
#include<vcg/simplex/vertex/base.h>
#include<vcg/simplex/face/base.h>
#include <vcg/complex/trimesh/base.h>
#include <vcg/complex/trimesh/update/topology.h>
#include <vcg/complex/trimesh/update/edges.h>
#include <vcg/complex/trimesh/update/bounding.h>
#include <vcg/complex/trimesh/update/quality.h>

#include <vcg/complex/trimesh/update/flag.h>
#include <vcg/complex/trimesh/clean.h>
#include <vcg/complex/intersection.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/space/index/spatial_hashing.h>
#include <vcg/complex/trimesh/closest.h>
#include <vcg/complex/trimesh/smooth.h>

// VCG File Format Importer/Exporter
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <wrap/io_trimesh/export_ply.h>
#include <vcg/complex/trimesh/update/color.h>

using namespace vcg;

class MyFace;
class MyEdge;
class MyVertex;
struct MyUsedTypes : public UsedTypes<	Use<MyVertex>		::AsVertexType,
Use<MyEdge>			::AsEdgeType,
Use<MyFace>			::AsFaceType>{};
class MyVertex  : public Vertex< MyUsedTypes, vertex::Coord3f, vertex::BitFlags, vertex::Normal3f, vertex::Mark,vertex::Color4b, vertex::Qualityf>{};
class MyFace    : public Face  <MyUsedTypes, face::VertexRef,face::BitFlags,face::Mark, face::Normal3f> {};
class MyMesh : public tri::TriMesh< vector<MyVertex>, vector<MyFace > >{};
typedef MyMesh::ScalarType ScalarType;

// Uncomment only one of the two following lines to test different data structures
typedef vcg::GridStaticPtr<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid;
//typedef vcg::SpatialHashTable<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid;

int main(int argc,char ** argv){
        char filename[256];
        int iteration;
        float lambda;
        bool lcheck=false;
        bool icheck=false;
        bool noout=false;
  if (argc< 3){
          printf("\n");
          printf("    Smooth a mesh using a taubin smooth\n");
          printf("    Usage: trismooth <input.mesh> <output.ply>\n");
          printf("       <input>        any common mesh file (any common mesh file).\n");
          printf("       [-it <integer>]     smoothing iterations (default is 10).\n");
          printf("       [-l <float>]        lambda value (default is 0.5).\n");
          printf("       <output>            smoothed mesh with taubin smooth applied (PLY Format).\n");

   
		return 0;
	}
        for (int i = 1; i < argc; i++) {

            if (i + 1 != argc) {// Check that we haven't finished parsing already

                      if (strcmp("-l", argv[i]) == 0) {
                      lambda = atof(argv[i + 1]);
                      lcheck=true;
                      if (i==argc-2)
                      {noout=true;}
                  }


                  if (strcmp("-it", argv[i]) == 0) {
                        iteration=atoi(argv[i + 1]);
                        icheck=true;
                        if (i== (argc-2))
                        {noout=true;
                        //printf("please specify output file");
                        }

                  }
                }
             }
      //  printf("%i\n",argc);
        if (lcheck == false)
        {lambda=0.5;
        }
        if (icheck == false)
        {iteration = 10;
        }
        if (noout==true)
        {printf("Error: please specify output file\n");
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

   //--------------------------------------------------------------------------------------//
  //
  //                                   PREPROCESS
  //
  // Update the bounding box and initialize max search distance
  // Remove duplicates and update mesh properties
  //--------------------------------------------------------------------------------------//

  int dup = tri::Clean<MyMesh>::RemoveDuplicateVertex(mesh);
        int unref =  tri::Clean<MyMesh>::RemoveUnreferencedVertex(mesh);
  if (dup > 0 || unref > 0)
                printf("Removed %i duplicate and %i unreferenced vertices from mesh %s\n",dup,unref,argv[1]);
  tri::UpdateBounding<MyMesh>::Box(mesh);
  tri::UpdateNormals<MyMesh>::PerFaceNormalized(mesh);
  //tri::UpdateNormals<MyMesh>::PerVertexAngleWeighted(mesh);
  tri::UpdateNormals<MyMesh>::NormalizeVertex(mesh);
  tri::Smooth<MyMesh>::VertexCoordTaubin(mesh,iteration,lambda,-0.53);
  tri::UpdateNormals<MyMesh>::PerFaceNormalized(mesh);
  tri::UpdateNormals<MyMesh>::PerVertexAngleWeighted(mesh);
  tri::UpdateNormals<MyMesh>::NormalizeVertex(mesh);



  //--------------------------------------------------------------------------------------//
 

  

  tri::io::ExporterPLY<MyMesh>::Save(mesh,filename,tri::io::Mask::IOM_VERTNORMAL, false); // in ASCII



  return 0;
}
