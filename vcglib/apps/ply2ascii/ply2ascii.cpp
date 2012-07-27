// Author: Stefan Schlager
// Date: 22 Juli 2010
#
#include <string.h>
#include <vector>
#include <cstddef>
using namespace std;

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
typedef vcg::GridStaticPtr<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid;
//typedef vcg::SpatialHashTable<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid;

int main(int argc,char ** argv){
	char filename[256];
        bool noclean = false, col =false;
  if (argc< 3){
		printf("\n");
    printf("    Convert a mesh to ascii ply format\n");
    printf("    Usage: ply2ascii <input.mesh> <output.ply> [option]\n");
    printf("       <mesh>        any common mesh file (any common mesh file).\n");
    printf("       <output.ply>      cleaned mesh with updated vertex normals saved in ascii PLY format.\n");
    printf("       --noclean      skip cleaning duplicate and unreferenced vertices.\n");
    printf("       --color        export vertex colors - if none are present, vertices will be colored white.\n");


    
   
		return 0;
	}
  for (int i = 1; i < argc; i++) {


                  if (strcmp("--noclean", argv[i]) == 0)
                  {
                    noclean = true;
                  }
 if (strcmp("--color", argv[i]) == 0)
                  {
                    col = true;
		  }
  }



    strcpy(filename, argv[2]);

	
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
if (noclean == false)
{
  int dup = tri::Clean<MyMesh>::RemoveDuplicateVertex(mesh);
  int unref =  tri::Clean<MyMesh>::RemoveUnreferencedVertex(mesh);
  if (dup > 0 || unref > 0)
    {
                printf("Removed %i duplicate and %i unreferenced vertices from mesh %s\n",dup,unref,argv[2]);
    }
}
/* if (col == true && noclean ==false)
   {
     int unref =  tri::Clean<MyMesh>::RemoveUnreferencedVertex(mesh);
     if (unref > 0)
       printf("Removed %i unreferenced vertices from mesh %s\n",unref,argv[2]);
       }*/
 
  tri::UpdateBounding<MyMesh>::Box(mesh);
  tri::UpdateNormals<MyMesh>::PerFaceNormalized(mesh);
  tri::UpdateNormals<MyMesh>::PerVertexAngleWeighted(mesh);
  tri::UpdateNormals<MyMesh>::NormalizeVertex(mesh);


  //--------------------------------------------------------------------------------------//
  int mask0 = tri::io::Mask::IOM_VERTNORMAL;
  if (col == true)
    {
      mask0  = mask0+tri::io::Mask::IOM_VERTCOLOR;
    }
  
  tri::io::ExporterPLY<MyMesh>::Save(mesh,filename,mask0, false); // in ASCII
  
 

  return 0;
}
