// Author: Andrea Tagliasacchi
// Date: 22 May 2010
// modified in the follwing months by Stefan Schlager

#include <vector>
using namespace std;

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
#include <vcg/complex/trimesh/smooth.h>

#include <vcg/complex/intersection.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/space/index/spatial_hashing.h>
#include <vcg/complex/trimesh/closest.h>

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
    if (argc<3){

        printf("\n");
        printf("    Compute a projection of the vertices of a reference mesh onto a target mesh and saves in ASCII PLY Format\n");
        printf("    Usage: mesh_closest <mesh> <mesh> <output.ply> \n");
        printf("       <mesh>        Reference mesh\n");
        printf("       <mesh>        Target mesh");
        printf("       <output.ply>      Optional: name of output (standard is project.mesh.ply)");
        printf("\n");
        printf("Saves the projected mesh cloud in when terminated.\n");
        printf("The vertex coordinates represent the projected samples (in same order as input.\n");
        printf("The vertex quality (requires ply files) represents the projection distance.\n");
		return 0;
	}

    else if (argc == 3)
            {
            strcpy(filename, "out_cloud.ply");
            }
            else
            {
            strcpy(filename, argv[3]);
            }
  
	MyMesh mesh;
  MyMesh in_cloud;
  MyMesh out_cloud;
	
  //--------------------------------------------------------------------------------------//
  //
  //                                   OPENING THE FILES  
  //
  //--------------------------------------------------------------------------------------//
  int err = tri::io::Importer<MyMesh>::Open(in_cloud,argv[1]);
  if(err) {
		printf("Error in reading %s: '%s'\n",argv[1],tri::io::Importer<MyMesh>::ErrorMsg(err));
		exit(-1);  
	}
  int err2 = tri::io::Importer<MyMesh>::Open(mesh,argv[2]);
  if(err2) {
		printf("Error in reading %s: '%s'\n",argv[2],tri::io::Importer<MyMesh>::ErrorMsg(err2));
		exit(-1);  
	}
  // Allocate space for projected cloud
  tri::Allocator<MyMesh>::AddVertices(out_cloud,in_cloud.vn);

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
                printf("Removed %i duplicate and %i unreferenced vertices from mesh %s\n",dup,unref,argv[2]);

  tri::UpdateBounding<MyMesh>::Box(mesh);
  tri::UpdateNormals<MyMesh>::PerFaceNormalized(mesh);
  tri::UpdateNormals<MyMesh>::PerVertexAngleWeighted(mesh);
  tri::Smooth<MyMesh>::VertexNormalLaplacian(mesh,2,false);

  tri::UpdateNormals<MyMesh>::NormalizeVertex(mesh);
  tri::UpdateQuality<MyMesh>::VertexConstant(out_cloud, 0);  
  float maxDist = mesh.bbox.Diag();
  float minDist = 1e-10;


  //--------------------------------------------------------------------------------------//
  //
  //                              INITIALIZE SEARCH STRUCTURES
  //
  // Update the FaceProjection flags needed for projection/distance queries
  // Create a static grid (for fast indexing) and fill it 
  //--------------------------------------------------------------------------------------//
  vcg::tri::FaceTmark<MyMesh> mf; 
  mf.SetMesh( &mesh );
  vcg::face::PointDistanceBaseFunctor<float> PDistFunct;
  tri::UpdateFlags<MyMesh>::FaceProjection(mesh);
	TriMeshGrid static_grid;
  printf("preprocessing mesh with %d faces\n", mesh.fn);
	static_grid.Set(mesh.face.begin(), mesh.face.end());

  
  //--------------------------------------------------------------------------------------//
  //
  //                                PERFORM DISTANCE QUERIES
  //
  //--------------------------------------------------------------------------------------//
  int t1=clock();
  for(int i=0; i<in_cloud.vn; i++){
    Point3f& currp = in_cloud.vert[i].P();
    Point3f& clost = out_cloud.vert[i].P();
    MyFace* f_ptr= GridClosest(static_grid, PDistFunct, mf, currp, maxDist, minDist, clost);
    in_cloud.vert[i].Q() = minDist;
    in_cloud.vert[i].P()=out_cloud.vert[i].P();
  }
  tri::UpdateBounding<MyMesh>::Box(in_cloud);
  tri::UpdateNormals<MyMesh>::PerFaceNormalized(in_cloud);
  tri::Smooth<MyMesh>::VertexNormalLaplacian(in_cloud,2,false);
  tri::UpdateNormals<MyMesh>::NormalizeVertex(in_cloud);

  int t2 = clock();
  tri::io::ExporterPLY<MyMesh>::Save(in_cloud,filename,tri::io::Mask::IOM_VERTNORMAL + tri::io::Mask::IOM_VERTQUALITY, false); // in ASCII
  printf("Completed projection of %d sample in %i msec\n", in_cloud.vn, t2-t1);


  //printf("%i ",tt);
  return 0;
}
