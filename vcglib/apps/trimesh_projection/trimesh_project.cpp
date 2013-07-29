// Author: Andrea Tagliasacchi
// Date: 22 May 2010
// modified in the follwing months by Stefan Schlager

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
#include <vcg/complex/algorithms/smooth.h>
#include <vcg/complex/append.h>

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
typedef vcg::SpatialHashTable<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid;

int main(int argc,char ** argv){
  bool nosmooth = false, sign = false;
         char filename[256] = "out_cloud.ply";
    if (argc<3){

        printf("\n");
        printf("    Compute a projection of a point cloud onto a mesh\n");
        printf("    Usage: trimesh_project <cloud> <mesh> [options]\n");

        printf("       <cloud>        Point cloud over which to project (PLY|STL|OBJ format).\n");

        printf("       <mesh>         Mesh model over which to project (PLY|STL|OBJ format).\n");
        printf("       --nosmooth     don't smooth vertex normals.\n");
	printf("       --sign         return signed distances.\n");
        printf("       -o             specify output file (PLY format).\n");
        printf("\n");
        printf("Saves the projected point cloud in the file out_cloud.ply when terminated.\n");
        printf("The vertex coordinates represent the projected samples (in same order as input.\n");
        printf("The vertex quality (requires ply files) represents the projection distance.\n");
		return 0;
	}


    for (int i = 1; i < argc; i++) {
  if (i + 1 != argc) {// Check simple switches

                        if (strcmp("-o", argv[i]) == 0)
                        {strcpy(filename, argv[i+1]);;
                        }
                    }
            }
    for (int i = 1; i < argc; i++) {
      if (i  != argc) {// Check options with arguments

                        if (strcmp("--nosmooth", argv[i]) == 0)
                        {nosmooth=true;
                        }
			if (strcmp("--sign", argv[i]) == 0)
                        {sign=true;
                        }
            }
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
  //tri::Allocator<MyMesh>::AddVertices(out_cloud,in_cloud.vn)
  vcg::tri::Append<MyMesh,MyMesh>::Mesh(out_cloud,in_cloud);

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

  //if (sign == true)
  //  {
  //   tri::UpdateNormals<MyMesh>::PerFaceNormalized(in_cloud);
  //    tri::UpdateNormals<MyMesh>::PerVertexAngleWeighted(in_cloud);
  //    tri::UpdateNormals<MyMesh>::NormalizeVertex(in_cloud);
  // 
  // }
  if (nosmooth == false)
  {
  tri::Smooth<MyMesh>::VertexNormalLaplacian(mesh,2,false);
  }
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
  //tri::UpdateFlags<MyMesh>::FaceProjection(mesh);
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
    int f_i = vcg::tri::Index(mesh, f_ptr);
    out_cloud.vert[i].Q() = minDist;
    MyMesh::CoordType tt = currp*0;
     for (int j=0; j <3;j++)
	      {
		if (&(mesh.face[f_i].V(j)->N()))
		  {
		    Point3f vdist = mesh.face[f_i].V(j)->P() - clost;
		    float weight = sqrt(vdist.dot(vdist));
		    if (weight > 0)
		      weight = 1/weight;
		    else 
		      weight = 1e12;
		      
		    tt +=(mesh.face[f_i].V(j)->N()*weight);
		  }
	      }
    float vl = sqrt(tt.dot(tt));
	    if (vl > 0 && &vl)//check for zero length normals
	      {
		tt=tt/vl;
	      }   
    
    
    if (sign == true)
      {
	Point3f dif = clost - currp;
	float signo = dif.dot(tt);	
	if (signo < 0)
	  { out_cloud.vert[i].Q() = -out_cloud.vert[i].Q() ;
	  }
	
      }
    out_cloud.vert[i].N() = tt;
  }

  int t2 = clock();
  tri::io::ExporterPLY<MyMesh>::Save(out_cloud,filename,tri::io::Mask::IOM_VERTNORMAL + tri::io::Mask::IOM_VERTQUALITY, false); // in ASCII
  printf("Completed projection of %d sample in %i msec\n", in_cloud.vn, t2-t1);


  //printf("%i ",tt);
  return 0;
}
