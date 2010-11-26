// Author: Andrea Tagliasacchi
// Date: 22 May 2010
#include <string.h>
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
#include <vcg/complex/trimesh/update/color.h>
#include <vcg/complex/trimesh/update/flag.h>
#include <vcg/complex/trimesh/clean.h>
#include <vcg/complex/intersection.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/space/index/spatial_hashing.h>
#include <vcg/complex/trimesh/closest.h>

// VCG File Format Importer/Exporter
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <wrap/io_trimesh/export_ply.h>

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
    printf("    Compute a projection of a point cloud onto a mesh\n");
    printf("    Usage: trimesh_project <cloud> <mesh>\n");
    printf("       <cloud>        Point cloud over which to project (PLY format).\n");
    printf("       <mesh>         Mesh model over which to project (PLY format).\n");
    printf("\n");    
    printf("Saves the projected point cloud in the file out_cloud.ply when terminated.\n");
    printf("The vertex coordinates represent the projected samples (in same order as input.\n");
    printf("The vertex quality (requires ply files) represents the projection distance.\n");
		return 0;
	}
    strcpy(filename, argv[2]);
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
  tri::UpdateNormals<MyMesh>::NormalizeVertex(mesh);
  //tri::UpdateBounding<MyMesh>::Box(mesh);
  tri::UpdateNormals<MyMesh>::NormalizeVertex(in_cloud);
  tri::UpdateQuality<MyMesh>::VertexConstant(out_cloud, 0);
  float maxDist = mesh.bbox.Diag();
  float minDist = 1e-10;
  float t;

  //--------------------------------------------------------------------------------------//
  //
  //                              INITIALIZE SEARCH STRUCTURES
  //
  // Update the FaceProjection flags needed for projection/distance queries
  // Create a static grid (for fast indexing) and fill it 
  //--------------------------------------------------------------------------------------//
    vcg::tri::FaceTmark<MyMesh> mf;
    mf.SetMesh( &mesh );
    vcg::RayTriangleIntersectionFunctor<true> FintFunct;
    vcg::face::PointDistanceBaseFunctor<float> PDistFunct;
    tri::UpdateFlags<MyMesh>::FaceProjection(mesh);
          TriMeshGrid static_grid;
    printf("preprocessing mesh with %d faces\n", mesh.fn);
          static_grid.Set(mesh.face.begin(), mesh.face.end());
    tri::UpdateFlags<MyMesh>::FaceProjection(mesh);
  
  //--------------------------------------------------------------------------------------//
  //
  //                                PERFORM PROJECTIONS AND UPDATE NORMALS
  //
  //--------------------------------------------------------------------------------------//
  int t1=clock();
  for(int i=0; i<in_cloud.vn; i++){

    vcg::Ray3f ray;
    Point3f orig = in_cloud.vert[i].P();
    Point3f dir = in_cloud.vert[i].N();

    ray.SetOrigin(orig);
    ray.SetDirection(dir);


    MyFace* f_ptr=GridDoRay(static_grid,FintFunct, mf, ray, maxDist, t);

        if (f_ptr)
        {   MyMesh::CoordType tt = in_cloud.vert[i].P()+in_cloud.vert[i].N()*t;
            int f_i = vcg::tri::Index(mesh, f_ptr);
            MyMesh::CoordType ti = (mesh.face[f_i].V(0)->N()+mesh.face[f_i].V(1)->N()+mesh.face[f_i].V(2)->N())/3;
            double t0;
            t0 = sqrt(ti[0]*ti[0]+ti[1]*ti[1]+ti[2]*ti[2]);
            out_cloud.vert[i].N() = ti/t0;
            out_cloud.vert[i].Q() = t;
            out_cloud.vert[i].P() = tt;

        }
        else
        {   printf("Couldn't trace landmark %d along ray on %s: closest point on target will be sought\n",i+1,filename);
            Point3f& currp = in_cloud.vert[i].P();
            Point3f& clost = out_cloud.vert[i].P();
            MyFace* f_ptr=GridClosest(static_grid, PDistFunct, mf, currp, maxDist, minDist, clost);
            int f_i = vcg::tri::Index(mesh, f_ptr);
            MyMesh::CoordType ti = (mesh.face[f_i].V(0)->N()+mesh.face[f_i].V(1)->N()+mesh.face[f_i].V(2)->N())/3;
            double t0;
            t0 = sqrt(ti[0]*ti[0]+ti[1]*ti[1]+ti[2]*ti[2]);
            out_cloud.vert[i].N() = ti/t0;
            out_cloud.vert[i].Q() = minDist;
        }

  }

  int t2 = clock();
  tri::io::ExporterPLY<MyMesh>::Save(out_cloud,"out_cloud.ply",tri::io::Mask::IOM_VERTNORMAL +tri::io::Mask::IOM_VERTQUALITY, false); // in ASCII
  printf("Completed projection of %d sample in %i msec\n", in_cloud.vn, t2-t1);

  return 0;
}
