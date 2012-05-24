// Author: Stefan Schlager
// Date: 26 November 2010
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
typedef vcg::GridStaticPtr<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid;
//typedef vcg::SpatialHashTable<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid;

int main(int argc,char ** argv){
    float thresh = 10;
    bool cloud = false;
    bool noout = true;
    bool inbound = false;
    bool strict = false;
    bool minray = false;
    char filename[256];

  if (argc < 3){
		printf("\n");
    printf("    Compute a projection of a point cloud onto a mesh along given rays\n");
    printf("    Usage: rayproject <reference mesh> < targetmesh> [-t | -cloud | --inbound | --strict | -o] \n");
    printf("       <reference mesh>        Point cloud over which to project (PLY format).\n");
    printf("       <target mesh>         Mesh model over which to project (PLY format).\n");
    printf("       -t <threshold> - delimits max distance to seek along ray. default is 10\n");
    printf("       -cloud - vertex normals have to be computed individually.\n");
    printf("       --inbound - search will be along oposite of normal first.\n");
    printf("       --strict - this options will write the value 1e12 into Vertex quality.\n");
    printf("       --minray - find the closest point looking in both directions");
    printf("                  if no face is hit by the ray.\n");
    printf("       -o <output.ply> - define output filename (and path) default is project.mesh.ply.\n");

    printf("\n");    
    printf("projects the vertices of the reference mesh onto the target mesh when terminated.\n");
    //printf("The vertex coordinates represent the projected samples (in same order as input.\n");
    printf("The vertex quality (requires ply files) represents the projection distance.\n");
		return 0;
	}


  for (int i = 1; i < argc; i++) {

      if (i + 1 != argc) {// Check that we haven't finished parsing already for option wtih value
	if (strcmp("-t", argv[i]) == 0)
	  {
	    thresh = atof(argv[i + 1]);
	  }
	if (strcmp("-o", argv[i]) == 0)
	  {
	    strcpy(filename, argv[i+1]);
	    if (i==argc-2)
	      {noout=false;}
	  }
      }
                if (strcmp("-cloud", argv[i]) == 0)
                {
                cloud = true;

            }
                
                if (strcmp("--inbound", argv[i]) == 0)
                {
                inbound = true;
                
            }
                if (strcmp("--minray", argv[i]) == 0)
                {
                minray = true;

            }
                if (strcmp("--strict", argv[i]) == 0)
                {
                strict = true;
            }
       }
  

    if (noout == true)
          {
          strcpy(filename, "project.mesh.ply");
          }

   // thresh=atof(argv[3]);
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
  //tri::Smooth<MyMesh>::VertexNormalLaplacian(mesh,2,false);

  tri::UpdateNormals<MyMesh>::NormalizeVertex(mesh);
  tri::UpdateBounding<MyMesh>::Box(mesh);
  //tri::UpdateNormals<MyMesh>::PerFaceNormalized(in_cloud);
  //tri::UpdateNormals<MyMesh>::PerVertexAngleWeighted(in_cloud);
  tri::UpdateBounding<MyMesh>::Box(in_cloud);
  if (cloud == false)
  {
  tri::UpdateNormals<MyMesh>::PerFaceNormalized(in_cloud);
  tri::UpdateNormals<MyMesh>::PerVertexAngleWeighted(in_cloud);
}
  tri::UpdateNormals<MyMesh>::NormalizeVertex(in_cloud);
  tri::UpdateQuality<MyMesh>::VertexConstant(in_cloud, 0);

  tri::UpdateQuality<MyMesh>::VertexConstant(out_cloud, 0);

  float maxDist = mesh.bbox.Diag();
  float minDist = 1e-10;
  float t;
  float tNew;

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
  //                                PERFORM DISTANCE QUERIES
  //
  //--------------------------------------------------------------------------------------//
  int t1=clock();
   //int count = 1;
  for(int i=0; i<in_cloud.vn; i++){

    vcg::Ray3f ray;
    Point3f orig = in_cloud.vert[i].P();
    Point3f dir = in_cloud.vert[i].N();

    if (inbound == true)
        {dir = -dir;
    }
    ray.SetOrigin(orig);
    ray.SetDirection(dir);


        MyFace* f_ptr=GridDoRay(static_grid,FintFunct, mf, ray, maxDist, t);

        if (f_ptr && t < thresh)
            {
                double tOut = t;

                if (minray == true)
                {
                    Point3f dir1 = -dir;
                    ray.SetDirection(dir1);
                    MyFace* f_ptr1=GridDoRay(static_grid,FintFunct, mf, ray, maxDist, tNew);
                // check if reverse ray finds a closer match
                    if (f_ptr1 && tNew < t)
                    {
                    dir = dir1;
                    t = tNew;
                    tOut = -t;
                    f_ptr = f_ptr1;
               }

            }

                   MyMesh::CoordType tt = in_cloud.vert[i].P()+dir*t;
                   int f_i = vcg::tri::Index(mesh, f_ptr);
                   MyMesh::CoordType ti = (mesh.face[f_i].V(0)->N()+mesh.face[f_i].V(1)->N()+mesh.face[f_i].V(2)->N())/3;
                   double t0;
                   t0 = sqrt(ti[0]*ti[0]+ti[1]*ti[1]+ti[2]*ti[2]);
                   out_cloud.vert[i].N() = ti/t0;
                   out_cloud.vert[i].Q() = tOut;
                   out_cloud.vert[i].P() = tt;
        }
        else
        {
            ray.SetDirection(-dir);
            MyFace* f_ptr=GridDoRay(static_grid,FintFunct, mf, ray, maxDist, t);

            if (f_ptr && t < thresh)

                {MyMesh::CoordType tt = in_cloud.vert[i].P()+dir*t;
                int f_i = vcg::tri::Index(mesh, f_ptr);
                MyMesh::CoordType ti = (mesh.face[f_i].V(0)->N()+mesh.face[f_i].V(1)->N()+mesh.face[f_i].V(2)->N())/3;
                double t0;
                t0 = sqrt(ti[0]*ti[0]+ti[1]*ti[1]+ti[2]*ti[2]);
                out_cloud.vert[i].N() = ti/t0;
                out_cloud.vert[i].Q() = -t;
                out_cloud.vert[i].P() = tt;

            }
            else
            {   //printf("Couldn't trace landmark %d along ray, closest point on target will be sought\n",i+1);
                Point3f& currp = in_cloud.vert[i].P();
                Point3f& clost = out_cloud.vert[i].P();
                MyFace* f_ptr=GridClosest(static_grid, PDistFunct, mf, currp, maxDist, minDist, clost);
                if (strict == false)
                    { out_cloud.vert[i].Q() = minDist;
                }
                else
                    {out_cloud.vert[i].Q() = 1e12;
                }

                int f_i = vcg::tri::Index(mesh, f_ptr);
                MyMesh::CoordType ti = (mesh.face[f_i].V(0)->N()+mesh.face[f_i].V(1)->N()+mesh.face[f_i].V(2)->N())/3;
                double t0;
                t0 = sqrt(ti[0]*ti[0]+ti[1]*ti[1]+ti[2]*ti[2]);
                out_cloud.vert[i].N() = ti/t0;
                out_cloud.vert[i].P() = clost;
            }
          }

    in_cloud.vert[i].P()=out_cloud.vert[i].P();
    in_cloud.vert[i].Q()=out_cloud.vert[i].Q();
    in_cloud.vert[i].N()=out_cloud.vert[i].N();
    }
  //--------------------------------------------------------------------------------------//
  //
  //                                UPDATE PROJECTED MESH
  //
  //--------------------------------------------------------------------------------------//




  tri::UpdateBounding<MyMesh>::Box(in_cloud);
  if (cloud == false)
    {
    tri::UpdateNormals<MyMesh>::PerFaceNormalized(in_cloud);
    tri::Smooth<MyMesh>::VertexNormalLaplacian(in_cloud,2,false);
    tri::UpdateNormals<MyMesh>::NormalizeVertex(in_cloud);

    }
  int t2 = clock();
  tri::io::ExporterPLY<MyMesh>::Save(in_cloud,filename,tri::io::Mask::IOM_VERTNORMAL +tri::io::Mask::IOM_VERTQUALITY, false); // in ASCII
  printf("Completed projection of %d sample in %i msec\n", in_cloud.vn, t2-t1);


  //printf("%i ",tt);
  return 0;
}

