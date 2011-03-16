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
     bool bcheck = false;
     bool scheck = false;
     bool noout = false;
     bool tolcheck = false;
     float tol = 1e3;
     float tol1;
  if (argc<3){
		printf("\n");
    printf("    Compute a projection of a point cloud onto a mesh along given normals\n    if no face is hit, the closest point on the target mesh is used.\n");
    printf("    Usage: trinorm_project <cloud> [options] <mesh>\n");
    printf("       <cloud>        Point cloud containing starting points, including normal\n");
    printf("                      information of vertices (PLY format).\n");
    printf("            -b = check both directions of normal.\n");
    printf("            --tolout <float> = max distance to look in normal direction.\n");
    printf("            --tolin <float> = max distance to look in opposite direction of the normal.\n");
    printf("       <mesh>         Mesh model over which to project (PLY format).\n");
    printf("\n");
    printf("Saves the projected point cloud in the file out_norm.ply when terminated.\n");
    printf("The vertex coordinates represent the projected samples (in same order as input.\n");
    printf("The vertex quality (requires ply files) represents the projection distance.\n");
		return 0;
	}
    for (int i = 1; i < argc; i++) {
  if (i + 1 != argc) {// Check that we haven't finished parsing already



                        if (strcmp("--tolout", argv[i]) == 0)
                        {tol = atof(argv[i+1]);
                            if (i==argc-2)
                            {noout=true;}
                        }

                        if (strcmp("--tolin", argv[i]) == 0)
                        {tol = atof(argv[i+1]);
                        tolcheck = true;
                        if (i==argc-2)
                            {noout=true;}
                        }
                    }
            }

    for (int i = 1; i < argc; i++) {
            if (strcmp("-b", argv[i]) == 0)
            {bcheck=true;
                if (i==argc-1)
                {noout=true;}
            }

        }



    if (noout==true)
           {printf("Error: please specify target mesh\n");
               return 0;
           }
    if (tolcheck == false)
        {tol1=tol;
        }
    strcpy(filename, argv[1]);
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
  int err2 = tri::io::Importer<MyMesh>::Open(mesh,argv[argc-1]);
  if(err2) {
                printf("Error in reading %s: '%s'\n",argv[argc-1],tri::io::Importer<MyMesh>::ErrorMsg(err2));
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

    if (bcheck == false)
    {
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

  else
  {MyFace* f_ptr=GridDoRay(static_grid,FintFunct, mf, ray, maxDist, t);

      if (f_ptr && t < tol) //check if face is hit and if distance is below threshold

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
      {
          MyFace* f_ptr=GridDoRay(static_grid,FintFunct, mf, ray, maxDist, t);
          ray.SetDirection(-dir);         // search the other direction of the normal
          if (f_ptr && t < tol1)  //check if face is hit and if distance is below threshold

              {MyMesh::CoordType tt = in_cloud.vert[i].P()+in_cloud.vert[i].N()*t;
              int f_i = vcg::tri::Index(mesh, f_ptr);
              MyMesh::CoordType ti = (mesh.face[f_i].V(0)->N()+mesh.face[f_i].V(1)->N()+mesh.face[f_i].V(2)->N())/3;
              double t0;
              t0 = sqrt(ti[0]*ti[0]+ti[1]*ti[1]+ti[2]*ti[2]);
              out_cloud.vert[i].N() = ti/t0;
              out_cloud.vert[i].Q() = t;
              out_cloud.vert[i].P() = tt;

          }




          else //find the closest point if theres nothing along the ray

          {  // printf("Couldn't trace landmark %d along ray on %s: closest point on target will be sought\n",i+1,target);
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

    }
  }

  int t2 = clock();
  tri::io::ExporterPLY<MyMesh>::Save(out_cloud,"out_norm.ply",tri::io::Mask::IOM_VERTNORMAL +tri::io::Mask::IOM_VERTQUALITY, false); // in ASCII
  printf("Completed projection of %d sample in %i msec\n", in_cloud.vn, t2-t1);

  return 0;
}
