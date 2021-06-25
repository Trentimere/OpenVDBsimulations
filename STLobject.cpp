#include <openvdb/openvdb.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/LevelSetPlatonic.h>
#include <vector>
#include <openvdb/points/PointConversion.h>
#include <openvdb/tools/ParticlesToLevelSet.h>
#include <openvdb/points/PointAdvect.h>
#include <openvdb/points/PointMove.h>
#include <openvdb/tools/PotentialFlow.h>
#include <openvdb/tools/ValueTransformer.h>
#include <openvdb/math/Math.h>
#include <openvdb/tools/MeshToVolume.h>

#include "ppm.h"
#include "STLParser.h"

// grid element operations
struct Local {
  static inline void intersection(const openvdb::Vec3s& a,
				  const bool& b,
				  openvdb::Vec3s& result) {
    result = b*a;
    }

  // used to set the color value of a color grid
  static inline void op(const openvdb::BoolGrid::ValueOnCIter& iter,
			openvdb::Vec3SGrid::Accessor &accessor)
  {
    accessor.setValue(iter.getCoord(), openvdb::Vec3s(100,0,0));
  }

  static inline void sum(const openvdb::Vec3s& a, const openvdb::Vec3s&b,
			 openvdb::Vec3s& result)
  {
    result = a+b;
  }

};

struct Scaled
{
  Scaled(float f, float dt): scale(f), dtheta(dt) {}
  inline void operator()(const openvdb::Vec3s& a, const openvdb::Vec3s& b,
			 openvdb::Vec3s& result) const
  {
    result = scale*10*a + dtheta*b;
  }
  float scale;
  float dtheta;
};

struct MultScale
{
  MultScale(float s, float r): scale(s), theta(r) {}
  inline void operator()(const openvdb::Vec3s& a, const openvdb::Vec3s& b,
			 openvdb::Vec3s& result) const
  {
    result = scale*(a*openvdb::math::Cos(theta) - b*openvdb::math::Sin(theta));
  }
  float scale;
  float theta;
};


struct MatMul {
    openvdb::math::Mat3s M;
    MatMul(const openvdb::math::Mat3s& mat): M(mat) {}
    inline void operator()(const openvdb::Vec3SGrid::ValueOnIter& iter) const {
        iter.setValue(M.transform(*iter));
    }
};

// This is taken from the example in the unit tests
class MyParticleList
{
protected:
    struct MyParticle {
        openvdb::Vec3R p, v;
        openvdb::Real  r;
    };
    openvdb::Real           mRadiusScale;
    openvdb::Real           mVelocityScale;
    std::vector<MyParticle> mParticleList;
public:

    typedef openvdb::Vec3R  PosType;

    MyParticleList(openvdb::Real rScale=1, openvdb::Real vScale=1)
        : mRadiusScale(rScale), mVelocityScale(vScale) {}
    void add(const openvdb::Vec3R &p, const openvdb::Real &r,
             const openvdb::Vec3R &v=openvdb::Vec3R(0,0,0))
    {
        MyParticle pa;
        pa.p = p;
        pa.r = r;
        pa.v = v;
        mParticleList.push_back(pa);
    }
    /// @return coordinate bbox in the space of the specified transfrom
    openvdb::CoordBBox getBBox(const openvdb::GridBase& grid) {
        openvdb::CoordBBox bbox;
        openvdb::Coord &min= bbox.min(), &max = bbox.max();
        openvdb::Vec3R pos;
        openvdb::Real rad, invDx = 1/grid.voxelSize()[0];
        for (size_t n=0, e=this->size(); n<e; ++n) {
            this->getPosRad(n, pos, rad);
            const openvdb::Vec3d xyz = grid.worldToIndex(pos);
            const openvdb::Real   r  = rad * invDx;
            for (int i=0; i<3; ++i) {
                min[i] = openvdb::math::Min(min[i], openvdb::math::Floor(xyz[i] - r));
                max[i] = openvdb::math::Max(max[i], openvdb::math::Ceil( xyz[i] + r));
            }
        }
        return bbox;
    }

    /// Return the total number of particles in list.
    ///  Always required!
    size_t size() const { return mParticleList.size(); }

    /// Get the world space position of n'th particle.
    /// Required by ParticledToLevelSet::rasterizeSphere(*this,radius).
    void getPos(size_t n,  openvdb::Vec3R&pos) const { pos = mParticleList[n].p; }


    void getPosRad(size_t n,  openvdb::Vec3R& pos, openvdb::Real& rad) const {
        pos = mParticleList[n].p;
        rad = mRadiusScale*mParticleList[n].r;
    }
    void getPosRadVel(size_t n,  openvdb::Vec3R& pos, openvdb::Real& rad, openvdb::Vec3R& vel) const {
        pos = mParticleList[n].p;
        rad = mRadiusScale*mParticleList[n].r;
        vel = mVelocityScale*mParticleList[n].v;
    }
    // The method below is only required for attribute transfer
    void getAtt(size_t n, openvdb::Index32& att) const { att = openvdb::Index32(n); }
};


int main(int argc, char* argv[])
{
    openvdb::initialize();

    // make sure only have 2 arguments
    std::string filename;
    if (argc > 2)
    {
	return 1;
    }
    // extract stl file as a mesh from the given file name
    filename = argv[1];
    struct Mesh object = parse_stl(filename);
    // vectors to store point and triangle attributes
    std::vector<openvdb::Vec3s> points;
    std::vector<openvdb::Vec3I> triangles;

    // extract the points and triangles from the mesh object
    for (size_t i = 0; i < object.sortedVertices.size(); i++)
    {
      struct Point p = object.sortedVertices[i];
      points.push_back(openvdb::Vec3s(p.x, p.y, p.z));
    }
    for (size_t i = 0; i < object.triangles.size(); i++)
    {
      struct TrianglePoints t = object.triangles[i];
      triangles.push_back(openvdb::Vec3I(t.v1, t.v2, t.v3));
    }

    // initialize the transform for the object
    float voxelSize = 0.1, halfWidth = 10.0f;
    openvdb::math::Transform::Ptr tform =
      openvdb::math::Transform::createLinearTransform(voxelSize);

    // create the level-set object from the points and triangles
    openvdb::FloatGrid::Ptr grid =
      openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(*tform, points,
					      triangles, halfWidth);

    grid->setName("Box");

    openvdb::math::Transform xform = grid->transform();

    // define a color grid
    std::string colorGridName = "Color";
    // create a mask of the level-set object
    // in order to assign color at the surface
    openvdb::BoolGrid::Ptr surfaceMask =
      openvdb::tools::interiorMask(*grid);
    openvdb::tools::dilateVoxels(surfaceMask->tree(), 1);
    // create a color grid with a black default background
    openvdb::Vec3SGrid::Ptr colorGrid =
      openvdb::Vec3SGrid::create(openvdb::Vec3s(255,255,255));
    colorGrid->setTransform(xform.copy());
    // update the voxels at the surface of the object
    openvdb::tools::transformValues(surfaceMask->cbeginValueOn(),
				    *colorGrid, Local::op);
    colorGrid->setName(colorGridName);

    int numPoints = 1000; // number of particles to simulate
    float cube_size = 15.0; // size of half of cube size that particles are in
    std::vector<openvdb::Vec3R> positions(numPoints, openvdb::Vec3R(0, 0, 0));
// Randomize the point positions.
    std::mt19937 generator(/*seed=*/0);
    std::uniform_real_distribution<> distribution(-0.5, 0.5);
    // randomly distribute particles within the cube
    double i, j, k;
      for (int n = 0; n < numPoints; n++)
      {
	i = distribution(generator)*cube_size + 25.0;
	j = distribution(generator)*cube_size;
	k = 0*distribution(generator)*cube_size;
	positions[n] = openvdb::Vec3R(i,j,k);
      }
	  
    
    openvdb::points::PointAttributeVector<openvdb::Vec3R> positionsWrapper(positions);
    // Create a transform using this voxel-size.
    openvdb::math::Transform::Ptr transform =
        openvdb::math::Transform::createLinearTransform(voxelSize);
    // Create a PointDataGrid containing these four points and using the
    // transform given. This function has two template parameters, (1) the codec
    // to use for storing the position, (2) the grid we want to create
    // (ie a PointDataGrid).
    // We use no compression here for the positions.
    openvdb::points::PointDataGrid::Ptr pointGrid =
        openvdb::points::createPointDataGrid<openvdb::points::NullCodec,
                        openvdb::points::PointDataGrid>(positions, *transform);
    // Set the name of the grid
    pointGrid->setName("Points");

    std::string vdbFile = "grids.vdb";
    std::string ppmFile = "out.ppm";

    // size of window showing animation
    int width = 1240;
    int height = 720;
    
    // create an SDL context
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* window = SDL_CreateWindow("PPM", 100, 100, width, height, SDL_WINDOW_OPENGL);
    SDL_GLContext context = SDL_GL_CreateContext(window);
    SDL_Event event;

    // openvdb::Vec3d offset(-10, 0, 0);
    // event loop
    std::cout << "Creating Neumann velocites" << std::endl;
    const int dilation = 3;
    openvdb::MaskGrid::Ptr domain =
      openvdb::tools::createPotentialFlowMask(*grid, dilation);

    // compute potential flow for a global wind velocity around a object

    openvdb::Vec3SGrid::Ptr z_rotation =
      openvdb::Vec3SGrid::create(openvdb::Vec3s(0));
    z_rotation->setTransform(xform.copy());

    openvdb::CoordBBox activeBbox(grid->evalActiveVoxelBoundingBox());
    activeBbox.expand(5);

    z_rotation->denseFill(activeBbox, openvdb::Vec3s(0));

    for (auto leaf = z_rotation->tree().beginLeaf(); leaf; ++leaf) {
      for (auto iter = leaf->beginValueOn(); iter; ++iter) {
	openvdb::Vec3s position = xform.indexToWorld(iter.getCoord().asVec3d());
	openvdb::Vec3s vel = openvdb::Vec3s(0, 0, 1).cross(position);
	iter.setValue(vel);
      }
    }
    openvdb::Vec3fGrid::Ptr neumannZrot =
      openvdb::tools::createPotentialFlowNeumannVelocities(*grid,
            *domain, z_rotation, openvdb::Vec3s(0));


    openvdb::math::pcg::State state =
      openvdb::math::pcg::terminationDefaults<float>();

    state.iterations = 500;
    state.absoluteError = 1e-3;
    
    openvdb::FloatGrid::Ptr potentialZrot =
      openvdb::tools::computeScalarPotential(*domain, *neumannZrot, state);
    openvdb::Vec3SGrid::Ptr flowVelZrot =
      openvdb::tools::computePotentialFlow(*potentialZrot, *neumannZrot);


    
    std::cout << "Simulating for 5 seconds" << std::endl;
    
    int timeSteps = 5;
    float t = 0.0;
    float dt = 0.005;

    float theta = M_PI_2;
    std::cout << theta << std::endl;

    openvdb::FloatGrid::Ptr objectGrid = grid->deepCopy();
    openvdb::tools::erodeVoxels(objectGrid->tree(), int(halfWidth/2)+1);
    //openvdb::tools::signedFloodFill(objectGrid->tree());
    openvdb::Vec3SGrid::Ptr surfaceGradient =
      openvdb::tools::gradient(*objectGrid);

    openvdb::FloatGrid::Ptr gridCopy = grid->deepCopy();
    openvdb::math::Transform &gridXForm = gridCopy->transform();
    openvdb::math::Transform &colorXForm = colorGrid->transform();
    
    const openvdb::math::Axis zaxis = openvdb::math::Z_AXIS;

    gridXForm.preRotate(theta, zaxis);
    colorXForm.preRotate(theta, zaxis);
    
    
    float speed = 10.0;
    openvdb::Vec3s objectVelocity = openvdb::Vec3s(speed,0,0);

    float dtheta = 0;
    for (int frame = 0; t < 5.0; frame++) {
        if (SDL_PollEvent(&event)) {
            // break if esc key is pressed or window is closed
            if (event.type == SDL_QUIT) 
                break;
            if (
                event.type == SDL_KEYUP && 
                event.key.keysym.sym == SDLK_ESCAPE
            )
                break;
        }
    
	t += dt*timeSteps;
	theta += dt*timeSteps*dtheta;

	
	gridXForm.preRotate(dt*timeSteps*dtheta, zaxis);
	gridXForm.postTranslate(dt*timeSteps*objectVelocity);
	colorXForm.preRotate(dt*timeSteps*dtheta, zaxis);
	colorXForm.postTranslate(dt*timeSteps*objectVelocity);

	/*
	openvdb::Vec3s velocity(openvdb::math::Cos(theta),
				-openvdb::math::Sin(theta), 0);
	openvdb::Vec3fGrid::Ptr neumann =
	openvdb::tools::createPotentialFlowNeumannVelocities(*grid,
	  *domain, openvdb::Vec3fGrid::Ptr(), velocity);
	*/
	openvdb::Vec3SGrid::Ptr temp = openvdb::Vec3SGrid::create();
	temp->tree().combine2(surfaceGradient->tree(), flowVelZrot->tree(),
			    Scaled(speed, 0));
	  
        temp->setTransform(gridXForm.copy());
        
	openvdb::math::Mat3s rotation =
	  openvdb::math::rotation<openvdb::math::Mat3s>(zaxis, theta);
	openvdb::tools::foreach(temp->beginValueOn(), MatMul(rotation));
	
	openvdb::points::advectPoints(*pointGrid, *temp,
				      4, dt, timeSteps);
	
	if (frame % 5 == 0)
	{
	  MyParticleList pList;
	  float radius = 0.2;
	  for (auto leafIter = pointGrid->tree().cbeginLeaf(); leafIter; ++leafIter)
	    {
	      // Extract the position attribute from the leaf by name (P is position).
	      const openvdb::points::AttributeArray& positionArray =
		leafIter->constAttributeArray("P");
	      // Create read-only handles for position and radius.
		openvdb::points::AttributeHandle<openvdb::Vec3f> positionHandle(positionArray);
	      // Iterate over the point indices in the leaf.
	      for (auto indexIter = leafIter->beginIndexOn(); indexIter; ++indexIter)
		{
		  // Extract the voxel-space position of the point.
		  openvdb::Vec3f voxelPosition = positionHandle.get(*indexIter);
		  // Extract the world-space position of the voxel.
		  openvdb::Vec3d xyz = indexIter.getCoord().asVec3d();
		  // Compute the world-space position of the point.
		  openvdb::Vec3f worldPosition =
		    grid->transform().indexToWorld(voxelPosition + xyz);
		  pList.add(worldPosition, openvdb::Real(radius));
		}
	    }
	  
	  openvdb::FloatGrid::Ptr renderGrid = 
	    openvdb::createLevelSet<openvdb::FloatGrid>(voxelSize, halfWidth);
	  
	  openvdb::tools::ParticlesToLevelSet<openvdb::FloatGrid> raster(*renderGrid);
  
	  raster.setGrainSize(1);//a value of zero disables threading
	  raster.rasterizeSpheres(pList);
	  raster.finalize();
  
	  openvdb::FloatGrid::Ptr temp =
	    openvdb::createLevelSet<openvdb::FloatGrid>(voxelSize, halfWidth);
  
	  const openvdb::math::Transform &targetXform = temp->transform();
	  // Compute a source grid to target grid transform.
	  // (For this example, we assume that both grids' transforms are linear,
	  // so that they can be represented as 4 x 4 matrices.)
	  openvdb::Mat4R xform =
	    gridXForm.baseMap()->getAffineMap()->getMat4() *
	    targetXform.baseMap()->getAffineMap()->getMat4().inverse();
	  // Create the transformer.
	  openvdb::tools::GridTransformer transformer(xform);
	  // Resample using nearest-neighbor interpolation.
	  transformer.transformGrid<openvdb::tools::PointSampler,
	  openvdb::FloatGrid>(*gridCopy, *temp);
	  temp->tree().prune();
	
	  openvdb::tools::csgUnion(*renderGrid, *temp);
	  // Create a VDB file object and write out the grid.
	  openvdb::io::File("grids.vdb").write({renderGrid, colorGrid});
	
	  std::string str = "vdb_render ";
	  str = str + vdbFile + " " + ppmFile + " -translate 30,0,100"
	    + " -lookat 30,0,0"
	    + " -res " + std::to_string(width) + "x" + std::to_string(height)
	    + " -color " + colorGridName
	    //+ " -shader normal"
	    + " -camera ortho"
	    + " -frame 60";
	  const char* command = str.c_str();
	  system(command);

	  // open file stream in binary mode
	  FILE* stream = fopen(ppmFile.c_str(), "rb");
  
	  // pixel map containing data from ppm
	  int depth;
	  unsigned int* pixmap = decodePPM(stream, &width, &height, &depth);
    
	  // clean up file stream
	  fclose(stream);

	  // draw calls here
	  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	  glRasterPos2i(-1.0f, -1.0f);
	  glDrawPixels(width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixmap);
	  SDL_GL_SwapWindow(window);
	  usleep(10000);
	}
	
    }
    

    // clean up SDL context
    SDL_GL_DeleteContext(context);
    SDL_Quit();
}
