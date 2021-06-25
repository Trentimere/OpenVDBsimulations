#include <openvdb/openvdb.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <vector>
#include <openvdb/points/PointConversion.h>
#include <openvdb/tools/ParticlesToLevelSet.h>
#include <openvdb/points/PointAdvect.h>
#include <openvdb/points/PointMove.h>
#include <openvdb/tools/PotentialFlow.h>
#include "ppm.h"

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

struct OffsetDeformer
{
  OffsetDeformer(const openvdb::Vec3d _offset, const float _dt)
    : offset(_offset), dt(_dt){ }
    template <typename LeafIterT>
    void reset(const LeafIterT&, size_t idx) { }
    template <typename IndexIterT>
    void apply(openvdb::Vec3d& position, const IndexIterT&) const
    {
      position += offset*dt;
    }

  openvdb::Vec3d offset;
  float dt;
};

int main()
{
    openvdb::initialize();
    float voxelSize = 0.1, halfWidth = 2.0f;
    // Create a FloatGrid and populate it with a narrow-band
    // signed distance field of a sphere.
    openvdb::FloatGrid::Ptr grid =
        openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
            /*radius=*/10.0, /*center=*/openvdb::Vec3f(-10, 0, 0),
            /*voxel size=*/voxelSize, /*width=*/halfWidth);
    // Associate some metadata with the grid.
    grid->insertMeta("radius", openvdb::FloatMetadata(50.0));
    // Name the grid "LevelSetSphere".
    grid->setName("Sphere");

    openvdb::Vec3SGrid::Ptr colorGrid =
      openvdb::Vec3SGrid::create(openvdb::Vec3s(100,0,0));
    colorGrid->setName("Color");

    int numPoints = 1000;
    float cube_size = 15.0;
    std::vector<openvdb::Vec3R> positions(numPoints, openvdb::Vec3R(0, 0, 0));
// Randomize the point positions.
    std::mt19937 generator(/*seed=*/0);
    std::uniform_real_distribution<> distribution(-0.5, 0.5);
    double i, j, k;
      for (int n = 0; n < numPoints; n++)
      {
	i = distribution(generator)*cube_size + 20.0;
	j = distribution(generator)*cube_size;
	k = distribution(generator)*cube_size;
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

    // compute potential flow for a global wind velocity around a sphere

    openvdb::Vec3s windVelocity(-20, 0, 0);
    openvdb::Vec3fGrid::Ptr neumann =
      openvdb::tools::createPotentialFlowNeumannVelocities(*grid,
            *domain, openvdb::Vec3fGrid::Ptr(), -windVelocity);

    openvdb::math::pcg::State state =
      openvdb::math::pcg::terminationDefaults<float>();

    state.iterations = 500;
    state.absoluteError = 1e-3;

    std::cout << "Calculating Potential" << std::endl;
    openvdb::FloatGrid::Ptr potential =
      openvdb::tools::computeScalarPotential(*domain, *neumann, state);

    std::cout << "Computing flow velocities" << std::endl;
    openvdb::Vec3SGrid::Ptr flowVel =
      openvdb::tools::computePotentialFlow(*potential, *neumann);
    openvdb::tools::changeBackground(flowVel->tree(), windVelocity);
    openvdb::io::File("flowGrid.vdb").write({flowVel});
    
    std::cout << "Simulating for 5 seconds" << std::endl;
    /*float dt = 0.002;
    openvdb::points::advectPoints(*pointGrid, *flowVel, 4, dt, 2500);
    std::cout << "Done simulating" << std::endl;
    
    */
    float t = 0.0;
    for (; t < 5.0; ) {
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
    
	float dt = 0.002;
	t += dt*50;
	//OffsetDeformer deformer(windVelocity, dt);
	// Move the points using this deformer
	//openvdb::points::movePoints(*pointGrid, deformer);
	openvdb::points::advectPoints(*pointGrid, *flowVel, 4, dt, 50);
	
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
	//grid->deepCopy();
	openvdb::tools::ParticlesToLevelSet<openvdb::FloatGrid> raster(*renderGrid);

	raster.setGrainSize(1);//a value of zero disables threading
	raster.rasterizeSpheres(pList);
	raster.finalize();
    
	// Create a VDB file object and write out the grid.
	openvdb::io::File("grids.vdb").write({renderGrid});

	std::string str = "vdb_render ";
	str = str + vdbFile + " " + ppmFile + " -translate -50,50,100"
	  + " -lookat 0,0,0"
	  + " -res " + std::to_string(width) + "x" + std::to_string(height);
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
    

    // clean up SDL context
    SDL_GL_DeleteContext(context);
    SDL_Quit();
}
