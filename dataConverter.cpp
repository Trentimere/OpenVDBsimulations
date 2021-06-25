#include <openvdb/openvdb.h>
#include <openvdb/points/PointConversion.h>


int main(int argc, char* argv[])
{
    // Initialize grid types and point attributes types.
    openvdb::initialize();

    // Create a VDB file object.
    if (argc != 2)
      exit(1);
    openvdb::io::File file(argv[1]);
    // Open the file.  This reads the file header, but not any grids.
    file.open();
    // Loop over all grids in the file and retrieve a shared pointer
    // to the one named "LevelSetSphere".  (This can also be done
    // more simply by calling file.readGrid("LevelSetSphere").)
    openvdb::GridBase::Ptr baseGrid;
    for (openvdb::io::File::NameIterator nameIter = file.beginName();
	 nameIter != file.endName(); ++nameIter)
      {
	// Read in only the grid we are interested in.
	if (nameIter.gridName() == "Points") {
	  baseGrid = file.readGrid(nameIter.gridName());
	} else {
	  std::cout << "skipping grid " << nameIter.gridName() << std::endl;
	}
      }
    file.close();
    // From the example above, "LevelSetSphere" is known to be a FloatGrid,
    // so cast the generic grid pointer to a FloatGrid pointer.
    openvdb::points::PointDataGrid::Ptr grid =
      openvdb::gridPtrCast<openvdb::points::PointDataGrid>(baseGrid);
    
    // Iterate over all the leaf nodes in the grid.
	std::cout << "x,y,z,diameter" << std::endl;
    for (auto leafIter = grid->tree().cbeginLeaf(); leafIter; ++leafIter) {
        // Verify the leaf origin.
        // Extract the position attribute from the leaf by name (P is position).
        const openvdb::points::AttributeArray& positionArray =
            leafIter->constAttributeArray("P");
        // Extract the radius attribute from the leaf by name (pscale is radius).
        const openvdb::points::AttributeArray& darray =
            leafIter->constAttributeArray("diameter");
        // Create read-only handles for position and radius.
        openvdb::points::AttributeHandle<openvdb::Vec3f> positionHandle(positionArray);
        openvdb::points::AttributeHandle<float> diameterHandle(darray);
        // Iterate over the point indices in the leaf.
        for (auto indexIter = leafIter->beginIndexOn(); indexIter; ++indexIter) {
            // Extract the voxel-space position of the point.
            openvdb::Vec3f voxelPosition = positionHandle.get(*indexIter);
            // Extract the world-space position of the voxel.
            openvdb::Vec3d xyz = indexIter.getCoord().asVec3d();
            // Compute the world-space position of the point.
            openvdb::Vec3f worldPosition =
                grid->transform().indexToWorld(voxelPosition + xyz);
            // Extract the radius of the point.
            float diameter = diameterHandle.get(*indexIter);
            // Verify the index, world-space position and radius of the point.
            std::cout << worldPosition(0) << "," << worldPosition(1) << ","
		      << worldPosition(2);
            std::cout << "," << diameter << std::endl;
        }
    }
}
