#ifndef TIL_WORLD_IMAGE_H
#define TIL_WORLD_IMAGE_H

// include from TIL library
namespace til {


	template < class Image >
	class WorldImage : public Image
	{
	public:

		

	private:
		
	private:
		
		numeric_array<double,3> m_voxelSize;
		numeric_array<double,3> m_origin;

		// I have no idea whether this should be double or float...
		Affine<double> m_world2coord;
	};

} // namespace til

#endif

