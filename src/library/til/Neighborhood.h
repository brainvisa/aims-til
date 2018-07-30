#ifndef TIL_NEIGHBORHOOD_H
#define TIL_NEIGHBORHOOD_H

// include from STL library
#include <cassert>
#include <vector>

// includes from TIL library
#include "til/labels.h"
#include "til/numeric_array.h"
#include "til/Range.h"

// Ignore specific warnings
#ifdef _MSC_VER
#pragma warning( push )
// X needs to have dll-interface to be used by clients of class Z
#pragma warning( disable : 4251 )
#endif

/// A macro to quickly index into a 3^3 cube.
/// Coordinates vary in [-1;1].
#define COFFSET(i,j,k) ((i)+3*(j)+9*(k)+13)

// To be undefined at the end of this file
#define FOR_ALL_NEIGHBORS		\
for (i=-1; i<=1; ++i)			\
	for (j=-1; j<=1; ++j)		\
		for (k=-1; k<=1; ++k)	\

#define APPY_ALL_NEIGHBORS(macro)	\
	macro(-1,-1,-1)	\
	macro( 0,-1,-1)	\
	macro(+1,-1,-1)	\
	macro(-1, 0,-1)	\
	macro( 0, 0,-1)	\
	macro(+1, 0,-1)	\
	macro(-1,+1,-1)	\
	macro( 0,+1,-1)	\
	macro(+1,+1,-1)	\
	macro(-1,-1, 0)	\
	macro( 0,-1, 0)	\
	macro(+1,-1, 0)	\
	macro(-1, 0, 0)	\
	macro( 0, 0, 0)	\
	macro(+1, 0, 0)	\
	macro(-1,+1, 0)	\
	macro( 0,+1, 0)	\
	macro(+1,+1, 0)	\
	macro(-1,-1,+1)	\
	macro( 0,-1,+1)	\
	macro(+1,-1,+1)	\
	macro(-1, 0,+1)	\
	macro( 0, 0,+1)	\
	macro(+1, 0,+1)	\
	macro(-1,+1,+1)	\
	macro( 0,+1,+1)	\
	macro(+1,+1,+1)	\


// Namespace
namespace til {


// The type encoding neighborhood internally
// Could be int, bool, char...
//typedef bool t_nh;



// The neighborhood class represents a discrete, 3x3 neighborhood
// used e.g. for 3d image topology operations

// The neighborhood is center on (0,0,0), i.e. its range is [-1,1]^3


/// Neighborhood as a list of vectors.
// NB: I have chosen to use a std::vector container rather than a std::list. This might
// seem awkward as all the operations used in this class are rather costly when using a vector.
// But the bottleneck here is when we iterate, not when we actually manipulate the neighborhood.
class TIL_API Neighborhood_list : public Neighborhood_label
{
public: // typedefs

  typedef Neighborhood_list                     Self;
  typedef std::vector<numeric_array<int,3> >    NeighborList;
  typedef NeighborList::iterator                iterator;
  typedef NeighborList::const_iterator          const_iterator;

public: // constructors & destructor

	/// Default constructor, empty neighborhood
	Neighborhood_list() : m_neighbors() { }
	
	/// Copy constructor from any type of neighborhood
	template < class TNeighborhood >
	Neighborhood_list(const TNeighborhood &nh) : m_neighbors() { this->init(nh); }

	/// From C array
	Neighborhood_list(const bool array[3*3*3]) : m_neighbors() { this->init(array); }

public: // initialization

	template < class TNeighborhood >
	typename enable_if<is_Neighborhood<TNeighborhood> >::type
	init(const TNeighborhood &nh)
	{
		// We first clear the neighborhood, in case their were something else before
		this->reset();
		//nh.for_all_neighbors(this->set(_1));
		// We add a neighbor whenever there is one in the old neighborhood
		typename TNeighborhood::const_iterator iNeighbor = nh.begin();
		for (; iNeighbor < nh.end(); ++iNeighbor)
		{
			if (nh.isNeighbor(*iNeighbor))
			{
				this->set(*iNeighbor);
			}
		}
	}

	/// From C array
	void init(const bool array[3*3*3])
	{
		this->reset();
		int i, j, k;
		FOR_ALL_NEIGHBORS
		{
			if (array[COFFSET(i,j,k)]) this->set(numeric_array<int,3>(i,j,k));
		}
	}

	/// A loop to execute something per neighbor
	template < typename TFunctor >
	void for_all_neighbors(const TFunctor &functor)
	{
		for (
			std::vector<numeric_array<int,3> >::const_iterator iNeighbor = m_neighbors.begin();
			iNeighbor != m_neighbors.end();
			++iNeighbor)
		{
			functor(*iNeighbor);
		}
	}

public: // iterators

	Self::iterator begin() { return m_neighbors.begin(); }
	Self::iterator end() { return m_neighbors.end(); }
	Self::const_iterator begin() const { return m_neighbors.begin(); }
	Self::const_iterator end() const { return m_neighbors.end(); }

public: // set & get

	void set(const numeric_array<int,3> &v)
	{
		(void)(v);
		// Check that input vector verifies basic assumptions
    assert(all_less_equal(v, numeric_array<int,3>(1,1,1)));
		assert(all_greater_equal(v, numeric_array<int,3>(-1,-1,-1)));
		// Check that the input neighbor is not already inside the list
		// If not, add the neighbor
	}
	void reset(const numeric_array<int,3> &v)
	{
		(void)(v);
		// Check that input vector verifies basic assumptions
		assert(all_less_equal(v, numeric_array<int,3>(1,1,1)));
		assert(all_greater_equal(v, numeric_array<int,3>(-1,-1,-1)));
		// Remove input vector
	}

	/*
	bool get()(const numeric_array<int,3> &v)
	{
	}
	*/

public: // functions
	
	void reset() { m_neighbors.clear(); }

	/*
	template < class TFunctor >
	INLINE void forAllNeighbors(const TFunctor &functor) const
	{
		for (
			NeighborList::const_iterator iNeighbor = m_neighbors.begin();
			iNeighbor < m_neighbors.end();
			++iNeighbor)
		{
			functor(*iNeighbor);
		}
	}
	*/

private: // typedefs

private: // classes

	struct Assign
	{
		void operator()(const numeric_array<int,3> &)
		{
		}
	};

private: // data
	NeighborList m_neighbors;
};


class TIL_API Neighborhood
{
public: // typedefs

	typedef Neighborhood Self;

public:	// constructors &destructor

	/// No neighbors
	Neighborhood() { this->reset(); }

	/// Copy constructor
	Neighborhood(const Neighborhood & nh) { this->init(nh); }

	/// from C array
	Neighborhood(const bool array[3*3*3]) { this->init(array); }

	/// manual
	Neighborhood(
		bool n000, bool n001, bool n002,
		bool n010, bool n011, bool n012,
		bool n020, bool n021, bool n022,
		bool n100, bool n101, bool n102,
		bool n110, bool n111, bool n112,
		bool n120, bool n121, bool n122,
		bool n200, bool n201, bool n202,
		bool n210, bool n211, bool n212,
		bool n220, bool n221, bool n222)
	{
		m_cube[0]= n000; m_cube[1]= n001; m_cube[2]= n002;
		m_cube[3]= n010; m_cube[4]= n011; m_cube[5]= n012;
		m_cube[6]= n020; m_cube[7]= n021; m_cube[8]= n022;
		m_cube[9]= n100; m_cube[10]= n101; m_cube[11]= n102;
		m_cube[12]= n110; m_cube[13]= n111; m_cube[14]= n112;
		m_cube[15]= n120; m_cube[16]= n121; m_cube[17]= n122;
		m_cube[18]= n200; m_cube[19]= n201; m_cube[20]= n202;
		m_cube[21]= n210; m_cube[22]= n211; m_cube[23]= n212;
		m_cube[24]= n220; m_cube[25]= n221; m_cube[26]= n222;
	}

public: // initialization

	// from another neighborhood
	void init(const Neighborhood &nh) { this->init(nh.m_cube); }

	// from C array
	void init(const bool array[3*3*3])
	{
		for (int i = 0; i < 3*3*3; ++i) m_cube[i] = array[i];
	}

public: // functions

	/// Remove all neighbors from neighborhood
	void reset() { for (int i = 0; i <= 3*3*3; ++i) m_cube[i] = false; }

	/// Put a neighbor at (i,j,k) (do nothing if already there)
	void set(int i, int j, int k)
	{
		this->setValue(i,j,k,true);
	}

	/// Put a neighbor at position v (do nothing if already there)
	void set(const numeric_array<int,3> &v)
	{
		this->set(EXPAND_VECTOR(v));
	}


	/// Remove the neighbor (i,j,k) (do nothing if it was not there in the first place)
	void remove(int i, int j, int k)
	{
		this->setValue(i,j,k,false);
	}

	/// Remove the neighbor v (do nothing if it was not there in the first place)
	void remove(const numeric_array<int,3> &v)
	{
		this->remove(EXPAND_VECTOR(v));
	}

	// The returned value evalutes to true iff there is a neighbor
	// at point (i,j,k)
	bool isNeighbor(int i, int j, int k) const
	{
		assert(contains(Range<int,3>(numeric_array<int,3>(-1,-1,-1),numeric_array<int,3>(1,1,1)), numeric_array<int,3>(i, j, k)));
		return m_cube[Self::offset(i,j,k)];
	}
	bool isNeighbor(const numeric_array<int,3> &v)
	{
		return this->isNeighbor(EXPAND_VECTOR(v));
	}

	template < typename TFunctor >
	void for_all_neighbors(const TFunctor &functor)
	{
		int i, j, k;
		FOR_ALL_NEIGHBORS
		{
			if (m_cube[COFFSET(i,j,k)])
			{
				functor(i,j,k);
			}
		}
	}

public: // interface to use Neigbors as template

	template <int i, int j, int k>
	INLINE bool isNeighbor() const
	{
		// Check that we are not out of range
		// NB: this could/should be compile-time checking
		assert(i >=-1 && i <= 1);
		assert(j >=-1 && j <= 1);
		assert(k >=-1 && k <= 1);
		return m_cube[Self::offset(i,j,k)];
	}


private: // functions

	void setValue(int i, int j, int k, bool value)
	{
		assert(contains(Range<int,3>(numeric_array<int,3>(-1,-1,-1),numeric_array<int,3>(1,1,1)), numeric_array<int,3>(i, j, k)));
		m_cube[Self::offset(i,j,k)] = value;
	}

	INLINE static int offset(int i, int j, int k) { return i+3*j+9*k+13; }

private: // data

	bool m_cube[3*3*3];
};




const Neighborhood N26(
1,1,1,
1,1,1,
1,1,1,

1,1,1,
1,0,1,
1,1,1,

1,1,1,
1,1,1,
1,1,1);

const Neighborhood N6(
0,0,0,
0,1,0,
0,0,0,

0,1,0,
1,0,1,
0,1,0,

0,0,0,
0,1,0,
0,0,0);

const Neighborhood N18(
0,1,0,
1,1,1,
0,1,0,

1,1,1,
1,0,1,
1,1,1,

0,1,0,
1,1,1,
0,1,0);

#define TNDECARG                \
bool n000, bool n001, bool n002,\
bool n010, bool n011, bool n012,\
bool n020, bool n021, bool n022,\
bool n100, bool n101, bool n102,\
bool n110, bool n111, bool n112,\
bool n120, bool n121, bool n122,\
bool n200, bool n201, bool n202,\
bool n210, bool n211, bool n212,\
bool n220, bool n221, bool n222 
#define TNARG    \
n000, n001, n002,\
n010, n011, n012,\
n020, n021, n022,\
n100, n101, n102,\
n110, n111, n112,\
n120, n121, n122,\
n200, n201, n202,\
n210, n211, n212,\
n220, n221, n222 


/// Template neighborhoods.
/// These are used to inline neighborhood operations
template <
bool n000, bool n001, bool n002,
bool n010, bool n011, bool n012,
bool n020, bool n021, bool n022,
bool n100, bool n101, bool n102,
bool n110, bool n111, bool n112,
bool n120, bool n121, bool n122,
bool n200, bool n201, bool n202,
bool n210, bool n211, bool n212,
bool n220, bool n221, bool n222 >
class TArgNeighborhood
{
public:
	/*
	template <int i, int j, int k>
	struct isNeighbor;
	
	template <>
	struct isNeighbor<-1,-1,-1>
	{
		static INLINE bool operator()() { return n000; }
	};
	*/

	INLINE bool isNeighbor(int i, int j, int k) const
	{
	if (i==-1 && j==-1 && k==-1) return n000;
	if (i==-1 && j==-1 && k== 0) return n001;
	if (i==-1 && j==-1 && k==+1) return n002;

	if (i==-1 && j== 0 && k==-1) return n010;
	if (i==-1 && j== 0 && k== 0) return n011;
	if (i==-1 && j== 0 && k==+1) return n012;

	if (i==-1 && j==+1 && k==-1) return n020;
	if (i==-1 && j==+1 && k== 0) return n021;
	if (i==-1 && j==+1 && k==+1) return n022;

	if (i== 0 && j==-1 && k==-1) return n100;
	if (i== 0 && j==-1 && k== 0) return n101;
	if (i== 0 && j==-1 && k==+1) return n102;

	if (i== 0 && j== 0 && k==-1) return n110;
	if (i== 0 && j== 0 && k== 0) return n111;
	if (i== 0 && j== 0 && k==+1) return n112;

	if (i== 0 && j==+1 && k==-1) return n120;
	if (i== 0 && j==+1 && k== 0) return n121;
	if (i== 0 && j==+1 && k==+1) return n122;

	if (i==+1 && j==-1 && k==-1) return n200;
	if (i==+1 && j==-1 && k== 0) return n201;
	if (i==+1 && j==-1 && k==+1) return n202;

	if (i==+1 && j== 0 && k==-1) return n210;
	if (i==+1 && j== 0 && k== 0) return n211;
	if (i==+1 && j== 0 && k==+1) return n212;

	if (i==+1 && j==+1 && k==-1) return n220;
	if (i==+1 && j==+1 && k== 0) return n221;
	if (i==+1 && j==+1 && k==+1) return n222;
	return false;
	}

	template <int i, int j, int k>
	INLINE bool isNeighbor() const { return this->isNeighbor(i,j,k); }

	/*
	template <> INLINE bool isNeighbor<-1,-1,-1>() const { return n000; }
	template <> INLINE bool isNeighbor<-1,-1, 0>() const { return n001; }
	template <> INLINE bool isNeighbor<-1,-1,+1>() const { return n002; }

	template <> INLINE bool isNeighbor<-1, 0,-1>() const { return n010; }
	template <> INLINE bool isNeighbor<-1, 0, 0>() const { return n011; }
	template <> INLINE bool isNeighbor<-1, 0,+1>() const { return n012; }

	template <> INLINE bool isNeighbor<-1,+1,-1>() const { return n020; }
	template <> INLINE bool isNeighbor<-1,+1, 0>() const { return n021; }
	template <> INLINE bool isNeighbor<-1,+1,+1>() const { return n022; }

	template <> INLINE bool isNeighbor< 0,-1,-1>() const { return n100; }
	template <> INLINE bool isNeighbor< 0,-1, 0>() const { return n101; }
	template <> INLINE bool isNeighbor< 0,-1,+1>() const { return n102; }

	template <> INLINE bool isNeighbor< 0, 0,-1>() const { return n110; }
	template <> INLINE bool isNeighbor< 0, 0, 0>() const { return n111; }
	template <> INLINE bool isNeighbor< 0, 0,+1>() const { return n112; }

	template <> INLINE bool isNeighbor< 0,+1,-1>() const { return n120; }
	template <> INLINE bool isNeighbor< 0,+1, 0>() const { return n121; }
	template <> INLINE bool isNeighbor< 0,+1,+1>() const { return n122; }

	template <> INLINE bool isNeighbor<+1,-1,-1>() const { return n200; }
	template <> INLINE bool isNeighbor<+1,-1, 0>() const { return n201; }
	template <> INLINE bool isNeighbor<+1,-1,+1>() const { return n202; }

	template <> INLINE bool isNeighbor<+1, 0,-1>() const { return n210; }
	template <> INLINE bool isNeighbor<+1, 0, 0>() const { return n211; }
	template <> INLINE bool isNeighbor<+1, 0,+1>() const { return n212; }

	template <> INLINE bool isNeighbor<+1,+1,-1>() const { return n220; }
	template <> INLINE bool isNeighbor<+1,+1, 0>() const { return n221; }
	template <> INLINE bool isNeighbor<+1,+1,+1>() const { return n222; }
*/
};

/*
template <int i, int j, int k, class TNeighborhood>
struct IsNeighbor;

template <class TNeighborhood>
struct IsNeighbor<-1,-1,-1,TNeighborhood>
{
	static const bool value = TNeighborhood::n000;
};

template <int i, int j, int k, class TNeighborhood>
bool isNeighbor(const TNeighborhood&);

template <int i, int j, int k, TNDECARG>
bool isNeighbor(const TArgNeighborhood<TNARG> &)
{
  return n000;
}

template <int i, int j, int k>
struct RPos {};

template <TNDECARG>
bool isNeighbor<TNARG>(const TArgNeighborhood<TNARG> &)
{
  return n000;
}
*/

// Below are some predefined template neighborhoods corresponding to standard 3D
// connectivities

/// Template neighborhood for standard 3D 6-connectivity.
/// Note that the center is not included.
typedef TArgNeighborhood<
0,0,0,
0,1,0,
0,0,0,
0,1,0,
1,0,1,
0,1,0,
0,0,0,
0,1,0,
0,0,0> TN6;

/// Template neighborhood for standard 3D 18-connectivity
/// Note that the center is not included.
typedef TArgNeighborhood<
0,1,0,
1,1,1,
0,1,0,
1,1,1,
1,0,1,
1,1,1,
0,1,0,
1,1,1,
0,1,0> TN18;

/// Template neighborhood for standard 3D 26-connectivity
/// Note that the center is not included.
typedef TArgNeighborhood<
1,1,1,
1,1,1,
1,1,1,
1,1,1,
1,0,1,
1,1,1,
1,1,1,
1,1,1,
1,1,1> TN26;

} // namespace

#undef FOR_ALL_NEIGHBORS

#ifdef _MSC_VER
#pragma warning ( pop )
#endif

#endif

