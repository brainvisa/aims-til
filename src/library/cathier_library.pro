TEMPLATE = lib
TARGET = guevara_cathier${BUILDMODEEXT}

INCBDIR =

#!include ../../config-cpp-library

HEADERS = \
# required for meshCleaner
   til/numeric_array.h \
   til/numeric_array_operators.h \
   til/numeric_array_tools.h \
   til/numeric_array.tpp \
# required for constellation
   til/sparse_vector.h \
   til/sparse_vector_operators.h \
   til/sparse_vector_policies.h \
   til/sparse_vector_tools.h \
   cathier/poly_solver.h \
   cathier/poly_solver_policies.h \
   cathier/poly_solver.tpp \
   cathier/geometrics.tpp \
   cathier/kdtree.tpp \
# useless ?
   til/Accumulator.h \
   til/basicFunctors.h \
   til/basic_iterator.h \
   til/BasicPixelActions.h \
   til/basic_range.h \
   til/boost_wrap.h \
   til/cat2type.h \
   til/ConditionalIterator.h \
   til/ext_declarations.h \
   til/functors.h \
   til/ghost_array.h \
   til/if_then_else.h \
   til/is_traits.h \
   til/KeysInterpolation.h \
   til/labelClasses.h \
   til/labels.h \
   til/loop.h \
   til/Matrix3.h \
   til/matrix3Tools.h \
   til/meta.h \
   til/misc_scalar_functions.h \
   til/miscTools.h \
   til/multi_array.h \
   til/NeighborhoodConfigurations.h \
   til/Neighborhood.h \
   til/neighborhoodTools.h \
   til/NormalizedCorrelation.h \
   til/PointerTraits.h \
   til/PointList.h \
   til/pointListTools.h \
   til/point_operators.h \
   til/point_tools.h \
   til/proba_distributions.h \
   til/Ptr.h \
   til/range_tools.h \
   til/SmartObject.h \
   til/stditerator.h \
   til/std_wrap.h \
   til/stl_declarations.h \
   til/StructureTensor.h \
   til/SymMatrix3.h \
   til/symmatrix3tools.h \
   til/templateTools.h \
   til/TensorProductInterpolation.h \
   til/TExprBasicFunctors.h \
   til/TExprConcatenation.h \
   til/TExprFunctions.h \
   til/TExpr.h \
   til/TExprMacros.h \
   til/TExprOperators.h \
   til/TExprPlaceHolders.h \
   til/til_common.h \
   til/til_declarations.h \
   til/til.h \
   til/til_math.h \
   til/traits.h \
   til/value_proxy.h \
   til/value_proxy_policies.h \
   til/misc_tools.tpp \
   cathier/accessors.h \
   cathier/aims_wrap.h \
   cathier/binary_tree.h \
   cathier/cat2type.h \
   cathier/convert.h \
   cathier/cyclic_iterator.h \
   cathier/declarations_external.h \
   cathier/fraction.h \
   cathier/fraction_policies.h \
   cathier/func_iterator.h \
   cathier/fuzzy_logic.h \
   cathier/geometrics.h \
   cathier/globalTraits.h \
   cathier/graph_accessors.h \
   cathier/graph.h \
   cathier/index_collection.h \
   cathier/kdtree.h \
   cathier/mesh_conversion.h \
   cathier/mesh_decimation.h \
   cathier/Mesh.h \
   cathier/MeshTraits.h \
   cathier/meshUtils.h \
   cathier/minTools.h \
   cathier/misc_sort.h \
   cathier/miscUtils.h \
   cathier/ordered_iterator.h \
   cathier/scalar_matrix.h \
   cathier/SparseVector.h \
   cathier/stditerator.h \
   cathier/std_wrap.h \
   cathier/triangle_mesh_geodesic_map.h \
   cathier/binary_tree.tpp \
   cathier/fraction_policies.tpp \
   cathier/fraction.tpp \
   cathier/fuzzy_logic.tpp \
   cathier/math_functions.tpp \
   cathier/mesh_conversion.tpp \
   cathier/mesh_decimation.tpp \
   cathier/mesh_utils.tpp \
   cathier/mics_tools.tpp \
   cathier/minTools.tpp \
   cathier/misc_sort.tpp \
   cathier/misc_utils.tpp \
   cathier/triangle_mesh_geodesic_map.tpp

SOURCES = \
   cathier/aims_wrap.cpp

