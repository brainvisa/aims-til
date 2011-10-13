
namespace til
{

  //---------------------------------------------------------------------------

  template < typename T, std::size_t N >
  NaryTree< T, N >::
  ~NaryTree()
  {
    if (&*(root()))
    {
      post_order_scan(root(), Destructor());
    }
  }
  
  //---------------------------------------------------------------------------

  template < typename T, std::size_t N >
  typename NaryTree< T, N >::iterator
  NaryTree< T, N >::
  addChild(Node * parent, std::size_t childNumber, T value)
  {
    // Create new node
    Node * newNode = new Node(value);
    ++m_size;
    // Check that child has a parent
    if (parent)
    {
      // Warn parent it has a new child with given number
      parent->child(childNumber) = newNode;
    }
    // If child has no parent, it is a root:
    // Check that there is no root yet
    else if (m_root)
    {
      // Can't have two roots: throw
      throw TreeRootAlreadyDefined();
    }
    else
    {
      // Set new root
      m_root = newNode;
    }
    return iterator(newNode);
  }

  //---------------------------------------------------------------------------

} // namespace til

