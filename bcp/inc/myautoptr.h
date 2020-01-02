#ifndef __MY_AUTO_PTR
#define __MY_AUTO_PTR

// #include <memory> -- in the main header

BEGIN_COMMON_NAMESPACE__

  /**
   *  A wrapper class to provide auto_ptr with reference semantics.  For
   *  example, an auto_ptr can be assigned (or constructed from) the result of
   *  a function which returns an auto_ptr by value.
   *
   *  All the auto_ptr_ref stuff should happen behind the scenes.
   */
/*
  template<typename _Tp1>
    struct auto_ptr_ref
    {
      _Tp1* _M_ptr;

      explicit
      auto_ptr_ref(_Tp1* __p): _M_ptr(__p) { }
    };
*/

/// An analog to auto_ptr allowing inclusion in containers
/// because it meets the CopyConstructible and Assignable
///   *  requirements for Standard Library.
/// This is possible because assigning is allowed only
/// from empty

  template<typename _Tp>
    class my_auto_ptr : public auto_ptr <_Tp>
    {
    public:
      typedef typename auto_ptr<_Tp>::element_type element_type;
      explicit
      my_auto_ptr(element_type* __p = 0) throw()
        : auto_ptr<_Tp>(__p) { }
      my_auto_ptr(const my_auto_ptr& __a) throw()
        : auto_ptr<_Tp>((element_type*)(NULL))
      { assert(NULL == __a.get()); }
      template<typename _Tp1>
        my_auto_ptr(const my_auto_ptr& __a) throw()
	: auto_ptr<_Tp>((element_type*)(NULL))
      { assert(NULL == __a.get()); }
      my_auto_ptr&
      operator=(const my_auto_ptr& __a) throw()
      {
	assert(NULL == __a.get());
	reset(0);
	return *this;
      }
      template<typename _Tp1>
        my_auto_ptr&
        operator=(const my_auto_ptr& __a) throw()
        {
	  assert(NULL == __a.get());
	  reset(0);
	  return *this;
	}

      template<typename _Tp1>
        operator auto_ptr_ref<_Tp1>() throw()
        { return auto_ptr_ref<_Tp1>(this->release()); }

      template<typename _Tp1>
        operator my_auto_ptr<_Tp1>() throw()
        { return my_auto_ptr<_Tp1>(this->release()); }

/// Once more for auto_ptr:
///////////////////////////
      my_auto_ptr(const auto_ptr<_Tp>& __a) throw()
        : auto_ptr<_Tp>((element_type*)(NULL))
      { assert(NULL == __a.get()); }
      template<typename _Tp1>
        my_auto_ptr(const auto_ptr<_Tp1>& __a) throw()
	: auto_ptr<_Tp>((element_type*)(NULL))
      { assert(NULL == __a.get()); }
      my_auto_ptr&
      operator=(const auto_ptr<_Tp>& __a) throw()
      {
	assert(NULL == __a.get());
	reset(0);
	return *this;
      }
      template<typename _Tp1>
        my_auto_ptr&
        operator=(const auto_ptr<_Tp1>& __a) throw()
        {
	  assert(NULL == __a.get());
	  reset(0);
	  return *this;
	}

      template<typename _Tp1>
        operator auto_ptr<_Tp1>() throw()
        { return auto_ptr<_Tp1>(this->release()); }

  };

END_COMMON_NAMESPACE__

#endif // __MY_AUTO_PTR
