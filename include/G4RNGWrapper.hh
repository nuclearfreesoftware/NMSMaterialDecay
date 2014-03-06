#ifndef G4RNGWrapper_h
#define G4RNGWrapper_h 1

template<class T>
class G4RNGWrapper { 
  public:
    static void set(T* object, double (T::*func)(void));
    static double rng(void);
  private:
    static T* m_obj;
    static double (T::*m_func)(void);
};

template<class T> T* G4RNGWrapper<T>::m_obj;

template<class T> double (T::*G4RNGWrapper<T>::m_func)(void);

template<class T> void G4RNGWrapper<T>::set(T* object, double (T::*func)(void)) {
  m_obj = object; m_func = func;
}

template<class T> double G4RNGWrapper<T>::rng(void) { return (m_obj->*m_func)(); }

#endif
