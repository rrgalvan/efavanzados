#include <iostream>
#include <vector>

class X { // Clase muy complicada, con vectores enormes, etc...
public:
  X(int _x) : x(_x) {} // Constructor
  void print() const { std::cout << x << std::endl; }
private:
  int x;
};


void print_x(X x) {
  x.print();
}

void print_x_ptr(const X* x)  {
  x->print();
}

void print_x_ref(const X& x) {
  x.print();
}

template<class T> class MiVector {
public:
  MiVector(int _n) : n(_n) { x = new T[n]; } // Constructor
  ~MiVector() { delete[] x; } // Destructor
  int get_n() const { return n; }
  T operator[](int i) const { return x[i]; } // Leer elemento i
  T& operator[](int i) { return x[i]; } // Modificar el elemento i
private:
  int n; // Número de elementos
  T* x;  // x es un array de n elementos de tipo T
};

int main() {
  using namespace std;

  X x(2);
  print_x(x);
  print_x_ptr(&x);
  print_x_ref(x);

  vector<int> v; // Vector de la biblioteca estándar

  MiVector<int> w(20);
  w[3] = 77;
  cout << w[3] << endl; // w.operator[](3)
}
