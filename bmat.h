#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <initializer_list>
#include <memory>
#include <thread>
#include <vector>
#include <utility>
#include <stdexcept>
#include <type_traits>
#include <tuple>

namespace blib
{
  template <typename T>
  class bmatrix;

  template <typename T>
  bmatrix<T> identity(size_t N);
  /*
   * Definition of Matrix
   * The Data Structure of Matrix
   * row 1: [1,2,3,4,5]
   * row 2: [1,2,3,4,5]
   * row 3: [1,2,3,4,5]*/
  template <typename T>
  class bmatrix
  {
  private:
    T **data_block;
    struct S
    {
      size_t row;
      size_t column;
    } size;

  public:
    ~bmatrix()
    {
      for (size_t i = 0; i < size.row; i++)
      {
        delete[] data_block[i];
      }
      delete[] data_block;
    }
    // Construction Functions
    bmatrix(size_t row,
            size_t
                column); // A new matrix of all zero with several rows and columns
    bmatrix(size_t row, size_t column,
            T **src); // Copy from a 2D array
    template <size_t row_number, size_t column_number>
    bmatrix(std::array<std::array<T, column_number>,
                       row_number> &&); // Start from an array
    bmatrix(size_t row, size_t column,
            std::initializer_list<std::initializer_list<T>>
                &&); // Start from an initializer_list
    template <typename P>
    bmatrix(bmatrix<P> &);
    // Friends about Construction
    // friend bmatrix<T> &make_row( size_t len,  T *src);
    // friend bmatrix<T> &make_column( size_t len,  T *src);
    // Change the value of matrix
    T &at(size_t row, size_t column) const; // Change one value of the matrix
    // template<typename P> friend double operator/( bmatrix<P> &,  bmatrix<P> &);
    //  Functions
    bmatrix<double> rref();
    double norm(size_t);
    T det();
    std::tuple<size_t, size_t> get_size() const;
    friend bmatrix<T> identity<T>(size_t N);
  };

  template<typename T>
  bmatrix<T> identity(size_t N);
} // namespace blib

template <typename T1, typename T2>
auto operator*(const blib::bmatrix<T1> &a, const blib::bmatrix<T2> &b) -> blib::bmatrix<decltype(T1{} * T2{})>;

template <typename T1, typename T2>
auto operator-(const blib::bmatrix<T1> &a, const blib::bmatrix<T2> &b) -> blib::bmatrix<decltype(T1{} - T2{})>;

template <typename T1, typename T2>
auto operator+(const blib::bmatrix<T1> &a, const blib::bmatrix<T2> &b) -> blib::bmatrix<decltype(T1{} + T2{})>;

