#include <algorithm>
#include <array>
#include <cstddef>
#include <exception>
#include <initializer_list>
#include <memory>
#include <thread>

namespace blib {
/*
 * The Data Structure of Matrix
 * row 1: [1,2,3,4,5]
 * row 2: [1,2,3,4,5]
 * row 3: [1,2,3,4,5]*/
template <typename T> class bmatrix {
private:
  T **data_block;
  struct {
    size_t row;
    size_t column;
  } size;

public:
  ~bmatrix() {
    for (size_t i = 0; i < size.row; i++) {
      delete[] data_block[i];
    }
    delete[] data_block;
  }
  // Construction Functions
  bmatrix(const size_t row,
          const size_t
              column); // A new matrix of all zero with several rows and columns
  bmatrix(const size_t row, const size_t column,
          const T **src); // Copy from a 2D array
  template <size_t row_number, size_t column_number>
  bmatrix(const std::array<std::array<T, column_number>,
                           row_number> &&); // Start from an array
  bmatrix();                                // An empty matrix
  bmatrix(const size_t row, const size_t column,
          const std::initializer_list<std::initializer_list<T>>
              &&); // Start from an initializer_list
  // Friends about Construction
  friend bmatrix<T> &make_row(const size_t len, const T *src);
  friend bmatrix<T> &make_column(const size_t len, const T *src);
  // Change the value of matrix
  T &at(size_t row, size_t column); // Change one value of the matrix
  bmatrix<T> &operator=(const bmatrix<T> &&);
  // Basic Operators of Matrices
  friend bmatrix<T> &operator*(const bmatrix<T> &&, const bmatrix<T> &&);
  friend bmatrix<T> &operator+(const bmatrix<T> &&, const bmatrix<T> &&);
  friend bmatrix<T> &operator-(const bmatrix<T> &&, const bmatrix<T> &&);
  friend double operator/(const bmatrix<T> &&, const bmatrix<T> &&);
};
} // namespace blib

template <typename T>
blib::bmatrix<T>::bmatrix(const size_t row, const size_t column) {
  this->size.row = row;
  this->size.column = column;
  this->data_block = new T *[row];
  for (size_t i = 0; i < row; i++) {
    data_block[i] = new T[column]();
  }
}

template <typename T>
blib::bmatrix<T>::bmatrix(const size_t row, const size_t column,
                          const T **src) {
  this->size.row = row;
  this->size.column = column;
  this->data_block = new T *[row];
  for (size_t i = 0; i < row; i++) {
    data_block[i] = new int[column];
    std::copy(*(src + i), *(src + i) + column, this->data_block[i]);
  }
}

template <typename T>
template <size_t row, size_t column>
blib::bmatrix<T>::bmatrix(const std::array<std::array<T, column>, row> &&src) {
  size.row = row, size.column = column;
  data_block = new T *[row];
  for (auto it : src) {
    *data_block = new T[column];
    std::copy(it.begin(), it.end(), *(data_block++));
    data_block -= row;
  }
}

template <typename T> T &blib::bmatrix<T>::at(size_t row, size_t column) {
  if (row < this->size.row && column < this->size.column) {
    return data_block[row][column];
  } else {
    struct : std::exception {
      const char *what() const noexcept override {
        return "Array out of range";
      }
    } e;
    throw e;
  }
}

template <typename T>
blib::bmatrix<T> &operator*(blib::bmatrix<T> &&a, blib::bmatrix<T> &&b) {
  if (a.size.column == b.size.row) {
    auto res = *new blib::bmatrix<T>;
    res.size.row = a.size.row, res.size.column = b.size.column;
    size_t dimension = a.size.column;
    std::thread *th = new std::thread[res.size.column * res.size.row];
    for (size_t i = 0; i < res.size.row; i++) {
      for (size_t j = 0; j < res.size.column; j++) {
        th[i * res.size.column + j] = std::thread(
            [&](size_t row, size_t column) {
              T r = 0;
              for (size_t i = 0; i < dimension; i++) {
                r += a.at(row, i) * b.at(i, column);
              }
              res.data_block[row][column] = r;
            },
            i, j);
      }
    }
    for (size_t i = 0; i < res.size.row; i++) {
      for (size_t j = 0; j < res.size.column; j++) {
        th[i * res.size.column + j].join();
      }
    }
	delete[] th;
  } else {
    struct : std::exception {
      const char *what() const noexcept override {
        return "Two matrices cannot be multiplied.";
      }
    } e;
    throw e;
  }
}
