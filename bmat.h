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
#include <type_traits>

namespace blib
{
  template <typename T>
  class bmatrix;

  template <typename T>
  bmatrix<T> operator*(bmatrix<T> &, bmatrix<T> &);
  template <typename T>
  bmatrix<T> operator+(bmatrix<T> &, bmatrix<T> &);
  template <typename T>
  bmatrix<T> operator-(bmatrix<T> &, bmatrix<T> &);
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
    T &at(size_t row, size_t column); // Change one value of the matrix
    // bmatrix<T> &operator=( bmatrix<T> &&);
    //  Basic Operators of Matrices
    friend bmatrix<T> operator* <>(bmatrix<T> &, bmatrix<T> &);
    friend bmatrix<T> operator+ <>(bmatrix<T> &, bmatrix<T> &);
    friend bmatrix<T> operator- <>(bmatrix<T> &, bmatrix<T> &);
    // template<typename P> friend double operator/( bmatrix<P> &,  bmatrix<P> &);
    //  Functions
    bmatrix<double> rref();
    double norm(size_t);
    T det();
    struct S get_size()
    {
      return this->size;
    }
    friend bmatrix<T> identity<>(size_t N);
  };
} // namespace blib

// Realize the function of matrix
template <typename T>
blib::bmatrix<T>::bmatrix(size_t row, size_t column)
{
  this->size.row = row;
  this->size.column = column;
  this->data_block = new T *[row];
  for (size_t i = 0; i < row; i++)
  {
    data_block[i] = new T[column]();
  }
}

template <typename T>
blib::bmatrix<T>::bmatrix(size_t row, size_t column,
                          T **src)
{
  this->size.row = row;
  this->size.column = column;
  this->data_block = new T *[row];
  for (size_t i = 0; i < row; i++)
  {
    data_block[i] = new int[column];
    std::copy(*(src + i), *(src + i) + column, this->data_block[i]);
  }
}

template <typename T>
template <size_t row, size_t column>
blib::bmatrix<T>::bmatrix(std::array<std::array<T, column>, row> &&src)
{
  size.row = row, size.column = column;
  data_block = new T *[row];
  for (auto it : src)
  {
    *data_block = new T[column];
    std::copy(it.begin(), it.end(), *(data_block++));
    data_block -= row;
  }
}

template <typename T>
T &blib::bmatrix<T>::at(size_t row, size_t column)
{
  if (row < this->size.row && column < this->size.column)
  {
    return data_block[row][column];
  }
  else
  {
    struct : std::exception
    {
      const char *what() const noexcept override
      {
        return "Array out of range";
      }
    } e;
    throw e;
  }
}

template <typename T>
blib::bmatrix<T> blib::operator*(blib::bmatrix<T> &a, blib::bmatrix<T> &b)
{
  if (a.size.column == b.size.row)
  {
    auto &res = *new blib::bmatrix<T>(a.size.row, b.size.column);
    res.size.row = a.size.row, res.size.column = b.size.column;
    size_t dimension = a.size.column, r_row = res.size.row,
           r_column = res.size.column;
    auto th = std::make_unique<std::thread[]>(res.size.row);
    for (size_t i = 0; i < r_row; i++)
    {
      th[i] = std::thread(
          [r_column, &a, &b, &res, dimension](size_t current_row)
          {
            auto th = std::make_unique<std::thread[]>(res.size.row);
            for (size_t j = 0; j < r_column; j++)
            {
              th[j] = std::thread(
                  [&a, &b, &res, current_row,
                   dimension](size_t current_column)
                  {
                    T result = 0;
                    for (int i = 0; i < dimension; i++)
                    {
                      result = result +
                               a.at(current_row, i) * b.at(i, current_column);
                    }
                    res.at(current_row, current_column) = result;
                  },
                  j);
            }
            for (size_t j = 0; j < r_column; j++)
            {
              th[j].join();
            }
          },
          i);
    }
    for (size_t i = 0; i < res.size.row; i++)
    {
      th[i].join();
    }
    return std::move(res);
  }
  else
  {
    struct : std::exception
    {
      const char *what() const noexcept override
      {
        return "Two matrices cannot be multiplied.";
      }
    } e;
    throw e;
  }
}

template <typename T>
blib::bmatrix<T> blib::operator+(blib::bmatrix<T> &a, blib::bmatrix<T> &b)
{
  if (a.size.row == b.size.row && a.size.column == b.size.column)
  {
    auto &res = *new blib::bmatrix<T>(a.size.row, a.size.column);
    res.size.row = a.size.row, res.size.column = a.size.column;
    auto th = std::make_unique<std::thread[]>(res.size.row);
    size_t r_row = a.size.row;
    for (size_t i = 0; i < r_row; i++)
    {
      th[i] = std::thread(
          [&a, &b, &res](size_t current_row)
          {
            for (size_t i = 0; i < a.size.column; i++)
            {
              res.at(current_row, i) =
                  a.at(current_row, i) + b.at(current_row, i);
            }
          },
          i);
    }
    for (size_t i = 0; i < r_row; i++)
    {
      th[i].join();
    }
    return std::move(res);
  }
  else
  {
    struct : std::exception
    {
      const char *what() const noexcept override
      {
        return "Two matrices should have the same row/column";
      }
    } e;
    throw e;
  }
}

template <typename T>
blib::bmatrix<T> blib::operator-(blib::bmatrix<T> &a, blib::bmatrix<T> &b)
{
  if (a.size.row == b.size.row && a.size.column == b.size.column)
  {
    auto &res = *new blib::bmatrix<T>(a.size.row, a.size.column);
    res.size.row = a.size.row, res.size.column = a.size.column;
    auto th = std::make_unique<std::thread[]>(res.size.row);
    size_t r_row = a.size.row;
    for (size_t i = 0; i < r_row; i++)
    {
      th[i] = std::thread(
          [&a, &b, &res](size_t current_row)
          {
            for (size_t i = 0; i < a.size.column; i++)
            {
              res.at(current_row, i) =
                  a.at(current_row, i) - b.at(current_row, i);
            }
          },
          i);
    }
    for (size_t i = 0; i < r_row; i++)
    {
      th[i].join();
    }
    return std::move(res);
  }
  else
  {
    struct : std::exception
    {
      const char *what() const noexcept override
      {
        return "Two matrices should have the same row/column";
      }
    } e;
    throw e;
  }
}

template <typename T>
blib::bmatrix<double> blib::bmatrix<T>::rref()
{
  blib::bmatrix<double> &res = *new blib::bmatrix<double>(size.row, size.column);

  // Step 1: Copy elements to res using threads
  for (int i = 0; i < size.row; ++i)
  {
    auto th = std::make_unique<std::thread[]>(size.column);
    for (int j = 0; j < size.column; ++j)
    {
      th[j] = std::thread([&res, i, j, this]()
                          { res.at(i, j) = static_cast<double>(this->at(i, j)); });
    }
    for (int j = 0; j < size.column; ++j)
    {
      if (th[j].joinable())
        th[j].join();
    }
  }

  size_t rows = res.get_size().row;
  size_t cols = res.get_size().column;

  for (size_t row = 0; row < rows; ++row)
  {
    // Step 2: Find the pivot row and swap with current row if needed
    if (std::abs(res.at(row, row)) < 1e-10)
    {
      for (size_t t = row + 1; t < rows; ++t)
      {
        if (std::abs(res.at(t, row)) >= 1e-10)
        {
          for (size_t p = 0; p < cols; ++p)
          {
            std::swap(res.at(row, p), res.at(t, p));
          }
          break;
        }
      }
    }

    // Step 3: Normalize the pivot row
    double pivot = res.at(row, row);
    if (std::abs(pivot) >= 1e-10)
    {
      for (size_t j = 0; j < cols; ++j)
      {
        res.at(row, j) /= pivot;
      }
    }

    // Step 4: Eliminate other rows using threads
    auto th = std::make_unique<std::thread[]>(rows);
    for (size_t i = 0; i < rows; ++i)
    {
      if (i != row)
      {
        double factor = res.at(i, row);
        th[i] = std::thread([&, i, row, factor]()
                            {
                    for (size_t j = 0; j < cols; ++j) {
                        res.at(i, j) -= factor * res.at(row, j);
                    } });
      }
    }
    for (size_t i = 0; i < rows; ++i)
    {
      if (i != row && th[i].joinable())
        th[i].join();
    }
  }

  return std::move(res);
}

template <typename T>
double blib::bmatrix<T>::norm(size_t p)
{
  if (this->size.row == 1)
  {
    double res;
    for (int i = 0; i < this->size.column; i++)
    {
      res += std::pow(std::abs(this->at(0, i)), p);
    }
    res = std::pow(res, 1 / p);
    return res;
  }
  else if (this->size.column == 1)
  {
    double res;
    for (int i = 0; i < this->size.row; i++)
    {
      res += std::pow(std::abs(this->at(i, 0)), p);
    }
    res = std::pow(res, 1 / p);
    return res;
  }
  else
  {
    struct : std::exception
    {
      const char *what() const noexcept override
      {
        return "The matrix should be a row or a column";
      }
    } e;
    throw e;
  }
}

template <typename T>
T blib::bmatrix<T>::det()
{
  if (this->size.row == this->size.column)
  {
    if (this->size.row == 1)
    {
      return this->at(0, 0);
    }
    else if (this->size.row == 2)
    {
      return this->at(0, 0) * this->at(1, 1) - this->at(0, 1) * this->at(1, 0);
    }
    else
    {
      T res = 0;
      for (int i = 0; i < this->size.column; i++)
      {
        blib::bmatrix<T> tmp(this->size.row - 1, this->size.column - 1);
        for (int j = 1; j < this->size.row; j++)
        {
          for (int k = 0; k < this->size.column; k++)
          {
            if (k < i)
            {
              tmp.at(j - 1, k) = this->at(j, k);
            }
            else
            {
              tmp.at(j - 1, k) = this->at(j, k + 1);
            }
          }
        }
        res += std::pow(-1, i) * this->at(0, i) * tmp.det();
      }
      return res;
    }
  }
  else
  {
    struct : std::exception
    {
      const char *what() const noexcept override
      {
        return "The matrix should be a square matrix";
      }
    } e;
    throw e;
  }
}

#include <iostream>

template <typename T>
blib::bmatrix<T> blib::identity(size_t N)
{
  auto &res = *new blib::bmatrix<T>(N, N);
  for (int i = 0; i < N; i++)
  {
    res.at(i, i) = 1;
  }
  return std::move(res);
}

template <typename T>
template <typename P>
blib::bmatrix<T>::bmatrix(blib::bmatrix<P> &src)
{
  this->size.row = src.get_size().row;
  this->size.column = src.get_size().column;
  this->data_block = new T *[size.row];
  for (size_t i = 0; i < size.row; i++)
  {
    data_block[i] = new T[size.column];
  }
  auto th = std::make_unique<std::thread[]>(size.row);
  for (size_t i = 0; i < size.row; i++)
  {
    th[i] = std::thread([&, i]()
                        {
            for (size_t j = 0; j < size.column; j++)
            {
                this->at(i, j) = static_cast<T>(src.at(i, j));
            } });
  }
  for (size_t i = 0; i < size.row; i++)
  {
    th[i].join();
  }
}
