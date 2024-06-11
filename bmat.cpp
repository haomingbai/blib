#include "bmat.h"

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
T &blib::bmatrix<T>::at(size_t row, size_t column) const
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

  //size_t rows = res.get_size().row;
  //size_t cols = res.get_size().column;
  auto [rows,cols] = res.get_size();

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
  //this->size.row = src.get_size().row;
  //this->size.column = src.get_size().column;
  auto [row,column] = src.get_size();
  this->size.row = row;
  this->size.column = column;
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

template <typename T1, typename T2>
auto operator+(const blib::bmatrix<T1> &a, const blib::bmatrix<T2> &b) -> blib::bmatrix<decltype(T1{} + T2{})>
{
  //if (a.get_size().row == b.get_size().row && a.get_size().column == b.get_size().column)
  if(a.get_size() == b.get_size())
  {
    using ResultType = decltype(T1{} + T2{});
    auto [r,c] = a.get_size();
    auto res = blib::bmatrix<ResultType>(r,c);
    auto th = std::make_unique<std::thread[]>(r);
    for (size_t i = 0; i < r; i++)
    {
      th[i] = std::thread([&, i]()
                          {
                for (size_t j = 0; j < c; j++)
                {
                    res.at(i, j) = a.at(i, j) + b.at(i, j);
                } });
    }
    for (size_t i = 0; i < r; i++)
    {
      th[i].join();
    }
    return res;
  }
  else
  {
    throw std::invalid_argument("Matrices must have the same dimensions for addition.");
  }
}

template <typename T1, typename T2>
auto operator-(const blib::bmatrix<T1> &a, const blib::bmatrix<T2> &b) -> blib::bmatrix<decltype(T1{} - T2{})>
{
  if (a.get_size() = b.get_size())
  {
    auto [r,c] = a.get_size();
    using ResultType = decltype(T1{} - T2{});
    auto res = blib::bmatrix<ResultType>(r,c);
    auto th = std::make_unique<std::thread[]>(r);
    for (size_t i = 0; i < r; i++)
    {
      th[i] = std::thread([&, i]()
                          {
                for (size_t j = 0; j < c; j++)
                {
                    res.at(i, j) = a.at(i, j) - b.at(i, j);
                } });
    }
    for (size_t i = 0; i < r; i++)
    {
      th[i].join();
    }
    return res;
  }
  else
  {
    throw std::invalid_argument("Matrices must have the same dimensions for subtraction.");
  }
}
template <typename T1, typename T2>
auto operator*(const blib::bmatrix<T1> &a, const blib::bmatrix<T2> &b) -> blib::bmatrix<decltype(T1{} * T2{})>
{
  //if (a.get_size().column == b.get_size().row)
  auto [r0,c0] = a.get_size();
  auto [r1,c1] = b.get_size();
  if(c0==r1)
  {
    using ResultType = decltype(T1{} * T2{});
    auto res = blib::bmatrix<ResultType>(r0, c1);
    auto th = std::make_unique<std::thread[]>(r0);
    for (size_t i = 0; i < r0; i++)
    {
      th[i] = std::thread([&, i]()
                          {
                for (size_t j = 0; j < c1; j++)
                {
                    ResultType sum = 0;
                    for (size_t k = 0; k < c0; k++)
                    {
                        sum += a.at(i, k) * b.at(k, j);
                    }
                    res.at(i, j) = sum;
                } });
    }
    for (size_t i = 0; i < r0; i++)
    {
      th[i].join();
    }
    return res;
  }
  else
  {
    throw std::invalid_argument("Number of columns in the first matrix must equal the number of rows in the second matrix for multiplication.");
  }
}

template <typename T>
std::tuple<size_t,size_t> blib::bmatrix<T>::get_size() const
{
  return std::make_tuple(size.row, size.column);
}