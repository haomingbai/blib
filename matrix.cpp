#include "bmat.h"
#include "bmat.cpp"

// 显式实例化 bmatrix<double> 类型
template class blib::bmatrix<double>;

namespace
{
    typedef blib::bmatrix<double> matrix;
}

// 显式实例化相关的非成员模板函数
template blib::bmatrix<double> blib::identity<double>(size_t);

// 显式实例化运算符重载
template blib::bmatrix<double> operator*(const blib::bmatrix<double>& a, const blib::bmatrix<double>& b);
template blib::bmatrix<double> operator-(const blib::bmatrix<double>& a, const blib::bmatrix<double>& b);
template blib::bmatrix<double> operator+(const blib::bmatrix<double>& a, const blib::bmatrix<double>& b);
