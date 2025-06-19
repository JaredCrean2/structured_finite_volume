#include "bla_wrapper.h"

#include <string>
#include <stdexcept>
#include <cassert>
#include <iostream>
#include "Kokkos_Core.hpp"

namespace structured_fv {
namespace linear_system {

namespace {

void check_error(const char* fname, int ret_val)
{
  if (ret_val != 0)
    throw std::runtime_error(std::string("Function ") + fname + " return value " + std::to_string(ret_val));
}

template <typename T>
lapack_int size(const std::vector<T>& v)
{
  return v.size();
}

enum CBLAS_TRANSPOSE get_cblas_trans(char trans)
{
  switch (trans)
  {
    case 'N': return CblasNoTrans;
    case 'n': return CblasNoTrans;
    case 'T': return CblasTrans;
    case 't': return CblasTrans;
    case 'C': return CblasConjTrans;
    case 'c': return CblasConjTrans;
    default:
      throw std::runtime_error(std::string("unrecognized transpose specifier: ") + trans);
  };
}

enum CBLAS_UPLO get_cblas_uplo(char uplo)
{
  switch(uplo)
  {
    case 'U': return CblasUpper;
    case 'u': return CblasUpper;
    case 'L': return CblasLower;
    case 'l': return CblasLower;
    default:
      throw std::runtime_error(std::string("unrecognized upper/lower specifier: ") + uplo);
  };
}

}

//-----------------------------------------------------------------------------
// GEMV
void gemv(char trans, blasint m, blasint n, double alpha, const double* A, blasint lda, 
          const double* x, blasint incx, double beta, double* y, blasint incy)
{
  cblas_dgemv(CblasColMajor, get_cblas_trans(trans), m, n, alpha, A, lda, x, incx, beta, y, incy);
}

void gemv(char trans, blasint m, blasint n, float alpha, const float* A, blasint lda, 
          const float* x, blasint incx, float beta, float* y, blasint incy)
{
  cblas_sgemv(CblasColMajor, get_cblas_trans(trans), m, n, alpha, A, lda, x, incx, beta, y, incy);
}

template <typename T>
void gemv(char trans, blasint m, blasint n, T alpha, const std::vector<T>& A, const std::vector<T>& x, T beta, std::vector<T>& y)
{
  assert(size(A) == m * n);
  assert(size(x) == n);
  assert(size(y) == m);

  gemv(trans, m, n, alpha, A.data(), m, x.data(), 1, beta, y.data(), 1);
}

template
void gemv<double>(char trans, blasint m, blasint n, double alpha, const std::vector<double>& A, const std::vector<double>& x, double beta, std::vector<double>& y);

template
void gemv<float>(char trans, blasint m, blasint n, float alpha, const std::vector<float>& A, const std::vector<float>& x, float beta, std::vector<float>& y);

template <typename T>
void gemv(char trans,  T alpha, ConstKokkosMat<T> A, ConstKokkosVec<T> x, T beta, KokkosVec<T> y)
{
  assert(A.extent(1) == x.extent(0));
  assert(A.extent(0) == y.extent(0));
  blasint lda = A.extent(0);
  blasint inc = 1;
  return gemv(trans, A.extent(0), A.extent(1), alpha, A.data(), lda, x.data(), inc, beta, y.data(), inc);
}

void gemv(char trans, double alpha, ConstKokkosMat<double> A, ConstKokkosVec<double> x, double beta, KokkosVec<double> y)
{
  return gemv<double>(trans, alpha, A, x, beta, y);
}

void gemv(char trans, float alpha, ConstKokkosMat<float> A, ConstKokkosVec<float> x, float beta, KokkosVec<float> y)
{
  return gemv<float>(trans, alpha, A, x, beta, y);
}

//-----------------------------------------------------------------------------
// SYMV
void symv(char uplo, blasint n, double alpha, const double* A, blasint lda, 
          const double* x, blasint incx, double beta, double* y, blasint incy)
{
  cblas_dsymv(CblasColMajor, get_cblas_uplo(uplo), n, alpha, A, lda, x, incx, beta, y, incy);
}

void symv(char uplo, blasint n, float alpha, const float* A, blasint lda, 
          const float* x, blasint incx, float beta, float* y, blasint incy)
{
  cblas_ssymv(CblasColMajor, get_cblas_uplo(uplo), n, alpha, A, lda, x, incx, beta, y, incy);
}

template <typename T>
void symv(char uplo, blasint n, T alpha, const std::vector<T>& A, const std::vector<T>& x, T beta, std::vector<T>& y)
{
  assert(size(A) == n * n);
  assert(size(x) == n);
  assert(size(y) == n);

  symv(uplo, n, alpha, A.data(), n, x.data(), 1, beta, y.data(), 1);
}

template
void symv<double>(char uplo, blasint n, double alpha, const std::vector<double>& A, const std::vector<double>& x, double beta, std::vector<double>& y);

template
void symv<float>(char uplo, blasint n, float alpha, const std::vector<float>& A, const std::vector<float>& x, float beta, std::vector<float>& y);

template <typename T>
void symv(char uplo, T alpha, ConstKokkosMat<T> A, ConstKokkosVec<T> x, T beta, KokkosVec<T> y)
{
  assert(A.extent(1) == x.extent(0));
  assert(A.extent(0) == y.extent(0));
  blasint lda = A.extent(0);
  blasint inc = 1;
  return symv(uplo, A.extent(0), alpha, A.data(), lda, x.data(), inc, beta, y.data(), inc);
}

void symv(char uplo, double alpha, ConstKokkosMat<double> A, ConstKokkosVec<double> x, double beta, KokkosVec<double> y)
{
  return symv<double>(uplo, alpha, A, x, beta, y);
}

void symv(char uplo, float alpha, ConstKokkosMat<float> A, ConstKokkosVec<float> x, float beta, KokkosVec<float> y)
{
  return symv<float>(uplo, alpha, A, x, beta, y);
}

//-----------------------------------------------------------------------------
// GESV

lapack_int gesv(lapack_int n, lapack_int nrhs, double* a, lapack_int lda, lapack_int* ipiv, double* b, lapack_int ldb)
{
  auto ret_val = LAPACKE_dgesv_work(LAPACK_COL_MAJOR, n, nrhs, a, lda, ipiv, b, ldb);
  check_error("dgesv", ret_val);

  return ret_val;
}

lapack_int gesv(lapack_int n, lapack_int nrhs, float* a, lapack_int lda, lapack_int* ipiv, float* b, lapack_int ldb)
{
  auto ret_val = LAPACKE_sgesv_work(LAPACK_COL_MAJOR, n, nrhs, a, lda, ipiv, b, ldb);
  check_error("sgesv", ret_val);

  return ret_val;
}

template <typename T>
lapack_int gesv(lapack_int n, std::vector<T>& a, std::vector<lapack_int>& ipiv, std::vector<T>& b)
{
  assert(n > 0);
  assert(n * n == size(a));
  assert(size(ipiv) == n);
  assert(size(b) % n == 0);
  lapack_int lda = n;
  lapack_int ldb = n;
  lapack_int nrhs = size(b) / n;

  return gesv(n, nrhs, a.data(), lda, ipiv.data(), b.data(), ldb);
}

template
lapack_int gesv<double>(lapack_int n, std::vector<double>& a, std::vector<lapack_int>& ipiv, std::vector<double>& b);

template
lapack_int gesv<float>(lapack_int n, std::vector<float>& a, std::vector<lapack_int>& ipiv, std::vector<float>& b);

template <typename T>
lapack_int gesv(KokkosMat<T> a, KokkosVec<lapack_int> ipiv, KokkosVec<T> b)
{
  assert(a.extent(0) == a.extent(1));
  assert(a.extent(0) == b.extent(0));
  assert(a.extent(0) == ipiv.extent(0));
  lapack_int lda = a.extent(0);
  lapack_int ldb = b.extent(0);

  return gesv(a.extent(0), 1, a.data(), lda, ipiv.data(), b.data(), ldb);
}

lapack_int gesv(KokkosMat<double> a, KokkosVec<lapack_int> ipiv, KokkosVec<double> b)
{
  return gesv<double>(a, ipiv, b);
}

 
lapack_int gesv(KokkosMat<float> a, KokkosVec<lapack_int> ipiv, KokkosVec<float> b)
{
  return gesv<float>(a, ipiv, b);
}

//-----------------------------------------------------------------------------
// GETRF

lapack_int getrf(lapack_int m, lapack_int n, double* a, lapack_int lda, lapack_int* ipiv)
{
  auto ret_val = LAPACKE_dgetrf_work(LAPACK_COL_MAJOR, m, n, a, lda, ipiv);
  check_error("dgetrf", ret_val);

  return ret_val;
}

lapack_int getrf(lapack_int m, lapack_int n, float* a, lapack_int lda, lapack_int* ipiv)
{
  auto ret_val = LAPACKE_sgetrf_work(LAPACK_COL_MAJOR, m, n, a, lda, ipiv);
  check_error("sgetrf", ret_val);

  return ret_val;
}

template <typename T>
lapack_int getrf(lapack_int m, lapack_int n, std::vector<T>& a, std::vector<lapack_int>& ipiv)
{
  assert(m > 0);
  assert(n > 0);
  assert(m * n == size(a));
  assert(size(ipiv) == std::min(m, n));
  lapack_int lda = m;

  return getrf(m, n, a.data(), lda, ipiv.data());
}

template
lapack_int getrf<double>(lapack_int m, lapack_int n, std::vector<double>& a, std::vector<lapack_int>& ipiv);

template
lapack_int getrf<float>(lapack_int m, lapack_int n, std::vector<float>& a, std::vector<lapack_int>& ipiv);

template <typename T>
lapack_int getrf(KokkosMat<T> a, KokkosVec<lapack_int> ipiv)
{
  assert(ipiv.extent(0) == std::min(a.extent(0), a.extent(1)));
  lapack_int lda = a.extent(0);
  return getrf(a.extent(0), a.extent(1), a.data(), lda, ipiv.data());
}


lapack_int getrf(KokkosMat<double> a, KokkosVec<lapack_int> ipiv)
{
  return getrf<double>(a, ipiv);
}

lapack_int getrf(KokkosMat<float> a, KokkosVec<lapack_int> ipiv)
{
  return getrf<float>(a, ipiv);
}

//-----------------------------------------------------------------------------
// GETRS

lapack_int getrs(char trans, lapack_int n, lapack_int nrhs, double* a, lapack_int lda, const lapack_int* ipiv, double* b, lapack_int ldb)
{
  auto ret_val = LAPACKE_dgetrs_work(LAPACK_COL_MAJOR, trans, n, nrhs, a, lda, ipiv, b, ldb);
  check_error("dgetrs", ret_val);

  return ret_val;
}

lapack_int getrs(char trans, lapack_int n, lapack_int nrhs, float* a, lapack_int lda, const lapack_int* ipiv, float* b, lapack_int ldb)
{
  auto ret_val = LAPACKE_sgetrs_work(LAPACK_COL_MAJOR, trans, n, nrhs, a, lda, ipiv, b, ldb);
  check_error("sgetrs", ret_val);

  return ret_val;
}

template <typename T>
lapack_int getrs(char trans, lapack_int n, std::vector<T>& a, const std::vector<lapack_int>& ipiv, std::vector<T>& b)
{
  assert(n > 0);
  assert(n * n == size(a));
  assert(size(ipiv) == n);
  assert(size(b) % n == 0);
  lapack_int lda = n;
  lapack_int ldb = n;
  lapack_int nrhs = size(b) / n;

  return getrs(trans, n, nrhs, a.data(), lda, ipiv.data(), b.data(), ldb);
}

template
lapack_int getrs<double>(char trans, lapack_int n, std::vector<double>& a, const std::vector<lapack_int>& ipiv, std::vector<double>& b);

template
lapack_int getrs<float>(char trans, lapack_int n, std::vector<float>& a, const std::vector<lapack_int>& ipiv, std::vector<float>& b);

template <typename T>
lapack_int getrs(char trans, KokkosMat<T> a, ConstKokkosVec<lapack_int> ipiv, KokkosVec<T> b)
{
  assert(a.extent(0) == a.extent(1));
  assert(ipiv.extent(0) == a.extent(0));
  assert(b.extent(0) == a.extent(0));
  lapack_int lda = a.extent(0);
  lapack_int ldb = b.extent(0);

  return getrs(trans, a.extent(0), 1, a.data(), lda, ipiv.data(), b.data(), ldb);
}


lapack_int getrs(char trans, KokkosMat<double> a, ConstKokkosVec<lapack_int> ipiv, KokkosVec<double> b)
{
  return getrs<double>(trans, a, ipiv, b);
}


lapack_int getrs(char trans, KokkosMat<float> a, ConstKokkosVec<lapack_int> ipiv, KokkosVec<float> b)
{
  return getrs<float>(trans, a, ipiv, b);
}

//-----------------------------------------------------------------------------
// POTRF

lapack_int potrf(char uplo, lapack_int n, double* a, lapack_int lda)
{
  auto ret_val = LAPACKE_dpotrf_work(LAPACK_COL_MAJOR, uplo, n, a, lda);
  check_error("dpotrf", ret_val);

  return ret_val;
}

lapack_int potrf(char uplo, lapack_int n, float* a, lapack_int lda)
{
  auto ret_val = LAPACKE_spotrf_work(LAPACK_COL_MAJOR, uplo, n, a, lda);
  check_error("spotrf", ret_val);

  return ret_val;
}

template <typename T>
lapack_int potrf(char uplo, lapack_int n, std::vector<T>& a)
{
  assert(n > 0);
  assert(size(a) == n*n);
  lapack_int lda = n;

  return potrf(uplo, n, a.data(), lda);
}

template
lapack_int potrf<double>(char uplo, lapack_int n, std::vector<double>& a);

template
lapack_int potrf<float>(char uplo, lapack_int n, std::vector<float>& a);

template <typename T>
lapack_int potrf(char uplo, KokkosMat<T> a)
{
  lapack_int lda = a.extent(0);
  return potrf(uplo, a.extent(0), a.data(), lda);
}


lapack_int potrf(char uplo, KokkosMat<double> a)
{
  return potrf<double>(uplo, a);
}

lapack_int potrf(char uplo, KokkosMat<float> a)
{
  return potrf<float>(uplo, a);
}

//-----------------------------------------------------------------------------
// POTRS

lapack_int potrs(char uplo, lapack_int n, lapack_int nrhs, const double* a, lapack_int lda, double* b, lapack_int ldb)
{
  auto ret_val = LAPACKE_dpotrs_work(LAPACK_COL_MAJOR, uplo, n, nrhs, a, lda, b, ldb);
  check_error("dpotrs", ret_val);

  return ret_val;
}

lapack_int potrs(char uplo, lapack_int n, lapack_int nrhs, const float* a, lapack_int lda, float* b, lapack_int ldb)
{
  auto ret_val = LAPACKE_spotrs_work(LAPACK_COL_MAJOR, uplo, n, nrhs, a, lda, b, ldb);
  check_error("spotrs", ret_val);

  return ret_val;
}

template <typename T>
lapack_int potrs(char uplo, lapack_int n, const std::vector<T>& a, std::vector<T>& b)
{
  assert( n > 0);
  assert(size(a) == n*n);
  assert(size(b) % n == 0);
  lapack_int lda = n;
  lapack_int ldb = n;
  lapack_int nrhs = size(b) / n;


  return potrs(uplo, n, nrhs, a.data(), lda, b.data(), ldb);
}

template
lapack_int potrs<double>(char uplo, lapack_int n, const std::vector<double>& a, std::vector<double>& b);

template
lapack_int potrs<float>(char uplo, lapack_int n, const std::vector<float>& a, std::vector<float>& b);

template <typename T>
lapack_int potrs(char uplo, ConstKokkosMat<T> a, KokkosVec<T> b)
{
  assert(a.extent(0) == a.extent(1));
  assert(a.extent(0) == b.extent(0));
  lapack_int lda = a.extent(0);
  lapack_int ldb = b.extent(0);

  return potrs(uplo, a.extent(0), 1, a.data(), lda, b.data(), ldb);
}

 
lapack_int potrs(char uplo, ConstKokkosMat<double> a, KokkosVec<double> b)
{
  return potrs<double>(uplo, a, b);
}
 
lapack_int potrs(char uplo, ConstKokkosMat<float> a, KokkosVec<float> b)
{
  return potrs<float>(uplo, a, b);
}

} // namespace
}