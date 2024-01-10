//
//  mat_mkl.cpp
//  MSDLimit
//
//  Created by Yi Hu on 10/1/18.
//
#include <stdlib.h>
#include <string.h>
#include <mkl.h>
#include <mkl_spblas.h>
#include <mkl_solvers_ee.h>
#include <vector>
#include <set>

#include "ant_hbd.hpp"

#if __INTEL_MKL__ >= 2017
    #define MKL_CAT 1
#else
    #error MKL library not found or version under 2017
    #define MKL_CAT 0
#endif

// Transfer matrix method of solving dynamics
int HybridAntWalker::run_matrix(const MapType& siterec, std::vector<LeathSiteNode> &sitelist)
{
  MKL_INT csize = (MKL_INT)sitelist.size(); // conversion from size_t to int
  MKL_INT trisize, count_all_links, count_all_nzs, rval;
  mat_limit_type rp, rq, row_count;
  int retval;
  int rdim, rtime;
  double dtmp, eigi;
  //double *mat_sd = NULL; // square distance
  // Matrix storage
  double *mat_vals = NULL;
  MKL_INT *mat_cols = NULL, *mat_rows = NULL;
#if MKL_CAT == 1
  sparse_matrix_t mat_handle;
  const struct matrix_descr descs = {SPARSE_MATRIX_TYPE_SYMMETRIC, // .type
    SPARSE_FILL_MODE_UPPER, // .mode
    SPARSE_DIAG_NON_UNIT, // .diag
  };
#else
  double emin = 0.0, emax = 1.0;
  MKL_INT ret_loop = 0;
#endif
  // Eigen solver
  char ls_which = 'L';
  MKL_INT pm[129] = {0};
  
  MKL_INT k0 = (csize <= 10? csize: (csize<=30? 5+csize/2: 20)), ret_k;
  double *ret_E, *ret_X, *ret_res;
  double *cterms;
  unsigned short *count_link = NULL;
  LeathSiteNode tmpsite;
  MapType::const_iterator itbuff;
  std::set<mat_limit_type> links[csize];
  count_link = (unsigned short*)malloc(csize*sizeof(unsigned short));
  memset(count_link, 0, csize*sizeof(unsigned short));
  count_all_links = 0;
  
  // Obtain link information
  for(rp=0; rp<csize; ++rp)
  {
    for(rdim = 0; rdim<DIM; ++rdim)
    {
      tmpsite = sitelist[rp] + unitvc[rdim];
      itbuff = siterec.find(tmpsite);
      if(itbuff != siterec.end())
      {
        rq = itbuff->second;
        if(rq != 0)
        {
          --rq; // indices in siterec start with 1
          ++count_link[rp];
          ++count_link[rq];
          ++count_all_links;
          if(rp < rq)
          {
            links[rp].insert(rq);
          }
          else
          {
            links[rq].insert(rp);
          }
        }
      }
    }
  }
  // Construct transfer matrix
  count_all_nzs = count_all_links + csize;
  mat_vals = (double*)mkl_malloc(count_all_nzs*sizeof(double), 64);
  mat_cols = (MKL_INT*)mkl_malloc(count_all_nzs*sizeof(MKL_INT), 64);
  mat_rows = (MKL_INT*)mkl_malloc((1+csize)*sizeof(MKL_INT), 64);
  rval = 0;
  for(rp=0; rp<csize; ++rp)
  {
    
    mat_vals[rval] = 2*DIM-count_link[rp];
    mat_cols[rval] = rp;
    mat_rows[rp] = rval;
    ++rval;
    for(std::set<mat_limit_type>::const_iterator sit=links[rp].begin(); sit!=links[rp].end(); ++sit)
    {
      mat_vals[rval] = 1;
      mat_cols[rval] = *sit;
      ++rval;
    }
    links[rp].clear();
  }
  mat_rows[csize] = rval;
  // Create sparse matrix object

#ifdef DEBUG
  /*
  for(rp=0; rp<count_all_nzs; ++rp)
  {
    std::cout << mat_vals[rp] << " ";
  }
  std::cout << "\n";
  for(rp=0; rp<count_all_nzs; ++rp)
  {
    std::cout << mat_cols[rp] << " ";
  }
  std::cout << "\n";
  for(rp=0; rp<=csize; ++rp)
  {
    std::cout << mat_rows[rp] << " ";
  }
  std::cout << std::endl;
   */
#endif
#if MKL_CAT == 1
  sparse_status_t spflag = mkl_sparse_d_create_csr(&mat_handle, SPARSE_INDEX_BASE_ZERO, csize, csize,
                           mat_rows, mat_rows+1, mat_cols, mat_vals);
  if(spflag != SPARSE_STATUS_SUCCESS)
  {
    std::cerr << "Sparse matrix creation failed" <<std::endl;
    retval = -1;
    goto lbl_ret_A;
  }
#endif
  // Solve eigenvalue problem
  // spflag = mkl_sparse_ee_init(pm); // extended input parameter
  feastinit(pm);
  ret_E = (double*)mkl_malloc(k0*sizeof(double), 64);
  ret_X = (double*)mkl_malloc(k0*csize*sizeof(double), 64);
  ret_res = (double*)mkl_malloc(k0*sizeof(double), 64);
#if MKL_CAT == 1
  // version 2017 and later
  spflag = mkl_sparse_d_ev(&ls_which, pm, mat_handle, descs, k0, &ret_k, ret_E, ret_X, ret_res);
  if(spflag != SPARSE_STATUS_SUCCESS)
  {
    std::cerr << "Solve matrix of size " << csize << " failed with code " << spflag <<std::endl;
    retval = -2;
    goto lbl_ret_B;
  }
#else
  // version before 2017
  // TODO
#endif
  // Calculate const terms
  cterms = (double*)mkl_malloc(ret_k*sizeof(double), 64);
  memset(cterms, 0, ret_k*sizeof(double));
  // constant terms.
  for(rp=0; rp<csize; ++rp)
  {
    // skip A[rp][rp] which = 0
    for(rq=rp+1; rq<csize; ++rq)
    {
      // D^2
      dtmp = (sitelist[rq] - sitelist[rp]).norm_squared();
      for(row_count=0; row_count<ret_k; ++row_count)
      {
        cterms[row_count] += dtmp*ret_X[row_count*csize+rp]*ret_X[row_count*csize+rq];
      }
    }
  }
  // rescaling eigenvalues
  for(rp=0; rp<ret_k; ++rp)
  {
    ret_E[rp] /= (2*DIM);

  }
  eigi = ret_E[0];
#ifdef DEBUG
  std::cout << csize << " " << cterms[ret_k-1]*2/csize << " " << ret_E[0] << " ";
  /*for(rp=0; rp<ret_k; rp+=2)
  {
      std::cout << ret_E[rp] << " ";
  }
  std::cout << "\n";*/
#endif
  for(row_count=0; row_count<ret_k; ++row_count)
    cterms[row_count] *= 2.0/csize;
  // write to tmparray
  for(rtime=0; rtime<arraylen; ++rtime)
  {
    tmparray[rtime] = cterms[ret_k-1];
    for(row_count=0; row_count<ret_k-1; ++row_count)
    {
      tmparray[rtime] += ret_E[row_count]*cterms[row_count];
      ret_E[row_count] *= ret_E[row_count];
    }
  }
#ifdef DEBUG
/*for(rp=0; rp<ret_k; rp+=2)
{
  std::cout << cterms[rp] << " ";
}*/
  std::cout << cterms[0] << "\n";
#endif
  retval = (int)(-log(1-eigi)/log(2.0));
  mkl_free(cterms);
lbl_ret_B:
  mkl_free(ret_res);
  mkl_free(ret_X);
  mkl_free(ret_E);
lbl_ret_A:
  mkl_sparse_destroy(mat_handle);
  //mkl_free(mat_sd);
  mkl_free(mat_rows);
  mkl_free(mat_cols);
  mkl_free(mat_vals);

  free(count_link);
  return retval;
}
