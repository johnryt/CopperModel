import numpy as np
#from diptest import _diptest
import warnings
import os

import numpy as np
import numpy as np

EMSGS = {1:'n must be >= 1',
         2:'x must be sorted in ascending order'}

# defined in _dip.c
#def externnogil():

#    void diptst(const double* x,
 #               const int* n_,
  #              double* dip,
   #             int* lo_hi,
    #            int* ifault,
     #           int* gcm,
      #          int* lcm,
       #         int* mn,
        #        int* mj,
         #       const int* min_is_0,
          #      const int* debug)

# wrapper to expose 'diptst' C function to Python space
def _dip(x, full_output, min_is_0, x_is_sorted, debug):

#    def runner():
 #       double[:] x_
  #      int n = x.shape[0]
   #     double dip = np.nan
    #    int[:] lo_hi = np.empty(4, dtype=np.int32)
     #   int ifault = 0
      #  int[:] gcm = np.empty(n, dtype=np.int32)
       # int[:] lcm = np.empty(n, dtype=np.int32)
        #int[:] mn = np.empty(n, dtype=np.int32)
        #int[:] mj = np.empty(n, dtype=np.int32)
        #int min_is_0_ = min_is_0
        #int debug_ = debug

    # input needs to be sorted in ascending order
    if not x_is_sorted:

        # force a copy to prevent inplace modification of input
        x = x.copy()

        # sort inplace
        x.sort()

    # cast to double
    #x_ = x.astype(np.double)

    #diptst(&x_[0], &n, &dip, &lo_hi[0], &ifault, &gcm[0], &lcm[0], &mn[0],
    #       &mj[0], &min_is_0_, &debug_)

   # if ifault:
    #    raise ValueError(EMSGS[ifault])

    if full_output:
        res_dict = {
            'xs':np.array(x_),
            'n':n,
            'dip':dip,
            'lo':lo_hi[0],
            'hi':lo_hi[1],
            'xl':x_[lo_hi[0]],
            'xu':x_[lo_hi[1]],
            'gcm':np.array(gcm[:lo_hi[2]]),
            'lcm':np.array(lcm[:lo_hi[3]]),
            'mn':np.array(mn),
            'mj':np.array(mj),
        }
        return dip, res_dict

    else:
        return dip


def dip(x, full_output=False, min_is_0=True, x_is_sorted=False, debug=0):
    """
    Hartigan & Hartigan's dip statistic

    The dip statistic measures multimodality in a sample by the maximum
    difference, over all sample points, between the empirical distribution
    function, and the unimodal distribution function that minimizes that
    maximum difference.

    Arguments:
    -----------
    x:              [n,] array  containing the input data

    full_output:    boolean, see below

    min_is_0:       boolean, if True the minimum value of the test statistic is
                    allowed to be zero in cases where n <= 3 or all values in x
                    are identical.

    x_is_sorted:    boolean, if True x is assumed to already be sorted in
                    ascending order

    debug:          int, 0 <= debug <= 3, print debugging messages

    Returns:
    -----------
    dip:    double, the dip statistic

    [res]:  dict, returned if full_output == True. contains the following
            fields:

            xs:     sorted input data as doubles
            n:      len(x)
            dip:    dip statistic
            lo:     indices of lower end of modal interval
            hi:     indices of upper end of modal interval
            xl:     lower end of modal interval
            xu:     upper end of modal interval
            gcm:    (last-used) indices of the greatest concave majorant
            lcm:    (last-used) indices of the least concave majorant

    Reference:
    -----------
        Hartigan, J. A., & Hartigan, P. M. (1985). The Dip Test of Unimodality.
        The Annals of Statistics.
    """

    return _dip(x, full_output, min_is_0, x_is_sorted, debug)


def diptest(x, min_is_0=True, boot_pval=False, n_boot=2000):
    """
    Hartigan & Hartigan's dip test for unimodality.

    For X ~ F i.i.d., the null hypothesis is that F is a unimodal distribution.
    The alternative hypothesis is that F is multimodal (i.e. at least bimodal).
    Other than unimodality, the dip test does not assume any particular null
    distribution.

    Arguments:
    -----------
    x:          [n,] array  containing the input data

    min_is_0:   boolean, see docstring for dip()

    boot_pval:  if True the p-value is computed using bootstrap samples from a
                uniform distribution, otherwise it is computed via linear
                interpolation of the tabulated critical values in dip_crit.txt.

    n_boot:     if boot_pval=True, this sets the number of bootstrap samples to
                use for computing the p-value.

    Returns:
    -----------
    dip:    double, the dip statistic

    pval:   double, the p-value for the test

    Reference:
    -----------
        Hartigan, J. A., & Hartigan, P. M. (1985). The Dip Test of Unimodality.
        The Annals of Statistics.

    """
    n = x.shape[0]
    D = dip(x, full_output=False, min_is_0=min_is_0)

    if n <= 3:
        warnings.warn('Dip test is not valid for n <= 3')
        pval = 1.0

    elif boot_pval:

        # random uniform vectors
        boot_x = np.random.rand(n_boot, n)

        # faster to pre-sort
        boot_x.sort(axis=1)
        boot_D = np.empty(n_boot)

        for ii in xrange(n_boot):
            boot_D[ii] = dip(boot_x[ii], full_output=False,
                             min_is_0=min_is_0, x_is_sorted=True)

        pval = np.mean(D <= boot_D)

    else:

        i1 = N.searchsorted(n, side='left')
        i0 = i1 - 1

        # if n falls outside the range of tabulated sample sizes, use the
        # critical values for the nearest tabulated n (i.e. treat them as
        # 'asymptotic')
        i0 = max(0, i0)
        i1 = min(N.shape[0] - 1, i1)

        # interpolate on sqrt(n)
        n0, n1 = N[[i0, i1]]
        fn = float(n - n0) / (n1 - n0)
        y0 = np.sqrt(n0) * CV[i0]
        y1 = np.sqrt(n1) * CV[i1]
        sD = np.sqrt(n) * D

        pval = 1. - np.interp(sD, y0 + fn * (y1 - y0), SIG)

    return D, pval


# [len(N), len(SIG)] table of critical values
curdir = os.path.dirname(os.path.realpath(__file__))
CV = np.loadtxt(os.path.join(curdir, 'dip_crit.txt'))

N = np.array([    4,     5,     6,     7,     8,     9,    10,    15,    20,
                 30,    50,   100,   200,   500,  1000,  2000,  5000, 10000,
              20000, 40000, 72000])

SIG = np.array([ 0.     ,  0.01   ,  0.02   ,  0.05   ,  0.1    ,  0.2    ,
                 0.3    ,  0.4    ,  0.5    ,  0.6    ,  0.7    ,  0.8    ,
                 0.9    ,  0.95   ,  0.98   ,  0.99   ,  0.995  ,  0.998  ,
                 0.999  ,  0.9995 ,  0.9998 ,  0.9999 ,  0.99995,  0.99998,
                 0.99999,  1.     ])