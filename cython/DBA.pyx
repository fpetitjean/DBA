'''
/*******************************************************************************
 * Copyright (C) 2018 Francois Petitjean
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
'''
from __future__ import division
import numpy as np
cimport numpy as np
from cpython cimport array
import array

cdef inline int min_c_int(int a, int b): return a if a <= b else b
cdef inline int max_c_int(int a, int b): return a if a >= b else b

__author__ ="Francois Petitjean"


def performDBA(double[:,:] series,unsigned char n_iterations=10, unsigned int w = -1):
    cdef:
        int n_series = len(series)
        int max_length = 0,i,length_padded,medoid_ind
        double [:,:] cost_mat,cost_mat_dtw,delta_mat,delta_map_padded
        signed char[:,:] path_mat
        double [:] center

    for i in range(n_series):
        max_length = max_c_int(max_length,series[i].shape[0])
    if w == -1:
        w = max_length

    cost_mat = np.zeros((max_length, max_length),dtype=np.float64)
    delta_mat = np.zeros((max_length, max_length),dtype=np.float64)
    path_mat = np.zeros((max_length, max_length), dtype=np.int8)

    medoid_ind = approximate_medoid_index(series,w,cost_mat,delta_mat)

    center = series[medoid_ind]

    for i in range(n_iterations):
        center = DBA_update(center, series, w, cost_mat, path_mat, delta_mat)

    return center

cdef int approximate_medoid_index(double[:,:] series, int w, double[:,:] cost_mat, double[:,:] delta_mat):
    cdef:
        int n_series = len(series)
        int n_samples_candidates = min_c_int(10,n_series)
        int n_samples_ss = min_c_int(20,n_series)
        int index_candidate,medoid_ind = -1,i
        double best_ss,ss
        long[:] indices_candidates = np.random.choice(range(n_series),n_samples_candidates,replace=False)
        long[:] indices_ss = np.random.choice(range(n_series),n_samples_ss,replace=False)
        double[:] candidate

    for i in range(indices_candidates.shape[0]):
        index_candidate = indices_candidates[i]
        candidate = series[index_candidate]
        ss = sum_of_squares_partial(candidate,series,indices_ss,w,cost_mat,delta_mat)
        if(medoid_ind==-1 or ss<best_ss):
            best_ss = ss
            medoid_ind = index_candidate
    return medoid_ind

cdef double sum_of_squares_partial(double[:] s,double[:,:] series,long[:] indices_for_sum,int w,double[:,:] cost_mat,double[:,:] delta_mat):
    cdef:
        double sum = 0
        int i
    for i in indices_for_sum:
        sum += squared_DTW(s,series[i],w,cost_mat,delta_mat)
    return sum

cdef double DTW(double[:] s,double[:] t,int w,double[:,:] cost_mat,double[:,:] delta_mat):
    return np.sqrt(squared_DTW(s,t,w,cost_mat,delta_mat))

cdef double squared_DTW(double[:] s, double[:] t, int w, double[:,:] cost_mat, double[:,:] delta_mat):
    cdef:
        int s_len = len(s)
        int t_len = len(t)
        int length = len(s)
        int i,j,stop_j,start_j,idx_infty_left
        double res
        double[:,:] delta_mat_dtw
        double diag,left,top

    '''
    Currently calculating the whole delta mat, which is not necessary, but just faster to program.
    This should be optimized later if speed required.
    '''
    fill_delta_mat_dtw(s, t, delta_mat)

    cost_mat[0, 0] = delta_mat[0, 0]
    for i in range(1, min_c_int(s_len,w+1)):
        cost_mat[i, 0] = cost_mat[i-1, 0]+delta_mat[i, 0]

    stop_j = min_c_int(t_len,w+1)
    for j in range(1, stop_j):
        cost_mat[0, j] = cost_mat[0, j-1]+delta_mat[0, j]
    if stop_j < t_len:
        cost_mat[0, stop_j] = np.inf

    for i in range(1, s_len):
        start_j = max_c_int(1, i-w)
        stop_j = min_c_int(t_len, i+w+1)
        idx_infty_left = i-w-1
        if idx_infty_left >=0 :
            cost_mat[i][idx_infty_left] = np.inf
        for j in range(start_j, stop_j):
            diag,left,top =cost_mat[i-1, j-1], cost_mat[i, j-1], cost_mat[i-1, j]

            if diag <= left :
                if diag <= top :
                    res = cost_mat[i-1, j-1]
                else:
                    res = cost_mat[i-1, j]
            else:
                if (left<=top):
                    res = cost_mat[i, j-1]
                else:
                    res = cost_mat[i-1, j]

            cost_mat[i, j] = res + delta_mat[i, j]
        if stop_j < t_len :
            cost_mat[i][stop_j] = np.inf
    return cost_mat[s_len-1,t_len-1]

cdef void fill_delta_mat_dtw(double[:]center, double[:]s, double[:,:] delta_mat):
    delta_numpy = np.asarray(delta_mat)
    np.subtract.outer(np.asarray(center), np.asarray(s), out=delta_numpy)
    np.square(delta_numpy, out=delta_numpy)


cdef double[:] DBA_update(double[:] center,double[:,:] series, unsigned int w, double[:,:] cost_mat,signed char[:,:] path_mat, double[:,:] delta_mat):
    cdef:
        signed char[:,:] options_argmin
        signed char[:] move
        int i,j,stop_j,start_j,idx_infty_left,s_len,center_length
        double diag,left,top,res
        double[:] updated_center,s
        int[:] n_elements

    options_argmin = np.array([[-1, -1], [0, -1], [-1, 0]],dtype=np.int8)
    center_length = center.shape[0]

    updated_center = np.zeros(center_length,dtype = np.float64)
    n_elements = np.array(np.zeros(center_length), dtype=np.intc)

    for s in series:
        s_len = s.shape[0]
        fill_delta_mat_dtw(center, s, delta_mat)
        cost_mat[0, 0] = delta_mat[0, 0]
        path_mat[0, 0] = -1

        for i in range(1, min_c_int(center_length,w+1)):
            cost_mat[i, 0] = cost_mat[i-1, 0]+delta_mat[i, 0]
            path_mat[i, 0] = 2

        stop_j = min_c_int(s_len,w+1)
        for j in range(1, stop_j):
            cost_mat[0, j] = cost_mat[0, j-1]+delta_mat[0, j]
            path_mat[0, j] = 1
        if stop_j < s_len:
            cost_mat[0, stop_j] = np.inf

        for i in range(1, center_length):
            start_j = max_c_int(1, i-w)
            stop_j = min_c_int(s_len, i+w+1)
            idx_infty_left = i-w-1
            if idx_infty_left >=0 :
                cost_mat[i][idx_infty_left] = np.inf
            for j in range(start_j, stop_j):
                diag,left,top =cost_mat[i-1, j-1], cost_mat[i, j-1], cost_mat[i-1, j]
                if(diag <=left):
                    if(diag<=top):
                        res = diag
                        path_mat[i,j] = 0
                    else:
                        res = top
                        path_mat[i,j] = 2
                else:
                    if(left<=top):
                        res = left
                        path_mat[i,j] = 1
                    else:
                        res = top
                        path_mat[i,j] = 2


                cost_mat[i, j] = res+delta_mat[i, j]
            if stop_j < s_len :
                cost_mat[i][stop_j] = np.inf

        i = center_length-1
        j = s_len-1

        while(path_mat[i, j] != -1):
            updated_center[i] += s[j]
            n_elements[i] += 1
            move = options_argmin[path_mat[i, j]]
            i += move[0]
            j += move[1]
        assert(i == 0 and j == 0)
        updated_center[0] += s[0]
        n_elements[0] += 1

    for i in range(center_length): #normalizing
        updated_center[i] /= n_elements[i]

    return updated_center
