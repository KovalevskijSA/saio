import math
import numpy as np
from numpy.linalg import linalg
import itertools as it

eps = 1e-10


def inv_matrix(A_inv, i, x, n):
    l = np.dot(A_inv, x)
    li = l[i]
    if l[i] == 0:
        return
    l[i] = -1
    l_ = [(-1 / li) * x for x in l]
    Q = np.eye(n)
    for j in range(0, n):
        Q[j][i] = l_[j]
    ans = np.dot(Q, A_inv)
    return ans


def dual_simplex_method(a, b, c, jb, d_inf, d_sup):
    m, n = len(a), len(c)
    ab = (np.array([a[:,j] for j in jb])).transpose()
    ab_inv = linalg.inv(ab)
    y = np.dot([c[i] for i in jb], ab_inv)
    koplan = np.dot(y, a) - c
    if any([i < -eps for i in koplan]):
        return None
    kapa = [0]*n
    j_n = [i for i in range(n) if i not in jb]
    j_plus = []
    j_minus = []
    for i in j_n:
        if koplan[i] >= 0:
            j_plus.append(i)
        else:
            j_minus.append(i)
    iter = 0
    while True:
        for i in j_minus:
            kapa[i] = d_sup[i]
        for i in j_plus:
            kapa[i] = d_inf[i]
        a_t = (np.array([a[:,i] for i in j_n])).transpose()
        kapa_b = np.dot(ab_inv, b - np.dot(a_t, [kapa[i] for i in j_n]))
        for i, j in enumerate(jb):
            kapa[j] = kapa_b[i]
        j_mins = [j for j in jb if kapa[j] < d_inf[j] or kapa[j] > d_sup[j]]
        if not j_mins:
            return kapa, np.dot(c, kapa), jb
        j_min = min(j_mins)
        mu_jk = -1
        if kapa[j_min] < d_inf[j_min]:
            mu_jk = 1
        e = np.zeros(m)
        k = jb.index(j_min)
        e[k] = 1
        dy = np.dot(e.dot(mu_jk), ab_inv)
        mu_j = np.dot(dy, a)
        s = [np.inf] * n
        for j in j_n:
            if (mu_j[j] < 0 and j in j_plus) or (
                    mu_j[j] > 0 and j in j_minus):
                s[j] = -koplan[j] / mu_j[j]
        j0 = np.argmin(s)
        if s[j0] == np.inf:
            return None
        koplan = koplan + np.dot(s[j0], mu_j)
        if j_plus.count(j0) == 0:
            if mu_jk == 1:
                j_plus.append(jb[k])
        else:
            index = j_plus.index(j0)
            if mu_jk == 1:
                j_plus[index] = jb[k]
            elif mu_jk == -1:
                j_plus.pop(index)
        jb[k] = j0
        j_n = [j for j in range(n) if j not in jb]
        j_minus = [j for j in j_n if j not in j_plus]
        ab_inv = inv_matrix(ab_inv, k, a[:,j0], m)
        iter += 1


def get_basis_index(a, n):
    m = len(a)
    j = list(range(n))
    for i in it.combinations(j, m):
        ii = list(i)
        ab = (np.array([a[:, j] for j in ii])).transpose()
        if np.linalg.det(ab) != 0:
            yield ii


def solve_lp(a, b, c, d_inf, d_sup):
    for jb in get_basis_index(a, len(c)):
        sol = dual_simplex_method(a, b, c, jb, d_inf, d_sup)
        if sol is not None:
            return sol
    return None


def method_branch_border(a, b, c, d_inf, d_sup, r):
    task_list = []
    res_list = []
    res_item = None

    task_list.append([d_inf, d_sup])
    br = 0
    iter = 0
    while(task_list):
        next = task_list.pop(0)
        sol = solve_lp(a, b, c, next[0], next[1])
        if not sol or sol[1] < r:
            continue
        x = sol[0]
        # j = get_index_of_fraction_item(x)
        # if j == -1:
        j = [i for i in range(len(x)) if abs(x[i] - math.floor(x[i] + 0.5)) > eps]
        if not j:
            res_item = sol
            res_list.append(sol)
            r = sol[1]
        else:
            j = j[0]
            l = math.floor(x[j])
            d_inf_t, d_sup_t = next[0][:], next[1][:]
            d_sup_t[j] = l
            print("iter{0} task1 {1} {2}".format(iter, d_sup_t, d_inf_t))
            task_list.append([d_inf_t, d_sup_t])
            d_inf_t, d_sup_t = next[0][:], next[1][:]
            d_inf_t[j] = l + 1
            task_list.append([d_inf_t, d_sup_t])
            print("iter{0} task1 {1} {2}".format(iter, d_sup_t, d_inf_t))
            iter += 1
            br += 2
    # for i in res_list:
    #     print('f= {0}\nx= {1}\n'.format(i[0], i[3]))
    print('ver1')
    print(br)
    if res_item:
        print('f= {0}\nx= {1}\n'.format(res_item[1], res_item[0]))
    else:
        print('Нет решений')


def main():
    # prim 1
    # c = [7, -2, 6, 0, 5, 2]
    # b = [-8, 22, 30]
    # a = [[1, -5, 3, 1, 0, 0],
    #      [4, -1, 1, 0, 1, 0],
    #      [2, 4,  2, 0, 0, 1]]
    # d_inf = [2, 1, 0, 0, 1, 1]
    # d_sup = [6, 6, 5, 2, 4, 6]
    #
    # # prim2
    # c = [7, -2, 6, 0, 5, 2]
    # b = [10, 8, 10]
    # a = [[1, 0, 3, 1, 0, 0],
    #      [0, -1, 1, 1, 1, 2],
    #      [-2, 4,  2, 0, 0, 1]]
    # d_inf = [0, 1, -1, 0, -2, 1]
    # d_sup = [3, 3, 6, 2, 4, 6]
    #
    # # task 1
    # c = [2, 1, -2, -1, 4, -5, 5, 5]
    # b = [40, 107, 61]
    # a = [[1, 0, 0, 12, 1, -3, 4, -1],
    #      [0, 1, 0, 11, 12, 3, 5, 3],
    #      [0, 0, 1, 1, 0, 22, -2, 1]]
    # d_inf = [0, 0, 0, 0, 0, 0, 0, 0]
    # d_sup = [3, 5, 5, 3, 4, 5, 6, 3]

    # # task 2
    c = [-1, 5, -2, 4, 3, 1, 2, 8, 3]
    b = [3, 9, 9, 5, 9]
    a = [[1, -3, 2, 0, 1, -1, 4, -1, 0],
         [1, -1, 6, 1, 0, -2, 2,  2, 0],
         [2,  2,-1, 1, 0, -3, 8, -1, 1],
         [4,  1, 0, 0, 1, -1, 0, -1, 1],
         [1,  1, 1, 1, 1,  1, 1,  1, 1]]
    d_inf = [0, 0, 0, 0, 0, 0, 0, 0, 0]
    d_sup = [5, 5, 5, 5, 5, 5, 5, 5, 5]

    # task 3
    # c = [2, 1, -2, -1, 4, -5, 5, 5, 1, 2]
    # b = [43.5, 107.3, 106.3]
    # a = [[1, 0, 0, 12, 1, -3, 4, -1, 2.5, 3],
    #      [0, 1, 0, 11, 12, 3, 5, 3, 4, 5.1],
    #      [0, 0, 1, 1,  0, 22, -2, 1, 6.1, 7]]
    # d_inf = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    # d_sup = [2, 4, 5, 3, 4, 5, 4, 4, 5, 6]

    # # task 4
    # c = [2, 1, -2, -1, 4, -5, 5, 5, 1, 2]
    # b = [8, 5, 4, 7, 8]
    # a = [[4, 0, 0, 0, 0, -3, 4, -1, 2, 3],
    #      [0, 1, 0, 0, 0,  3, 5, 3, 4,  5],
    #      [0, 0, 1, 0, 0, 22,-2, 1, 6,  7],
    #      [0, 0, 0, 1, 0, 6, -2, 7, 5, 6],
    #      [0, 0, 0, 0, 1, 5, 5, 1, 6, 7]]
    # d_inf = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    # d_sup = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10]
    #
    # task 5
    # c = [7, -2, 6, 0, 5, 2]
    # b = [-8, 22, 30]
    # a = [[1, -5, 3, 1, 0, 0],
    #      [4, -1, 1, 0, 1, 0],
    #      [2, 4, 2, 0, 0, 1]]
    # d_inf = [2, 1, 0, 0, 1, 1]
    # d_sup = [6, 6, 5, 2, 4, 6]
    #
    # # task 6
    # c = [2, 1, -2, -1, 4, -5, 5, 5]
    # b = [30, 78, 5]
    # a = [[1, 0, 0, 3, 1, -3, 4, -1],
    #      [0, 1, 0, 4, -3, 3, 5, 3],
    #      [0, 0, 1, 1, 0, 2, -2, 1]]
    # d_inf = [0, 0, 0, 0, 0, 0, 0, 0]
    # d_sup = [5, 5, 3, 4, 5, 6, 6, 8]

    # task 7
    # c = [7, 5, -2, 4, 3, 1, 2, 8, 3]
    # b = [18, 18, 30, 15, 18]
    # a = [[1, -3, 2, 0, 1, -1, 4, -1, 0],
    #      [1, -1, 6, 1, 0, -2, 2, 2, 0],
    #      [2, 2, -1, 1, 0, -3, 2, -1, 1],
    #      [4, 1, 0, 0, 1, -1, 0, -1, 1],
    #      [1, 1, 1, 1, 1,  1, 1,  1, 1]]
    # d_inf = [0, 0, 0, 0, 0, 0, 0, 0, 0]
    # d_sup = [8, 8, 8, 8, 8, 8, 8, 8, 8]

    # task 8
    # c = [1, 2, 3, -1, 4, -5, 6]
    # b = [26, 185, 32.5]
    # a = [[1, 0, 1, 0, 4, 3, 4],
    #      [0, 1, 2, 0, 55, 3.5, 5],
    #      [0, 0, 3, 1, 6, 2, -2.5]]
    # d_inf = [0, 1, 0, 0, 0, 0, 0]
    # d_sup = [1, 2, 5, 7, 8, 4, 2]

    # task 9
    # c = [1, 2, 3, 1, 2, 3, 4]
    # b = [58, 66.3, 36.7, 13.5]
    # a = [[2, 0, 1, 0, 0, 3, 5],
    #      [0, 2, 2.1, 0, 0, 3.5, 5],
    #      [0, 0, 3, 2, 0, 2, 1.1],
    #      [0, 0, 3, 0, 2, 2, -2.5]]
    # d_inf = [1, 1, 1, 1, 1, 1, 1]
    # d_sup = [2, 4, 4, 5, 8, 7, 7]

    # task 10
    # c = [-2, 1, -2, -1, 8, -5, 3, 5, 1, 2]
    # b = [27, 6, 18]
    # a = [[1, 0, 0, 1, 1, -3, 4, -1, 3, 3],
    #      [0, 1, 0, -2, 1, 1, 7, 3, 4, 5],
    #      [0, 0, 1, 1, 0, 2, -2, 1, -4, 7]]
    # d_inf = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    # d_sup = [8, 7, 6, 7, 8, 5, 6, 7, 8, 5]

    a = np.array(a)
    method_branch_border(a, b, c, d_inf, d_sup, -np.inf)

if __name__ == '__main__':
    main()
