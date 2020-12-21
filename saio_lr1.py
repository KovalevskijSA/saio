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
    ab = a[:, jb]
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
        if koplan[i] > -eps:
            j_plus.append(i)
        else:
            j_minus.append(i)
    iter = 0
    while True:
        print('iter ', iter)
        for i in j_minus:
            kapa[i] = d_sup[i]
        for i in j_plus:
            kapa[i] = d_inf[i]
        kapa_b = np.dot(ab_inv, b - np.dot(a[:,j_n], [kapa[i] for i in j_n]))
        for i, j in enumerate(jb):
            kapa[j] = kapa_b[i]

        print('kapa_b', kapa_b)
        print('kapa', kapa)
        j_mins = [j for j in jb if kapa[j] < d_inf[j] or kapa[j] > d_sup[j]]
        if not j_mins:
            print("\nОптимальный план: \n", kapa)
            print('значение целевой функции: ', np.dot(c, kapa))
            return kapa, np.dot(c, kapa)
        j_k = j_mins[0]
        mu_jk = -1
        if kapa[j_k] < d_inf[j_k]:
            mu_jk = 1
        e = np.zeros(m)
        k = jb.index(j_k)
        e[k] = 1
        dy = np.dot(e.dot(mu_jk), ab_inv)
        mu_j = np.dot(dy, a)
        s = [np.inf] * n
        for j in j_n:
            if (mu_j[j] < 0 and j in j_plus) or (
                    mu_j[j] > 0 and j in j_minus):
                s[j] = -koplan[j] / mu_j[j]
        j0 = np.argmin(s)
        print('step=', s)
        if s[j0] == np.inf:
            print("задача не имеет решения")
            return None
        koplan = koplan + np.dot(s[j0], mu_j)
        if j0 in j_plus:
            j_plus.remove(j0)
        if mu_jk == 1:
            j_plus.append(jb[k])

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


def main():
    # prim 1
    c = [3, 2, 0, 3, -2, -4]
    b = [2, 5, 0]
    a = [[2, 1, -1, 0, 0, 1],
         [1, 0,  1, 1, 0, 0],
         [0, 1,  0, 0, 1, 0]]
    jb = [3, 4, 5]
    d_inf = [0, -1, 2, 1, -1, 0]
    d_sup = [2, 4, 4, 3, 3, 5]
    #
    # # task 1
    # c = [7, -2, 6, 0, 5, 2]
    # b = [-7, 22, 30]
    # a = [[1, -5, 3, 1, 0, 0],
    #      [4, -1, 1, 0, 1, 0],
    #      [2, 4, 2, 0, 0, 1]]
    # jb = [3, 4, 5]
    # d_inf = [2, 1, 0, 0, 1, 1]
    # d_sup = [6, 6, 5, 2, 4, 6]

    # task 2
    # c = [3, 0.5, 4, 4, 1, 5]
    # b = [15, 0, 13]
    # a = [[1, 0, 2, 2, -3, 3],
    #      [0, 1, 0, -1, 0, 1],
    #      [1, 0, 1, 3, 2, 1]]
    # jb = [3, 4, 5]
    # d_inf = [0, 0, 0, 0, 0, 0]
    # d_sup = [3, 5, 4, 3, 3, 4]

    # task 3
    # c = [2, 1, -2, -1, 4, -5, 5, 5]
    # b = [40, 107, 61]
    # a = [[1, 0, 0, 12, 1, -3, 4, -1],
    #      [0, 1, 0, 11, 12, 3, 5, 3],
    #      [0, 0, 1, 1,  0, 22, -2, 1]]
    # jb = [3, 4, 5]
    # d_inf = [0, 0, 0, 0, 0, 0, 0, 0]
    # d_sup = [3, 5, 5, 3, 4, 5, 6, 3]

    # task 4
    # c = [-1, 5, -2, 4, 3, 1, 2, 8, 3]
    # b = [3, 9, 9, 5, 9]
    # a = [[1, -3, 2, 0, 1, -1, 4, -1, 0],
    #      [1, -1, 6, 1, 0, -2, 2, 2, 0],
    #
    #      [2, -2, -1, 1, 0, -3, 8, -1, 1],
    #      [4, 1, 0, 0, 1, -1, 0, -1, 1],
    #      [1, 1, 1, 1, 1, 1, 1, 1, 1]]
    # d_inf = [0, 0, 0, 0, 0, 0, 0, 0, 0]
    # d_sup = [5, 5, 5, 5, 5, 5, 5, 5, 5]

    # task 5
    # c = [1, 2, 1, -3, 3, 1, 0]
    # b = [1, 4, 7]
    # a = [[1, 7, 2, 0, 1, -1, 4],
    #      [0, 5, 6, 1, 0, -3, -2],
    #      [3, 2, 2, 1, 1, 1, 5]]
    # jb = [1, 3, 5]
    # d_inf = [-1, 1, -2, 0, 1, 2, 4]
    # d_sup = [3, 2, 2, 5, 3, 4, 5]

    # # task 6
    # c = [0, 1, 2, 1, -3, 4, 7]
    # b = [1.5, 9, 2]
    # a = [[2, -1, 1, 0, 0, -1, 3],
    #      [0, 4, -1, 2, 3, -2, 2],
    #      [3, 1, 0, 1, 0, 1, 4]]
    # jb = [3, 4, 5]
    # d_inf = [0, 0, -3, 0, -1, 1, 0]
    # d_sup = [3, 3, 4, 7, 5, 3, 2]

    # # task 7
    # c = [0, -1, 1, 0, 4, 3]
    # b = [2, 2, 5]
    # a = [[2, 1, 0, 3, -1, 1],
    #      [0, 1, -2, 1, 0, 3],
    #      [3, 0, 1, 1, 1, 1 ]]
    # d_inf = [2, 0, -1, -3, 2, 1]
    # d_sup = [7, 3, 2, 3, 4, 5]

    # # task 8
    # c = [2, -1, 2, 3, -2, 3, 4, 1]
    # b = [4, 12, 4]
    # a = [[1, 3, 1, -1, 0, -3, 2, 1],
    #      [2, 1, 3, -1, 1, 4, 1, 1],
    #      [-1, 0, 2, -2, 2, 1, 1, 1]]
    # d_inf = [-1, -1, -1, -1, -1, -1, -1, -1]
    # d_sup = [ 2,  3,  1,  4,  3,  2,  4,  4]

    solve_lp(np.array(a), b, c, d_inf, d_sup)


if __name__ == '__main__':
    main()
