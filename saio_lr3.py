import numpy as np
from numpy.linalg import linalg
import itertools as it
import math

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


def dual_simplex_method2(a, b, c, jb):
    m, n = len(a), len(c)
    ab = a[:, jb]
    ab_inv = linalg.inv(ab)
    y = np.dot([c[i] for i in jb], ab_inv)
    koplan = np.dot(y, a) - c
    if any([i < -eps for i in koplan]):
        return None

    while True:
        kapa_b = np.dot(ab_inv, b)
        if min(kapa_b) > -eps:
            kapa = [0] * n
            for j, i in zip(kapa_b, jb):
                kapa[i] = j
            return kapa, jb, np.dot(c, kapa)
        k = np.argmin(kapa_b)
        j_n = [j for j in range(n) if j not in jb]
        e = np.zeros(m)
        e[k] = 1
        mu = e.dot(ab_inv.dot(a))
        if all(mu[i] >= 0 for i in j_n):
            return None
        s = [np.inf] * n
        for j in j_n:
            if mu[j] < 0:
                s[j] = -koplan[j] / mu[j]
        j0 = np.argmin(s)
        y = y + (s[j0] * e).dot(b)
        koplan = koplan + s[j0] * mu
        jb[k] = j0
        ab_inv = inv_matrix(ab_inv, k, a[:, j0], m)


def get_basis_index(a, n):
    m = len(a)
    j = list(range(n))
    for i in it.combinations(j, m):
        ii = list(i)
        ab = (np.array([a[:, j] for j in ii])).transpose()
        if np.linalg.det(ab) != 0:
            yield ii


def solve_lp(a, b, c):
    for jb in get_basis_index(a, len(c)):
        sol = dual_simplex_method2(a, b, c, jb)  # return (x, jb, f) or None
        if sol is not None:
            return sol
    return None


def is_integer(i):
    return abs(i - math.floor(i + 0.5)) < eps


def correct_task(a, b, c, jb, j_r):
    n, m = len(c), len(b)
    p = n - m
    j_r.sort(reverse=True)
    print(j_r)
    j_rm = [j - p for j in j_r]
    for index in j_r:
        for i in range(index - p + 1, m):
            if i in j_rm or a[i][index] == 0:
                continue
            k = a[i][index] / a[index - p][index]
            b[i] -= k * b[index - p]
            a[i, :] -= k * a[index - p, :]

    for i in j_rm:
        del b[i]
    for i in j_r:
        jb.remove(i)
        del c[i]
    a = np.delete(a, j_r, axis=1)
    a = np.delete(a, j_rm, axis=0)
    return a, b, c, jb


def get_j_k(f_prev, f_cur, x, jb):
    if abs(f_prev - f_cur) > eps:
        j_k = max(
            [(j, x[j] - math.floor(x[j])) for j in jb if not is_integer(x[j])],
            key=lambda p: p[1],
            default=(-1, 0))[0]
    else:
        j_k = min([j for j in jb if not is_integer(x[j])], default=-1)
    return j_k


def method_gomory(a, b, c):
    f_prev = np.inf
    iter = 0
    n_p = len(c)
    sol = solve_lp(a, b, c)
    while(sol):
        m, n = len(b), len(c)
        x, jb, f = sol
        print(f'\niter {iter}\nx {x}\njb {jb}\nf {f}')
        if n > n_p + 3:
            j_r = [j for j in jb if j >= n_p]
            if j_r:
                print(f'delete, j {j_r}')
                a, b, c, jb = correct_task(a, b, c, jb, j_r)
                m, n = len(b), len(c)

        j_k = get_j_k(f_prev, f, x, [j for j in jb if j < n_p])
        print('j_k', j_k)
        if j_k == -1:
            print(f'\nsol\nx {x[:n_p]}\nf {f}')
            return

        k = jb.index(j_k)
        ab_inv = linalg.inv(a[:, jb])
        print('ab_inv', ab_inv)

        e = np.zeros(m)
        e[k] = 1
        y = np.dot(e, ab_inv)
        print('y', y)
        a_j = np.dot(y, a)
        print('a_j', a_j)
        betta = np.dot(y, b)
        print('betta', betta)

        f_j = [-(i - math.floor(i)) if not is_integer(i) else 0 for i in a_j]
        # f_j = [-(i % 1) if not is_integer(i) else 0 for i in a_j]
        if all([f_j[i] > -eps for i in range(n) if i not in jb]):
            print('нет решений')
            return
        e = np.zeros((m + 1, 1))
        e[m] = 1
        a = np.vstack([a, f_j])
        a = np.hstack([a, e])
        b.append(-(betta - math.floor(betta)))
        c.append(0)
        f_prev = f
        jb.append(n)
        sol = dual_simplex_method2(a, b, c, jb)
        iter += 1


def main():
    # # prim 1
    # c = [-3.5, 1, 0, 0, 0]
    # b = [15, 6, 0]
    # a = [[5, -1, 1, 0, 0],
    #      [-1, 2, 0, 1, 0],
    #      [-7, 2, 0, 0, 1]]
    #
    # # prim 2
    # c = [1, -1, 0, 0, 0]
    # b = [4, 3, 7]
    # a = [[5, 3, 1, 0, 0],
    #      [-1, 2, 0, 1, 0],
    #      [1, -2, 0, 0, 1]]
    #
    #
    # # task 1
    # c = [7, -2, 6, 0, 5, 2]
    # b = [-8, 22, 30]
    # a = [[1, -5, 3, 1, 0, 0],
    #      [4, -1, 1, 0, 1, 0],
    #      [2,  4, 2, 0, 0, 1]]
    #
    # #
    # # # task 2
    # c = [-1, 5, -2, 4, 3, 1, 2, 8, 3]
    # b = [3, 9, 9, 5, 9]
    # a = [[1, -3, 2, 0, 1, -1, 4, -1, 0],
    #      [1, -1, 6, 1, 0, -2, 2,  2, 0],
    #      [2,  2,-1, 1, 0, -3, 8, -1, 1],
    #      [4,  1, 0, 0, 1, -1, 0, -1, 1],
    #      [1,  1, 1, 1, 1,  1, 1,  1, 1]]
    #
    # # #
    # # # task 3
    # c = [2, 1, -2, -1, 4, -5, 5, 5]
    # b = [40, 107, 61]
    # a = [[1, 0, 0, 12, 1, -3, 4, -1],
    #      [0, 1, 0, 11, 12, 3, 5, 3],
    #      [0, 0, 1, 2, 0, 22, -2, 1]]
    #
    # # # task 4
    # c = [2, 1, -2, -1, 4, -5, 5, 5, 1, 2]
    # b = [153, 123, 112]
    # a = [[1, 2, 3, 12, 1, -3, 4, -1, 2, 3],
    #      [0, 2, 0, 11, 12, 3, 5, 3, 4, 5],
    #      [0, 0, 2, 1,  0, 22, -2, 1, 6, 7]]
    #
    # # task 4
    # c = [1, 2, 1, -1, 2, 3]
    # b = [7, 16, 6]
    # a = [[2, 1, -1, -3, 4, 7],
    #      [0, 1,  1,  1, 2, 4],
    #      [6, -3, -2, 1, 1, 1]]
    #
    # # #
    # # # # task 6
    # c = [10, 2, 1, 7, 6, 3, 1]
    # b = [12, 27, 19]
    # a = [[0, 7, 1, -1, -4, 2, 4],
    #      [5, 1, 4, 3, -5, 2, 1],
    #      [2, 0, 3, 1, 0, 1, 5]]
    #
    # # # task 7
    # c = [2, 9, 3, 5, 1, 2, 4]
    # b = [6, 3, 7, 7]
    # a = [[0, 7, -8, -1, 5, 2, 1],
    #      [3, 2, 1, -3, -1, 1, 0],
    #      [1, 5, 3, -1, -2, 1, 0],
    #      [1, 1, 1,  1,  1, 1, 1]]
    # #
    # # task 8
    # c = [-1, -3, -7, 0, -4, 0, -1]
    # b = [4, 8, 24]
    # a = [[1, 0, -1, 3, -2, 0, 1],
    #      [0, 2, 1, -1, 0, 3, -1],
    #      [1, 2, 1, 4, 2, 1, 1]]
    # #-16
    #
    # # # task 9
    # c = [-1, 5, -2, 4, 3, 1, 2, 8, 3]
    # b = [3, 9, 9, 5, 9]
    # a = [[1, -3, 2, 0, 1, -1, 4, -1, 0],
    #      [1, -1, 6, 1, 0, -2, 2, 2, 0],
    #      [2, 2, -1, 1, 0, -3, 2, -1, 1],
    #      [4, 1, 0, 0, 1, -1, 0, -1, 1],
    #      [1, 1, 1, 1, 1, 1,  1,  1, 1]]
    # 25
    # a = np.array(a)
    # method_gomory(a, b, c)

    tasks = {
        'prim 1': {
            'c': [-3.5, 1, 0, 0, 0],
            'b': [15, 6, 0],
            'a': [[5, -1, 1, 0, 0],
                  [-1, 2, 0, 1, 0],
                  [-7, 2, 0, 0, 1]]
        },
        # 'prim 2': {
        #     'c': [1, -1, 0, 0, 0],
        #     'b': [4, 3, 7],
        #     'a': [[5, 3, 1, 0, 0],
        #           [-1, 2, 0, 1, 0],
        #           [1, -2, 0, 0, 1]]
        # },

        # 'task 1': {
        #     'c': [7, -2, 6, 0, 5, 2],
        #     'b': [-8, 22, 30],
        #     'a': [[1, -5, 3, 1, 0, 0],
        #           [4, -1, 1, 0, 1, 0],
        #           [2, 4, 2, 0, 0, 1]]
        # },
        # 'task 2': {
        #     'c': [-1, 5, -2, 4, 3, 1, 2, 8, 3],
        #     'b': [3, 9, 9, 5, 9],
        #     'a': [[1, -3, 2, 0, 1, -1, 4, -1, 0],
        #           [1, -1, 6, 1, 0, -2, 2, 2, 0],
        #           [2, 2, -1, 1, 0, -3, 8, -1, 1],
        #           [4, 1, 0, 0, 1, -1, 0, -1, 1],
        #           [1, 1, 1, 1, 1, 1, 1, 1, 1]]
        # },
        # 'task 3': {
        #     'c': [2, 1, -2, -1, 4, -5, 5, 5],
        #     'b': [40, 107, 61],
        #     'a': [[1, 0, 0, 12, 1, -3, 4, -1],
        #           [0, 1, 0, 11, 12, 3, 5, 3],
        #           [0, 0, 1, 2, 0, 22, -2, 1]]
        # },
        # 'task 4': {
        #     'c': [2, 1, -2, -1, 4, -5, 5, 5, 1, 2],
        #     'b': [153, 123, 112],
        #     'a': [[1, 2, 3, 12, 1, -3, 4, -1, 2, 3],
        #           [0, 2, 0, 11, 12, 3, 5, 3, 4, 5],
        #           [0, 0, 2, 1, 0, 22, -2, 1, 6, 7]]
        # },
        # 'task 5': {
        #     'c': [1, 2, 1, -1, 2, 3],
        #     'b': [7, 16, 6],
        #     'a': [[2, 1, -1, -3, 4, 7],
        #           [0, 1, 1, 1, 2, 4],
        #           [6, -3, -2, 1, 1, 1]]
        # },

        # 'task 6': {
        #     'c': [10, 2, 1, 7, 6, 3, 1],
        #     'b': [12, 27, 19],
        #     'a': [[0, 7, 1, -1, -4, 2, 4],
        #           [5, 1, 4, 3, -5, 2, 1],
        #           [2, 0, 3, 1, 0, 1, 5]]
        # },
        # 'task 7': {
        #     'c': [2, 9, 3, 5, 1, 2, 4],
        #     'b': [6, 3, 7, 7],
        #     'a': [[0, 7, -8, -1, 5, 2, 1],
        #           [3, 2, 1, -3, -1, 1, 0],
        #           [1, 5, 3, -1, -2, 1, 0],
        #           [1, 1, 1, 1, 1, 1, 1]]
        # },
        # 'task 8': {
        #     'c': [-1, -3, -7, 0, -4, 0, -1],
        #     'b': [4, 8, 24],
        #     'a': [[1, 0, -1, 3, -2, 0, 1],
        #           [0, 2, 1, -1, 0, 3, -1],
        #           [1, 2, 1, 4, 2, 1, 1]]
        # },

        # 'task 9': {
        #     'c': [-1, 5, -2, 4, 3, 1, 2, 8, 3],
        #     'b': [3, 9, 9, 5, 9],
        #     'a': [[1, -3, 2, 0, 1, -1, 4, -1, 0],
        #           [1, -1, 6, 1, 0, -2, 2, 2, 0],
        #           [2, 2, -1, 1, 0, -3, 2, -1, 1],
        #           [4, 1, 0, 0, 1, -1, 0, -1, 1],
        #           [1, 1, 1, 1, 1, 1, 1, 1, 1]]
        # }
    }

    for i, item in tasks.items():
        print('-'*20 + i)
        method_gomory(np.array(item['a']), item['b'], item['c'])


if __name__ == '__main__':
    main()
