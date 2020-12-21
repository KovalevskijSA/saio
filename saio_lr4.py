from prettytable import PrettyTable
# # n - resurs, m - level
def get_plans(b, x, n, m, plan, plans):
    if m < 0:
        plan.reverse()
        plans.append(plan)
        return
    max_sum = max(b[m][:n])

    i_max_l = []
    for i in range(n):
        if b[m][i] == max_sum:
            i_max_l.append(i)
    for j in i_max_l:
        for i in x[m][j]:
            new_plan = plan.copy()
            new_plan.append(i)
            get_plans(b, x, n - i, m - 1, new_plan, plans)


def get_plan2(b, x, n, m):
    plan = []
    k = n
    for i in range(m - 1, -1, -1):
        max_sum = max(b[i][:k])
        j = b[i][:k].index(max_sum)
        plan.append(x[i][j][0])
        k -= x[i][j][0]
    return list(reversed(plan))


def dp(fx):
    m = len(fx)
    n = len(fx[0]) if m else 0
    b = [[0] * n for i in range(m)]
    x = [[[0] for j in range(n)] for i in range(m)]
    for i in range(n):
        b[0][i] = fx[0][i]
        x[0][i] = [i]

    for i in range(1, m):
        for j in range(n):
            bj = [fx[i][k] + b[i - 1][j - k] for k in range(j + 1)]
            max_sum = max(bj)
            k_list = [i for i in range(len(bj)) if bj[i] == max_sum]
            b[i][j] = max_sum
            x[i][j] = k_list

    table = PrettyTable([str(i) for i in range(n)])
    for i in range(m):
        table.add_row(['{0}{1},'.format(b[i][j], x[i][j]) for j in range(n)])
    print(table)
    plan = []
    plans = []
    get_plans(b, x, n, m - 1, plan, plans)
    print(plans)
    print(get_plan2(b, x, n, m))



def main():
    # fx = [[0, 3, 3, 3, 3, 3, 3],
    #       [0, 3, 3, 3, 3, 3, 3],
    #       [0, 3, 3, 3, 3, 3, 3]]
    # task 1
    # fx = [[0, 1, 2, 2, 4, 5, 6],
    #       [0, 2, 3, 5, 7, 7, 8],
    #       [0, 2, 4, 5, 6, 7, 7]]
    # # task 2
    # fx = [[0, 1, 2, 2, 4, 5, 6],
    #       [0, 2, 3, 5, 7, 7, 8],
    #       [0, 2, 4, 5, 6, 7, 7]]
    # # task 3
    # fx = [[0, 1, 1, 3, 6, 10, 11],
    #       [0, 2, 3, 5, 6, 7, 13],
    #       [0, 1, 4, 4, 7, 8, 7]]
    #
    # # task 4
    # fx = [[0, 1, 2, 4, 8, 9, 9, 23],
    #       [0, 2, 4, 6, 6, 8, 10, 11],
    #       [0, 3, 4, 7, 7, 8, 8, 24]]
    # # task 5
    # fx = [[0, 3, 3, 6, 7, 8, 9, 14],
    #       [0, 2, 4, 4, 5, 6, 8, 13],
    #       [0, 1, 1, 2, 3, 3, 10, 11]]
    # # task 6
    # fx = [[0, 2, 2, 3, 5, 8, 8, 10, 17],
    #       [0, 1, 2, 5, 8, 10, 11, 13, 15],
    #       [0, 4, 4, 5, 6, 7, 13, 14, 14],
    #       [0, 1, 3, 6, 9, 10, 11, 14, 16]]
    #
    # # task 6
    # fx = [[0, 1, 3, 4, 5, 8, 9, 9, 11, 12, 12, 14],
    #       [0, 1, 2, 3, 3, 3, 7, 12, 13, 14, 17, 19],
    #       [0, 4, 4, 7, 7, 8, 12, 14, 14, 16, 18, 22],
    #       [0, 5, 5, 5, 7, 9, 13, 13, 15, 15, 19, 24]]
    #
    # # task 7
    fx = [[0, 4, 4, 6, 5, 12, 12, 15, 16, 19, 19, 19],
          [0, 1, 1, 1, 4, 7, 8, 8, 13, 13, 19, 20],
          [0, 2, 5, 6, 7, 8, 9, 11, 11, 13, 13, 18],
          [0, 1, 2, 4, 5, 7, 8, 8, 9, 9, 15, 19],
          [0, 2, 5, 7, 8, 9, 10, 10, 11, 14, 17, 21]]

    # # task 8
    # fx = [[0, 1, 2, 2, 2, 3, 5, 8, 9, 13, 14],
    #       [0, 1, 3, 4, 5, 5, 7, 7, 10, 12, 12],
    #       [0, 2, 2, 3, 4, 6, 6, 8, 9, 11, 17],
    #       [0, 1, 1, 1, 2, 3, 9, 9, 11, 12, 15],
    #       [0, 2, 7, 7, 7, 9, 9, 10, 11, 12, 13],
    #       [0, 2, 5, 5, 5, 6, 6, 7, 12, 18, 22]]
    dp(fx)

if __name__ == '__main__':
    main()