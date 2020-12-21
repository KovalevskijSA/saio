inf = float('inf')


def floyd(d, n):
    r = [list(range(n)) for _ in range(n)]
    for j in range(n):
        for i in range(n):
            for k in range(n):
                if d[i][k] > d[i][j] + d[j][k]:
                    d[i][k] = d[i][j] + d[j][k]
                    r[i][k] = r[i][j]
    return d, r


def main():
    # d = [[0, 9, inf, 3, inf, inf, inf, inf],
    #      [9, 0,   2, inf, 7, inf, inf, inf],
    #      [inf, 2, 0,   2, 4,   8,   6, inf],
    #      [3, inf, 2,   0, inf, inf, 5, inf],
    #      [inf, 7, 4, inf,  0, 10, inf, inf],
    #      [inf, inf, 8, inf, 10, 0, 7, inf],
    #      [inf, inf, 6, 5, inf, 7, 0, inf],
    #      [inf, inf, inf, inf, 9, 12, 10, 0],
    #      ]
    #
    # dr, r = floyd(d, len(d))
    # print(dr)
    # print(r)

    d = [[0, 3, 2, 6, inf, inf, inf, inf, inf],
         [inf, 0, inf, 2, inf, inf, inf, inf, inf],
         [inf, inf, 0, inf, inf, 4, inf, inf, inf],
         [inf, inf, 3, 0, 1, inf, 6, inf, inf],
         [inf, inf, inf, inf, 0, inf, 7, 5, inf],
         [inf, inf, inf, inf, 5, 0, inf, 4, inf],
         [inf, inf, inf, inf, inf, inf, 0, 2, 4],
         [inf, inf, inf, inf, inf, inf, inf, 0, 4],
         [inf, inf, inf, inf, inf, inf, inf, inf, 0],
         ]

    dr, r = floyd(d, len(d))
    print(dr)
    print(r)

if __name__ == '__main__':
    main()