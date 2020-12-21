from queue import PriorityQueue


def dijkstra_algorithm(g, s):
    n = len(g)
    d = [float('inf')] * n
    p = [-1] * n
    d[s] = 0
    u = [False] * n
    for i in range(n):
        v = None
        for j in range(n):
            if not u[j] and (v is None or d[j] < d[v]):
                v = j
        u[v] = True
        for k, j in g[v]:
            if d[v] + j < d[k]:
                d[k] = d[v] + j
                p[k] = v
    return d, p


def dijkstra_algorithm2(g, s):
    n = len(g)
    d = [float('inf')] * n
    p = [-1] * n
    d[s] = 0
    q = PriorityQueue()
    q.put((d[s], s))
    while not q.empty():
        dis, v = q.get()
        if dis == d[v]:
            for i, l in g[v]:
                if d[v] + l < d[i]:
                    d[i] = d[v] + l
                    p[i] = v
                    q.put((d[i], i))
    return d, p


def main():
    # g = [(index, wight), ()...]
    tasks = {
        'task 1': {
            'g': [[(1, 5), (7, 3)], # 0
                  [(2, 2), (6, 3)], # 1
                  [(4, 5)],         # 2
                  [(2, 2)],         # 3
                  [(3, 1), (9, 2)],
                  [(2, 4), (4, 1), (8, 2), (7, 6)], # 5
                  [(0, 2), (2, 2), (5, 5)],
                  [(1, 1), (6, 4), (8, 1)], # 7
                  [(9, 5)], # 8
                  [(3, 6), (5, 3)]]
        },
        'task 2': {
                'g': [[(1, 6), (2, 2), (6, 2)], # 0
                      [(2, 5), (5, 6)], # 1
                      [(5, 1)],         # 2
                      [(4, 3), (2, 2)],         # 3
                      [(7, 4)],
                      [(0, 4), (4, 6), (6, 3), (7, 7)], # 5
                      [(7, 4)],
                      [(1, 1), (8, 1)], # 7
                      [(3, 1), (4, 5), (5, 2)],
                      [(10, 1)]] # 8
            }
    }

    for i, j in tasks.items():
        print('-'*20 + i)
        d, p = dijkstra_algorithm2(j['g'], 0)
        print(f'd {d}\np {p}')


if __name__ == '__main__':
    main()