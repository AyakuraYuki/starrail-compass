# -*- coding: utf-8 -*-

import copy


def gcd(a: int, b: int) -> int:
    """
    返回 a、b 两数的最大公约数
    """
    while b:
        a, b = b, a % b
    return a


def egcd(a: int, b: int) -> tuple[int, int, int]:
    """
    扩展欧几里得算法
    返回 a、b 两数的最大公约数 g 同时，找到 x、y，使他们满足贝祖等式 ax + by = gcd(a, b)

    :returns g, x, y
    """
    if a == 0:
        return b, 0, 1
    else:
        g, y, x = egcd(b % a, a)
        return g, x - (b // a) * y, y


def mod_inv(a: int, m: int) -> int:
    g, x, y = egcd(a, m)
    assert g == 1
    return x % m


class GaussMatrix:
    def __init__(self, m_matrix: list[list[int]], m_mod: int):
        self.matrix: list[list[int]] = copy.deepcopy(m_matrix)
        self.d: list[list[int]] = []
        self.r = len(m_matrix)
        self.c = len(m_matrix[0])
        self.N = len(m_matrix[0]) - 1
        self.mod = m_mod
        self.count = 1
        self.error_str = "unknown error"

    def swap_row(self, row_a, row_b):
        (self.d[row_a], self.d[row_b]) = (self.d[row_b], self.d[row_a])

    def swap_col(self, column_a, column_b):
        for j in range(self.r):
            (self.d[j][column_a], self.d[j][column_b]) = (self.d[j][column_b], self.d[j][column_a])

    def inv_result(self, r, n) -> list[int] | None:
        b = self.d[n][self.N]
        a = self.d[n][n]
        m = self.mod
        k = gcd(a, m)
        for j in range(n + 1, self.N):
            b = (b - (self.d[n][j] * r[j] % m)) % m

        if k == 1:
            return [mod_inv(a, m) * b % m]
        else:
            if k == gcd(k, b):
                a /= k
                b /= k
                m /= k
                x0 = mod_inv(a, m) * b % m
                x = []
                for i in range(k):
                    x.append(x0 + m * i)
                return x
        return None

    def find_min_gcd_row_col(self, i, j) -> list[int]:
        for k in range(i, self.r):
            for _l in range(j, self.c - 1):
                if gcd(self.d[k][_l], self.mod) == 1:
                    return [k, _l]

        def add_min_gcd(a, b, m) -> list[int]:
            _r = [m, 1]
            _g = gcd(a, b)
            if _g:
                _i = a / _g
                for _j in range(int(_i)):
                    _g = gcd((a + _j * b) % m, m)
                    if _g < _r[0]:
                        _r[0] = _g
                        _r[1] = _j
                    if _g == 1:
                        break
            return _r

        r = [self.mod, 1, i, i + 1, j]
        for k in range(i, self.r):
            for kk in range(k + 1, self.r):
                for _l in range(j, self.c - 1):
                    rr = add_min_gcd(self.d[k][_l], self.d[kk][_l], self.mod)
                    if rr[0] < r[0]:
                        r[0] = rr[0]
                        r[1] = rr[1]
                        r[2] = k
                        r[3] = kk
                        r[4] = _l
                        pass
                    if 1 == rr[0]:
                        break
        g = r[0]
        n = r[1]
        k = r[2]
        kk = r[3]
        _l = r[4]

        if n and g < self.mod:
            self.d[k] = list(map(lambda x, y: (x + n * y) % self.mod, self.d[k], self.d[kk]))
        return [k, _l]

    def mul_row(self, i, k, j):
        a = self.d[k][j]
        b = self.d[i][j]

        def get_mul(_a: int, _b: int, m: int) -> int | float | None:
            _k = gcd(_a, m)
            if _k == 1:
                return mod_inv(_a, m) * _b % m
            else:
                if _k == gcd(_k, _b):
                    return mod_inv(int(_a / _k), int(m / _k)) * (_b / _k) % (m / _k)
            return None

        if b:
            mul = get_mul(a, b, self.mod)
            if mul is None:
                print_matrix(self.d)
                assert (mul is not None)
            self.d[i] = list(map(lambda x, y: (y - x * mul) % self.mod, self.d[k], self.d[i]))

    def guess(self):
        self.d = copy.deepcopy(self.matrix)
        for i in range(self.r):
            for j in range(self.c):
                self.d[i][j] = self.matrix[i][j] % self.mod

        if self.r < self.N:
            self.d.extend([[0] * self.c] * (self.N - self.r))

        index = [x for x in range(self.N)]
        for i in range(self.N):
            tmp = self.find_min_gcd_row_col(i, i)
            if tmp:
                self.swap_row(i, tmp[0])
                (index[i], index[tmp[1]]) = (index[tmp[1]], index[i])
                self.swap_col(i, tmp[1])
            else:
                self.error_str = "no min"
                return None

            for k in range(i + 1, self.r):
                self.mul_row(k, i, i)

        if self.r > self.N:
            for i in range(self.N, self.r):
                for j in range(self.c):
                    if self.d[i][j]:
                        self.error_str = "r(A) != r(A~)"
                        return None

        for i in range(self.N):
            self.count *= gcd(self.d[i][i], self.mod)

        if self.count > 100:
            self.error_str = "solution too more:%d" % self.count
            return None

        result: list[list[int]] = [[0] * self.N]
        for col in range(self.N - 1, -1, -1):  # reverse iterate
            new_result = []
            for row in result:
                inv_ret = self.inv_result(row, col)
                if inv_ret:
                    for column in inv_ret:
                        copied_r = row[:]  # deep copy
                        copied_r[col] = column
                        new_result.append(copied_r)

                else:
                    self.error_str = "no inv:col=%d" % col
                    return None

            result = new_result

        for i in range(len(result)):
            def xchg(a, b):
                result[i][b] = int(a)

            list(map(xchg, result[i][:], index))

        return result


################################################################################
# test
################################################################################


def print_array(x: list[int]):
    prn = "\t["
    for j in x:
        if j:
            prn += "%3d, " % j
        else:
            prn += "  0, "

    print(prn[:-2] + "],")


def print_matrix(x: list[list[int]]):
    print("[")
    for i in x:
        print_array(i)
    print("]")


def solve_matrix(_mod: int, _matrix: list[list[int]]):
    g = GaussMatrix(_matrix, _mod)
    ret = g.guess()
    if not ret:
        print("error:")
        print_matrix(g.d)
        print("error_str:", g.error_str)
    else:
        print(ret[0])


if __name__ == "__main__":
    mod = 6
    matrix = [
        [-1, 0, -1, mod - 4],
        [-3, -3, 0, mod - 0],
        [0, 3, 3, mod - 0]
    ]
    solve_matrix(mod, matrix)
