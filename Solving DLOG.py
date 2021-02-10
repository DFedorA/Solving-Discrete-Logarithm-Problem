# coding: utf8
import random
import numpy as np
import time
from sympy import n_order, factorint
from sympy.ntheory.modular import crt
import gmpy2
from gmpy2 import mpz, isqrt


def algorithm_shanks(g, p, h, o=None):
    if o is None:
        o = n_order(g, p)
    n = 1 + isqrt(o)
    array = {}
    last_big_pow = 1
    array[1] = 0
    for j in range(1, n):
        y = (last_big_pow * g) % p
        last_big_pow = y
        if y in array:
            continue
        array[y] = j
    if h in array:
        return 0 * n + array[h]
    big_pow = gmpy2.powmod(g, -n, p)
    last_big_pow = (h * big_pow) % p
    for j in range(1, n):
        y = last_big_pow
        last_big_pow = (last_big_pow * big_pow) % p
        if y in array:
            return j * n + array[y]
    raise Exception('Нет решения')


def prime_factorization_in_degree(n): 
    array_final = np.array([], int)
    dict_prime = factorint(n)
    for k, p in dict_prime.items():
        array_final = np.append(array_final, mpz(pow(k, p)))
    return array_final


def algorithm_polig_hellman(g, p, h):
    N = gmpy2.mpz(n_order(g, p))
    array_y = np.array([], int)
    dict_count_multipliers = prime_factorization_in_degree(N)
    for i in range(0, len(dict_count_multipliers)):
        extent = gmpy2.t_div(N, dict_count_multipliers[i])
        array_y = np.append(array_y, algorithm_shanks(gmpy2.powmod(g, extent, p), p, gmpy2.powmod(h, extent, p)))
    return crt(dict_count_multipliers, array_y)[0]


def calc_degree(g, array_x, p, e, m):
    x = 1
    for j in range(0, len(array_x)):
        x = x * (pow(g, int(-array_x[j] * int(pow(p, e + j))), m))
    return x


def reduction_dlp(g, m, h):
    ord = n_order(g, m)
    N = factorint(ord)
    x_f = []
    for p, e in N.items():
        extent = gmpy2.t_div(ord, pow(p, e))
        gi = gmpy2.powmod(g, extent, m)
        hi = gmpy2.powmod(h, extent, m)
        array_x = []
        h1 = pow(hi, int(pow(p, e - 1)), m)
        g1 = pow(gi, int(pow(p, e - 1)), m)
        array_x.append(algorithm_shanks(g1, m, h1, p) % pow(p, e))
        for i in range(1, e):
            h2 = ((pow(hi, int(pow(p, e - i - 1)))) * calc_degree(gi, array_x, p, e - i - 1, m)) % m
            array_x.append(algorithm_shanks(g1, m, h2, p) % pow(p, e))
        x_f_value = 0
        for i in range(0, e):
            x_f_value = (x_f_value + int(array_x[i]) * (int(pow(p, i)))) % pow(p, e)
        x_f.append(x_f_value)
    return crt(prime_factorization_in_degree(ord), x_f)[0]


def generate_prime(bits):
    rand_state = gmpy2.random_state(random.randint(20, 300))
    temp = gmpy2.mpz_rrandomb(rand_state, bits)
    return gmpy2.next_prime(temp)


def check_for_compare_numbers(a, m):
    return gmpy2.gcd(a, m) == 1


if __name__ == '__main__':
    check = False
    p = g = 0
    while check is False:
        p = generate_prime(40)
        g = generate_prime(38)
        check = check_for_compare_numbers(p, g)
    x = generate_prime(32)
    h = gmpy2.powmod(g, x, p)
    print(f"p: {p} g:{g} h:{h}")

    start_time = time.time()
    print(x,reduction_dlp(g, p, h))
    print("Результат теста алгоритма подгрупп: ", reduction_dlp(g, p, h) == x)
    print("--- Время выполнения: %s  ---" % (time.time() - start_time))

    start_time = time.time()
    print("Результат теста алгоритма Полига-Хеллмана: ", algorithm_polig_hellman(g, p, h) == x)
    print("--- Время выполнения: %s  ---" % (time.time() - start_time))

    start_time = time.time()
    print("Результат теста алгоритма Шенкса: ", algorithm_shanks(g, p, h) == x)
    print("--- Время выполнения: %s  ---" % (time.time() - start_time))




# g^x = h (mod p)


# Пример 8 (65 бит)
# p = mpz(1152903929600671559)
# g = mpz(1152851135879315581)
# h = mpz(575993103851904534)
# 188497496704843170
# --- 0.05800485610961914 seconds --- algorithm_polig_hellman
