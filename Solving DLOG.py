import math
import numpy as np
import time
from sympy import n_order
from sympy.ntheory.modular import crt
import gmpy2
from gmpy2 import mpz,isqrt


def algorithm_shanks(g, p, h):
    n = 1 + isqrt(n_order(g,p))
    array = {}
    last_big_pow = 1
    array[1] =0
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

def calc_degree(g,array_x,p,e,m):
    x = 1
    for j in range(0,len(array_x)):
        x = x * (pow(g,int(-array_x[j] * int(pow(p,e+j))),m))
    return x

def prime_factorization_in_dagree(n):
    array_final = np.array([], int)
    unique, counts = np.unique(prime_factorization(n), return_counts=True)
    for i in range(0, len(counts)):
        if counts[i] != 1:
            array_final = np.append(array_final, pow(unique[i], counts[i]))
        else:
            array_final = np.append(array_final, unique[i])
    return array_final

def prime_factorization(n): 
    array = np.array([], int)
    while n > 1:
        i = 2
        f = 0
        while 1:
            if n % i == 0:
                # n = n // i
                n = gmpy2.t_div(n,i)
                array = np.append(array, gmpy2.mpz(i))
                f = 1
                break
            else:
                i += 1
        if f == 1:
            continue
    return array

def prime_factorization_and_dagree(n): 
    unique, counts = np.unique(prime_factorization(n), return_counts=True)
    return dict(zip(unique, counts))

def algorithm_polig_hellman(g, p, h):
    N = gmpy2.mpz(n_order(g,p))
    array_y = np.array([], int)
    dict_count_multipliers = prime_factorization_in_dagree(N)
    for i in range(0, len(dict_count_multipliers)):
        extent = gmpy2.t_div(N, dict_count_multipliers[i])
        array_y = np.append(array_y, algorithm_shanks(gmpy2.powmod(g, extent, p), p, gmpy2.powmod(h, extent, p)))
    return (crt(dict_count_multipliers,array_y)[0])

def calc_degree(g,array_x,p,e,m):
    x = 1
    for j in range(0,len(array_x)):
        x = x * (pow(g,int(-array_x[j] * int(pow(p,e+j))),m))
    return x

def reduction_dlp(g, m, h):
    ord = n_order(g, m)
    N = prime_factorization_and_dagree(ord)
    now_modal = prime_factorization_in_dagree(ord)
    x_f=[]
    for y in range(len(N)):
        p, e = list(N.items())[y]
        array_x = []
        h1 = pow(h, int(pow(p, e - 1)), m)
        g1 = pow(g, int(pow(p, e - 1)), m)
        array_x.append(algorithm_shanks(g1, m, h1) % now_modal[y])
        for i in range(1,e):
            h2 = ((pow(h, int(pow(p, e - i - 1)))) * calc_degree(g,array_x,p, e - i - 1, m)) % m
            array_x.append(algorithm_shanks(g1, m, h2) % now_modal[y])
        x_f_value = 0
        for i in range(0, e):
            x_f_value = (x_f_value + int(array_x[i]) * (int(pow(p, i))))% now_modal[y]
        x_f.append(x_f_value)
    return (crt(prime_factorization_in_dagree(ord),x_f)[0])


# g^x = h (mod p)

# Пример  (65 бит)
# p = mpz(1152903929600671559)
# g = mpz(1152851135879315581)
# h = mpz(575993103851904534)


# start_time = time.time()
# print(reduction_dlp(g, p, h))
# print("--- %s seconds ---" % (time.time() - start_time))

start_time = time.time()
print(algorithm_polig_hellman(g, p, h))
print("--- %s seconds ---" % (time.time() - start_time))


start_time = time.time()
print(algorithm_shanks(g, p, h))
print("--- %s seconds ---" % (time.time() - start_time))

