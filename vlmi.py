import multiprocessing
import timeit


def simulate(a,b,c):
    result = a + b + c
    return result


def run_parallel(data):
    pool = multiprocessing.Pool()
    results = pool.starmap(simulate, data)
    pool.close()
    pool.join()
    print(len(results))


def run_sync(data):
    results = list(map(lambda l: simulate(*l), data))
    print(len(results))


def main(is_parallel, data):
    if is_parallel:
        runtime = timeit.timeit(
            lambda:run_parallel(data),
            number=1
        )
    else:
        runtime = timeit.timeit(
            lambda:run_sync(data),
            number=1
        )
    print(runtime)

def say_hello(): print("hello")
