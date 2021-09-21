from multiprocessing import Pool
from pprint import PrettyPrinter

def compute(row, all):
    # row, all = arg
    line = all[row]
    ret = [line + ele for ele in all]
    return ret


if __name__ == "__main__":
    all = [i for i in range(10)]

    pool = Pool(6)

    inputs = [(i, all) for i in all]

    res = pool.starmap(compute, inputs)
    pp = PrettyPrinter()
    pp.pprint(res)

