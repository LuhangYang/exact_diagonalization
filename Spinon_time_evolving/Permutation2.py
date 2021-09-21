from pprint import PrettyPrinter
#import MyTimer

#timer = MyTimer.MyTimer()
pp = PrettyPrinter()


def permutation(length: int, one: int, zero: int):
    # NOTE: dp table
    #  { len1:
    #       {
    #       "one1": {"zero1":{ "1":[],"0":[]}}
    #       "one2": {"zero1":{ "1":[],"0":[]}}
    #       }
    #    len2:
    #       {
    #       "one1": {"zero1":{ "1":[],"0":[]}}
    #       "one2": {"zero1":{ "1":[],"0":[]}}
    #       }
    #  }

    table = dict()
    for i in range(length):
        table[i + 1] = dict()
        for j in range(i + 1):
            table[i + 1][j] = dict()
            for k in range(i + 1):
                # NOTE start with 1 and start with 0 in string of length i+1
                #   the string has j ones and k zeros
                table[i + 1][j][k] = {"1": [], "0": []}

    # NOTE len=1 one = 0, zero = 0
    table[1][0][0] = {"1": ["1"], "0": ["0"]}
    # # NOTE len=1 one = 1, zero = 0
    # table[1][1][0] = {"1": [], "0": []}
    # # NOTE len=1 one = 0, zero = 1
    # table[1][0][1] = {"1": [], "0": []}
    # # NOTE len=1 one = 1 zero = 1
    # table[1][1][1] = {"1": [], "0": []}

    # NOTE len=2 one = 0, zero = 0
    table[2][0][0] = {"1": ["10"], "0": ["01"]}
    # NOTE len=2 one = 1, zero = 0
    table[2][1][0] = {"1": ["11"], "0": []}
    # NOTE len=2 one = 0, zero = 1
    table[2][0][1] = {"1": [], "0": ["00"]}
    # NOTE len=2 one = 1 zero = 1
    table[2][1][1] = {"1": [], "0": []}

    for l in range(2, length + 1):
        one_range = l if l < one else one
        for ones in range(one_range + 1):
            upper = l - ones
            zero_range = zero if zero < upper else upper
            for zeros in range(zero_range + 1):
                # NOTE generate current from previous smaller string
                # print(l, ones, zeros)
                if ones == 0 and zeros == 0:
                    table[l][ones][zeros]["1"] = add_one(
                        table[l - 1][ones][zeros]["0"])
                    table[l][ones][zeros]["0"] = add_zero(
                        table[l - 1][ones][zeros]["1"])
                # table[l][ones][zeros]["1"] = table[l - 1][ones - 1]["zeros"]["1"]
                elif ones == 0:

                    if l - 1 < zeros:
                        continue

                    if l - 1 == zeros:
                        table[l][ones][zeros]["0"] = add_zero(
                            table[l - 1][ones][zeros - 1]["0"])
                        continue

                    # NOTE zeros are not 0
                    table[l][ones][zeros]["0"] = add_zero(
                        table[l - 1][ones][zeros]["1"]) + add_zero(
                        table[l - 1][ones][zeros - 1]["0"])

                    table[l][ones][zeros]["1"] = add_one(
                        table[l - 1][ones][zeros]["0"])

                    pass
                elif zeros == 0:
                    # NOTE ones are not 0
                    # TODO
                    if l - 1 < ones:
                        continue
                    if l - 1 == ones:
                        table[l][ones][zeros]["1"] = add_one(
                            table[l - 1][ones - 1][zeros]["1"])
                        continue
                    table[l][ones][zeros]["1"] = add_one(
                        table[l - 1][ones][zeros]["0"]) + add_one(
                        table[l - 1][ones - 1][zeros]["1"])

                    table[l][ones][zeros]["0"] = add_zero(
                        table[l - 1][ones][zeros]["1"])
                    pass
                else:
                    if ones + zeros > l - 1:
                        # NOTE impossible
                        # print(l, ones, zeros, "impossible")
                        continue
                    if ones + zeros == l - 1:
                        table[l][ones][zeros]["0"] = add_zero(
                            table[l - 1][ones][zeros - 1]["0"])
                        table[l][ones][zeros]["1"] = add_one(
                            table[l - 1][ones - 1][zeros]["1"])

                        continue
                    # NOTE zeros are not 0
                    table[l][ones][zeros]["0"] = add_zero(
                        table[l - 1][ones][zeros]["1"]) + add_zero(
                        table[l - 1][ones][zeros - 1]["0"])

                    # NOTE ones are not 0
                    table[l][ones][zeros]["1"] = add_one(
                        table[l - 1][ones][zeros]["0"]) + add_one(
                        table[l - 1][ones - 1][zeros]["1"])
                    pass

    res = table[length][one][zero]["1"] + table[length][one][zero]["0"]

    for bin_str in res:
        if check(bin_str, one, zero) == False:
            raise ValueError("Error: " + bin_str)

    return res


def add_one(lst: list):
    res = ["1" + item for item in lst]
    return res


def add_zero(lst: list):
    res = ["0" + item for item in lst]
    return res



def val_to_bin_str(value: int, length: int):
    f = "{:0" + str(length) + "b}"
    return f.format(value)


def check(bin_str: str, one: int, zero: int) -> bool:
    ones = 0
    zeros = 0
    for i in range(len(bin_str) - 1):
        if bin_str[i] == "1" and bin_str[i + 1] == "1":
            ones += 1
            if ones > one:
                return False
        elif bin_str[i] == "0" and bin_str[i + 1] == "0":
            zeros += 1
            if zeros > zero:
                return False
    return ones == one and zero == zeros



#timer.start()

if __name__ == "__main__":

    res = permutation(8, 2, 0)
    print(res)
    print(len(res))

    res = permutation(8, 4, 2)
    print(res)
    print(len(res))


