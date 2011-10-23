from collections import defaultdict

def AfterBreak(expr, sigs):
    breakpoints = []
    # Pad by 1 so end caes don't automatically pass
    for break_sample in range(1, len(expr) - 1):
        sign = (expr[0] > expr[break_sample])
        for i1 in range(break_sample):
            for i2 in range(break_sample, len(expr)):
                if sigs[i1][i2] == False or (expr[i1] > expr[i2] != sign):
                    break
            # "else continue break" construct breaks out of outer loop as well
            else:
                continue
            break
        else:
            # Monotonic split means we made it through all pairs of samples
            breakpoints.append(break_sample)
    return breakpoints

def FindPeak(expr, sigs):
    peak_points = []

    for peak_point in range(1, len(expr) - 1):
        if ((expr[peak_point - 1] < expr[peak_point] < expr[peak_point +1]) 
            or (expr [peak_point - 1] > expr[peak_point] > expr[peak_point + 1])):
            continue
            # This isn't the actual peak
                                                           

        sign = 1 if expr[peak_point] > expr[peak_point - 1] else -1

        prepost = [False, False]
        for other in range(len(expr)):
            if other == peak_point:
                continue
            if sign * expr[peak_point] > sign * expr[other] and sigs[peak_point][other]:
                prepost[peak_point < other] = True
                # If peak_point < other , then update prepost[1] (post),
                # otherwise update prepost[0] (pre)
        if prepost == [True, True]:
            peak_points.append(peak_point)
    return peak_points


def test_case_1():
    "A monotonic increasing series of significance"

    expr = [0, 1, 2, 3, 4, 5]
    sigs = defaultdict(lambda: defaultdict(lambda: True))

    print "Case 1", AfterBreak(expr, sigs), FindPeak(expr, sigs)

def test_case_2():
    "A single, obvious breakpoint"

    expr = [0, 0, 0, 1, 1, 1]
    sigs = {0:{}, 1:{}, 2:{}, 3:{}, 4:{}, 5:{}}
    sigs[0][0] = sigs[0][1] = sigs[0][2] = False
    sigs[0][3] = sigs[0][4] = sigs[0][5] = True
    sigs[1] = sigs[2] = sigs[0]
    sigs[3][0] = sigs[3][1] = sigs[3][2] = True
    sigs[3][3] = sigs[3][4] = sigs[3][5] = False
    sigs[4] = sigs[5] = sigs[3]

    print "Case 2:", AfterBreak(expr, sigs), FindPeak(expr, sigs)


def test_case_3():
    "A monotonic increasing series of no significance"

    expr = [0, 1, 2, 3, 4, 5]
    sigs = defaultdict(lambda: defaultdict(lambda: False))

    print "Case 3:", AfterBreak(expr, sigs), FindPeak(expr, sigs)


def test_case_4():
    "A single, obvious breakpoint, with noise"

    expr = [1, 3, 2, 15, 14, 17]
    sigs = {0:{}, 1:{}, 2:{}, 3:{}, 4:{}, 5:{}}
    sigs[0][0] = sigs[0][1] = sigs[0][2] = False
    sigs[0][3] = sigs[0][4] = sigs[0][5] = True
    sigs[1] = sigs[2] = sigs[0]
    sigs[3][0] = sigs[3][1] = sigs[3][2] = True
    sigs[3][3] = sigs[3][4] = sigs[3][5] = False
    sigs[4] = sigs[5] = sigs[3]

    print "Case 4:", AfterBreak(expr, sigs), FindPeak(expr, sigs)

def test_case_5():
    "A rising peak in the middle, difference with neighbors not significant"
    expr = [2,3,4,3,2]
    sigs = {}
    sigs[0] = {0:False, 1:False, 2:True, 3:False, 4:False}
    sigs[1] = {0:False, 1:False, 2:False, 3:False, 4:False}
    sigs[2] = {0:True, 1:False, 2:False, 3:False, 4:True,}
    sigs[3] = sigs[1]
    sigs[4] = sigs[0]

    print "Case 5:", AfterBreak(expr, sigs), FindPeak(expr, sigs)
    

if __name__ == "__main__":
    test_case_1()
    test_case_2()
    test_case_3()
    test_case_4()
    test_case_5()

