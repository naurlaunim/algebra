from first_mod__update import *

def lnz(*args):
    a = matrix(*args)
    if a.len > a.dim:
        return False
    else:
        if a.rank() < a.len:
            return False
        else:
            return True
    ''' check is some vectors lnz '''

def mlnz(*args):
    assert not lnz(*args)
    a = matrix(*args) # make mx from vectors
    r = a.rank() # check it's rank
    b = [matrix(e) for e in args] # make one-vec mx from any vec
    c = [] # list for all lnz
    for k in range(r-1): # ?
        for i in range(a.len):
            for mx in b: # i don't know why, but it works
                x = [e for e in mx]
                vec = a.m[i]
                x.append(vec)

                if lnz(*x):
                    x = matrix(*x)
                    c.append(x)
        b = c
        c = []

    for e in b:
        # remove duplicates
        e = {tuple(el) for el in e}
        if e not in c:
            c.append(e)
    return c
    ''' find all mlnz systems '''


def lcv(x, a):
    assert type(a) in [list, tuple]
    assert lnz(*a), 'please, use some mlnz of this'
    assert not lnz(x, *a), 'can not use vector lnz with base'
    assert len(a) == len(x) == len(a[0]), 'dimension error'
    A = [e for e in a]
    A.append(x)
    m = matrix(*A)
    m = m.tran()
    m = m.gaustrt()
    m = m.gausrtr()
    return [e[-1] for e in m]
    ''' express vector x as a linear combination of basis vectors

    please, make list or tuple from basis vectors for correct work
    if you don't sure your vectors will set the base, use other functions as a preparation '''

def change_base(e1, e2):
    assert type(e1) == type(e2) == list and len(e1) == len(e2) == len(e1[0]) == len(e2[0])
    assert lnz(*e1) and lnz(*e2)
    mx = []
    for vec in e2:
        x = lcv(vec, e1)
        mx.append(x)
    m = matrix(*mx)
    return m.tran()

    '''find the transformation matrix from e1 to e2

    please, use list for both of them'''

def lm_change_base(e1, e2, f1, f2, A):
    assert all(type(e) == list for e in [e1, e2, f1, f2]) and type(A) == matrix
    assert all(lnz(*vec) for vec in [e1, e2, f1, f2])
    assert len(e1) == len(e1[0]) == len(e2) == len(e2[0]) == A.dim
    assert len(f1) == len(f1[0]) == len(f2) == len(f2[0]) == A.len
    dimV = A.dim
    dimW = A.len
    if e1 != e2:
        S = change_base(e1, e2)
    else: # can change only one base
        S = matrix.e_mx(dimV)
    if e1 == f1 and e2 == f2: # if f: V -> V
        Tinv = S.invmx()
    elif f1 != f2:
        T = change_base(f1, f2)
        Tinv = T.invmx()
    else:
        T = matrix.e_mx(dimW)
        Tinv = T
    B = Tinv*A*S
    return B
    '''find matrix of a linear map (f: V -> W) after changing the basis

    to input: lm type - matrix, basises - list of vectors
    [e] - base of vector space V, [f] - base W, A - linear map matrix in [e1] and [f1]'''

def bf_change_base(e1, e2, A):
    assert type(e1) == type(e2) == list and type(A) == matrix
    assert len(e1) == len(e1[0]) == len(e2) == len(e2[0]) == A.dim
    assert  lnz(*e1) and lnz(*e2)
    S = change_base(e1, e2)
    B = S.tran()*A*S
    return B
    '''find matrix of a bilinear function (f: VxV -> P) in new base

    [e1], [e2] - old and new bases, type - list; A - bilinear form in [e1], type - matrix'''





if __name__ == '__main__':
    # theme 2, ex. 29
    print(*mlnz([1,1,1,1,0], [1,1,-1,-1,-1], [2,2,0,0,-1], [1,1,5,5,2], [1,-1,-1,0,0]), sep = '\n')
    # 2, 30
    print(lnz([1,1,1], [1,1,2], [1,2,3]))
    print(lcv([6, 9, 14], [[1,1,1], [1,1,2], [1,2,3]]))
    # 2, 32
    a = [[1,1,1,1], [1,2,1,1], [1,1,2,1], [1,3,2,3]]
    b = [[1,0,3,3], [-2,-3,-5,-4], [2,2,5,4],[-2,-3,-4,-4]]
    print(lnz(*a) and lnz(*b))
    print(change_base(a, b))
    # 4, 38
    a = [[4,2,1], [5,3,2], [3,2,1]]
    b = [[1,4,0], [1,3,1], [1,2,3]]
    m = matrix([-41,-51,-30], [-95,-112,-67], [219,262,156])
    print(lm_change_base(e1 = a, f1 = a, f2 = b, e2 = b, A = m))
