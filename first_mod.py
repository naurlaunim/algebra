def ismx(a):
    A = [k for k in a]
    l = len(A[0])
    for i in range(len(A)):
        if len(A[i]) != l:
            return False
        else:
            for j in A[i]:
                if j != float(j):
                    return False
    return True
    ''' checks whether the object is a matrix of real numbers '''

def e_mx(n):
    a = []
    for i in range(n):
        l = [0 for j in range(n)]
        l[i] = 1
        a.append(l)
    return a
    ''' creates an identity matrix of rank n '''

def tran(a):
    return list(map(list, zip(*a.copy())))
    ''' returns the transposed matrix '''

def add(a, b):
    A = [k for k in a]
    B = [l for l in b]
    if len(A) == len(B):
        row = len(A)
        col = len(A[0])
        c = []
        for i in range(row):
            if len(A[i]) == len(B[i]) == col:
                d = []
                for j in range(col):
                    p = A[i][j] + B[i][j]
                    d.append(p)
            else:
                return None
            c.append(d)
        return(c)
    else:
        return None
    ''' addition of matrices '''

def mul(a, b):
    A = [k for k in a]
    B = [l for l in b]
    if ismx(A) and ismx(B) and len(A[0]) == len(B):
        m = len(A)
        n = len(B)
        k = len(B[0])
    else:
        return None
    C = []
    for i in range(m):
        p = []
        for j in range(k):
            q = 0
            for f in range(n):
                q += A[i][f]*B[f][j]
            q = round(q, 3)
            p.append(q)
        C.append(p)
    return C
    '''  matrix multiplication '''

def fullred(a):
    if ismx(a):
        A = [h for h in a]
    else: return None
    C = []
    while C != A:
        C = A.copy()
        for k in [2, 3, 5, 7, 11, 13, 17, 19]:
            for i in range(len(A)):
                for j in range(len(A[i])):
                    if A[i][j]%k == 0:
                        t = 1
                    else:
                        t = 0
                        break
                if t == 1:
                    l = [f/k for f in A[i]]
                    A.pop(i)
                    A.insert(i, l)
        for k in [2, 3, 5, 7, 11, 13, 17, 19]:
            for j in range(len(A[0])):
                for i in range(len(A)):
                    if A[i][j]%k == 0:
                        t = 1
                    else:
                        t = 0
                        break
                if t == 1:
                    for i in range(len(A)):
                        A[i][j] = A[i][j]/k
    return A
    ''' reducing the matrix in rows and columns '''
def rowred(a):
    if ismx(a):
        A = [h for h in a]
    else: return None
    C = []
    while C != A:
        C = A.copy()
        for k in [2, 3, 5, 7, 11, 13, 17, 19]:
            for i in range(len(A)):
                for j in range(len(A[i])):
                    if A[i][j]%k == 0:
                        t = 1
                    else:
                        t = 0
                        break
                if t == 1:
                    l = [f/k for f in A[i]]
                    A.pop(i)
                    A.insert(i, l)
    return(A)
    ''' reducing the matrix in rows

    it may be userful for the system of equations '''

def gaustrt(a):
    if ismx(a):
        A = [u for u in a]
        m = len(A)
        n = len(A[0])
    else:
        return None
    for j in range(min(n, m)):
        if A[j][j] != 1:
            for i in range(j+1, m):
                if A[i][j] == 1:
                    A[i], A[j] = A[j], A[i]
                    break
        if A[j][j] == 0:
            for i in range(j+1, m):
                if A[i][j] != 0:
                    A[i], A[j] = A[j], A[i]
                    break
        if A[j][j] == 0:
            return A
        for i in range(j+1, m):
            if A[i][j] != 0:
                l = [(o/A[j][j])*A[i][j] for o in A[j]]
                A[i] = [A[i][p] - l[p] for p in range(n)]
                A[i] = [round(p, 3) for p in A[i]]
        A = rowred(A)
    return A
    ''' straight move of Gauss's algorithm

    it does not touch the columns and stop if find empty one'''
def gausrtr(a):
    if gaustrt(a) == a:
        A = [el for el in a]
    else:
        return 'invalid move'

    r = rank(A)

    for h in range(r):
        A = rowred(A)
        j = r - h - 1 # the row which we subtract
        A[j] = [o/A[j][j] for o in A[j]] # turn the diagonal element in '1'
        for i in range(0, j): # a row from which subtract
            l = [o*A[i][j] for o in A[j]] # multiply each j-row item
            A[i] = [A[i][p] - l[p] for p in range(len(l))] # subtract from the i-row j-row multiplied by Aij
            A[i] = [round(p, 3) for p in A[i]]
    return A
    ''' returning move of Gauss's algorithm

    (!) do not use with zero-columns '''

def rank(a):
    copy = [tuple(el) for el in a] # some another thing that doesn't work without tuples
    A = [list(el) for el in copy]
    A = fullred(A)
    A = gaustrt(A)
    for j in range(min(len(A), len(A[0]))):
        while A[j][j] == 0:
            for i in A:
                i.append(i[j])
                del i[j]
                #translocation columns
            B = gaustrt(A)
            while B[j][j] == 0:
                # print('1')
                # print(*A, sep = '\n')
                # print()
                # print(*B, sep = '\n') #
                #is the row empty?
                rank = j
                for i in A[j]:
                    if i != 0:
                        #no
                        rank = j+1
                        break
                if rank > j:
                    break
                else:
                    #is there another rows to check?
                    if j != len(A)-1:
                        #yes
                        del A[j]
                        del B[j]
                        break
                    else:
                        #no
                        return j
            if B[j][j] != 0:
                if B == A:
                    # print('3')
                    # print(*A, sep = '\n')
                    # print()
                    # print(*B, sep = '\n') #
                    return j+1
                else:
                    # print('2')
                    # print(*A, sep = '\n')
                    # print()
                    # print(*B, sep = '\n') #
                    A = B
    return min(len(A), len(A[0]))

def invmx(a):
    if ismx(a) and len(a) == len(a[0]) == rank(a):
        A = [el for el in a]
    else:
        return 'invalid move'

    B = e_mx(len(A))
    for i in range(len(A)):
        A[i] = A[i] + B[i]

    m = len(A) # code from gaustrt
    n = len(A[0])
    for j in range(min(n, m)):
        if A[j][j] != 1:
            for i in range(j+1, m):
                if A[i][j] == 1:
                    A[i], A[j] = A[j], A[i]
                    break
        for i in range(j+1, m):
            if A[i][j] != 0:
                l = [(o/A[j][j])*A[i][j] for o in A[j]]
                A[i] = [A[i][p] - l[p] for p in range(n)]

    for h in range(m): # code from gausrtr
        j = m - h - 1
        A[j] = [o/A[j][j] for o in A[j]]
        for i in range(0, j):
            l = [o*A[i][j] for o in A[j]]
            A[i] = [A[i][p] - l[p] for p in range(len(l))]
    for i in range(len(A)):
        B[i] = [A[i][j] for j in range(len(A[i])) if j >= len(A)]
        for j in range(m):
            B[i][j] = round(B[i][j] , 3)
    return B
    ''' get inverse matrix '''

def el_eij(i, j, l, n):
    A = e_mx(n)
    A[i][j] = l
    return A
def el_ei(i, l, n):
    A = e_mx(n)
    A[i][i] = l
    return A
def el_p(i, j, n):
    A = e_mx(n)
    A[i][j], A[i][i] = A[i][i], A[i][j]
    A[j][i], A[j][j] = A[j][j], A[j][i]
    return A

def sys(a, b):
    if ismx(a) and len(b) == len(a):
        A = [el for el in a]
        B = [el for el in b]
        l = len(b)
    else:
        return 'invalid move'
    C = A.copy()
    A = [tuple(el) for el in A] # without this it doesn't work, nobody knows why
    for i in range(l):
        C[i].append(b[i])
    A = [list(el) for el in A]
    # print(*A, sep = '\n')
    # print()
    # print(*C, sep = '\n')
    # print() # (1)

    if rank(C) == rank(A):
        type = 'sum'
        if rank(C) == rank(A) == len(A[0]):
            type = 'vyz'
    else:
        type = 'nes'
    # print(type)
    # print(*A, sep = '\n')
    # print()
    # print(*C, sep = '\n') # (1) special strings for check about the thing that doesn't work without tuples

    if type == 'nes':
        return type

    if type == 'vyz':
        C = gaustrt(C)
        C = gausrtr(C) # this matrix can have more rows than its rank, but can not have more columns, so we use it safely
        listx = [C[i][-1] for i in range(len(A[0]))]
        return type, listx

    if type == 'sum':
        r = rank(C)
        n = len(A[0])
        xchanges = []
        C = gaustrt(C)
        while any(C[o][o] == 0 for o in range(r)):
            for j in range(r):
                if C[j][j] == 0: # diagonal element
                    for e in range(j+1, n):
                        if C[j][e] != 0:
                        # find not empty row
                            for i in C: # rows
                                i[e], i[j] = i[j], i[e]
                                # translocation columns
                            xchanges.append((j, e))
                    if all(C[j][e] == 0 for e in range(len(C[j]))):
                    # find empty row
                        del C[j]
                    break # stop check diagonal to make gaustrt
            C = gaustrt(C)

        C = gausrtr(C)
        print(*C, sep = '\n')
        print()
        fsr = []
        for j in range(r, n): # columns of free variables
            print(j)
            vec = [(1*(-C[i][j]) + C[i][-1]) for i in range(r)]
            print(vec)
            if r != j:
                print('1', j, r)
                print(vec)
                vec += [0 for o in range(r,j)]
            vec += [1]
            print(vec)
            if j != n-1:
                print('2', j, n-1)
                print(vec)
                vec += [0 for o in range(j+1, n)]
            print(vec)
            fsr.append(vec)
        # we don't going to forget positions of variables
        if xchanges != []:
            xchanges = xchanges[::-1]
            fsr = list(zip(*fsr))
            for i in xchanges:
                fsr[i[0]], fsr[i[1]] = fsr[i[1]], fsr[i[0]]
            fsr = list(zip(*fsr))
        return type, fsr

