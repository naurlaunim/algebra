class matrix:
    def __init__(self, *args):
        assert all(type(e) in [list, tuple] for e in args)
        assert all(len(args[i]) == len(args[i-1]) for i in range(1, len(args)))
        for vec in args:
            assert all(type(e) in [float, int, complex] for e in vec)
        self.m = list(args)
        self.len = len(args)
        self.dim = len(args[0]) # matrix consists of vectors from arithmetic vector space dimension n
        self.changelist = []

    def __iter__(self):
        return (e for e in self.m)


    def e_mx(n):
        a = []
        for i in range(n):
            l = [0 for j in range(n)]
            l[i] = 1
            a.append(list(l))
        return matrix(*a)

    def __str__(self):
        # return str(*self.m, sep = '\n')
        s = '[' + str(self.m[0])
        for i in range(1, len(self.m)):
            s += '\n ' + str(self.m[i])
        return s + ']'

    def tran(self):
        return matrix(*map(list, zip(*self.m)))

    def __add__(self, other):
        assert self.len == other.len and self.dim == other.dim
        row = self.len
        col = self.dim
        c = []
        for i in range(row):

            d = [self.m[i][j] + other.m[i][j] for j in range(col)]
            c.append(d)
        return matrix(*c)

    def __mul__(self, other):
        assert self.dim == other.len
        m = self.len
        n = other.len
        k = other.dim
        C = []
        for row in range(m):
            p = []
            for col in range(k):
                q = 0
                for el in range(n):
                    q += self.m[row][el]*other.m[el][col]
                q = round(q, 3)
                p.append(q)
            C.append(p)
        return matrix(*C)

    def rowred(self):
        A = [list(e) for e in self.m]
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
                        l = [e/k for e in A[i]]
                        A.pop(i)
                        A.insert(i, l)
        # A = [tuple(e) for e in A]
        return matrix(*A)
        ''' reducing the matrix in rows

        it may be userful for the system of equations '''

    def fullred(self):
        A = [list(e) for e in self.m]
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
        # A = [tuple(e) for e in A]
        return matrix(*A)
        ''' reducing the matrix in rows and columns '''

    def is_square(self):
        return self.len == self.dim

    def el_eij(i, j, l, n):
        A = matrix.e_mx(n)
        x = list(A.m[i])
        x[j] = l
        A.m[i] = x
        return A
    def el_ei(i, l, n):
        A = matrix.e_mx(n)
        x = [e*l for e in A.m[i]]
        A.m[i] = tuple(x)
        return A
    def el_p(i, j, n):
        A = matrix.e_mx(n)
        A = [list(e) for e in A.m]
        A[i][j], A[i][i] = A[i][i], A[i][j]
        A[j][i], A[j][j] = A[j][j], A[j][i]
        A = [tuple(e) for e in A]
        return matrix(*A)

    def gaustrt(self):
        a = self
        m = self.len
        n = self.dim
        l = []

        for j in range(min(n, m)):
            A = [list(e) for e in a.m]
            if A[j][j] != 1:
                for i in range(j+1, m):
                    if A[i][j] == 1:
                        A[i], A[j] = A[j], A[i]
                        l.append((i, j))
                        break
            if A[j][j] == 0:
                for i in range(j+1, m):
                    if A[i][j] != 0:
                        A[i], A[j] = A[j], A[i]
                        l.append((i, j))
                        break
            if A[j][j] == 0:
                # A = [tuple(e) for e in A]
                a = matrix(*A)
                a = a.rowred()
                a.changelist = l
                return a
            for i in range(j+1, m):
                if A[i][j] != 0:
                    l = [(o/A[j][j])*A[i][j] for o in A[j]]
                    A[i] = [A[i][p] - l[p] for p in range(n)]
                    A[i] = [round(p, 3) for p in A[i]]
            # A = [tuple(e) for e in A]
            a = matrix(*A)
            a = a.rowred()
        a.changelist = l
        return a
        ''' straight move of Gauss's algorithm

        it does not touch the columns and stop if find empty one'''

    def gausrtr(self):
        a = self.gaustrt()
        assert a.m == self.m, 'sequence error'
        r = self.rank()
        a = self

        for h in range(r):
            a = a.rowred()
            A = [list(e) for e in a.m]

            j = r - h - 1 # the row which we subtract
            A[j] = [o/A[j][j] for o in A[j]] # turn the diagonal element in '1'
            for i in range(0, j): # a row from which subtract
                l = [o*A[i][j] for o in A[j]] # multiply each j-row item
                A[i] = [A[i][p] - l[p] for p in range(len(l))] # subtract from the i-row j-row multiplied by Aij
                A[i] = [round(p, 3) for p in A[i]]
            a = matrix(*A)
        return a
        ''' returning move of Gauss's algorithm

        (!) do not use with zero-columns '''

    def rank(self):
        a = self.fullred()
        a = a.gaustrt()

        A = [list(e) for e in a.m]
        for j in range(min(len(A), len(A[0]))):
            while A[j][j] == 0:
                for i in A:
                    i.append(i[j])
                    del i[j]
                    #translocation columns

                B = [tuple(e) for e in A]
                b = matrix(*B)
                b = b.gaustrt()
                B = [list(e) for e in b.m]

                while B[j][j] == 0:
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
                        return j+1
                    else:
                        A = B
        return min(len(A), len(A[0]))

    def invmx(self):
        assert self.is_square() and self.len == self.rank()
        A = [list(e) for e in self.m]
        l = self.len
        b = matrix.e_mx(l)
        B = [list(e) for e in b]
        for i in range(l):
            A[i] = A[i] + B[i]

        m = len(A) # code from gaustrt
        n = len(A[0])
        for j in range(m):
            # if A[j][j] != 1:
            #     for i in range(j+1, m):
            #         if A[i][j] == 1:
            #             A[i], A[j] = A[j], A[i]
            #             break
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
        # B = [tuple(e) for e in B]
        return matrix(*B)
        ''' get inverse matrix '''



