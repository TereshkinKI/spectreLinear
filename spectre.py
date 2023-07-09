def xxx(a, b):
    c = [0] * n
    for i in range(n):
        if a[i] == b[i]:
            c[i] = '0'
        else:
            c[i] = '1'
    return ''.join(c)


def n_ones(a):
    cnt = 0
    for bit in a:
        if bit =='1':
            cnt += 1
    return cnt


file_in = open('in.txt')
A = []

k = 0
for x in file_in:
    A.append(x.rstrip('\n'))
    k += 1

n = len(A[0])
res = set()
zero = '0' * n

for outer in range(2**k):
    multiplier = outer % 2
    outer //= 2
    tmp = A[0] if multiplier else zero

    for inner in range(1, k):
        multiplier = outer % 2
        outer //= 2
        if multiplier:
            tmp = xxx(tmp, A[inner])
    res.add(tmp)

gist = {}

for num in res:
    tmp = n_ones(num)
    if tmp in gist:
        gist[tmp] += 1
    else:
        gist[tmp] = 1

file_out = open("out.txt", mode="w")
for i in range(n+1):
    if i in gist:
        print(i, '\t', gist[i], file=file_out)