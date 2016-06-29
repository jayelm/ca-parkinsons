from collections import defaultdict
NCLUS = 6


def init():
    init = [range(1, NCLUS + 1) for x in xrange(NCLUS)]
    for index, val in enumerate(init):
        val.pop(index)
    return init


with open('./diffs.txt', 'r') as fin:
    lines = fin.readlines()

for line in lines:
    variable, rest = line.split(':')
    others = [map(int, pair.split('-')) for pair in rest.split()]
    nonsig = defaultdict(list)
    for a, b in others:
        nonsig[a].append(b)
        nonsig[b].append(a)
    lst = init()
    for a, bs in nonsig.iteritems():
        # Since 0 index
        for b in bs:
            lst[a - 1].remove(b)
    print '{}\t{}'.format(
        variable, '\t'.join([''.join(map(str, l)) for l in lst])
    )
