# An implementation of the widespread Huffman Algorithm

def huffman(pmf, d):
    '''Return a D-ary Huffman code.'''
    if (d < 2):
        raise Exception('d must be greater than or equal to 2!')    

    # recursion base case
    if (len(pmf) <= d):
        ensemble = {}
        for i, key in enumerate(pmf):
         ensemble[key] = str(i)
 
        return ensemble

    # get 'd' lowest probablities and their symbols
    min_x_set = []
    px_set = []

    for i in range(d):
        min_x_set.append(min(pmf, key=pmf.get))
        px_set.append(pmf.pop(min_x_set[i]))

    pmf[''.join(min_x_set)] = sum(px_set)
    
    # recursive step
    h = huffman(pmf, d)

    # update new ensemble
    node = h.pop(''.join(min_x_set))   

    for i in range(d):
        h[min_x_set[i]] = str(node)+str(i)

    return h
