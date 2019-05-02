import huffman
import filecmp
import math

x_pmf = {
    'A': 0,
    'C': 0,
    'T': 0,
    'G': 0
}

x_entropy = 0
k_set = []
code_lengths = {}

def resolve_huffman_input(ensemble, d):
    if not ((len(ensemble) - 1) / (d - 1) == int((len(ensemble) - 1) / (d - 1))):
        resolved = False
        i = 1
        while(not resolved):
            ensemble['k{}'.format(i)] = 0
            k_set.append('k{}'.format(i))
            if (len(ensemble) - 1) / (d - 1) == int((len(ensemble) - 1) / (d - 1)):
                resolved = True
                return ensemble
    else:
        return ensemble          

# Read in Bacillus Subtilis bacterium genome and establish x_pmf
with open('NC_000964.3.seq') as bsb_genome:
    dna_length = 0
    for line in bsb_genome:
       for c in line:
           dna_length += 1
           x_pmf[c] = x_pmf[c] + 1

    for key in x_pmf:
        x_pmf[key] = x_pmf[key] / dna_length

    print('X pmf = {}\n'.format(x_pmf))   
    code_lengths['dna_length'] = dna_length

# Encode for D = 2, 3, 4
for i in range(2, 5):
    print('Encoding for d = {}\n'.format(i))
    huffman_code = huffman.huffman(resolve_huffman_input(x_pmf.copy(), i), i)
    print(huffman_code)
    with open('Encoded-D{}-NC_000964.3.seq'.format(i), 'w') as encoding:
        code_length = 0
        with open('NC_000964.3.seq') as bsb_genome:
            for line in bsb_genome:
                  for c in line:
                      code = str(huffman_code[c])
                      code_length += len(code)
                      encoding.write(code)

    code_lengths['d{}'.format(i)] = code_length     
    print('done\n')

# decode huffman codes
for i in range(2, 5):
    print('Decoding for d =', str(i))
    huffman_code = huffman.huffman(resolve_huffman_input(x_pmf.copy(), i), i)
    huffman_decoder = dict(zip(huffman_code.values(), huffman_code.keys()))
    with open('Decoded-D{}-NC_000964.3.seq'.format(i), 'w') as decoding:
        with open('Encoded-D{}-NC_000964.3.seq'.format(i)) as encoding:
            buffer = []
            for line in encoding:
                for c in line:
                    buffer.append(c)
                    if (''.join(buffer) in huffman_decoder):
                        decoding.write(huffman_decoder[''.join(buffer)])
                        buffer = []
    print('done\n')

# verify the correctness the decoded huffman codes
for i in range(2, 5):
    print('Checking correctness of decoded huffman code for d = {}'.format(i))
    if(filecmp.cmp('NC_000964.3.seq', 'Decoded-D{}-NC_000964.3.seq'.format(i))):
        print('done\n')
    else:
        raise Exception('Incorrect decoding of huffman code for d ={}'.format(i))    

# Compute entropies, empirical average length    h and average codeword lengths
for i in range(2, 5):
    for key in x_pmf:
         x_entropy += x_pmf[key] * (math.log(x_pmf[key]) / math.log(i))

    x_entropy = -1 * x_entropy    
    print('H(X) = {} for d = {}'.format(x_entropy, i))
    print('Calculating empirical average length of the huffman encoding for d = {}'.format(i))
    huffman_code = huffman.huffman(resolve_huffman_input(x_pmf.copy(), i), i)

    for k in k_set:
        if k in huffman_code.keys():
             huffman_code.pop(k)

    empirical_average_length = code_lengths['d{}'.format(i)] / code_lengths['dna_length']
    print('Empirical average length = {}\ndone'.format(empirical_average_length))
    print('Determining average codeword length for the huffman encoding for d = {}'.format(i))
    print('Average codeword length = {}'.format(len(''.join(huffman_code.values())) / len(huffman_code)))
    print('done\n')
    x_entropy = 0
