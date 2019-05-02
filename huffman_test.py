import huffman
import filecmp

x_pmf = {
    'A': 0,
    'C': 0,
    'T': 0,
    'G': 0
}

# Read in Bacillus Subtilis bacterium genome and establish x_pmf
with open('NC_000964.3.seq') as bsb_genome:
    dna_length = 0
    for line in bsb_genome:
       for c in line:
           dna_length += 1
           x_pmf[c] = x_pmf[c] + 1

    for key in x_pmf:
        x_pmf[key] = x_pmf[key] / dna_length

# Encode for D = 2, 3, 4
for i in range(2, 5):
    print('Encoding for D =', str(i))
    huffman_code = huffman.huffman(x_pmf.copy(), i)
    print(huffman_code)
    with open('Encoded-D' + str(i) + '-NC_000964.3.seq', 'w') as encoding:
        with open('NC_000964.3.seq') as bsb_genome:
            for line in bsb_genome:
                  for c in line:
                      encoding.write(str(huffman_code[c]))
    print('done')

# decode huffman codes
for i in range(2, 5):
    print('Decoding for D =', str(i))
    huffman_code = huffman.huffman(x_pmf.copy(), i)
    huffman_decoder = dict(zip(huffman_code.values(), huffman_code.keys()))
    with open('Decoded-D' + str(i) + '-NC_000964.3.seq', 'w') as decoding:
        with open('Encoded-D' + str(i) + '-NC_000964.3.seq') as encoding:
            buffer = []
            for line in encoding:
                for c in line:
                    buffer.append(c)
                    if (''.join(buffer) in huffman_decoder.keys()):
                        decoding.write(huffman_decoder[''.join(buffer)])
                        buffer = []
    print('done')

# verify the correctness the decoded huffman codes
for i in range(2, 5):
    print('Checking correctness of decoded huffman code for D = ', str(i))
    if(filecmp.cmp('NC_000964.3.seq', 'Decoded-D' + str(i) + '-NC_000964.3.seq')):
        print('done')
    else:
        raise Exception('Incorrect decoding of huffman code for D = ', str(i))    
