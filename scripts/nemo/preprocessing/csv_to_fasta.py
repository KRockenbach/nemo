#!/usr/bin/python3

'''
BSD 3-Clause License

Copyright (c) 2023, F. Teufel and M. H. Gíslason

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

# for source refer to https://github.com/graph-part/graph-part


import pandas as pd
import argparse



def convert_to_fasta(csv_file: str, 
                    sequence_col:int, 
                    identifier_col: int, 
                    label_col: int = None, 
                    priority_col: int = None, 
                    label_name: str = 'label',
                    priority_name: str='priority',
                    delimiter='|'):

    df = pd.read_csv(csv_file)

    cols = df.columns

    seqs = df[cols[sequence_col]]

    # make fasta headers
    headers = '>' + df[cols[identifier_col]]
    if label_col is not None:
        headers = headers + delimiter + label_name + '=' + df[cols[label_col]].astype('string')
    if priority_col is not None:
        headers = headers + delimiter + priority_name + '=' + df[cols[priority_col]].astype('string')

    target_path = csv_file[:-3] + 'fasta'
    with open(target_path, 'w') as f:
            for head, seq in zip(headers, seqs):
                f.write(head+'\n')
                f.write(seq+'\n')

    print(f'Created {target_path}')



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', '-f', type=str, help='Path the the .csv file to be converted.')
    parser.add_argument('--sequence_col', '-s', type=int, help='Index of the column that contains the amino acid sequence.')
    parser.add_argument('--identifier_col', '-id', type=int, help='Index of the column that contains the sequence identifier')
    parser.add_argument('--label_col', '-la', type=int, default=None, help='Index of the column that contains the label of the sequence.')
    parser.add_argument('--priority_col', '-pr', type=int, default=None, help='Index of the column that contains the priority group of the sequence.')
    parser.add_argument('--delimiter', '-d', type=str, default='|', help='Delimiter to use in the generated fasta headers.')
    parser.add_argument('--label_name', '-ln', type=str, default='label', help='Name of the label in the generated headers.')
    parser.add_argument('--priority_name', '-pn', type=str, default='priority', help='Name of the priority in the generated headers.')
    args = parser.parse_args()

    convert_to_fasta(args.file, 
                     args.sequence_col, 
                     args.identifier_col, 
                     args.label_col, 
                     args.priority_col, 
                     args.label_name, 
                     args.priority_name, 
                     args.delimiter
                     )

if __name__ == '__main__':
    main()
