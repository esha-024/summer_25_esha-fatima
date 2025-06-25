#functions
import sys

def seq_concat(a,b):
    seq1=a
    seq2=b
    concat= seq1 + seq2

    return print("concatenation is" ,concat)


if __name__=="__main__":

 if len(sys.argv) !=3:
    sys.exit('usage python file_name arg_1, arg_2, arg_3')

seq_1 = sys.argv[1]
seq_2 = sys.argv[2]
seq_concat(seq_1,seq_2)


    