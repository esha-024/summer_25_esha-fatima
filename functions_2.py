import sys

def gc_content(a):

    dna=a

    length = len(dna)
    gc_count = dna.count('G') + dna.count('C')
    gc_content = int((gc_count / length)*100)
    
    return print("gc content", gc_content)
    


if __name__=="__main__":

 if len(sys.argv) !=2:
    sys.exit('usage python file_name arg_1, arg_2, arg_3')

seq = sys.argv[1]

gc_content(seq)