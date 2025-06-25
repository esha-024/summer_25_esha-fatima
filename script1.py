import sys

i=sys.argv
print(i)
if len(sys.argv)!=3 :
    sys.exit() 
    dna=sys.argv[1]
    rna=sys.argv[2]
print('the dna seq',dna)
print('the rna seq',rna)


