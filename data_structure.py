#DATA STRUCTURE
# list

gene_list=["gene1","gene2","gene3"]
gene_list.append("gene4")  #insert 
gene_list.insert(2,"gene4") #insert at specific index
gene_list.remove("gene1")    #remove an element
remove=gene_list.pop()     #return the removed element
print(remove)
print(gene_list.index("gene3"))  #print the index of a specific element
print("gene list",gene_list)      #print whole gene list

#TUPLE

gene_1=(100,200),(500,600)
for start, end in gene_1:
    print("start position",start,"end position",end)

# dictionary
gene_dic={'brca':1,'p53':2,'gc35':3} #adding genes to dictionary
gene_dic['efgr']=5
print('gene dictionary' , gene_dic)
    

sequences = {" Seq1 ": " ATGC ", " Seq2 ": " GCTA "} #linking Sequences to IDs.
print (" Sequence Seq1 :", sequences [" Seq1 "])