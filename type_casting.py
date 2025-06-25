#-----------Type Casting----------

#  Converting Length to String:

length = len("ATGCATGC")
length_str = str(length)
print("Sequence length is " + length_str + " bases.")

#  Converting Input to Float:

expression = float(input("Enter expression level: "))
print(f"Expression as float: {expression}")

#  Converting Count to String:
gc_count = "ATGGCC".count('G') + "ATGGCC".count('C')
print("GC Count: " + str(gc_count))
