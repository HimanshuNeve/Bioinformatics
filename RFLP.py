#Restriction enzyme digestion of fASTA seq

'''def rflp (fasta_seq , restriction_site):
    seq_list = []
    for seq in fasta_seq.strip().split('\n'):
        if not seq.startswith('>'):
            seq_list.append(seq.strip())
    
    contineous_seq = "".join(seq_list)

    rflp_fragments = contineous_seq.split(restriction_site)

    return rflp_fragments

fasta_seq = input("Enter a FASTA sequence: ")
restriction_site = input("Enter the restriction site: ")

r_fragments = rflp(fasta_seq, restriction_site)

print(f"Original FASTA seq is: {fasta_seq}")
if len(r_fragments) > 1:
    print(f"the restriction site {restriction_site} was found and the sequence is digested")
    print(f"The restriction fragments: {r_fragments}")
else:
    print(f"the restriction site {restriction_site} was not found")
    print(f"the original segment is {r_fragments[0]}")'''



'''above code is with error and below code is debugged using gpt.
 The bug was showing error in reading the larger fasta sequence while imputng the sequence from the user,
thus we used open() to read the fasta file.'''

def rflp (fasta_file_path , restriction_site):
    seq_list = []
    with open (fasta_file_path ,"r" ) as f:
        for seq in f :
            seq = seq.strip()
            if not seq.startswith('>'):
                seq_list.append(seq)
    
    continuous_seq = "".join(seq_list)

    rflp_fragments = continuous_seq.split(restriction_site)

    return continuous_seq , rflp_fragments

fasta_file_path = input("Enter the FASTA file path: ")
restriction_site = input("Enter the restriction site: ")

continuous_seq, rflp_fragments = rflp(fasta_file_path, restriction_site)

print(f"Original FASTA seq is:\n {continuous_seq}")
if len(rflp_fragments) > 1:
    print(f"\nThe restriction site {restriction_site} was found and the sequence is digested")
    print(f"\nThe restriction fragments:\n {rflp_fragments}")
else:
    print(f"The restriction site {restriction_site} was not found")
    print(f"\nThe original segment is {rflp_fragments[0]}")
