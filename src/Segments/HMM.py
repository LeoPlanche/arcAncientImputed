import numpy as np

states = [0,1]

def initS(p) -> np.array:
    S = np.zeros(2)
    S[0]=1-p
    S[1]=p
    return S
    
#T = Time admixture
def initA(T,r,L,a) -> np.array:
    A = np.zeros((2,2))
    A[0][1]=T*r*L*a
    A[1][0]=T*r*L*(1-a)
    A[0][0]=1-A[0][1]
    A[1][1]=1-A[1][0]
    return A



def initB(m,L,Ti,Ta) -> np.array:
    B = np.zeros((2,75))
    meani = m*L*Ti
    meana = m*L*Ta
    B[0][0] = np.exp(-meani)
    B[1][0] = np.exp(-meana)
    sumi=B[0][0]
    suma=B[1][0]
    for i in range(1,len(B[0])):
        B[0][i]=B[0][i-1]*meani/i
        B[1][i]=B[1][i-1]*meana/i
        sumi=sumi+B[0][i]
        suma=suma+B[1][i]
    B[0][0]=B[0][0]+1-sumi
    B[1][0]=B[1][0]+1-suma
    return B



def viterbi(V, initial_distribution, a, b):
    T = len(V)
    M = a.shape[0]
 
    omega = np.zeros((T, M))
    omega[0, :] = np.log(initial_distribution * b[:, V[0]])
 
    prev = np.zeros((T - 1, M))
 
    for t in range(1, T):
        for j in range(M):
            # Same as Forward Probability
            probability = omega[t - 1] + np.log(a[:, j]) + np.log(b[j, V[t]])
            prev[t - 1, j] = np.argmax(probability)
            # This is the probability of the most probable state (2)
            omega[t, j] = np.max(probability)
 
    # Path Array
    S = np.zeros(T)
    # Find the most probable last hidden state
    last_state = np.argmax(omega[T - 1, :])
 
    S[0] = last_state
 
    backtrack_index = 1
    for i in range(T - 2, -1, -1):
        S[backtrack_index] = prev[i, int(last_state)]
        last_state = prev[i, int(last_state)]
        backtrack_index += 1
 
    # Flip the path array since we were backtracking
    S = np.flip(S, axis=0)
 
    # Convert numeric values to actual hidden states
 
    result = []
    for s in S:
        if s == 0:
            result.append(0)
        else:
            result.append(1)
 
    return result


def create_obs(seq1,seq2,cut):
    k=0
    maxpos = min(max(seq1[0][1:]),max(seq2[0][1:]))
    seq = np.zeros((int(int(maxpos)/cut)+1),dtype=int) #List of seq, a seq is a list with nb of mutations
    for i in range(0,len(seq1[0])):
        j=int(int(seq1[0][i])/1000)
        while k<len(seq2[0])-1 and seq2[0][k] < seq1[0][i]:
            k+=1
        if seq2[0][k] == seq1[0][i]:
            if not (seq2[1][k] in seq1[1][i]) or  not (seq2[2][k] in seq1[1][i]):
                seq[j]=seq[j]+1    
    return seq

def get_HMM_tracts(seq):
    migrating_tracts = []
    maxi = 0
    for e in seq:
        if e>maxi:
            maxi=e
            
    for i in range(maxi+1):
        migrating_tracts.append([])
    start=0
    for i in range(1,len(seq)):
        if seq[i]!=seq[i-1]:
            migrating_tracts[seq[i-1]].append([start,i-1])
            start=i
    migrating_tracts[seq[len(seq)-1]].append([start,len(seq)-1])
    return migrating_tracts

def write_tracts(file, tracts, chr, filename, similarity_neanderthal,similarity_altai,similarity_cha, similarity_denisova\
                ,snps_archaic, snps_neanderthal, snps_altai, snps_cha, snps_denisova):
    for tract, sim_nean, sim_altai, sim_cha, sim_deni, snp_arc, snp_nean, snp_altai, snp_cha, snp_deni\
          in zip(tracts,similarity_neanderthal,similarity_altai,similarity_cha,similarity_denisova,\
                 snps_archaic, snps_neanderthal, snps_altai, snps_cha, snps_denisova):
        file.write(f"chr{chr}\t{tract[0]*1000}\t{tract[1]*1000}\t{sim_nean}\t{sim_altai}\t{sim_cha}\t{sim_deni}\t{snp_arc}\t{snp_nean}\t{snp_altai}\t{snp_cha}\t{snp_deni}\t{filename}\n")
    file.flush()
    file.close() 