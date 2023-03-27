import numpy as np



def Pauli_mult(a,b):
    """a and b corresponds to Paulis. For now we don't care about the phase. This phase can be corrected 
    after the state preparation."""
    if a == 'X' and b == 'X':
        result = 'I'
        phase = 1
    if a == 'X' and b == 'Y':
        result = 'Z'
        phase = 1j
    if a == 'X' and b == 'Z':
        result = 'Y'
        phase = -1j
    if a == 'Y' and b == 'X':
        result = 'Z'
        phase = 1j
    if a == 'Y' and b == 'Y':
        result = 'I'
        phase = 1
    if a == 'Y' and b == 'Z':
        result = 'X'
        phase = 1j
    if a == 'Z' and b == 'X':
        result = 'Y'
        phase = 1j
    if a =='Z' and b =='Y':
        result = 'X'
        phase = -1j
    if a =='Z' and b =='Z':
        result = 'I'
        phase = 0
    if a =='I':
        result = b
        phase = 0
    if b =='I':
        result =a
        phase= 0
    return result




def stab_multi(stab_a,stab_b):
    """The two input corresponds to list of Paulis. The output gives the result after multiplication.
    Input form example: ['X','X'],['Y','Y'].
    Out form example :['Z','Z']
    """
    m = len(stab_a)
    stab_after = [Pauli_mult(stab_a[i],stab_b[i]) for i in range(0,m)]
    return stab_after


def echelon_gauge_column(stab,no,stab_new_final):
    '''Inputs:
    stab: the partial old stabilizer generator
    stab_new_final: the partial new stabilizer generator. 
    no: corresponds to the column number under consideration. 
    This is a subroutine used in the final echelon gauge function. 
    
    Algorithmic details: 
    1. we specify the column number under consideration. We then look for the Paulis in the column. 
    2. We rearrange the rows such that the rows holding the non-trivial Paulis (upto) are the top two rows. 
    3. We then make all the Paulis below the rows holding non trivial Paulis as I by row multiplication. 
    4. The rows which contain the non-trivial Paulis go to the stab_new_final and the left over ones go to stab_new1
    
    '''
    
    list_col_paulis = [stab[i][no] for i in range(0,len(stab))]# get the non-trivial Paulis
    uniq_elements = [i for i in list(set(list_col_paulis)) if i != 'I']#get the list of non-trivial Paulis in a column
    stab_0_element = [list_col_paulis.index(uniq_elements[i]) for i in range(0,len(uniq_elements))] # get the index of the first occurrence of Paulis

    if uniq_elements!=[]:
    #Below we have moved around the rows
        stab_new = []
        [stab_new.append(stab[j]) for j in stab_0_element] # add the rows with non-trivial Pauli to a new list
        for i in range(0,len(stab)):
            if i not in stab_0_element:
                stab_new.append(stab[i]) # add the rest of the rows

        # We now start with the multiplication of elements to make the columns identity.
        for i in range(len(uniq_elements),len(stab_new)):
#             print(i)
            if stab_new[i][no] in uniq_elements:
                stab_new11 = stab_new[uniq_elements.index(stab_new[i][no])] 
                stab_new[i] = stab_multi(stab_new[i],stab_new11)
            elif stab_new[i][no]== 'I':
                stab_new


            else:
                
                for k in range(0,len(uniq_elements)):
                    stab_new[i] = stab_multi(stab_new[k],stab_new[i])
                    
        stab_new_final = stab_new_final + stab_new[0:len(uniq_elements)]
        stab_new1 = stab_new.copy()
        [stab_new1.remove(i) for i in stab_new_final if i in stab_new1]
    if uniq_elements ==[]:

        if stab!=[]:
            stab_new_final.append(stab[0])
            stab_new1= stab.copy()
            stab_new1.pop(0)
        if stab == []:
            stab_new1 = []
    return stab_new_final,stab_new1


def echelon_gauge(stabilizer_generator):
    """The stabilizer_generator is the input and outputs the echelon transformation of the stabilizer_generator. """
    no_columns = len(stabilizer_generator[0])
    stab = stabilizer_generator
    stab_new_final = []
    for no in range(0,no_columns):
        stab_new_final, stab = echelon_gauge_column(stab,no,stab_new_final)
    return  stab_new_final



def entropy_inter(echelon_stab):
    """The input is the stabilizer generator in the echelon form of the generator.
    The output is list of occurences of non-trivial Paulis in each row"""
    index_pauli = []
    echelon_stab = np.array(echelon_stab)
    num_rows, num_cols = echelon_stab.shape
    
    index_pauli = [np.argmax(echelon_stab[j1,:]!= 'I') for j1 in range(num_rows)]
    
    return index_pauli
    
    
    
    
def entropy(echelon_stab,photon_number):
    """
    Inputs: 
    echelon_stab: the stabilizer generator in the echelon gauge
    photon_number: the total number of photons in the graph state
    Output: 
    This returns the list of entanglement entropy across the partition 012...x|x+1...photon_number, where x is the index number. """
    entropy1 = []
    for x in range(0,photon_number):
        
        index_pauli = entropy_inter(echelon_stab)
        indices = list(filter(lambda y: index_pauli[y] > x, range(len(index_pauli))))

        inter2 = len(indices)
    
        entropy1.append(photon_number-x-inter2-1) # The minus one is there because the numbering the x starts from 0
    return entropy1


def steps_algo_1(stabilizer_generator,number_emit,photon_number):
    """This returns the time_count for the first algorithm. The inputs are
    stabilizer_generator: stabilizer generators (not necessarily in the echelon gauge)
    number_emit: The total number of emitters
    photon_number: The total number of photons in the graph state. 
    The output gives us the time steps. 
    
    This algorithm depends on the the entropy function and doesn't make use of the internal structure of the stabilizer. 
    """
    echelon_stab = echelon_gauge(stabilizer_generator)
    entropy_all = entropy(echelon_stab,photon_number)
    
    count = 1
    number_left = photon_number-number_emit
    while number_left > 0:
        count = count +1 
        emit_available = number_emit - entropy_all[number_left-1]
        if entropy_all[number_left-1]>= entropy_all[number_left-2]:
            absorb_photon = emit_available + 1
            number_left = number_left - absorb_photon
        else: 
            absorb_photon = emit_available
            number_left = number_left - absorb_photon
    return count
        
            
        
def stab_with_emitter(stabilizer_generator,emitter_number):
    """Takes in as input the stabilizer generator and the number of emitters."""
    """As an output gives the stabilizer generator of the photons + emitters"""
    # Adding the extra columns
    photon_number = len(stabilizer_generator[0])
    emi_stab = ['I' for i in range(emitter_number)]
    emi_stab = np.array(emi_stab)
    emi_stab1 = np.tile(emi_stab, (photon_number, 1))
    stabilizer_generator1 = np.append(stabilizer_generator, emi_stab1, axis=1)
    
    
    # Adding extra rows
    stab_new1 = ['I' for i in range(len(stabilizer_generator1[0]))]
    stab_new1 = np.array([stab_new1])
    stab_new1 = np.tile(stab_new1 , (emitter_number, 1))
    stabilizer_generator1 = np.append(stabilizer_generator1, stab_new1, axis=0)
    
    list1 = [k+photon_number for k in range(0,emitter_number)]
    for l in list1:
        stabilizer_generator1[l,l] = 'Z'
        
    stabilizer_generator1 = stabilizer_generator1.tolist()
        

    return stabilizer_generator1



def apply_hadamard(m,stabilizer_generator):
    '''Applies Hadamard to the m qubit. The numbering starts from 0.'''
    stabilizer_generator1 = stabilizer_generator.copy()
    for stab in stabilizer_generator1:
        #print(stab[m])
        if stab[m] == 'X':
            stab[m] ='Z'
        elif stab[m] == 'Z':
            stab[m] = 'X'
    return stabilizer_generator1



def apply_phase(m,stabilizer_generator):
    '''Applies phase to the m qubit. The numbering starts from 0.'''
    stabilizer_generator1 = stabilizer_generator.copy()
    for stab in stabilizer_generator1:
        #print(stab[m])
        if stab[m] == 'X':
            stab[m] ='Y'
        elif stab[m] == 'Y':
            stab[m] = 'X'
    return stabilizer_generator1



def apply_CNOT(n1,n2,stabilizer_generator):
    stabilizer_generator1 = stabilizer_generator.copy()
    'Apply CNOT to the n1 and n2 qubit. n1 is the control gate. The numbering starts from 0.'
    for stab in stabilizer_generator1:
        m1 = stab[n1]
        m2 = stab[n2]

        if m1 == 'X':
            if m2 == 'I':
                stab[n1],stab[n2] = 'X','X'
            elif m2 == 'X':
                stab[n1],stab[n2] = 'X','I'
            elif m2 == 'Y':
                stab[n1],stab[n2] = 'Y','Z'
            elif m2 == 'Z':
                stab[n1],stab[n2] = 'Y','Y'
        elif m1 == 'Y':
            if m2 == 'I':
                stab[n1],stab[n2] = 'Y','X'
            elif m2 == 'X':
                stab[n1],stab[n2] = 'Y','I'
            elif m2 == 'Y':
                stab[n1],stab[n2] = 'X','Z'
            elif m2 == 'Z':
                stab[n1],stab[n2] = 'X','Y'
        elif m1 == 'Z':
            if m2 == 'I':
                stab[n1],stab[n2] = 'Z','I'
            elif m2 == 'X':
                stab[n1],stab[n2] = 'Z','X'
            elif m2 == 'Y':
                stab[n1],stab[n2] = 'I','Y'
            elif m2 == 'Z':
                stab[n1],stab[n2] = 'I','Z'
        elif m1 == 'I':
            if m2 == 'I':
                stab[n1],stab[n2] = 'I','I'
            elif m2 =='X':
                stab[n1],stab[n2] = 'I','X'
            elif m2 == 'Y':
                stab[n1],stab[n2] = 'Z','Y'
            elif m2 == 'Z':
                stab[n1],stab[n2] = 'Z','Z'
        
    return stabilizer_generator1



def time_reverse_measure(stabilizer_generator,p1,emitter,total_photon,total_emitter):
    '''stabilizer_generator : the set of generators
    p1: the photon being absorbed (the photon numbering starts from 0)
    emitter: the emitter used for absorbtion (the emitter numbering starts from 1)
    total_photon: the total number of photons 
    total_emitter: the total number_emitters'''
    
    """Algorithm: The action of the time-reverse measurement is to exchange the stabilizer of p1 and emitter 
    and apply a H on the both. Due to the structure of the stabilizer generator for using the time reversed meas
    all the stab of photon except 1 would be converted to I and the emitter stab would be H transform of the stab of p1  """
    stabilizer_generator1 = stabilizer_generator.copy()
    stabilizer_generator1 = echelon_gauge(stabilizer_generator1)
    gate_set = []
    for stab in stabilizer_generator1:
        stab[p1],stab[emitter+total_photon-1] = stab[emitter+total_photon-1],stab[p1] # exchange the emitter and p1 stabilizers
    new_stabilizer = apply_hadamard(emitter+total_photon-1,stabilizer_generator1) # apply hadamard on the emitter
    gate_set.append(['time_reverse',p1,emitter+total_photon-1])
    
    
    
    for stab in stabilizer_generator1: 
        if stab.count('I') != total_photon+total_emitter-1: # convert all except one stab of p1 to I. 
            stab[p1] = 'I'
    
    return stabilizer_generator1,gate_set



def emitter_disentangle(stabilizer_generator1,total_emitter,total_photon):
    """We disentangle the emitters. 
    Inputs: 
    stabilizer_generator : the stabilizer generator for photons + emitter
    total_emitter : the total number of emitters
    total_photon: the total number of photons.
    Output:
    The new stabilizer generator, gates applied during this step, the emitters that are free (the ordering is total photon+ free_emitter)."""
    stabilizer_generator = stabilizer_generator1.copy()
    emit_index_list = []
    free_emitter_list = []
    gate_set = []
    # In this for loop we search for stabilizers that are acting just on the emitters
    for stab in stabilizer_generator: 
        if stab[:-total_emitter].count('I') == total_photon:
            emit_index_list.append(stabilizer_generator.index(stab))
            #print(emit_index_list)
    
    # in this loop, we apply gates to disentangle the emitters.       
    for emit in emit_index_list:
        for j in range(len(stabilizer_generator[0])):
            if stabilizer_generator[emit][j] == 'X':
                stabilizer_generator = apply_hadamard(j,stabilizer_generator)
                gate_set.append(['h',j])
            elif stabilizer_generator[emit][j] =='Y':
                stabilizer_generator = apply_phase(j,stabilizer_generator)
                stabilizer_generator = apply_hadamard(j,stabilizer_generator)
                gate_set.append(['phase',j])
                gate_set.append(['h',j])
        #print(stabilizer_generator)
        zindex = [i1 for i1, x in enumerate(stabilizer_generator[emit]) if x == 'Z']
        #print(zindex)
        free_emitter = zindex[0]
        free_emitter_list.append(free_emitter)
        zindex.pop(0)
        for k1 in zindex:
            #print(k1)
            stabilizer_generator = apply_CNOT(k1,free_emitter,stabilizer_generator)
            gate_set.append(['CNOT',k1,free_emitter])
            #print(stabilizer_generator)
            
        # In this step, we make sure to remove redundant Zs present in the stabilizer generator. 
        for j in range(0,len(stabilizer_generator)): 
            
            if (stabilizer_generator[j][free_emitter]!= 'I' and stabilizer_generator[j] != stabilizer_generator[emit]):
                stabilizer_generator[j] = stab_multi(stabilizer_generator[j],stabilizer_generator[emit])

        
    return stabilizer_generator,gate_set,free_emitter_list




def normal_stab(stabilizer_generator1,photon_absorbtion_num,total_photon,total_emitter):
    """
    stabilizer_generator1: the stabilizer generators
    photon_absorbtion_num: the photon to be absorbed (the ordering starts from 0)
    total_photon: the total number of photons
    total_emitter: the total numer of emitters
    
    The function returns the stabilizer to use for photon absorbtion. In principle, there can be multiple stabilizer
    that we can use. This function returns the stabilizer that will require minimal gates
    
    """
    stabilizer_generator1 = echelon_gauge(stabilizer_generator1)
    stabilizer_generator = stabilizer_generator1.copy()
    total_systems = total_photon+total_emitter
    possible_stabili = []
    number_emitter = []
    gate_set = []
    #We get the desired stabilizer with minimum number of emitters in use
    
    
    for stab in stabilizer_generator:
        a3 = photon_absorbtion_num
        if (stab[0:a3].count('I') == photon_absorbtion_num and stab[a3]!='I'):
            possible_stabili.append(stabilizer_generator.index(stab))
            number_emitter.append(stab[total_photon:].count('I'))
    
    
    stabilizer_used_index = number_emitter.index(min(number_emitter))
    stabilizer_use = stabilizer_generator[possible_stabili[stabilizer_used_index]]

    return stabilizer_use
    
    
    
def photon_absorbtion(stabilizer_generator1,photon_absorbtion_num,total_photon,total_emitter):
    """In this step, we absorb the photon. 
    Inputs: the stabilizer generators of the state, the photon to be absorbed, total photon numbers and total number of emitters.
    Outputs: the stabilizer_generator of the state, the gate set, the leftover free emitters."""
    stabilizer_generator1 = echelon_gauge(stabilizer_generator1)
    stabilizer_generator = stabilizer_generator1.copy()
    total_systems = total_photon+total_emitter
    possible_stabili = []
    number_emitter = []
    gate_set = []
    
    stabilizer_use = normal_stab(stabilizer_generator1,photon_absorbtion_num,total_photon,total_emitter)

    
    #----- we disentangled the emitters and absorb the photons -----
    for i in range(0, total_systems):
        if stabilizer_use[i]=='X':
            apply_hadamard(i,stabilizer_generator)
            gate_set.append(['h',i])
        elif stabilizer_use[i]=='Y':
            apply_phase(i,stabilizer_generator)
            apply_hadamard(i,stabilizer_generator)
            gate_set.append(['phase',i])
            gate_set.append(['h',i])
    zindex = [i1 for i1, x in enumerate(stabilizer_use) if x == 'Z']

    
    free_emitter = zindex[-1]
    zindex.remove(free_emitter)

    for k1 in zindex:
        if k1 >= total_photon:
            stabilizer_generator = apply_CNOT(k1,free_emitter,stabilizer_generator)
            gate_set.append(['CNOT',k1,free_emitter])
    #print(stabilizer_generator)
    stabilizer_generator = apply_CNOT(free_emitter,photon_absorbtion_num,stabilizer_generator)
    gate_set.append(['CNOT',free_emitter,photon_absorbtion_num])
    
    

    #----take care of the extra 'Z' lying around----
    for stab in stabilizer_generator: 
        if (stab[photon_absorbtion_num]== 'Z' and stab != stabilizer_use):
            #print(stab)
            stab[photon_absorbtion_num]='I'
    return stabilizer_generator,gate_set,free_emitter



def photon_absorbtion_stab_specified(stabilizer_generator1,photon_absorbtion_num,total_photon,total_emitter,stab,emitter_to_use):
    """In this step, we absorb the photon. 
    Inputs: the stabilizer generators of the state, the photon to be absorbed, total photon numbers,total number of emitters, stabilizer to use
    the emitter to use for absorbtion (the count starts from total_photon).
    Output: stabilizer_generator of the state and the gates applied"""

    stabilizer_use1 = stabilizer_generator1.index(stab)
    #stabilizer_generator1 = echelon_gauge(stabilizer_generator1)
    #stabilizer_used_index = number_emitter.index(min(number_emitter))

    stabilizer_use = stabilizer_generator1[stabilizer_use1]
    gate_set = []
    total_systems = total_photon + total_emitter
    
    for i in range(0, total_systems):
        if stabilizer_use[i]=='X':
            apply_hadamard(i,stabilizer_generator1)
            gate_set.append(['h',i])
        elif stabilizer_use[i]=='Y':
            apply_phase(i,stabilizer_generator1)
            apply_hadamard(i,stabilizer_generator1)
            gate_set.append(['phase',i])
            gate_set.append(['h',i])
    zindex1 = [i1 for i1, x in enumerate(stabilizer_use) if x == 'Z']

    [zindex1.remove(emitter_to_use)]

    for k1 in zindex1:
        if k1 >= total_photon:
            stabilizer_generator = apply_CNOT(k1,emitter_to_use,stabilizer_generator1)
            gate_set.append(['CNOT',k1,emitter_to_use])
    stabilizer_generator = apply_CNOT(emitter_to_use,photon_absorbtion_num,stabilizer_generator1)
    gate_set.append(['CNOT',emitter_to_use,photon_absorbtion_num])
    
    

    for stab in stabilizer_generator: 
        if (stab[photon_absorbtion_num]== 'Z' and stab != stabilizer_use):
            stab[photon_absorbtion_num]='I'
    return stabilizer_generator,gate_set



def stabilizer_gen(n,node,connected_node):
    stab = ['I' for i in range(0,n)]
    stab[node] = 'X'
    for node1 in connected_node: 
        stab[node1] = 'Z'
    return stab


def first_step(stabilizer_generator, no_of_emitters, no_of_photons):
    """In the first step, we absorb as many photons as possible. The algorithm takes in as input: 
    the stabilizer generator, the numbe of emitters and number of photons"""
    
    """The output gives the stabilizer generator, the gates applied in this time step, the photons absorbed and the next photon to absorb"""
    
    stabilizer_generator = echelon_gauge(stabilizer_generator)
    stabilizer_generator1 = stab_with_emitter(stabilizer_generator,no_of_emitters)
    gate_applied = []
    photon_to_absorb = no_of_photons-1
    free_emitter_list1 = [i for i in range(1,no_of_emitters+1)]
    gate_111=[]
    total_systems = no_of_emitters + no_of_photons
    photon_absorb_time_count = []
    photon_list = []
    
    for j in free_emitter_list1:
        stabilizer_generator1,gate1 = time_reverse_measure(stabilizer_generator1,photon_to_absorb,j,no_of_photons,no_of_emitters)
        stabilizer_generator1 = echelon_gauge(stabilizer_generator1)
        gate_applied.append(gate1)
        photon_list.append(photon_to_absorb)
        photon_to_absorb = photon_to_absorb - 1
        
    photon_absorb_time_count.append(photon_list)
    return stabilizer_generator1, gate_applied, photon_absorb_time_count, photon_to_absorb 
    
    
def search_for_stab(stabilizer_generator1, no_of_emitters, no_of_photons, emitter_used_list, photon_to_absorb):
    """In this function, we try to find a stabilizers so that the photon can be absorbed in parrallel."""
    """Inputs: stabilizer generator, number of emitters, number of photons, list of emitters that have already been used, and the photon to absorb."""
    """Outputs: the stabilizer to use, and emitter_to_use_list"""
    
    
    
    """The function follows the following algorithm: 
    1. We first check if there exists a stabilizer with non-Pauli operator only on the photon to absorb --- apart from the emitters
    2. If we find such a stabilizer, we check if the emitters that have been used before consists of only I or Z. 
    3. If so, we check if there exists an emitter having a non-trivial Pauli on the unused emitters. 
    4. If we find such a stabilizer, we exit the loop.
    
    """
    
    stab = []
    number1 = []
    emitter_to_use_list = []
    for number in range(len(stabilizer_generator1)):
        stab_check_bool = []
        # checking for first condition
        if (stabilizer_generator1[number][:photon_to_absorb].count('I') == photon_to_absorb and stabilizer_generator1[number][photon_to_absorb]!='I'):
            # checking for second condition
            for j222 in emitter_used_list:
                if ((stabilizer_generator1[number][j222]=='I' or stabilizer_generator1[number][j222]=='Z')):
                    stab_check_bool.append(1)

                else:
                     stab_check_bool.append(0) # if one of the emitter used fails, we can't do parralellization




            if 0 not in stab_check_bool: 
                # checking for condition 3
                a = stabilizer_generator1[number].copy()
                a1 = [i for i,x in enumerate(a) if x!='I']

                a2 = [x1 for x1 in a1 if x1 not in range(no_of_photons)]# list of emitters with non-trivial Pauli on them
                emitter_to_use1 = list(set(a2) - set(emitter_used_list)) # changed ^->-

                if emitter_to_use1 != []:
                    stab.append(stabilizer_generator1[number].copy())
                    number1.append(number)
                    emitter_to_use_list.append(emitter_to_use1)
                    break;

    return stab,number1,emitter_to_use_list

def start_parrallelization(stabilizer_generator1, no_of_emitters, no_of_photons, emitter_used_list, photon_to_absorb,entropy_all,free_emitter_list11,photon_list,gate_applied):
    """In this function, we start the parallel absorbtion process. We need to specify the emitter_used_list, free_emitter_list11 as well
    
    
    The algorithm works as follows:
    
    1. There are two while loops, the first one to check if there exists more free emitters. Free emitters can be used for time-reversed measurement. 
    2. The second one exists for checking if absorbtion is possible in parrallel.
    
    
    """
    free_emitter_use1 = 0
    stab = []
    #photon_list = []
    stabilizer_generator1 = echelon_gauge(stabilizer_generator1)
    #gate_applied = []
    while free_emitter_use1 == 0:
        
        absorbtion_parrallel = 0
        if len(emitter_used_list) == no_of_emitters:
            #print('exiting loop')
            break

        while absorbtion_parrallel == 0:
            if len(emitter_used_list) == no_of_emitters:
                #print('exiting loop')
                break

            if entropy_all[photon_to_absorb] >= entropy_all[photon_to_absorb-1]:

                    stab,number1,emitter_to_use_list = search_for_stab(stabilizer_generator1, no_of_emitters, no_of_photons, emitter_used_list, photon_to_absorb)
            else: 
                stab = []



            if stab == []:
                absorbtion_parrallel = 1

            if absorbtion_parrallel ==0:
#                 print(stab)
                emitter_to_use = emitter_to_use_list[0][0]
                stab_consider = stab[0]
                #emitter_to_use = emitter_to_use_list2[0]
                stabilizer_generator1, gate3 = photon_absorbtion_stab_specified(stabilizer_generator1,photon_to_absorb,no_of_photons,no_of_emitters,stab_consider,emitter_to_use)
                gate_applied.append(gate3)
                #print(photon_to_absorb)
                photon_list.append(photon_to_absorb)
                photon_to_absorb = photon_to_absorb - 1
                stabilizer_generator1 = echelon_gauge(stabilizer_generator1)

                emitter_used_list.append(emitter_to_use)

        if (free_emitter_list11!=[] and photon_to_absorb>=0):

            emitter_num = free_emitter_list11[0]
            stabilizer_generator1,gate1 = time_reverse_measure(stabilizer_generator1,photon_to_absorb,emitter_num,no_of_photons,no_of_emitters)
            gate_applied.append(gate1)
            photon_list.append(photon_to_absorb)
            photon_to_absorb = photon_to_absorb - 1
            
            emitter_used_list.append(emitter_num+ no_of_photons-1) # make this consistent with the notation. 
            free_emitter_list11.pop(0)
            stabilizer_generator1 = echelon_gauge(stabilizer_generator1)

        else: 
            free_emitter_use1=1
            #print('exiting loop')
                    
    return stabilizer_generator1, gate_applied, photon_to_absorb,photon_list


def algorithm2(stabilizer_generator, no_of_emitters, no_of_photons):
    """In this function, obtain the circuits, and the count for the photon emission process."""
    
    """It takes as input the stabilizer generator, the number of emitters, and the number of photons. 
    
    The algorithm works as follows: 
    1. We first absorb as many photons as possible through time reverse measurement. 
    2. We then, free as many emitters as possible. These can be used later for photon absorbtion. 
    3. We next try to absorb the photon by checking the entropy conditions. 
    4. We then start the parralellization process. 
    5. process2-5 continue till no more photons are left to absorb, 
    
    Output: 
    stabilizer_generator1: the final stabilizer_generator1 for sanity purpose
    
    gate_111: A lists of list telling us the gate sequence.
    
    count: the total number of timesteps taken
    
    photon_absorb_time_count: the photons that are absorbed in each time step. It is a list of lists. 
    """
    
    
    gate_111 = []
    stabilizer_generator = echelon_gauge(stabilizer_generator)
    entropy_all = entropy(stabilizer_generator,no_of_photons)
    stabilizer_generator1, gate_applied,photon_absorb_time_count, photon_to_absorb = first_step(stabilizer_generator, no_of_emitters, no_of_photons)

    gate_111.append(gate_applied)
    count = 1
    free_emitter_list1 = []
    gate_applied = []

    
    
    while photon_to_absorb >= 0:
        print('*****')
        print(photon_to_absorb)
        #--------------------------- We obtain the list of free emitters
        
        photon_list = []
        emitter_used_list = []

        if entropy_all[photon_to_absorb]< no_of_emitters:

            stabilizer_generator1,gate2,free_emitter_list2 = emitter_disentangle(stabilizer_generator1,no_of_emitters,no_of_photons)
            gate_applied.append(gate2)
            free_emitter_list3 = [emit111- no_of_photons+1 for emit111 in free_emitter_list2]
            free_emitter_list1 = free_emitter_list1 + free_emitter_list3
        free_emitter_list11 = list(set(free_emitter_list1))
        if no_of_emitters-entropy_all[photon_to_absorb] != len(free_emitter_list11):
            print('going bad')
        #---------------------------------------- free_emitter_list11 is the list of free emitters
        # start the absorbtion process
        stabilizer_generator1 = echelon_gauge(stabilizer_generator1)
        
        
        if entropy_all[photon_to_absorb] >= entropy_all[photon_to_absorb-1]:
            print(photon_to_absorb)

            stabilizer_generator1, gate3, emitter_used= photon_absorbtion(stabilizer_generator1,photon_to_absorb,no_of_photons,no_of_emitters)
            gate_applied.append(gate3)
            photon_list.append(photon_to_absorb)
            photon_to_absorb = photon_to_absorb - 1
            stabilizer_generator1 = echelon_gauge(stabilizer_generator1)
            emitter_used_list.append(emitter_used)

        
        

        if photon_to_absorb >= 0:

            stabilizer_generator1, gate_applied, photon_to_absorb,photon_list = start_parrallelization(stabilizer_generator1, no_of_emitters, no_of_photons, emitter_used_list, photon_to_absorb,entropy_all,free_emitter_list11,photon_list,gate_applied)





        count = count +1
        gate_111.append(gate_applied)
        gate_applied = []
        free_emitter_list1=[]
        photon_absorb_time_count.append(photon_list)
    
    
    return stabilizer_generator1,gate_111,count, photon_absorb_time_count
