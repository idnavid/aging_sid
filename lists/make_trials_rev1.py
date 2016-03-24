#! /usr/bin/python
import sys

def store_files(age_dict,age, data_list, train_session, test_session):
    fin = open(data_list)
    for i in fin:
        session = i.strip().split('_')[4]
        if (int(session) in range(1,7)):
            age_groups[age]['train'].append(i.strip().split('/')[-1])
        elif (int(session) in range(11,14)):
            age_groups[age]['test'].append(i.strip().split('/')[-1])
    return age_dict


def check_spkrsession(basename1,basename2):
    """
        string format: Gender_AgeRange_SpeakerID_SessionID_ChunkID_ChunkID2
        e.g.: M_A_AAA_11_01_02.htk
        """
    list1 = basename1.split('_')
    list2 = basename2.split('_')
    if (list1[2]==list2[2] and list1[3]==list2[3]):
        return -1
    elif (list1[2]==list2[2] and list1[3]!=list2[3]):
        return 1
    elif list1[2]!=list2[2]:
        return 2

def calculate_diff(basename1,basename2):
    session1 = int(basename1.split('_')[3])
    session2 = int(basename2.split('_')[3])
    return abs(session1 - session2)

fin = open('all_marp_male.lst')
all_marp = []
for i in fin:
    all_marp.append(i.strip())

fin.close()

fout = open('marp_trials_diff7','w')
for i in all_marp:
    for j in all_marp:
        if i!=j:
            base1 = i.split('/')[-1]
            base2 = j.split('/')[-1]
            diff = calculate_diff(base1,base2)
            if (diff == 7):
                state = check_spkrsession(base1,base2)
                if state==-1:
                    pass
                elif state==1:
                    fout.write(base1.split('.')[0]+' '+base2.split('.')[0]+" target\n")
                elif state==2:
                    fout.write(base1.split('.')[0]+' '+base2.split('.')[0]+" nontarget\n")
fout.close()










