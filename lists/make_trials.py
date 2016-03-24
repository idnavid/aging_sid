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

age_groups = {'A':{'train':[],'test':[]}, 'B':{'train':[],'test':[]}, 'C':{'train':[],'test':[]}, 'D':{'train':[],'test':[]}}
train_sesh = int(sys.argv[1])
test_sesh = int(sys.argv[2])


for j in age_groups:
    store_files(age_groups,j,'trial_feat_'+j+'.lst',train_sesh, test_sesh)

#for i in age_groups:
#    for j in age_groups[i]:
#        print i, j, age_groups[i][j]


fout = open('speaker_models','w')
for i in ['A']:
    file_list = age_groups[i]['train']
    for j in file_list:
        spkr = j.split('_')[2]
        fileid = j.split('.')[0]
        fout.write(spkr+' '+fileid+'\n')
fout.close()

for i in age_groups['A']['train']:
    base1 = i.split('.')[0]
    for j in age_groups['A']['test']:
        base2 = j.split('.')[0]
        state = check_spkrsession(base1,base2)
        if state==-1:
            pass
        elif state==1:
            print base1.split('_')[2], base2, "target"
        elif state==2:
            print base1.split('_')[2], base2, "nontarget"









