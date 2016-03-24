#! /usr/bin/python 

import sys

target_age = int(sys.argv[1])

fcallinfo = open('swb2_callinfo.tbl')
call_dict = {}
for i in fcallinfo:
    line_list = i.strip().split(',')
    fileid = line_list[0]
    spkrid = line_list[2]
    spkrch = line_list[3]
    if not(spkrid in call_dict):
        call_dict[spkrid] = [fileid+spkrch]
    else:
        call_dict[spkrid].append(fileid+spkrch)

fcallinfo.close()

fsphlist = open('all_swb2.lst')
file_dict = {}
for i in fsphlist:
    line = i.strip()
    line_list = line.split('/')
    fileid = line_list[-1][:-4]
    file_dict[fileid] = line
fsphlist.close()

fspkrinfo = open('swb2_male_spkrs.tbl')
spkr_dict = {}
for i in fspkrinfo:
    line_list = i.strip().split(',')
    spkrid = line_list[0]
    spkrgender = line_list[1]
    spkrage = line_list[2]
    try:
        if int(spkrage) == target_age:
            for j in call_dict[spkrid]:
                print file_dict[j]+' '+spkrid
    except:
        pass





