#!/usr/bin/python

def collect_expected_seqlengths(inf_ribos, prot_dna):
    from Bio import SeqIO
    dicty = {}
    lens = []
    ribos = []
    for inf in inf_ribos:
        print("inf here", inf)
        inf_handle = open(inf, 'r')
        inf_rec = SeqIO.read(inf_handle, 'fasta')
        inf_id_split = inf_rec.id.split('_')
        inf_id = inf_id_split[1]
        print("inf_id", inf_id, inf_id_split)
        if prot_dna == 'protein':
            seqlen = len(inf_rec.seq)
        elif prot_dna == 'dna':
            seqlen = 3*len(inf_rec.seq)
        try:
            dicty[inf_id].append(int(seqlen))
        except:
            dicty[inf_id] = int(seqlen)
    print(dicty)
    return dicty
    

def concatenate_diamond_matches(infile, prot_dna):
    import collections
    import datetime
    import os
    from Bio import SeqIO

    handle = open(infile, 'r')
   # inf_two = open('atccs.txt', 'r')
   # inf_ribos = [lin.rstrip() for lin in inf_two]
   # dicty = collect_expected_seqlengths(inf_ribos, prot_dna)
   # ribo_names_field = 1 #ribo_names_field
   # print("ribo_names_field", ribo_names_field)
    d = collections.defaultdict(dict)
    records = list(SeqIO.parse(handle, "fasta"))
    outf_strains = set()
    #ribosomal_proteins = ["rpsB","rpsC","rpsD","rpsE","rpsF","rpsG","rpsH","rpsI","rpsJ","rpsK","rpsL","rpsM","rpsO","rpsP","rpsQ","rpsR","rpsS","rpsT","rplA","rplB","rplC","rplD","rplE","rplF","rplI", "rplJ", "rplK", "rplL", "rplM", "rplN", "rp10", "rplP", "rplQ", "rplR", "rplS", "rplT", "rplU", "rplV", "rplW", "rplX", "rpmA", "rpmC", "rpmI", "ftsZ","metK","uvrC","pyrH","purH","miaA","recN","recA","murD","mraY","truB","tsaD","tpiA","rsmH","dfp","gpsA","typA","rsmE","adk","rimM","gmk","ribF","secY","nusB","rsmG","rimP","recR","frr","pth","atpG","smpB","recO","infC","lepA","tsaE","rsmD","ybeY","birA","pyrB","coaD","rnc","grpE","atpH","rbfA","rsfS","infA","dnaA","aspS","purL","purF","serS","nusA","tyrS","atpD","rpoC","pgk","pheS","guaA","radA","eno","rpoA","proS","metG","aroA","atpA","murB","ruvB","alaS","purM","purD","tig","uvrB","aroC","era","der","purA","topA","ruvX","rpoB","rseP","mnmA"];
    ribosomal_proteins = ["rplN","rplP", "rplR", "rplB", "rplV", "rplX", "rplC", "rplD", "rplE", "rplF", "rpsJ", "rpsQ", "rpsS", "rpsC", "rpsH", "rp10"]
    outfile_strains = open("strains_missing_ribos.txt", 'a+')
##uncomment print statements for debugging, ID changes can cause concatenation to fail

#Loop through annotation and create a dictionary containing each ribosomal sequence per isolate genome

    for rec in records:
        id = rec.id.split('|')
        print("id to split on |", id)
        shortid = id[0]#.split('_')
        print("shortid to split on _", shortid)
        print("length of id", len(id))
        if len(id) > 1:
            newd = id[0].split('_')
            newid = "{}_{}".format(newd[0],newd[1])
        elif len(id) == 2:
            newid = id[0]
        else:
            print("id length problem", id)
        print("newid", newid)
        #rib = shortid #[0] #shortid[int(ribo_names_field)]
        rib = id[1]
        print("rib", rib)
        seqy = str(rec.seq)
        if seqy.endswith('*'):
             #remove the stop character sometimes included in protein translation files
             newseq = seqy[:-1]
        else:
             newseq = seqy
        for ribo in ribosomal_proteins:
            if rib == ribo:
                #check for exact match
               # if  len(newseq) - 30 < dicty[rib] < len(newseq) + 30:   
                    try:
                        d[newid][ribo].append(newseq)
                    except:
                        d[newid][ribo] = [newseq]
                    #    print(d[newid][ribo])
                        if len(d[newid][ribo]) > 1:
                            print("duplicate annotation, will search with diamond blast", ribo, newid)
                            del d[newid][ribo]
             #   else:
             #       print("length is incorrect for {}".format(rib), len(newseq), dicty[rib])
             #       outf_strains.add("{}".format(newid))
                   
    #create datestamp file
    #print("this is dictionary", d)
    outname = datetime.date.today().strftime("%d-%m-%y")+"concatenated_ribosomal_proteins_db.fasta"
    if os.path.exists(outname):
       #if files have been collected from the annotation and we are now collecting from diamond blast process
       outname = outname+"_2"
    else:
       #mark strains for diamond blast process if they were missing some sequences
       pass

    outfile = open(outname, 'a+')
    for key, value in d.items():
         #test we have all 16 sequences, otherwise mark genome for diamond blast searches or ignore
         #print("key {} length value {}".format(key, value))
         if len(value) == 16:
                print("ok in")
                # try:
                joinstring = "{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}".format(''.join(d[key]["rplN"]),''.join(d[key]["rplP"]),''.join(d[key]["rplR"]),''.join(d[key]["rplB"]),''.join(d[key]["rplV"]),''.join(d[key]["rplX"]),''.join(d[key]["rplC"]),''.join(d[key]["rplD"]),''.join(d[key]["rplE"]),''.join(d[key]["rplF"]),''.join(d[key]["rpsJ"]),''.join(d[key]["rpsS"]),''.join(d[key]["rpsC"]),''.join(d[key]["rpsH"]),''.join(d[key]["rpsH"]),''.join(d[key]["rp10"]))
                outfile.write(">{}\n".format(key))
                outfile.write(''.join(joinstring)+"\n")
                 #except:
                 #    print("excepted")
                 #    pass
                 
         else:
              print(key, "Missing some ribosomal proteins, will search with Diamond Blast", len(value))
              outf_strains.add("{}\t{}".format(key,value.keys()))
    for stray in outf_strains:
        outfile_strains.write('{}\n'.format(stray))

import sys
concatenate_diamond_matches(sys.argv[1], sys.argv[2])
