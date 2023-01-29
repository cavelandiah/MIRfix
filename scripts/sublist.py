def sublist(queue, configurer, level, filename, args):
    logid = scriptname+'.sublist: '
    configurer(queue, level)
    log.debug(logid+'Starting to process '+str(filename))
    try:
        filesdir=str(args.famdir)#dir for families
        command="list"
        mapfile=str(args.mapping)#mapping of mir to mirfam
        matfile=str(args.mature)#mature sequences
        matrdir=args.maturedir#directory for mature files
        genomes_file=args.genomes #file that contains all genomes
        tProcessed=0#all start with 't', are for the total of all families
        tRemoved=0
        tTotalnumberofSequences=0
        tpredicted=0
        tnoprediction=0
        tflippednotcorrected=0
        tflippedcorrected=0
        twith2mats=0
        twith0mats=0
        twith1mats=0
        toldshanon=0
        tnewshanon=0
        numberoffamilies=0
        tcountcorrected=0
        #moved from global:
        listofnew=[]
        listofnewloop=[]
        listoldstatus=[]
        listofoldloop=[]#add seq id and old seq
        listofold=[]#add seq id and old seq
        listofboth=[]#add the id, that shouldn't be added to the new file
        listremovedbroken=[]
        listremovedscore=[]
        listremovedN=[]
        listofmirstar=[]
        listnomat=[]
        list2mat=[]
        list2matcoor=[]
        list1matcoor=[]
        listmatcoor=[]
        listnogenomes=[]#list of IDs that no existing genome to search, for their species
        listnotingenome=[]#list of IDs that are not found in the genome searched, for their species
        listremovenoloop = [] #Collect candidates without loop after curation
        listhighmfe = [] #Collect all candidates with high MFE after all curation
        listlonghairpin = [] #Collect all candidates with a final reported hairpin
        listbadmature = [] #Collect all candidates that failed mature prediction
        nomats=0
        templong=[]
        userflanking=0
        listmisalignedcorr=[]
        listgoodnew=[]
        listcorrected=[]
        listcorrectedori=[]

        log.debug(logid+filename)
        numberoffamilies+=1
        countcorrected=0
        countcorrectedTonew=0
        Processed=0
        Removed=0
        TotalnumberofSequences=0
        predicted=0
        noprediction=0
        flippednotcorrected=0
        flippedcorrected=0
        with2mats=0
        with0mats=0
        with1mats=0
        oldshanon=0
        newshanon=0
        Totalproc=[]
        Totalrem=[]
        Nomatspred=[]
        Nomatsnotpred=[]
        flippednotcorr=[]
        flippedcorr=[]
        seq0mat=[]
        seq1mat=[]
        seq2mat=[]
        newshan=[]
        oldshan=[]
        suma=[]
        sumall=[]
        nomats=0#in case no mat exists and no mats at all
        listnomatremoved=[]
        flagnomatexists=False

        ## prepare output directory for current family
        filename = str( filename).strip()
        outdir = str( args.outdir) + filename + ".out/"
        makeoutdir( outdir, args.force)
        filen = filesdir + filename + ".fa"

        OldShanon=0
        NewShanon=0
        tempcountsucnomat=0
        infile=""
        outfile=""
        infile=filen[:]
        outfile=outdir+filename+"-tempshan.aln"
        clustaline = ClustalwCommandline("clustalw2", infile=infile, outfile=outfile)
        stdoutshan,stdershan=clustaline()
        alignTostock(outfile)
        OldShanon=CalShanon(outfile+'.stk')
        log.debug(logid+str(["OldShanon",OldShanon]))
        os.remove(outfile)
        os.remove(filesdir+filename+".dnd")
        os.remove(outfile+".stk")
        infile=""
        outfile=""
        userflanking=int(args.extension)
        indexed_genomes = index_genomes(genomes_file) # Create blastn databases

        with openfile(filen) as fl:
            for rec in SeqIO.parse(fl,'fasta'):
                pidsplit=rec.description.split()
                pid=str(pidsplit[1])
                mat2seq=str(rec.seq)
                coorflag=0
                smat=""
                emat=""
                with openfile(mapfile) as mf:
                    for line in mf:
                        linesplit=line.split()
                        if len(linesplit)>8 and pid in line:
                            list2mat.append(linesplit[2].strip())
                            numofmat=((len(linesplit)+1)//3)-1
                            smat=linesplit[4]#first mat ID
                            emat=linesplit[4+numofmat-1].strip()#last mat ID
                            mtf = openfile(matfile)
                            first = None
                            second = None
                            for k in SeqIO.parse(mtf,'fasta'):
                                if smat in k.description:
                                    first = str(k.seq)
                                if emat in k.description:
                                    second = str(k.seq)
                                if first and second:
                                    list2mat.append(first)#add the first mat seq
                                    list2mat.append(second)#add the second mat seq
                                    break
                            if first and not second:
                                list2mat.append(first)#add the first mat seq
                                list2mat.append('')#empty string for second
                            if second and not first:
                                list2mat.append('')#empty string for first
                                list2mat.append(second)#add the second mat seq

        log.debug(logid+str(["list 2 mat is here",list2mat]))

        if os.path.isfile(outdir+'nomat-'+filename.strip()+'.fa'):#before calling checknomat in Submit
            log.debug(logid+"The file "+filename.strip()+" already processed, will be done again")
            os.remove(outdir+'nomat-'+filename+'.fa')

        flagnomatexists,nomats,listremovedbroken,listremovedscore,listremovedN,listnomat=checknomat(filen,mapfile,matfile,outdir,filename,listremovedbroken,listremovedscore,listremovedN,nomats,listnomat)

        if flagnomatexists and nomats!=-1:
            log.debug(logid+"flagnomatexists"+str(flagnomatexists)+';'+str(nomats))
            if os.path.isfile(outdir+filename.strip()+"-new.fa"):
                os.remove(outdir+filename.strip()+"-new.fa")

            listnomatremoved=listremovedbroken+listremovedscore+listremovedN
            log.debug(logid+str(["all removed",listnomatremoved]))

            if len(listnomatremoved)>0:
                log.debug(logid+"listnomatremoved"+str(listnomatremoved))
                fl = openfile(filen)
                for record in SeqIO.parse(fl, 'fasta'):
                    tempdes=record.description
                    tempdeslst=record.description.split()
                    if tempdeslst[1].strip() not in listremovedscore and tempdeslst[1].strip() not in listremovedbroken and tempdeslst[1].strip() not in listremovedN:# and tempdeslst[1].strip() not in lst2mat:
                        newprecfile=open(outdir+filename+"-new.fa",'a')
                        newprecfile.write(">"+str(tempdes)+"\n"+str(record.seq)+"\n")
                        newprecfile.close()

            for f in [outdir+"tempfold.fa", outdir+"tempmat.fa", outdir+"temptoalign.aln",outdir+"temptoalign.dnd", outdir+"temptoalign.fa", outdir+"temptoalign.txt", outdir+"temptofold.fa", outdir+"temptopredict.aln", outdir+"temptopredict.dnd", outdir+"temptopredict.fa"]:
                if os.path.isfile(f):
                    os.remove(f)

            if len(listnomatremoved)>0:
                listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, listnogenomes, listnotingenome, templong, listgoodnew=flip(filename.strip(), outdir+filename+"-new.fa", outdir, mapfile, matfile, listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, list2mat, listnogenomes, listnotingenome, templong, listgoodnew, indexed_genomes, args)
            else:
                listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, listnogenomes, listnotingenome, templong, listgoodnew=flip(filename.strip(), filen, outdir, mapfile, matfile, listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, list2mat, listnogenomes, listnotingenome, templong, listgoodnew, indexed_genomes, args)

        elif flagnomatexists and nomats==-1:
            log.debug(logid+"flagnomatexists"+str(flagnomatexists)+';'+str(nomats))
            with open(outdir+filename.strip()+"-summ.txt","a") as summaryfile:
                summaryfile.write("no matures for the sequences and no related matures in the mapping file\n")

        elif not flagnomatexists:
            log.debug(logid+str(["this flip",flagnomatexists]))
            listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, listnogenomes, listnotingenome, templong, listgoodnew=flip(filename.strip(), filen, outdir, mapfile, matfile, listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, list2mat, listnogenomes, listnotingenome, templong, listgoodnew, indexed_genomes, args)#filename: filename/family, filen: the file itself(with the directory)

        log.debug(logid+"listofnew: "+str(listofnew))
        if os.path.isfile(outdir+filename+"-new.fa"):
            os.remove(outdir+filename+"-new.fa")

        if len(listofnew)>0:
            for i in range(0,len(listofnew)):
                if i%2==0:
                    log.debug(logid+"new"+str([">"+listofnew[i],listofnew[i+1]]))
                    with open(outdir+filename+"-res.fa","a") as familyfileres:
                        familyfileres.write(">"+str(listofnew[i])+"\n"+str(listofnew[i+1].replace('T','U'))+"\n")

        if len(listofnewloop)>0:
            for i in range(0,len(listofnewloop)):
                if i%2==0:
                    log.debug(logid+"new loop"+str([">"+listofnewloop[i],listofnewloop[i+1]]))
                    with open(outdir+filename+"-res.fa","a") as familyfileres:
                        familyfileres.write(">"+str(listofnewloop[i])+"\n"+str(listofnewloop[i+1].replace('T','U'))+"\n")

        if len(listofold)>0:
            for i in range(0,len(listofold)):
                if i%2==0:
                    log.debug(logid+"old"+str([">"+listofold[i],listofold[i+1]]))
                    with open(outdir+filename+"-res.fa","a") as familyfileres:
                        familyfileres.write(">"+str(listofold[i])+"\n"+str(listofold[i+1].replace('T','U'))+"\n")

        if len(listofoldloop)>0:
            for i in range(0,len(listofoldloop)):
                if i%2==0:
                    log.debug(logid+"old loop"+str([">"+listofoldloop[i],listofoldloop[i+1]]))
                    with open(outdir+filename+"-res.fa","a") as familyfileres:
                        familyfileres.write(">"+str(listofoldloop[i])+"\n"+str(listofoldloop[i+1].replace('T','U'))+"\n")

        listnomatbroken=listremovedbroken
        listnomatscore=listremovedscore
        listnomatN = listremovedN

        if len(list2mat)>0:
            log.debug(logid+"list2mat: "+str(list2mat))
            for i in range(0,len(list2mat),3):
                precID=list2mat[i]
                fmatseq=list2mat[i+1].strip()#first mat id
                ematseq=list2mat[i+2].strip()#last mat id
                fl = openfile(filen)
                for record in SeqIO.parse(fl,'fasta'):
                    pdescsplit=record.description.split()
                    pid=str(pdescsplit[1].strip())
                    pseq=str(record.seq)
                    xcut=userflanking
                    ycut=userflanking
                    precDes=record.description
                    precitem=precDes.split()
                    specie=precitem[2].strip()+" "+precitem[3]
                    mspos=None
                    mepos=None
                    mspos=pseq.find(fmatseq)
                    mepos=pseq.rfind(ematseq)  # rfind returns the last index of the match

                    xcutseq=len(pseq[:mspos])#used in case the seq not found in the genome
                    ycutseq=len(pseq[mepos+1:])

                    if precID==pid:
                        # On miRBase some mature were assigned for the same miRNAs.
                        if mspos == -1 or mepos == -1:
                            log.error(logid+'Referred mir or mir* from '+pid+' did not fit into the precursor, maybe some reference is repeated on mapping file?')
                            sys.exit()

                        long2matseq=""
                        long2matseq, listnogenomes, listnotingenome=getindex2mat(pseq.replace("U", "T"), specie, precID, precDes, listnogenomes, listnotingenome, args)  # from here on we have the reverse complement if seq on minus strand
                        log.debug(logid+'long2matseq: '+str(long2matseq))

                        if long2matseq!="":
                            (fspos, fepos, cutpseq) = find_positions(long2matseq, pseq.replace("U", "T"),
                                                                             fmatseq.replace("U", "T"), ematseq.replace("U", "T"),
                                                                             userflanking)

                            with open(outdir+filename+"-res.fa","a") as familyfileres:
                                familyfileres.write(">"+str(record.description)+"\n"+str(cutpseq).replace('T','U')+"\n")
                            with open(outdir+filename.strip()+"-Final.fasta","a") as familyfileresfinal:
                                familyfileresfinal.write(">"+str(record.description)+"\n"+str(cutpseq).replace('T','U')+"\n")

                        elif long2matseq=="":
                            if xcutseq>=xcut and xcut>=0 and xcut<=50:
                                fspos=mspos-xcut
                            elif xcutseq<xcut and xcut>=0 and xcut<=50:#if requested flankings more than available for this seq, because not found in genome,then take all available
                                fspos=0
                            elif xcutseq>10 and (xcut>50 or xcut<0):
                                fspos=mspos-10
                            elif xcutseq<=10 and (xcut>50 or xcut<0):
                                fspos=0

                            if ycutseq>=ycut and ycut>=0 and ycut<=50:
                                fepos=mepos+ycut+1
                            elif ycutseq<ycut and ycut>=0 and ycut<=50:
                                fepos=len(pseq)
                            elif ycutseq>10 and (ycut>50 or ycut<0):
                                fepos=mepos+11
                            elif ycutseq<=10 and (ycut>50 or ycut<0):
                                fepos=len(pseq)

                            cutpseq=pseq[fspos:fepos+len(ematseq)]
                            with open(outdir+filename+"-res.fa","a") as familyfileres:
                                familyfileres.write(">"+str(record.description)+"\n"+str(cutpseq).replace('T','U')+"\n")
                            with open(outdir+filename.strip()+"-Final.fasta","a") as familyfileresfinal:
                                familyfileresfinal.write(">"+str(record.description)+"\n"+str(cutpseq).replace('T','U')+"\n")


        if len(listofmirstar)>0:
            mirstarfile=open(outdir+filename.strip()+"-mirstar.fa","a")
            mirstarmapfile=open(outdir+filename.strip()+"-mirstar-map.txt","a")

            for star in range(0,len(listofmirstar),4):#devided into 4s..#1:star desc 2:mirstar seq 3:positions 4:related prec ID
                mirstarfile.write(">"+str(listofmirstar[star].strip())+" "+str(listofmirstar[star+3])+"\n"+str(listofmirstar[star+1])+"\n")
                tempsplit=listofmirstar[star].split()#to get only the ID of the mirstar
                starsplit=listofmirstar[star].split()
                mirstarmapfile.write(str(listofmirstar[star+3])+" "+str(starsplit[1])+" "+str(listofmirstar[star+2])+"\n")
            mirstarmapfile.close()
            mirstarfile.close()
            del listofmirstar[:]

        if ".fa" in filename:
            filename=filename[:filename.find(".fa")]
            resultfastafile=outdir+filename.strip()+"-res.fa"
        else:
            resultfastafile=outdir+filename.strip()+"-res.fa"

        rff = openfile(resultfastafile)
        for frec in SeqIO.parse(rff, 'fasta'):
            resprecdesc=str(frec.description)
            resfilesplit=(frec.description).split()
            resprecid=str(resfilesplit[1].strip())
            mat2seq=str(frec.seq).replace('T','U')
            mat1seq=str(frec.seq).replace('T','U')
            log.debug(["res",resprecid,mat2seq,mat1seq])

            mf = openfile(mapfile)
            for line in mf:
                coorflag =0   #
                linesplit=line.split()
                if len(linesplit)>8 and resprecid.strip() in linesplit:
                    numofmat=((len(linesplit)+1)//3)-1
                    list2matcoor.append(linesplit[2].strip())
                    firstmat=linesplit[4]
                    lastmat=linesplit[4+numofmat-1]

                    mtf = openfile(matfile)
                    first = None
                    second = None
                    curmatseq = None
                    curmatstar = None
                    # Here use find because both matures were reported at the
                    # Beginning and we are iterating over the complete mature file
                    for record in SeqIO.parse(mtf, 'fasta'):
                        sequence=str(record.seq)
                        if firstmat.strip() == record.description.split(" ")[1]: # specific label
                            first = 'Found'
                            curmatseq = sequence

                        if lastmat.strip() == record.description.split(" ")[1]: # specific label
                            second = 'Found'
                            curmatstar = sequence
                    (startmat, startmatstar, finalseq) = find_positions(mat1seq, mat1seq,curmatseq, curmatstar,userflanking)
                    log.debug(logid+str([curmatseq,startmat, startmatstar, finalseq]))
                    if startmat == None or startmatstar == None:
                        log.debug(logid+'Not possible to locate the mapping referred mir or mir* on the mature file for '+resprecid+' with '+mat2seq)
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        finalseq=curmatseq
                        list2matcoor.append("NULL")#matstarID/here no mir* found
                        list2matcoor.append("NULL")#matstarID/here no mir* found
                        list2matcoor.append(startmat)
                        list2matcoor.append(endmat)
                        list2matcoor.append(startmatstar)
                        list2matcoor.append(endmatstar)
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        listbadmature.append(resprecid)
                        break

                    endmat = startmat+len(curmatseq)-1
                    endmatstar = (startmatstar+len(curmatstar))-1
                    log.debug(logid+str(["heres new",mat2seq,finalseq,startmat,endmat,startmatstar,endmatstar,curmatseq,curmatstar]))
                    if first and second and startmat != None and startmatstar != None:
                        coorflag=1

                    if coorflag==1:
                        list2matcoor.append(firstmat.strip())
                        list2matcoor.append(lastmat.strip())
                        list2matcoor.append(startmat)
                        list2matcoor.append(endmat)
                        list2matcoor.append(startmatstar)
                        list2matcoor.append(endmatstar)

                elif len(linesplit)<8 and resprecid in linesplit and resprecid.strip() not in templong:
                    star=False
                    nstar=False
                    curmatID=linesplit[4]
                    startmat=0
                    endmat=0
                    startmatstar=0
                    endmatstar=0
                    list1matcoor.append(resprecid)
                    finalseq="A"
                    startfinalseq=0
                    endfinalseq=0
                    mtf = openfile(matfile)
                    for reco in SeqIO.parse(mtf, 'fasta'):
                        splitreco=reco.description.split()

                        if curmatID == splitreco[1]:
                            curmatseq=str(reco.seq)
                            nstar=True
                            break

                    mtfs = openfile(outdir+filename.strip()+"-mirstar.fa")
                    for starrec in SeqIO.parse(mtfs,'fasta'):
                        curmatsplit=(starrec.description).split()
                        curmatsplit1=(curmatsplit[1]).split('-')
                        starrecID=curmatsplit1[0]

                        if (curmatID).strip()==(starrecID).strip() and (resprecid.strip() in starrec.description):
                            curmatstar=str(starrec.seq)
                            star=True
                            break

                    # Here genome did not exists, then it is reeplaced by the same precursor seq to obtain
                    # the mature coordinates
                    (coortemp1, coortemp2, finalseq) = find_positions(mat1seq, mat1seq,
                                                 curmatseq, curmatstar,
                                                 userflanking)

                    if coortemp2<coortemp1:
                        tempseqex=curmatseq[:]
                        curmatseq=curmatstar[:]
                        curmatstar=tempseqex[:]

                    log.debug(logid+str(["no long",coortemp1,coortemp2,curmatseq,curmatstar]))
                    startmat=coortemp1
                    endmat=(startmat+len(curmatseq))-1

                    if star:
                        startmatstar=coortemp2
                        endmatstar=startmatstar+int(len(curmatstar)-1)
                        log.debug(logid+str(['coor star not',startmatstar,endmatstar,finalseq]))#,startmatstarlong,endmatstarlong)

                    if star and nstar and startmatstar!=-1 and endmatstar!=-1 and startmatstar>endmat:
                        list1matcoor.append(curmatID.strip())#mat original
                        list1matcoor.append(curmatsplit[1].strip())#matstarID
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        with open(outdir+filename.strip()+"-Final.fasta","a") as familyfileresfinal:
                            familyfileresfinal.write(">"+resprecdesc+"\n"+finalseq+"\n")
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        break
                    elif star and nstar and startmatstar!=-1 and endmatstar!=-1 and startmatstar<endmat:#to put them in order
                        list1matcoor.append(curmatsplit[1].strip())#matstarID
                        list1matcoor.append(curmatID.strip())#mat original
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        with open(outdir+filename.strip()+"-Final.fasta","a") as familyfileresfinal:
                            familyfileresfinal.write(">"+resprecdesc+"\n"+finalseq+"\n")
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        break

                    elif (not star and nstar) or startmatstar==-1 or endmatstar==-1:
                        startmat=int(mat1seq.find(curmatseq))
                        endmat=startmat+len(curmatseq)-1
                        finalseq=mat1seq
                        startmstar=0
                        endmatstar=len(mat1seq)-1
                        list1matcoor.append(curmatID.strip())#mat original
                        list1matcoor.append("NULL")#matstarID/here no mir* found
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        with open(outdir+filename.strip()+"-Final.fasta","a") as familyfileresfinal:
                            familyfileresfinal.write(">"+resprecdesc+"\n"+finalseq+"\n")
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        break

                    elif (not star and not nstar) or startmatstar==-1 or endmatstar==-1:#impossible but just in case
                        startmat=0
                        endmat=len(mat1seq)-1
                        startmstar=0
                        endmatstar=len(mat1seq)-1
                        finalseq=mat1seq
                        list1matcoor.append("NULL")#mat original/here no mir
                        list1matcoor.append("NULL")#matstarID/here no mir* found
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        break

                elif len(linesplit)<8 and resprecid in linesplit and resprecid.strip() in templong:
                    star=False
                    nstar=False
                    curmatID=linesplit[4]
                    curmatseq=None
                    curmatstar=None
                    startmat=0
                    endmat=0
                    startmatstar=0
                    endmatstar=0
                    longseq=""
                    startmat=0
                    finalseq="A"
                    startmatstar=0
                    endmatstar=0
                    list1matcoor.append(resprecid)
                    indexlongmat=int(templong.index(resprecid.strip()))
                    longseq=str(templong[indexlongmat+1]).replace('T','U')

                    with openfile(matfile) as mtf:
                        for reco in SeqIO.parse(mtf, 'fasta'):
                            splitreco=reco.description.split()
                            if curmatID == splitreco[1]:
                                curmatseq=str(reco.seq)
                                nstar=True
                                break

                    log.debug(logid+'Curmatseq: '+str(curmatseq))

                    with openfile(outdir+filename.strip()+"-mirstar.fa") as msfa:
                        for starrec in SeqIO.parse(msfa,'fasta'):
                            curmatsplit=(starrec.description).split()
                            curmatsplit1=(curmatsplit[1]).split('-')
                            starrecID=curmatsplit1[0]
                            if (curmatID).strip()==(starrecID).strip():
                                curmatstar=str(starrec.seq)
                                star=True
                                break

                    if not curmatstar:
                        log.error(logid+'Not possible to define curmatstar for '+ resprecid +' in '+str(matfile)+' and '+str(outdir+filename.strip()+"-mirstar.fa"))
                        sys.exit()

                    log.debug(["coor1temp",longseq,curmatseq,curmatstar,resprecid])
                    (coortemp1, coortemp2, finalseq) = find_positions(longseq, mat2seq,
                                                 curmatseq, curmatstar,
                                                 userflanking)

                    if coortemp1 == None or coortemp2 == None:
                        log.debug(logid+'Not possible to locate miR or miR* in ' + resprecid + " " + curmatseq)
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        finalseq=curmatseq
                        list1matcoor.append("NULL")#mat original/here no mir
                        list1matcoor.append("NULL")#matstarID/here no mir* found
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        listbadmature.append(resprecid)
                        break

                    if coortemp2<coortemp1:
                        tempseqex=curmatseq[:]
                        curmatseq=curmatstar[:]
                        curmatstar=tempseqex[:]

                    #Look for mir in final seq:
                    startmirfinal = int(finalseq.find(curmatseq))
                    endmirfinal = int(startmirfinal + len(curmatseq) - 1)
                    subsetnomir = finalseq[endmirfinal:]
                    indexupdate = finalseq.find(subsetnomir)
                    #Look mirstar in final seq
                    startmirstarfinal = subsetnomir.find(curmatstar)
                    endmirstarfinal = int(startmirstarfinal + len(curmatstar) - 1)
                    #update coordinates
                    startmat=startmirfinal
                    endmat=endmirfinal
                    startmatstar = startmirstarfinal + indexupdate
                    endmatstar = endmirstarfinal + indexupdate
                    log.debug(['coor',startmat, endmat, startmatstar,endmatstar])
                    # Previous evaluations:
                    (structure, mfe, loopsize, length_precursor) = evaluate_final_hairpin(finalseq, startmat, endmat, startmatstar, endmatstar, resprecid)

                    if mfe > -10:
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        finalseq=mat1seq
                        list1matcoor.append("NULL")#mat original/here no mir
                        list1matcoor.append("NULL")#matstarID/here no mir* found
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        listhighmfe.append(resprecid)
                        break
                    if length_precursor > 200:
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        finalseq=mat1seq
                        list1matcoor.append("NULL")#mat original/here no mir
                        list1matcoor.append("NULL")#matstarID/here no mir* found
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        listlonghairpin.append(resprecid)
                        break
                    if loopsize <= 1:
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        finalseq=mat1seq
                        list1matcoor.append("NULL")#mat original/here no mir
                        list1matcoor.append("NULL")#matstarID/here no mir* found
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        listremovenoloop.append(resprecid)
                        break

                    if star and nstar and startmatstar!=-1 and endmatstar!=-1 and startmatstar>endmat and loopsize > 1:
                        list1matcoor.append(curmatID.strip())#mat original
                        list1matcoor.append(curmatsplit[1].strip())#matstarID
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        with open(outdir+filename.strip()+"-Final.fasta","a") as familyfileresfinal:
                            familyfileresfinal.write(">"+resprecdesc+"\n"+finalseq+"\n")
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        break

                    elif star and nstar and startmatstar!=-1 and endmatstar!=-1 and startmatstar<endmat and loopsize > 1:#to put them in order
                        list1matcoor.append(curmatsplit[1].strip())#matstarID
                        list1matcoor.append(curmatID.strip())#mat original
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        with open(outdir+filename.strip()+"-Final.fasta","a") as familyfileresfinal:
                            familyfileresfinal.write(">"+resprecdesc+"\n"+finalseq+"\n")
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        break

                    # Here is not possible to calculate the loop
                    elif (not star and nstar) or startmatstar==-1 or endmatstar==-1:
                        startmat=int(mat1seq.find(curmatseq))
                        endmat=(startmat+len(curmatseq))-1
                        finalseq=mat1seq
                        startmstar=0
                        endmatstar=len(mat1seq)-1
                        list1matcoor.append(curmatID.strip())#mat original
                        list1matcoor.append("NULL")#matstarID/here no mir* found
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        with open(outdir+filename.strip()+"-Final.fasta","a") as familyfileresfinal:
                            familyfileresfinal.write(">"+resprecdesc+"\n"+finalseq+"\n")
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        break

                    elif (not star and not nstar) or startmatstar==-1 or endmatstar==-1:
                        startmat=0
                        endmat=len(mat1seq)-1
                        startmstar=0
                        endmatstar=len(mat1seq)-1
                        finalseq=mat1seq
                        list1matcoor.append("NULL")#mat original/here no mir
                        list1matcoor.append("NULL")#matstarID/here no mir* found
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        break

        log.debug([len(list2matcoor), len(list1matcoor), list2matcoor, list1matcoor])
        listmatcoor=list2matcoor+list1matcoor

        mi=0
        mk=0
        log.debug(len(listmatcoor))
        log.debug(len(listmatcoor)/7)

        while mi<=int(len(listmatcoor))-7:
            mk=mk+1
            r=0

            for n in range(mk+1,int(len(listmatcoor)/7)+1):
                with open(outdir+filename.strip()+"-Final.anc","a") as anchorcoorfile:
                    if listmatcoor[mi+1] == "NULL" or listmatcoor[mi+3] == -1 or listmatcoor[mi+10+r] == -1 or listmatcoor[mi+3] == 0 or listmatcoor[mi+10+r] == 0:
                        continue
                    else:
                        anchorcoorfile.write(str(mk)+" "+str(n)+" "+str(listmatcoor[mi+3]+1)+" "+str(listmatcoor[mi+10+r]+1)+" "+str(22)+" "+str(1)+"\n")
                    if listmatcoor[mi+1] == "NULL" or listmatcoor[mi+5] == -1 or listmatcoor[mi+12+r] == -1 or listmatcoor[mi+5] == 0 or listmatcoor[mi+12+r] == 0:
                        continue
                    else:
                        anchorcoorfile.write(str(mk)+" "+str(n)+" "+str(listmatcoor[mi+5]+1)+" "+str(listmatcoor[mi+12+r]+1)+" "+str(22)+" "+str(1)+"\n")
                r=r+7

            mi=mi+7

        log.debug(["here listmatcoor 1 ",listmatcoor])
        maxidesc=0

        finalstk=open(outdir+filename.strip()+'.stk','a')
        finalstk.write('# STOCKHOLM 1.0\n')
        if matrdir:
            fs=os.environ["DIALIGN2_DIR"]=matrdir
            log.debug(matrdir)
        f1=os.popen("dialign2-2 -n -anc -fa "+outdir+filename.strip()+'-Final.fasta')
        log.debug(f1)
        f1.close()

        if os.path.isfile(outdir+filename.strip()+'-Final.fa'):
            doalifold(outdir+filename.strip()+"-Final.fa",outdir)
            for rec in SeqIO.parse(openfile(outdir+filename.strip()+'-Final.fa'),'fasta'):
                if len(rec.description.strip())>maxidesc:
                    maxidesc=len(rec.description.strip())

                if maxidesc<len('#=GC SS_cons'):
                    maxidesc=len('#=GC SS_cons')

            for rec in SeqIO.parse(openfile(outdir+filename.strip()+'-Final.fa'),'fasta'):
                finalstk.write(str(rec.description.strip())+" "*(maxidesc-len(rec.description.strip())+2)+str(rec.seq)+"\n")

            struct=getstructure(outdir+'alifoldtemp.txt')
            os.remove(outdir+'alifoldtemp.txt')
            os.remove(outdir+filename.strip()+"-Final.ali")
            finalstk.write('#=GC SS_cons'+" "*(maxidesc-len('#=GC SS_cons')+2)+str(struct))
            finalstk.close()

            stki = openfile(outdir+filename.strip()+'.stk')
            alignment=AlignIO.read(stki,"stockholm")
            countcorrected=0
            countcorrectedTonew=0
            sthalfgaps=0
            ndhalfgaps=0
            sthalfsum=0
            ndhalfsum=0
            stnucnum=0
            ndnucnum=0
            numofseqs=0
            totalstnucnum=0
            totalndnucnum=0

            for record in alignment:
                lenseq=len(str(record.seq))
                seq=str(record.seq)
                sthalf=seq[0:int((lenseq/2))]
                ndhalf=seq[int((lenseq/2)):]
                sthalfgaps=sthalf.count('-')
                stnucnum=len(sthalf)-sthalfgaps
                totalstnucnum=totalstnucnum+stnucnum
                ndhalfgaps=ndhalf.count('-')
                ndnucnum=len(ndhalf)-ndhalfgaps
                totalndnucnum=totalndnucnum+ndnucnum
                sthalfsum=sthalfgaps+sthalfsum
                ndhalfsum=ndhalfgaps+ndhalfsum
                numofseqs=numofseqs+1
            log.debug([sthalfsum,ndhalfsum,totalstnucnum,totalndnucnum])
            # Here calculated the average of nt on left side of align (stnucavg) and on the right side (ndnucavg) on the alignment
            stnucavg=totalstnucnum/numofseqs
            ndnucavg=totalndnucnum/numofseqs
            alignment=AlignIO.read(outdir+filename.strip()+'.stk',"stockholm")
            # Here iterate again but with calculated average scores
            for record in alignment:
                lenseq=len(str(record.seq))
                seq=str(record.seq)
                corid=str(record.id)
                sthalf=seq[0:int((lenseq/2))]
                ndhalf=seq[int((lenseq/2)):]
                sthalfgaps=sthalf.count('-')
                stnucnum=len(sthalf)-sthalfgaps
                ndhalfgaps=ndhalf.count('-')
                ndnucnum=len(ndhalf)-ndhalfgaps
                nucnum=stnucnum+ndnucnum

                if (stnucnum>(stnucavg*1.5) and (ndnucnum<(ndnucavg/2) or ndnucnum==0) ) or (stnucnum>0.7*nucnum and (ndnucnum<(ndnucavg/2) or ndnucnum==0)):#-(stnucavg*0.1)):
                    log.debug([record.id,'st'])
                    log.debug([stnucnum,ndnucnum,nucnum])
                    log.debug(["list new good",corid,listgoodnew])
                    countcorrected,countcorrectedTonew,listmisalignedcorr,listcorrected,listcorrectedori=correct(corid.strip(),userflanking,countcorrected,countcorrectedTonew,listofnew,listofnewloop,listoldstatus,templong,listmisalignedcorr,listcorrected,listcorrectedori,listgoodnew)

                if (ndnucnum>(ndnucavg*1.5) and (stnucnum<(stnucavg/2) or stnucnum==0)) or (ndnucnum>0.7*nucnum and (stnucnum<(stnucavg/2) or stnucnum==0)):#-(ndnucavg*0.1)):
                    log.debug([record.id,'nd'])
                    log.debug([stnucnum,ndnucnum,nucnum])
                    log.debug(["list new good",corid,listgoodnew])
                    countcorrected,countcorrectedTonew,listmisalignedcorr,listcorrected,listcorrectedori=correct(corid.strip(),userflanking,countcorrected,countcorrectedTonew,listofnew,listofnewloop,listoldstatus,templong,listmisalignedcorr,listcorrected,listcorrectedori,listgoodnew)

        if len(listmisalignedcorr)>0:
            familyfilerescorrected=open(outdir+filename+"-corrected.fasta","a")
            mirstarcorrectedfile=open(outdir+filename+"-mirstar-corrected.fa","a")

            for corrrecord in SeqIO.parse(openfile(outdir+filename.strip()+'-Final.fasta'),'fasta'):
                splitdes=(corrrecord.description).split()
                splitID=splitdes[1].strip()

                if splitID.strip() in listmisalignedcorr:
                    coind=listmisalignedcorr.index(splitID.strip())
                    corcoorind=listmatcoor.index(splitID.strip())
                    listmatcoor[corcoorind+3]=listmisalignedcorr[coind+2]#set again the coordinates of the corrected sequence, in the liost related to build the anc file
                    listmatcoor[corcoorind+4]=listmisalignedcorr[coind+3]
                    listmatcoor[corcoorind+5]=listmisalignedcorr[coind+4]
                    listmatcoor[corcoorind+6]=listmisalignedcorr[coind+5]
                    familyfilerescorrected.write(">"+str(corrrecord.description)+"\n"+str(listmisalignedcorr[coind+1])+"\n")

                else:
                    familyfilerescorrected.write(">"+str(corrrecord.description)+"\n"+str(corrrecord.seq)+"\n")

            familyfilerescorrected.close()
            log.debug(["here listmatcoor 2 ",listmatcoor])

            for mstarrec in SeqIO.parse(openfile(outdir+filename.strip()+'-mirstar.fa'),'fasta'):
                splitdes=(mstarrec.description).split()
                descstar=str(mstarrec.description)
                oriseq=str(mstarrec.seq)
                splitID=str(splitdes[-1]).strip()

                if splitID.strip() in listmisalignedcorr:
                    newcorrind=listmisalignedcorr.index(splitID.strip())
                    orien=str(listmisalignedcorr[newcorrind+6])

                    if orien=='3p':
                        newsstar=listmisalignedcorr[newcorrind+2]
                        newestar=listmisalignedcorr[newcorrind+3]
                    else:
                        newsstar=listmisalignedcorr[newcorrind+4]
                        newestar=listmisalignedcorr[newcorrind+5]

                    correctseq=str(listmisalignedcorr[newcorrind+1])
                    newmatstar=correctseq[newsstar:newestar+1]
                    mirstarcorrectedfile.write(">"+descstar.strip()+"\n"+newmatstar+"\n")
                else:
                    mirstarcorrectedfile.write(">"+descstar.strip()+"\n"+oriseq+"\n")

            mirstarcorrectedfile.close()

            log.debug(["here listmatcoor 3 ",listmatcoor])
            mi=0
            mk=0
            while mi<=int(len(listmatcoor))-7:
                mk=mk+1
                r=0
                for n in range(mk+1,int(len(listmatcoor)/7)+1):
                    with open(outdir+filename.strip()+"-corrected.anc","a") as anchorcoorfilecorrected:
                        anchorcoorfilecorrected.write(str(mk)+" "+str(n)+" "+str(listmatcoor[mi+3])+" "+str(listmatcoor[mi+10+r])+" "+str(22)+" "+str(1)+"\n")
                        anchorcoorfilecorrected.write(str(mk)+" "+str(n)+" "+str(listmatcoor[mi+5])+" "+str(listmatcoor[mi+12+r])+" "+str(22)+" "+str(1)+"\n")
                    r=r+7

                mi=mi+7

            maxidesc=0
            finalstkcorrected=open(outdir+filename.strip()+'corrected.stk','a')
            finalstkcorrected.write('# STOCKHOLM 1.0\n')
            if matrdir:
                fe=os.environ["DIALIGN2_DIR"]=matrdir
            f11=os.popen("dialign2-2 -n -anc -fa  "+outdir+filename.strip()+'-corrected.fasta')
            f11.close()

            doalifold(outdir+filename.strip()+"-corrected.fa",outdir)

            for rec in SeqIO.parse(openfile(outdir+filename.strip()+'-corrected.fa'),'fasta'):
                if len(rec.description.strip())>maxidesc:
                    maxidesc=len(rec.description.strip())

                if maxidesc<len('#=GC SS_cons'):
                    maxidesc=len('#=GC SS_cons')

            for rec in SeqIO.parse(openfile(outdir+filename.strip()+'-corrected.fa'),'fasta'):
                finalstkcorrected.write(str(rec.description.strip())+" "*(maxidesc-len(rec.description.strip())+2)+str(rec.seq)+"\n")

            struct=getstructure(outdir+'alifoldtemp.txt')
            os.remove(outdir+'alifoldtemp.txt')
            os.remove(outdir+filename.strip()+"-corrected.fa")
            os.remove(outdir+filename.strip()+"-corrected.ali")
            finalstkcorrected.write('#=GC SS_cons'+" "*(maxidesc-len('#=GC SS_cons')+2)+str(struct))
            finalstkcorrected.close()
            NewShanon=CalShanon(outdir+filename.strip()+'corrected.stk')

        else:
            if os.path.isfile(outdir+filename.strip()+'-Final.fa'):
                os.remove(outdir+filename.strip()+"-Final.fa")
                NewShanon=CalShanon(outdir+filename.strip()+'.stk')
                log.debug(["stk file studied is: "+outdir+filename.strip()+'.stk'])
                log.debug(["final.fa exists",NewShanon])

            else:
                NewShanon=OldShanon
                log.debug(["final.fa DO NOT exists",NewShanon])

        log.debug(["new shanon",NewShanon])
        finalcoor=open(outdir+filename+"-FinalCoor.txt","a")
        log.debug(["here listmatcoor 4 ",len(listmatcoor),listmatcoor])

        for l in range(0,len(listmatcoor),7):
            finalcoor.write(str(listmatcoor[l]).strip()+" "+str(listmatcoor[l+1]).strip()+" "+str(listmatcoor[l+2]).strip()+" "+str(listmatcoor[l+3]).strip()+" "+str(listmatcoor[l+4]).strip()+" "+str(listmatcoor[l+5]).strip()+" "+str(listmatcoor[l+6]).strip()+"\n")

        finalcoor.close()

        with open(outdir+filename.strip()+"-summ.txt","a") as summaryfile:
            summaryfile.write("---------------------------Original precursors with bad positioned matures ----------------------------\n")
            if len(listofoldloop)>0:
                for k in range(0,len(listofoldloop)):
                    if k%2==0:
                        tempsplit=listofoldloop[k].split()
                        summaryfile.write(str(tempsplit[1].strip())+"\n")
            else:
                summaryfile.write("---> NO Original precursors with bad positioned matures\n")
            summaryfile.write("\n")
            summaryfile.write("---------------------------flipped/changed precursors ----------------------------\n")

            if len(listofnew)>0:
                for k in range(0,len(listofnew)):
                    if k%2==0:
                        tempsplit=listofnew[k].split()
                        summaryfile.write(str(tempsplit[1].strip())+"\n")
            else:
                summaryfile.write("--->NO flipped/changed precursors\n")
            summaryfile.write("\n")
            summaryfile.write("---------------------------flipped/changed precursors; inloop ----------------------------\n")

            if len(listofnewloop)>0:
                for k in range(0,len(listofnewloop)):
                    if k%2==0:
                        tempsplit=listofnewloop[k].split()
                        summaryfile.write(str(tempsplit[1].strip())+"\n")
            else:
                summaryfile.write("---> NO flipped/changed precursors; inloop\n")
            summaryfile.write("\n")

            summaryfile.write("---------------------------Flipped/changed precursors don't fit with the final alignment (changed back to original)----------------------------\n")
            if len(listcorrected)>0:
                for k in (listcorrected):
                    summaryfile.write(k.strip()+"\n")
            else:
                summaryfile.write("--->None of the changed/flipped precursors was changed back to original\n")

            summaryfile.write("\n")

            if len(listcorrectedori)>0:
                summaryfile.write("---------------------------Original precursors, Change at the Alignment level (to fit the alignment----------------------------\n")
                for k in (listcorrectedori):
                    summaryfile.write(k.strip()+"\n")

            summaryfile.write("\n")
            summaryfile.write("---------------------------Precursors without annotated mature----------------------------:\n")

            if len(listnomat)>0:
                summaryfile.write("#Precursors without annotated mature; Successfully predicted mature---------------------------\n")
                tempcountsucnomat=0
                for k in range(0,len(listnomat)):
                    if listnomat[k] not in listnomatbroken and listnomat[k] not in listnomatscore and listnomat[k] not in listnomatN:
                        tempcountsucnomat=tempcountsucnomat+1
                        summaryfile.write(str(listnomat[k].strip())+"\n")

                if tempcountsucnomat==0:
                    summaryfile.write("---> NO Precursors without annotated mature; Successfully predicted mature\n")

                summaryfile.write("\n")
                summaryfile.write("#Precursors without annotated mature; bad positioned predicted mature---------------------------\n")

                if len(listnomatbroken)>0:
                    for k in listnomatbroken:
                        summaryfile.write(k+"\n")

                else:
                    summaryfile.write("--->NO Precursors without annotated mature; bad positioned predicted mature:\n")

                summaryfile.write("\n")
                summaryfile.write("#Precursors without annotated mature; no similar mature---------------------------\n")

                if len(listnomatscore)>0:
                    for k in listnomatscore:
                        summaryfile.write(k+"\n")

                else:
                    summaryfile.write("--->NO Precursors without annotated mature; no similar mature\n")

            else:
                summaryfile.write("--->All precursors have annotated matures\n")

            summaryfile.write("\n")
            summaryfile.write("---------------------------Original precursors totally removed ----------------------------\n")

            if len(listofboth)>0:
                for i in range(0,len(listofboth)):
                    if i%2==0:
                        tempsplit=listofboth[i].split()
                        summaryfile.write(str(tempsplit[1].strip())+"\n")
            elif len(listnomatN)>0:
                for j in range(0, len(listnomatN)):
                    summaryfile.write(str(listnomatN[j].strip()) + "\n")
            elif len(listremovenoloop)>0:
                for k in range(0, len(listremovenoloop)):
                    summaryfile.write(str(listremovenoloop[k].strip()) + "\n")
            elif len(listhighmfe)>0:
                for k in range(0, len(listhighmfe)):
                    summaryfile.write(str(listhighmfe[k].strip()) + "\n")
            elif len(listbadmature)>0:
                for k in range(0, len(listbadmature)):
                    summaryfile.write(str(listbadmature[k].strip()) + "\n")
            elif len(listlonghairpin)>0:
                for k in range(0, len(listlonghairpin)):
                    summaryfile.write(str(listlonghairpin[k].strip()) + "\n")
            else:
                summaryfile.write("---> NO Original precursors totally removed ----------------------------\n")

            summaryfile.write("\n")
            summaryfile.write("#Precursors without loop regions (<= 1 nt)---------------------------\n")
            if len(listremovenoloop)>0:
                for j in range(0, len(listremovenoloop)):
                    summaryfile.write(str(listremovenoloop[j].strip()) + "\n")
            else:
                summaryfile.write("---> All precursors reported a valid loop region ----------------------------\n")

            summaryfile.write("\n")
            summaryfile.write("# Precursors with high N content---------------------------\n")
            if len(listnomatN)>0:
                for j in range(0, len(listnomatN)):
                    summaryfile.write(str(listnomatN[j].strip()) + "\n")
            else:
                summaryfile.write("---> All precursors reported valid nucleotide composition with N content > 60 % ----------------------------\n")

            summaryfile.write("\n")
            summaryfile.write("# Precursors with high MFE---------------------------\n")
            if len(listhighmfe)>0:
                for j in range(0, len(listhighmfe)):
                    summaryfile.write(str(listhighmfe[j].strip()) + "\n")
            else:
                summaryfile.write("---> All precursors reported a valid MFE range ----------------------------\n")

            summaryfile.write("\n")
            summaryfile.write("# Precursors with unsuccessful mature prediction---------------------------\n")
            if len(listbadmature)>0:
                for j in range(0, len(listbadmature)):
                    summaryfile.write(str(listbadmature[j].strip()) + "\n")
            else:
                summaryfile.write("---> All precursors reported canonical mature positions ----------------------------\n")

            summaryfile.write("\n")
            summaryfile.write("# Precursors with reported precursor size > 200 nt---------------------------\n")
            if len(listlonghairpin)>0:
                for j in range(0, len(listlonghairpin)):
                    summaryfile.write(str(listlonghairpin[j].strip()) + "\n")
            else:
                summaryfile.write("---> All precursors reported a valid size ----------------------------\n")

            summaryfile.write("\n")


            summaryfile.write("---------------------------Precursors without given genome file(s)----------------------------\n")

            if len(listnogenomes)>0:
                for i in listnogenomes:
                    summaryfile.write(i.strip()+"\n")

            else:
                summaryfile.write("---> All the precursors have given genomes to search----------------------------\n")

            summaryfile.write("\n")
            summaryfile.write("---------------------------Precursors NOT found in their given genomes file(s)----------------------------\n")

            if len(listnotingenome)>0:
                for i in listnotingenome:
                    summaryfile.write(i.strip()+"\n")

            else:
                summaryfile.write("--->All the precursors were found in their given genomes file(s)----------------------------\n")


        #Total
        Processed=(len(listofold)/2)+(len(list2mat)/3)+(len(listofoldloop)/2)+(len(listofnew)/2)+(len(listofnewloop)/2)#+len(listnomat)
        Removed=(len(listofboth)/2 + len(listnomatN) + len(listremovenoloop) + len(listhighmfe) + len(listbadmature) + len(listlonghairpin))
        TotalnumberofSequences=Processed+Removed
        #no mats/predicted-no prediction
        predicted=tempcountsucnomat
        noprediction=len(listnomat)-tempcountsucnomat-len(listnomatN)-len(listremovenoloop)-len(listhighmfe)-len(listbadmature)-len(listlonghairpin)

        #Corrected/notcorrected
        flippednotcorrected=int((len(listofnew)/2))+int((len(listofnewloop)/2))-int(countcorrected)
        flippedcorrected=int(countcorrected)+int(countcorrectedTonew)

        #sequences distribution
        with2mats=len(list2mat)/3
        with0mats=len(listnomat)
        with1mats=TotalnumberofSequences-with2mats-with0mats

        #Shanon
        oldshanon=OldShanon/10
        newshanon=NewShanon/10

        tProcessed+=Processed#all start with 't', are for the total of all families
        tRemoved+=Removed
        tTotalnumberofSequences+=TotalnumberofSequences
        tpredicted+=predicted
        tnoprediction+=noprediction
        tflippednotcorrected+=flippednotcorrected
        tflippedcorrected+=flippedcorrected
        twith2mats+=with2mats
        twith0mats+=with0mats
        twith1mats+=with1mats
        toldshanon+=oldshanon
        tnewshanon+=newshanon
        tcountcorrected=int(countcorrected)+int(countcorrectedTonew)+tcountcorrected
        #Total
        Totalproc = [Processed, 0, 0, 0, 0]
        Totalrem = [Removed, 0, 0, 0, 0]
        #nomats
        Nomatspred = [0,0,0,predicted,0]
        Nomatsnotpred = [0,0,0,noprediction,0]
        #flipped
        flippednotcorr=[0,0,flippednotcorrected,0,0]
        flippedcorr=[0,0,flippedcorrected,0,0]
        #sequences
        seq0mat=[0,with0mats,0,0,0]
        seq1mat=[0,with1mats,0,0,0]
        seq2mat=[0,with2mats,0,0,0]

        #Shanon
        newshan=[0,0,0,0,newshanon]
        oldshan=[0,0,0,0,oldshanon]
        suma=[]
        sumall=[]

        for i in range(0,len(Totalproc)):
            suma.append(seq0mat[i]+seq1mat[i])

        sumall.append(Totalproc[0])
        sumall.append(Totalrem[0])
        sumall.append(Nomatspred[3])
        sumall.append(Nomatsnotpred[3])
        sumall.append(flippednotcorr[2])
        sumall.append(flippedcorr[2])
        sumall.append(seq0mat[1])
        sumall.append(seq1mat[1])
        sumall.append(seq2mat[1])
        sumall.append(newshan[4])
        sumall.append(oldshan[4])

        with open(outdir+filename.strip()+"-summ.txt","a") as summaryfile:
            summaryfile.write("---------------------------Results In Numbers----------------------------\n")
            summaryfile.write("*Number of remained precursors= "+str(int((len(listofold)/2)+(len(list2mat)/3)))+"\n")
            summaryfile.write("*Number of remained precursors with bad positioned matures= "+str(int((len(listofoldloop)/2)))+"\n")
            summaryfile.write("*Number of flipped(changed) precursors= "+str(int((len(listofnew)/2)))+"\n")
            summaryfile.write("*Number of flipped(changed) precursors with bad positioned matures= "+str(int((len(listofnewloop)/2)))+"\n")
            summaryfile.write("*Number of removed precursors= "+str(int(Removed))+"\n")
            summaryfile.write("*Number of precursors without a given matures= "+str(int(len(listnomat)))+"\n")
            summaryfile.write("*Number of precursors with successfully predicted matures= "+str(int(len(listnomat)-len(listremovedbroken)-len(listremovedscore)-len(listremovedN)-len(listremovenoloop)-len(listhighmfe)-len(listbadmature)-len(listlonghairpin)))+"\n")
            summaryfile.write("*Number of precursors without a given genome file= "+str(int(len(listnogenomes)))+"\n")
            summaryfile.write("*Number of precursors not found in their given genomes= "+str(int(len(listnotingenome)))+"\n")
            summaryfile.write("---------------------------Numbers Used For The Graph----------------------------\n")
            summaryfile.write("Total number of sequences="+str(int(TotalnumberofSequences))+"\n")
            summaryfile.write("Processed="+str(int(Processed))+"\n")
            summaryfile.write("Removed="+str(int(Removed))+"\n")
            summaryfile.write("Predicted="+str(int(predicted))+"\n")
            summaryfile.write("without given mature/no predicted mature="+str(int(noprediction))+"\n")
            summaryfile.write("Changed(flipped)="+str(int(flippednotcorrected))+"\n")
            summaryfile.write("Misaligned (shifted) precursors corrected at the end="+str(int(flippedcorrected))+"\n")
            summaryfile.write("Number of precursors without annotated mature(s)="+str(int(with0mats))+"\n")
            summaryfile.write("Number of precursors with one annotated mature="+str(int(with1mats)-int(with0mats))+"\n")
            summaryfile.write("Number of precursors with two annotated matures="+str(int(with2mats))+"\n")
            summaryfile.write("Old Entropy="+str(oldshanon)+"\n")
            summaryfile.write("New Entropy="+str(newshanon)+"\n")
            summaryfile.close()

        listofoldloopjson=[]
        for k in range(0,len(listofoldloop)):
            if k%2==0:
                tempsplit=listofoldloop[k].split()
                listofoldloopjson.append(str(tempsplit[1].strip()))

        listofnewjson=[]
        for k in range(0,len(listofnew)):
            if k%2==0:
                tempsplit=listofnew[k].split()
                listofnewjson.append(str(tempsplit[1].strip()))

        listofnewloopjson=[]
        for k in range(0,len(listofnewloop)):
            if k%2==0:
                tempsplit=listofnewloop[k].split()
                listofnewloopjson.append(str(tempsplit[1].strip()))

        listsuccpred=[]
        for k in range(0,len(listnomat)):
            if listnomat[k] not in listnomatbroken and listnomat[k] not in listnomatscore and listnomat[k] not in listnomatN and listnomat[k] not in listremovenoloop and listnomat[k] not in listhighmfe and listnomat[k] not in listbadmature and listnomat[k] not in listlonghairpin:
                listsuccpred.append(str(listnomat[k].strip()))

        listremovedjson=[]
        for i in range(0,len(listofboth)):
            if i%2==0:
                tempsplit=listofboth[i].split()
                listremovedjson.append(str(tempsplit[1].strip()))
        for i in range(0,len(listnomatN)):
            listremovedjson.append(listnomatN[i].strip())
        for j in range(0, len(listremovenoloop)):
            listremovedjson.append(listremovenoloop[j].strip())
        for k in range(0, len(listhighmfe)):
            listremovedjson.append(listhighmfe[k].strip())
        for k in range(0, len(listbadmature)):
            listremovedjson.append(listbadmature[k].strip())
        for k in range(0, len(listlonghairpin)):
            listremovedjson.append(listlonghairpin[k].strip())
        data = {
            "Processed":int(Processed),
            "Precursors totally removed":{
                "Number":int(len(listremovedjson)),
                "IDs":listremovedjson,
            },
            "Remained precursors":int((len(listofold)/2)+(len(list2mat)/3)-len(listnomatN)-len(listremovenoloop)-len(listhighmfe)-len(listbadmature)-len(listlonghairpin)),
            "Remained precursors with bad positioned matures":{
                "Number":int(len(listofoldloopjson)),#int((len(listofoldloop)/2)),
                "IDs":listofoldloopjson,
            },
            "Flipped(changed) precursors":{
                "Number":int(len(listofnewjson)),
                "IDs":listofnewjson,
            },
            "Flipped(changed) precursors with bad positioned matures":{
                "Number":int(len(listofnewloopjson)),
                "IDs":listofnewloopjson,
            },
            "Flipped/changed precursors changed back to original":{
                "Number":int(len(listcorrected)),
                "IDs":listcorrected,
            },
            "Flipped/changed precursors changed back to original":{
                "Number":int(len(listcorrected)),
                "IDs":listcorrected,
            },
            "Original precursors, Change at the Alignment level":{
                "Number":int(len(listcorrectedori)),
                "IDs":listcorrectedori,
            },
            "Precursors without a given matures":{
                "Number":int(len(listnomat)-len(listremovedbroken)-len(listremovedscore)),
                "IDs":listnomat,
            },
            "Precursors with successfully predicted matures":{
                "Number":int(len(listsuccpred)),
                "IDs":listsuccpred,
            },
            "Precursors with bad positioned predicted matures":{
                "Number":int(len(listnomatbroken)),
                "IDs":listnomatbroken,
            },
            "Precursors without annotated and without predicted mature":{
                "Number":int(len(listnomatscore)),
                "IDs":listnomatscore,
            },
            "Precursors with high N content":{
                "Number": int(len(listnomatN)),
                "IDs": listnomatN,
            },
            "Precursors without loop":{
                "Number": int(len(listremovenoloop)),
                "IDs": listremovenoloop,
            },
            "Precursors with high MFE":{
                "Number": int(len(listhighmfe)),
                "IDs": listhighmfe,
            },
            "Precursors with unsuccessful mature prediction":{
                "Number": int(len(listbadmature)),
                "IDs": listbadmature,
            },
            "Precursors with long precursors":{
                "Number": int(len(listlonghairpin)),
                "IDs": listlonghairpin,
            },
            "Precursors without a given genome":{
                "Number":int(len(listnogenomes)),
                "IDs":listnogenomes,
            },
            "Precursors not found in their genomes":{
                "Number":int(len(listnotingenome)),
                "IDs":listnotingenome,
            },
            "Old Entropy":oldshanon,
            "New Entropy":newshanon,
        }
        jsonData = json.dumps(data)

        with open(outdir+filename.strip()+'.json', 'w') as f:
            json.dump(jsonData, f)

        Processed=0
        Removed=0
        TotalnumberofSequences=0
        predicted=0
        noprediction=0
        flippednotcorrected=0
        flippedcorrected=0
        with2mats=0
        with0mats=0
        with1mats=0
        oldshanon=0
        newshanon=0
        countcorrected=0
        countcorrectedTonew=0
        del Totalproc[:]
        del Totalrem[:]
        del Nomatspred[:]
        del Nomatsnotpred[:]
        del flippednotcorr[:]
        del flippedcorr[:]
        del seq0mat[:]
        del seq1mat[:]
        del seq2mat[:]
        del newshan[:]
        del oldshan[:]
        del suma[:]
        del sumall[:]
        del listofnew[:]#add seq id and new sequence
        del listofnewloop[:]#add seq id and new sequence
        del listofoldloop[:]#add seq id and old seq
        del listofold[:]#add seq id and old seq
        del listofboth[:]#add the id, that shouldn't be added to the new file
        del listremovedbroken[:]
        del listremovedscore[:]
        del listremovedN[:]
        del listnomat[:]
        del list2mat[:]
        del listmatcoor[:]
        del listnogenomes[:]#list of IDs that no existing genome to search, for their species
        del listnotingenome[:]

        familyfileres.close()
        log.debug("done")
        log.debug([listofold,listofnew,listofnewloop,listofoldloop,listremovedbroken,listremovedscore,listremovedN,listofmirstar])
        log.debug(list2mat)
        if os.path.isfile(outdir+filename.strip()+"-res.fa"):
            os.remove(outdir+filename.strip()+"-res.fa")

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

##############

def main(args):
    try:
        #  Logging configuration
        logdir = args.logdir
        logfile = str.join(os.sep,[os.path.abspath(logdir),scriptname+'.log'])
        makelogdir(logdir)
        makelogfile(logfile)
        #  Multiprocessing with spawn
        nthreads = args.cores or 2
        multiprocessing.set_start_method('spawn')  # set multiprocessing start method to safe spawn
        pool = multiprocessing.Pool(processes=nthreads, maxtasksperchild=1)
        queue = multiprocessing.Manager().Queue(-1)
        listener = multiprocessing.Process(target=listener_process, args=(queue, listener_configurer, logfile, args.loglevel))
        listener.start()
        worker_configurer(queue, args.loglevel)
        log.info(logid+'Running '+scriptname+' on '+str(args.cores)+' cores.')
        log.info(logid+'CLI: '+sys.argv[0]+' '+'{}'.format(' '.join( [shlex.quote(s) for s in sys.argv[1:]] )))
        lfams = []
        with openfile(args.families) as filelist:
            for line in filelist:
                line = line.strip()
                if line.endswith( ".fa"):
                    line = line[:-3]
                ## check if there are already output directories
                ## to prevent errors
                if os.path.exists( args.outdir + line + ".out/") and not args.force:
                    sys.exit( '\nError: At least one output directory already exists! Won\'t override!\nTo override the original output folder, please specify the \'--force\' option.\n')
                lfams.append( line)

        for fam in lfams:
            pool.apply_async(sublist, args=(queue, worker_configurer, args.loglevel, fam, args))
        pool.close()
        pool.join()
        queue.put_nowait(None)
        listener.join()
        sys.exit()

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

##############################MAIN##############################

if __name__ == '__main__':

    logid = scriptname+'.main: '
    try:
        args = parseargs()

        # find_executable('clustalw2') or sys.exit('Please install clustalw2 to run this')
        find_executable('dialign2-2') or sys.exit('Please install dialign2-2 to run this')

        main(args)
        sys.exit()

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))
