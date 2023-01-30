def checknomat(precfile,mapfile,matfile,directory,precfilename,listremovedbroken,listremovedscore,listremovedN,nomats,listnomat):
    logid = scriptname+'.checknomat: '
    try:
        flagnomatexists=False
        directory=os.path.abspath(directory)+"/"
        matIDs=[]
        precnomat=[]
        countnomat=0#count number of precs without mature

        for record in SeqIO.parse(openfile(precfile), 'fasta'):
            precAllID=record.description
            precID=record.id
            precseq=record.seq
            flag=0
            mf = openfile(mapfile)
            for line in mf:
                item=line.split()
                if precID in line:#and not 2 mat
                    flag=1
                    numofmat=((len(item)+1)//3)-1
                    for i in range (0,numofmat):
                        matIDs.append(item[4+i])
            if flag==0:
                subprecfile=open(directory+'nomat-'+precfilename+'.fa','a')
                countnomat=countnomat+1
                precnomat.append(precID)

                subprecfile.write(">"+precAllID+"\n"+str(precseq.strip())+"\n")
                subprecfile.close()

        if countnomat>0 and len(matIDs)>0:
            nomats=countnomat
            flagnomatexists=True
            log.debug("here is True")
            for i in matIDs:
                mtf = openfile(matfile)
                for record in SeqIO.parse(mtf, 'fasta'):
                    if i.strip() in record.description:
                        with open(directory+'tempmat.fa','a') as tempmaturefile:
                            tempmaturefile.write(">"+record.description+"\n"+str(record.seq)+"\n")

            for prec in SeqIO.parse(openfile(directory+'nomat-'+precfilename+'.fa'),'fasta'):
                listofmat=[]
                listofscore=[]
                splitprecdisc=prec.description.split()#because the clustal reads only id, so it will only write in thr result the id without MI00....
                prectempid=splitprecdisc[0]+"-"+splitprecdisc[1]
                listnomat.append((splitprecdisc[1]).strip())
                for mat in SeqIO.parse(openfile(directory+'tempmat.fa'),'fasta'):
                    splitmatdisc=mat.description.split()
                    mattempid=splitmatdisc[0]+"-"+splitmatdisc[1]#because the clustal reads only id, so it will only write in thr result the id without MIMAT...
                    with open(directory+'temptoalign.fa','w') as temptoalign:
                        temptoalign.write(">"+prectempid+"\n"+str(prec.seq)+"\n"+">"+mattempid+"\n"+str(mat.seq)+"\n")
                    infile=directory+'temptoalign.fa'
                    outfile=directory+'temptoalign.aln'
                    cline = ClustalwCommandline("clustalw", infile=infile, outfile=outfile)

                    stdout,stder=cline()
                    with open(directory+'temptoalign.txt','w') as fi:
                        fi.write(stdout)

                    for line in openfile(directory+'temptoalign.txt'):
                        if 'Sequence 2:' in line:
                            lineitem=line.split()
                            lineitem1=lineitem[2].split('-')
                            tempmatID=lineitem1[-1]
                            listofmat.append(tempmatID)

                        if 'Sequences (1:2)' in line:
                            lineitem=line.split()
                            score=float(lineitem[4])
                            listofscore.append(score)
                maxscore=max(listofscore)
                indexofbest=listofscore.index(maxscore)
                bestMatID=listofmat[indexofbest]
                prectempsplit=prec.description.split(" ")
                log.debug(["best for ",prectempsplit[1]," is ",bestMatID, "of ", listofscore, listofmat])

                if maxscore>=21:
                    for mattopredict in SeqIO.parse(openfile(directory+'tempmat.fa'),'fasta'):
                        if bestMatID in mattopredict.description:
                            matpredictid=mattopredict.id#the first part of the ID which be read and mentioned in the alignment file
                            matpredictID=str(mattopredict.description)#the full ID (with description to write ti full in the fasta)
                            matpredictseq=str(mattopredict.seq)
                            break
                    temptopredictfa=open(directory+'temptopredict.fa','w')
                    tempinfilepredict=directory+'temptopredict.fa'
                    tempoutfilepredict=directory+'temptopredict.aln'
                    temptopredictfa.write(">"+matpredictID+"\n"+matpredictseq+"\n"+">"+prec.description+"\n"+str(prec.seq))
                    temptopredictfa.close()
                    predictcline = ClustalwCommandline("clustalw", infile=tempinfilepredict, outfile=tempoutfilepredict)
                    predictcline()

                    matId=matpredictid
                    prectempsplit=prec.description.split()
                    mattempsplit=matpredictID.split()
                    newmatID=prectempsplit[0]+"-mat "+mattempsplit[1]+"/"+prectempsplit[1]+" "+prectempsplit[2]+" "+prectempsplit[3]
                    filename=getfilename(precfile)

                    listremovedbroken,listremovedscore,listremovedN=predict(tempoutfilepredict,matId,newmatID,matfile,filename,prec.description,mapfile,directory,listremovedbroken,listremovedscore,listremovedN)
                    stktemptopredict=tempoutfilepredict+".stk"
                    os.remove(str(stktemptopredict))

                elif maxscore<21:
                    prectempsplit=prec.description.split()
                    log.debug(["if not score",prectempsplit[1]])
                    if prectempsplit[1] not in listremovedscore:
                        listremovedscore.append(prectempsplit[1].strip())

        elif countnomat>0 and len(matIDs)==0:
            flagnomatexists=True
            nomats=-1
            return flagnomatexists,nomats,listremovedbroken,listremovedscore,listremovedN,listnomat
        elif countnomat<=0:
            log.debug("do the normal procedure")
            #flip here
            nomats=0
            flagnomatexists=False
            log.debug("all precursors has mature, at least one")
        return flagnomatexists,nomats,listremovedbroken,listremovedscore,listremovedN,listnomat
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))
