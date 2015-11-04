BWA_PATH="bwa"
SAMTOOLS_PATH="samtools"
REFINER_PATH="./TERefiner_1"
JELLYFISH_PATH=""
VELVET_PATH=""
THREADS=15
OUTPUT_FOLDER="./REPdenovo_Output/"
VERBOSE=True

def setParameters(pbwa,psamtools,prefiner,pjellyfish,pvelvet,t,pout,vbs):
    global  BWA_PATH
    BWA_PATH=pbwa
    global  SAMTOOLS_PATH
    SAMTOOLS_PATH=psamtools
    global REFINER_PATH
    REFINER_PATH=prefiner
    global  JELLYFISH_PATH
    JELLYFISH_PATH=pjellyfish
    global VELVET_PATH
    VELVET_PATH=pvelvet
    global THREADS
    THREADS=t
    global OUTPUT_FOLDER
    OUTPUT_FOLDER=pout
    global VERBOSE
    VERBOSE=vbs

def getBWAPath():
    global BWA_PATH
    return BWA_PATH

def getSamtoolsPath():
    global SAMTOOLS_PATH
    return SAMTOOLS_PATH

def getRefinerPath():
    global REFINER_PATH
    return REFINER_PATH

def getJellyfishPath():
    global JELLYFISH_PATH
    return JELLYFISH_PATH

def getVelvetPath():
    global VELVET_PATH
    return VELVET_PATH

def getThreadsNum():
    global THREADS
    return THREADS

def getOutputFolder():
    global OUTPUT_FOLDER
    return OUTPUT_FOLDER

def getVerbose():
    global VERBOSE
    return VERBOSE

def printCommand(cmd):
    if VERBOSE!=0:
        print "Running command: "+cmd+" ..."