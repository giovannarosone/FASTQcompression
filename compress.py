#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path, shutil, struct

Version = "v0.1"


Description = """Tool to compresse FASQ files

For example
  {exe} INPUT.fastq
will produce the output files XXX and YYY
 
--------------------------
Command line options:
--------------------------
""".format(exe=sys.argv[0])

gsufsort_exe = "external/gsufsort/gsufsort"
header_split = "sed -n 1~4p"
qs_split = "sed -n 4~4p"
dna_split = "sed -n 2~4p"
gzip_exe = "gzip -9 -k -f"
zip7_exe = "7z a -mx9 -mmt12"

smooth_exe = "src/fastqcompression"

def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input', help='input file name', type=str, nargs='+')
    parser.add_argument('-o', '--out', help='output base name (def. input base name)', default="", type=str)  
    #parser.add_argument('--delete', help='delete output files',action='store_true')
    parser.add_argument('-1', '--step1', help='stop after step 1 (debug only)',action='store_true')
    parser.add_argument('-2', '--step2', help='stop after step 2 (debug only)',action='store_true')
    parser.add_argument('-3', '--step3', help='stop after step 3 (debug only)',action='store_true')
    parser.add_argument('-4', '--step4', help='stop after step 4 (debug only)',action='store_true')
    parser.add_argument('--original', help='do not call step 2',action='store_true')
    parser.add_argument('--all', help='run all competitors',action='store_true')
    parser.add_argument('-v',  help='verbose: extra info in the log file',action='store_true')
    args = parser.parse_args()
    # ---- check number of input files and define basename
    check_input(args)
    # ---- create and open log file
    logfile_name = args.basename + ".log"
    # get main directory
    args.dir = os.path.split(sys.argv[0])[0]
    print("Sending logging messages to file:", logfile_name)

    with open(logfile_name,"w") as logfile:
    
        print(">>> fastq-bwt version " + Version,file=logfile)
        show_command_line(logfile)
        logfile.flush()

        if len(args.out)==0 : args.out=args.input[0]

        # temporary files
        args.tmp = []

        # --- step1: compute BWT and LCP array
        start = time.time()
        if(step1(args, logfile, logfile_name)!=True):
            sys.exit(1)
        print("Elapsed time: {0:.4f}".format(time.time()-start))
        if args.step1:
            print("Exiting after step 1 as requested")
            return
    
        # --- step2: smooth BWT and QS sequences 
        start = time.time()
        if(step2(args, logfile, logfile_name)!=True):
            sys.exit(1)
        print("Elapsed time: {0:.4f}".format(time.time()-start))
        if args.step2:
            print("Exiting after step 2 as requested")
            return

        args.stream = []

        # --- step3: 
        start = time.time()
        if(step3(args, logfile, logfile_name)!=True):
            sys.exit(1)
        print("Elapsed time: {0:.4f}".format(time.time()-start))
        if args.step3:
            print("Exiting after step 3 as requested")
            return

        # --- step4: 
        start = time.time()
        if(step4(args, logfile, logfile_name)!=True):
            sys.exit(1)
        print("Elapsed time: {0:.4f}".format(time.time()-start))
        if args.step4:
            print("Exiting after step 3 as requested")
            return

        # compressed files
        args.output = []

        # --- step5: compress QS, reads and headers separatedly 
        start = time.time()
        if(step5(args, logfile, logfile_name)!=True):
            sys.exit(1)
        print("Elapsed time: {0:.4f}".format(time.time()-start))

        # ---- final report
        insize = os.path.getsize(args.input[0])

        print("=== results ==="); 
        print("Original:\t{0:.2f} MB".format(insize/(1024*1024)))
        outsize = 0
        for f in args.output:
            size = os.path.getsize(f)
            if(args.v): print(f, " - ", size)
            outsize += size
        print("Compressed:\t{0:.2f} MB".format(outsize/(1024*1024)))
        print("Ratio = {0:.2f}".format(outsize/insize))

        # --- extra: compress INPUT.fastq with default method 
        if not args.all:
            return

        #gzip
        start = time.time()
        if(gzip(args, logfile, logfile_name)!=True):
            sys.exit(1)
        outsize = os.path.getsize(args.input[0]+".gz")
        print("Compressed:\t{0:.2f} MB".format(outsize/(1024*1024)))
        print("Ratio = {0:.2f}".format(outsize/insize))

        print("Elapsed time: {0:.4f}".format(time.time()-start))

        #7zip
        start = time.time()
        if(zip7(args, logfile, logfile_name)!=True):
            sys.exit(1)
        outsize = os.path.getsize(args.input[0]+".7z")
        print("Compressed:\t{0:.2f} MB".format(outsize/(1024*1024)))
        print("Ratio = {0:.2f}".format(outsize/insize))

        print("Elapsed time: {0:.4f}".format(time.time()-start))


    return True

##

def step1(args, logfile, logfile_name):
    print("--- Step 1 ---", file=logfile); logfile.flush()
    exe = os.path.join(args.dir, gsufsort_exe)
    options = ""
    if len(args.out)>0 : options+="-o "+args.out
    else : options+="-o "+args.input[0]
    command = "{exe} {ifile} --bwt --lcp --qs {opt}".format(exe=exe, ifile=args.input[0], opt=options)
    print("=== gsufsort ==="); print(command)
    # tmp files
    args.tmp.append(args.basename+".bwt")
    args.tmp.append(args.basename+".bwt.qs")
    return execute_command(command, logfile, logfile_name)

def step2(args, logfile, logfile_name):
    if args.original:
        print("--- Step 2 ---", file=logfile); logfile.flush()
        command = "cp "+ args.input[0] +" "+args.out+".fq" 
        print(command)
        os.system(command)
    else:
        print("--- Step 2 ---", file=logfile); logfile.flush()
        exe = os.path.join(args.dir, smooth_exe)
        options = "-e " + args.tmp[0] + " -q " + args.tmp[1] + " -f " + args.input[0]+" -o "+args.out+".fq"
        command = "{exe} {ifile} {opt}".format(exe=exe, ifile=args.input[0], opt=options)
        print("=== smooth-qs ===")
        print(command)
        return execute_command(command, logfile, logfile_name)
    return True

##
def step3(args, logfile, logfile_name):
    print("--- Step 3 ---", file=logfile); logfile.flush()
    ##
    exe = header_split 
    ifile = args.out+".fq"
    ofile = args.out+".h"
    command = "{exe} {ifile} > {ofile}".format(exe=exe, ifile=ifile, ofile=ofile)
    print("=== header ===")
    print(command)
    os.system(command)
    args.stream.append(args.basename+".h")
    return True
##

def step4(args, logfile, logfile_name):
    print("--- Step 4 ---", file=logfile); logfile.flush()
    ##
    exe = dna_split 
    ifile = args.out+".fq"
    ofile = args.out+".dna"
    command = "{exe} {ifile} > {ofile}".format(exe=exe, ifile=ifile, ofile=ofile)
    print("=== streams ===")
    print(command)
    os.system(command)
    args.stream.append(args.basename+".dna")
    ##
    exe = qs_split 
    ofile = args.out+".qs"
    command = "{exe} {ifile} > {ofile}".format(exe=exe, ifile=ifile, ofile=ofile)
    print(command)
    args.stream.append(args.basename+".qs")
    logfile.flush()
    os.system(command)
    return True

def step5(args, logfile, logfile_name):
    print("--- Step 5 ---", file=logfile); logfile.flush()
    #exe = gzip_exe
    exe = zip7_exe
    print("=== compression ===")
    for f in args.stream:
        #ofile = f+".gz"
        #command = "{exe} {ifile}".format(exe=exe, ifile=f)
        ofile = f+".7z"
        command = "{exe} {ofile} {ifile}".format(exe=exe, ifile=f, ofile=ofile)
        print(command)
        execute_command(command, logfile, logfile_name)
        args.output.append(ofile)
    return True

def gzip(args, logfile, logfile_name):
    print("--- gzip ---", file=logfile); logfile.flush()
    exe = gzip_exe
    print("=== gzip ===")
    command = "{exe} {ifile}".format(exe=exe, ifile=args.input[0])
    print(command)
    return execute_command(command, logfile, logfile_name)

def zip7(args, logfile, logfile_name):
    print("--- 7z ---", file=logfile); logfile.flush()
    exe = zip7_exe
    ofile = args.input[0]+".7z"
    print("=== 7z ===")
    command = "{exe} {ofile} {ifile}".format(exe=exe, ifile=args.input[0], ofile=ofile)
    print(command)
    return execute_command(command, logfile, logfile_name)

########

# check correctness of number of input file and define basename for output
def check_input(args):
    if len(args.out)==0:       # specify basename for input files gap+merge
        args.basename = args.input[0]
    else:
        args.basename = args.out
    return True

# compute hash digest for a file 
def file_digest(name,logfile):
    try:
        hash_command = "{exe} {infile}".format(exe=shasum_exe, infile=name)
        hashsum = subprocess.check_output(hash_command.split(),stderr=logfile)
        hashsum = hashsum.decode("utf-8").split()[0]
    except:
        hashsum = "Error!" 
    return hashsum  

# execute command: return True is everything OK, False otherwise
def execute_command(command,logfile,logfile_name):
    try:
        subprocess.check_call(command.split(),stdout=logfile,stderr=logfile)
    except subprocess.CalledProcessError:
        print("Error executing command line:")
        print("\t"+ command)
        print("Check log file: " + logfile_name)
        return False
    return True

def show_command_line(f):
    f.write("Python command line: ") 
    for x in sys.argv:
        f.write(x+" ")
    f.write("\n")   

if __name__ == '__main__':
    main()
