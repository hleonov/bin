def parse_args(options,argv):
    # manage input. This could also be done using pymacs, but
    # 1. the "fancy command line handling" is a bit to "talkative" for my case
    # 2. it would be the only reason to require pymacs which might not always be
    #    available
    #
    # arguments:
    #    options: list of options to expect. List of lists for each option:
    #             ['name','type','default']
    
    
    # prepare parsed args and option-sets
    # define dictionary for parsed arguments
    parsed_args = dict()

    file_options = set()                   # file options expect a single file
    multi_file_options = set()             # multi file options expect a list of files
    value_options = set()
    multi_value_options = set()
    flag_options = set()                   # flag options don't expect a value - they are true 
                                           # if set, false otherwise

    # sort options into the appropriate sets and initialise the parsed_args dictionary
    for option in options:
        if option[1]=='FILE':
            file_options.add(option[0])
            parsed_args[option[0]]=option[2]              # set default
        elif option[1]=='MULTI_FILE':
            multi_file_options.add(option[0])
            parsed_args[option[0]]=list()                 # default for lists:
            if type(option[2])==list:
                parsed_args[option[0]].extend(option[2])  # if default was list
            else:
                parsed_args[option[0]].append(option[2])  # if default was a value
        elif option[1]=='VALUE':
            value_options.add(option[0])
            parsed_args[option[0]]=option[2]
        elif option[1]=='MULTI_VALUE':
            multi_value_options.add(option[0])
            parsed_args[option[0]]=list()
            if type(option[2])==list:
                parsed_args[option[0]].extend(option[2])
            else:
                parsed_args[option[0]].append(option[2])
        elif option[1]=='FLAG':
            flag_options.add(option[0])
            parsed_args[option[0]]=option[2]
        else:
            print "Error - unknown option type: <%s>" % option
            exit(1)
    
    
    # define union sets of options
    single = file_options.union(value_options)
    multi = multi_file_options.union(multi_value_options)
    all_options = single.union(multi.union(flag_options))
    
    # the whole process is a finite automaton. This variable contains it's state
    mode = ''
    
    # parse all arguments but the first (which is the name of the script)
    for a in argv[1:]:
        # case 1: the argument is a defined option
        if a in all_options:
            if a in flag_options:
                # a flag option is just switched true. the mode is set to nothing
                parsed_args[a]=True
                mode = ''
            elif a in single:
                # for a single-valued option, this one is initialized and the mode 
                # is set
                parsed_args[a]=''
                mode = a
            else:
                # for a multi-valued option, the parsed result is a list which
                # is initialized
                parsed_args[a]=[]
                mode = a
        # case 2: the argument is not defined, hence a value
        else:
            # a single mode is saved, the mode is set to '' to denote this option
            # is set
            if mode in single:
                parsed_args[mode]=a
                mode = ''
            # multi-value lists are not reset as they can go on
            elif mode in multi:
                parsed_args[mode].append(a)
            # if the mode has been reset (or not set at all), the value can not
            # be assigned to an option => error.
            else:
                print "Error: unknown command line option: " + a
                exit(1)
    
    # TODO: add defaults, more complex error handling (empty values)
    return parsed_args

if __name__ == "__main__":
    import sys
    pa = parse_args([['-file1','FILE','file1.pdb'],['-file2','FILE','file2.pdb'],['-mfile1','MULTI_FILE','mfile1.pdb'],['-mfile2','MULTI_FILE','mfile2.pdb'],['-mvalue1','MULTI_VALUE',[]],['-flag1','FLAG',True],['-flag2','FLAG',False]],sys.argv)
    print pa
    
