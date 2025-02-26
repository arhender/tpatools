from pathlib import Path
import sys

def eV_to_nm(inputarray):
    return 4.135667516E-15 * 2.9979E8 * 1E9 / inputarray

def filepath_searcher(suppliedpath):
    """
    For integration in the command-line scripts.
    This function processes the filepath input supplied to the script (if any).

    If no filepath provided, the current directory is searched for any matching a list of fallback/default 2PA output names 

    If a path to an existing file is given, simply returns said filepath. 
    If invalid, exits and requests a restart

    If a path to a directory is given instead of a file, the given directory is searched for files matching the the names provided in the fallback list

    """
    fallback_names = ['escf.out', 'bse.out', 'tpa.out', 'tddft.out', 'td-dft.out']

    filepath = None
    if suppliedpath is None:
        for name in fallback_names:
            if Path(name).is_file():
                filepath = Path(name)
                break
        
        if filepath is None: 
            print(f'No output filename has been provided. As a fallback, the program has checked and found no files from the following default namelist in the directory:\n\n{"\n".join(fallback_names)}\n\nPlease retry and specify the name of your output file' )
            sys.exit()
    else:
        if Path(suppliedpath).is_file():
            filepath = Path(suppliedpath)
        else:
            for name in fallback_names:
                if (Path(suppliedpath) / name).is_file():
                    filepath = Path(suppliedpath) / name
                    break

            if filepath is None:
                print('The provided filepath does not exist, please check your input and try again')
                sys.exit()

    return filepath



