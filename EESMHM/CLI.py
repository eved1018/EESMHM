import argparse
import sys
import os
import urllib.request

def cli():

    files = [i for i in os.listdir() if i.endswith(".pdb")]
    parser = argparse.ArgumentParser()
    parser.add_argument('-pdb',
                        help=f"please choose an input file from {files} or add to input folder", default=None)
    parser.add_argument('-qc', '--query_chain', help="please choose a query chain", default=None)
    parser.add_argument('-ic', '--interacting_chain',
                        help="please choose the interacting chain", default=None)
    parser.add_argument('-config', '--config_file',
                        help="config file containg mutation data", default="config.txt")

    parser.add_argument('-foldx', '--foldx', default=False, help="foldx path")
    parser.add_argument('-evoef', '--evoef', default=False, help="evoef path")
    args = parser.parse_args()
    pdb = args.pdb
    qc = args.query_chain
    ic = args.interacting_chain
    config  = args.config_file
    foldx = args.foldx
    evoef = args.evoef
    pdb, pdb_file = pdbManager(pdb, files)
    if qc is None or ic is None:
        qc, ic, pdb_file = getChains(pdb_file)
    return pdb, pdb_file, qc, ic, config, foldx, evoef


def pdbManager(pdb, files):
    if pdb is None:
        pdb = input("Please input a pdb id (ex. 1i8l): ")
        pdb = pdb.lower()
    if not pdb.endswith(".pdb"):
        pdb_file = pdb + ".pdb"
    else:
        pdb_file = pdb
        pdb = pdb.replace(".pdb", "")
    if pdb_file not in files:
        download = input(f"Download {pdb} from the RCSB? [y,n]: ")
        if download == "y":
            download_pdb(pdb, "")
        else:
            print(f"{pdb} not downloaded please add {pdb} to input folder")
            sys.exit()
    return pdb, pdb_file


def download_pdb(pdbcode, datadir, downloadurl="https://files.rcsb.org/download/"):
    """
    Downloads a PDB file from the Internet and saves it in a data directory.
    :param pdbcode: The standard PDB ID e.g. '3ICB' or '3icb'
    :param datadir: The directory where the downloaded file will be saved
    :param downloadurl: The base PDB download URL, cf.
        `https://www.rcsb.org/pages/download/http#structures` for details
    :return: the full path to the downloaded PDB file or None if something went wrong
    """
    pdbfn = pdbcode + ".pdb"
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)
    try:
        urllib.request.urlretrieve(url, outfnm)
        return outfnm
    except Exception as err:
        print(str(err), file=sys.stderr)
        return None


def getChains(pdb_file):
    chains = set()
    has_icode = False
    with open(f"{pdb_file}", "r") as pdb_fh:
        for line in pdb_fh:
            if line.startswith("ATOM"):
                chains.add(line[21])
                if line[26] != ' ':
                    has_icode = True
    qc = input(f"Please select query [{chains}]: ")
    if qc not in chains:
        print("chain not found")
        sys.exit()

    chains.remove(qc)
    ic = input(f"Please select interacting chain [{chains}]: ")
    if ic not in chains:
        print("chain not found")
        sys.exit()
    if has_icode:
        pdb_file = fixInsert(pdb_file)
    return qc, ic ,pdb_file

def run(fhandle, option_list):
    """
    Delete insertion codes (at specific residues).
    By default, removes ALL insertion codes on ALL residues. Also bumps
    the residue numbering of residues downstream of each insertion.
    This function is a generator.
    Parameters
    ----------
    fhandle : a line-by-line iterator of the original PDB file.
    option_list : list
        List of insertion options to act on.
        Example ["A9", "B12"]. An empty list performs the default
        action.
    Yields
    ------
    str (line-by-line)
        The modified (or not) PDB line.
    """

    option_set = set(option_list)  # empty if option_list is empty

    # Keep track of residue numbering
    # Keep track of residues read (chain, resname, resid)
    offset = 0
    prev_resi = None
    seen_ids = set()
    clean_icode = False
    records = ('ATOM', 'HETATM', 'ANISOU', 'TER')
    for line in fhandle:
        if line.startswith(records):
            res_uid = line[17:27]  # resname, chain, resid, icode
            id_res = line[21] + line[22:26].strip()  # A99, B12
            has_icode = line[26].strip()  # ignore ' ' here

            # unfortunately, this is messy but not all PDB files follow a nice
            # order of ' ', 'A', 'B', ... when it comes to insertion codes..
            if prev_resi != res_uid:  # new residue
                # Does it have an insertion code
                # OR have we seen this chain + resid combination before?
                # #2 to catch insertions WITHOUT icode ('A' ... ' ' ... 'B')
                if (has_icode or id_res in seen_ids):
                    # Do something about it
                    # if the user provided options and this residue is in them
                    # OR if the user did not provide options
                    if (not option_set) or (id_res in option_set):
                        clean_icode = True
                    else:
                        clean_icode = False
                else:
                    clean_icode = False

                prev_resi = res_uid

                if id_res in seen_ids:  # offset only if we have seen this res.
                    offset += 1

            if clean_icode:  # remove icode
                line = line[:26] + ' ' + line[27:]

            # Modify resid if necessary
            resid = int(line[22:26]) + offset
            line = line[:22] + str(resid).rjust(4) + line[26:]
            seen_ids.add(id_res)

            # Reset offset on TER
            if line.startswith('TER'):
                offset = 0
        yield line


def fixInsert(pdb_file):
    has_icode = False
    with open(f"{pdb_file}", "r") as pdb_fh:
        for line in pdb_fh:
            if line.startswith("ATOM"):
                if line[26] != ' ':
                    has_icode = True

    if has_icode:
        # Check Input
        pdbfh = open(f"{pdb_file}", 'r')
        option_list = []
        # Do the job
        new_pdb = run(pdbfh, option_list)

        with open(f"fixed_{pdb_file}", "w+") as _file:
            for line in new_pdb:
                _file.write(line)
        pdbfh.close()
        return f"fixed_{pdb_file}"
    else:
        return pdb_file




