import os, hashlib, subprocess, requests
import tqdm
import argparse
import xml.etree.ElementTree as ET
import pandas as pd
import collections, logging, time

STUDY_ACCESSION_KEY = 'study_accession'
RUN_ACCESSION_KEY = 'run'
NCBI_API_KEY_ENV = 'NCBI_API_KEY'
BIOPROJECT_ACCESSION_KEY = 'bioproject'

md5_cache_dir = os.path.join(os.path.expanduser('~'), '.ena_fastq_validator')

def add_api_key(other_params):
        if NCBI_API_KEY_ENV in os.environ:
            other_params['api_key'] = os.environ[NCBI_API_KEY_ENV]
        return other_params

def _retry_request(description, func):
    '''Retry a reqests.post or requests.get 3 times, returning the request
    when OK, otherwise raising an Exception'''

    num_retries = 3
    sleep_time = 60
    def retrying(i, num_retries=3):
        if i < num_retries-1:
            logging.warning("Retrying request (retry {} of {})".format(i+1, num_retries-1))
    
    for i in range(num_retries):
        try:
            this_res = func()
            if not this_res.ok:
                logging.warning("Request not OK when {}: {}: {}".format(description, this_res, this_res.text))
                logging.warning("Sleeping for {} seconds before retrying".format(sleep_time))
                time.sleep(60)
                retrying(i)
            else:
                return this_res
        except Exception as e:
            logging.warning("Exception raised when {}: {}".format(description, e))
            logging.warning("Sleeping for {} seconds before retrying".format(sleep_time))
            time.sleep(60)
            retrying(i)
    raise Exception("Failed to {} after {} attempts".format(description, num_retries))

def efetch_accession_from_ids(webenv, accessions, num_ids):
    print("Fetching accession using SRA IDs from NCBI...")
    data_frames = []
    retmax = num_ids + 10
    res = _retry_request(
        'efetch_from_ids',
        lambda: requests.get(
            url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            params = add_api_key({
                "db": "sra",
                "tool": "Budding",
                "email": "lirui@big.ac.cn",
                "webenv": webenv,
                "query_key": 1
                }),
            ))
    if not res.ok:
        raise Exception("HTTP Failure when requesting efetch from IDs: {}: {}".format(res, res.text))
    root = ET.fromstring(res.text)
    print("Got accession using SRA IDs!")

    def try_get(func):
        try:
            return func()
        except AttributeError:
            return ''
        except KeyError:
            return None

    if root.find("ERROR") is not None:
        raise Exception("Error when fetching metadata: {}".format(root.find("ERROR").text)) # type: ignore

    # Some samples such as SAMN13241871 are linked to multiple runs e.g. SRR10489833
    accessions_set = None if accessions is None else set(accessions)

    for pkg in root.findall('EXPERIMENT_PACKAGE'):
        d = collections.OrderedDict()
        for run in pkg.findall('./RUN_SET/RUN'):
            accession_here = run.attrib['accession']
            if accessions_set is None or accession_here in accessions_set:
                d['run'] = try_get(lambda: run.attrib['accession'])
                data_frames.append(d)

    return pd.DataFrame(data_frames)

def fetch_runs_from_bioproject(bioproject_accession):
    retmax = 10000
    res = requests.get(
        url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        params = {
            "db": "sra",
            "term": f"{bioproject_accession}[BioProject]",
            "tool": "Budding",
            "email": "lirui@big.ac.cn",
            "retmax": retmax,
            "usehistory": "y",
            },
        )
    if not res.ok:
        raise Exception("HTTP Failure when requesting search from bioproject: {}: {}".format(res, res.text))
    root = ET.fromstring(res.text)
    print(f"Got SRA IDs using BioProject: {bioproject_accession}")
    
    # sra_ids = list([c.text for c in root.find('IdList')])
    id_list_element = root.find('IdList')
    if id_list_element is not None:
        sra_ids = list([c.text for c in id_list_element])
    else:
        raise Exception("No IdList element found in response: {}".format(res.text))
    
    if len(sra_ids) == retmax:
        print("Unexpectedly found the maximum number of results for this query, possibly some results will be missing")
    
    #webenv = root.find('WebEnv').text
    webenv_element = root.find('WebEnv')
    if webenv_element is not None:
        webenv = webenv_element.text
    else:
        raise Exception("No WebEnv element found in response: {}".format(res.text))

    # Now convert the IDs into runs
    metadata = efetch_accession_from_ids(webenv, None, len(sra_ids))
    return metadata['run'].to_list()

def check_md5sum(file_path):
        hash_md5 = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()

def get_md5_from_ena(accession, cache = True):
    md5_dic = {}
    print("Start generating md5 dictionary...")
    
    # This function get a full list of fastq accessions from NCBI, instead of using filenames
    fastq_lis = fetch_runs_from_bioproject(accession)
    if cache:
        print("Using cache to search md5...")
        full_fastq_set = set(fastq_lis)
        cache_fastq_lis = []
        with open(md5_cache_dir) as f:
            lines = f.readlines()
            for line in lines:
                cache_fastq_lis.append(line.split('\t')[0])
                md5_dic[line.split('\t')[0]] = line.split('\t')[1].split(';')
        cache_fastq_set = set(cache_fastq_lis)
        need_query_lis = list(full_fastq_set - cache_fastq_set)
        if len(need_query_lis) == 0:
            print("All runs are cached, no need to fetch md5 from ENA")
            return md5_dic
        else:
            print("Some runs are not cached, start fetching md5 from ENA, this may take a while...")
            fq_bar = tqdm.tqdm(need_query_lis)
            res_lis = [] # list of cached md5s
            for fastq in fq_bar:
                fastq_ass = fastq.split("_")[0]
                query_url = f"http://www.ebi.ac.uk/ena/portal/api/filereport?accession={fastq_ass}&result=read_run&fields=fastq_md5"
                query_p = subprocess.run(["wget", "-qO-", query_url], stdout = subprocess.PIPE, encoding = "utf-8")
                
                res_lis.append(str(query_p.stdout).split('\n')[1])
                res = str(query_p.stdout).split('\n')[1].split('\t')
                ass, md5s = res[0], res[1]
                md5_dic[ass] = md5s.split(';')
                
            with open(md5_cache_dir, 'a') as out:
                for md5 in res_lis:
                    out.write(md5 + '\n')
            
    else:       
        print("Start fetching md5 from ENA..., this may take a while")
        fq_bar = tqdm.tqdm(fastq_lis)
        
        res_lis = [] # list of cached md5s
        for fastq in fq_bar:
            fastq_ass = fastq.split("_")[0]
            query_url = f"http://www.ebi.ac.uk/ena/portal/api/filereport?accession={fastq_ass}&result=read_run&fields=fastq_md5"
            query_p = subprocess.run(["wget", "-qO-", query_url], stdout = subprocess.PIPE, encoding = "utf-8")
            
            res_lis.append(str(query_p.stdout).split('\n')[1])
            res = str(query_p.stdout).split('\n')[1].split('\t')
            ass, md5s = res[0], res[1]
            md5_dic[ass] = md5s.split(';')
    
        with open(md5_cache_dir, 'a') as out:
            for md5 in res_lis:
                out.write(md5 + '\n')
    
            
    print("md5 dictionary generated!")
    return md5_dic

def cache_optimize():
    # TODO : remove duplicate in .ena_fastq_validator
    pass
    
def report():
    # TODO: repeort md5 failed files verbosely
    pass


def validate(accession, dir):
    fastq_dir = dir
    fastq_lis = os.listdir(fastq_dir)
    fail_cnt = 0
    md5_dic = get_md5_from_ena(accession)

    print("Start calculating md5")
    for fastq in fastq_lis:
        ass = fastq.split("_")[0]
        dig = int(fastq.split(".")[0].split("_")[1])
        full_path = os.path.join(fastq_dir, fastq)
    
        actual_md5 = check_md5sum(full_path)
        expected_md5 = str(md5_dic[ass][dig - 1]).replace("\n", "")
        if actual_md5 == expected_md5:
            # print(f"{ass}_{dig} success!")
            pass
        else:
            print(f"{ass}_{dig} failed!")
            fail_cnt += 1
    print(f'md5 failed count: {fail_cnt}')
    if fail_cnt == 0:
        print("project download success!")
        
    


def main():
    parser = argparse.ArgumentParser(
    prog = 'Budding',
    description = 'Budding is an ENA fastq validator')
    parser.add_argument('-i', '--accession', type = str, help = 'bioproject accession', action = "store", nargs = '*')
    parser.add_argument('-d', '--dir', type = str, help = 'fastq directory', action = "store", nargs = '*')
    # parser.add_argument('-c', '--cache', type = bool, help = 'use cache', action = "store", nargs = '*')
    args = parser.parse_args()
    for study in zip(args.accession, args.dir):
        validate(study[0], study[1])
    
    # validate("PRJNA498125", "./sc_mouse/PRJNA498125/")

if __name__ == '__main__':
    main()