import pandas as pd
import sqlalchemy as sa
import pandas as pd
import numpy as np

def mysql_pull(mysql_usr:str,mysql_pswd:str,db,query,mysql_ip:str='127.0.0.1',mysql_port:int = 3306):

    # Set up the database connection
    database_username = mysql_usr
    database_password = mysql_pswd
    database_hostname = f'{mysql_ip}:{mysql_port}'
    database_name = db

    # Create the connection string
    connection_string = f'mysql+pymysql://{database_username}:{database_password}@{database_hostname}/{database_name}'

    # Create the engine
    engine = sa.create_engine(connection_string)

    dataDF = pd.read_sql_query(query,engine)

    return dataDF

def author_json_extract(json_path:str,target_id:str,mysql_pass:str,db:str):
    df = pd.read_json(json_path, orient='index')
    df = df.transpose()
    df = pd.concat([df.training,df.test,df.validation])
    df = df.dropna()
    df = pd.DataFrame(df.to_list())
    df = df.rename(columns={0:'comp_id',1:target_id})
    comp_id_tuple = tuple(df.comp_id.to_list())
    df_smiles = get_smiles_chembl(comp_id_tuple,mysql_pass=mysql_pass,db=db)
    df = pd.merge(df,df_smiles,how='left',left_on='comp_id',right_on='comp_id')
    return df

def get_sequence_chembl(comp_id_tuple:tuple, mysql_pass:str, db:str = 'chembl_27'):
    '''given a tuple with comp_id it returns a pandas dataframe with all the compounds chembl id and its canonical smiles'''
    query_smiles = f'''
    SELECT trgd.chembl_id AS target_chemb_id,
    bcs.sequence AS target_sequence



SELECT md.chembl_id, cmpstc.canonical_smiles, cmpstc.standard_inchi, cmpstc.standard_inchi_key
FROM compound_structures AS cmpstc 
JOIN molecule_dictionary AS md ON md.molregno = cmpstc.molregno
WHERE md.chembl_id IN {tuple(comp_id_tuple)}'''
    df_smiles = mysql_pull('sjinich',mysql_pass,db,query=query_smiles)
    df_smiles = df_smiles.rename(columns={'chembl_id':'comp_id','canonical_smiles':'smiles'})
    return df_smiles

authors_zip_file = 'target_training_datasets.zip'

import zipfile
from io import BytesIO
from os import path

import getpass
pass_mysql = getpass.getpass('mysql_pass:')

folder = 'data_autores_inchi/'

with zipfile.ZipFile(authors_zip_file,mode='r') as archive:
    files = archive.namelist()
    del files[0]
    for file in files:
        comp_id = file[:file.find('.zip')]
        filedata = BytesIO(archive.read(file))
        with zipfile.ZipFile(filedata) as archive_in:
            filedata_json = BytesIO(archive_in.read(path.join(comp_id,'train_val_test_dict.json')))
            df = author_json_extract(filedata_json,comp_id,pass_mysql,'chembl_27')
            print(f'{file} extracted')
            print(f'missing smiles: {len(df[df.smiles.isna()])}')
            df.to_pickle(path.join(folder,comp_id+'.pickle'))