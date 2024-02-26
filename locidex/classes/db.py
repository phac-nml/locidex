import json
import os
import pathlib


class db_config:
    config = {}
    input_file = None
    required_fields = []
    status = True
    messages = []

    def __init__(self,input_file,required_fields):
        self.input_file = input_file
        if not os.path.isfile(self.input_file):
            self.messages.append(f'Error {self.input_file} input file does not exist')
            self.status = False
        self.parse_file()
        self.required_fields = required_fields
        (self.status,k) = self.validate_dict_keys(self.config, self.required_fields)
        if k != '':
            self.messages.append(f'Error required field missing from config file key:{k}')
            self.status = False

        return

    def parse_file(self):
        with open(self.input_file) as file:
            self.config = json.load(file)


    def validate_dict_keys(self,data_dict, required_keys):
        for k in required_keys:
            if not k in data_dict:
                return [False, k]
        return [True, '']


class search_db_conf:
    input_dir = None
    files = {
        'dir':[],
        'file':[]
    }
    config_file_path = ''
    config_obj = None
    meta_file_path = ''
    blast_paths = {

    }

    status = True
    messages = []

    def __init__(self, input_dir,db_basenames,config_required_fields,):
        self.input_dir = input_dir
        if not os.path.isdir(self.input_dir):
            self.messages.append(f'Error {self.input_dir} input file does not exist')
            self.status = False
            return
        self.get_dir_files()

        p = self.get_file_path_by_name(db_basenames['config'])

        if p == '':
            conf_name = db_basenames['config']
            self.messages.append(f'Error missing expected config file {conf_name} in db {input_dir}')
            self.status = False
            return

        self.config_file_path = p
        self.config_required_fields = config_required_fields
        self.config_obj = db_config(self.config_file_path,self.config_required_fields)
        self.status = self.config_obj.status
        self.messages = self.messages + self.config_obj.messages

        p = self.get_file_path_by_name(db_basenames['meta'])
        if p == '':
            conf_name = db_basenames['meta']
            self.messages.append(f'Error missing expected meta file {conf_name} in db {input_dir}')
            self.status = False
            return

        self.meta_file_path = p

        p = self.get_file_path_by_name('blast','dir')

        self.blast_paths['nucleotide'] = os.path.join(p,"nucleotide/nucleotide")
        self.blast_paths['protein'] = os.path.join(p, "protein/protein")


        if len(self.blast_paths) == 0:
            self.messages.append(f'Error could not find any expected blast dbs in path: {input_dir}, expected {db_basenames}')
            self.status = False
            return


    def get_dir_files(self):
        d = pathlib.Path(self.input_dir).rglob('*')
        for item in d:
            type = 'file'
            if item.is_dir():
                type = 'dir'
            self.files[type].append([f"{item.resolve()}",os.path.basename(item)])

    def get_file_path_by_name(self,filename,t='file'):
        for f in self.files[t]:
            if filename == f[1]:
                return f[0]
        return ''





