import json,os

class metadata:
    input_file = None
    required_fields = []
    data = {}
    status = True
    messages = []

    def __init__(self,input_file):
        self.input_file = input_file
        if not os.path.isfile(self.input_file):
            self.messages.append(f'Error {self.input_file} input file does not exist')
            self.status = False
            return self.status
        self.data = self.parse_file()
        return self.status

    def parse_file(self):
        with open(self.input_file) as file:
            self.config = json.load(file)

    def validate_dict_keys(self,data_dict, required_keys):
        for k in required_keys:
            if not k in data_dict:
                return [False, k]
        return [True, k]
