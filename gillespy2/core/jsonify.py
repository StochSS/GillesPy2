class Jsonify:
    def to_json(self):
        return self.public_vars()

    def from_json(self, json_object):
        self.__dict__ = json_object

    def public_vars(self):
        return { k:v for k, v in vars(self).items() if not k.startswith("_") }

    def encode_dict(self, dict):
        return list(map(lambda x: x.to_json(), dict.values()))