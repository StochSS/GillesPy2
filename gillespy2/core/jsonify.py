class Jsonify:
    def to_json(self):
        return self.public_vars()

    def from_json(json_object):
        pass

    def public_vars(self):
        return { k:v for k, v in vars(self).items() if not k.startswith("_") }

    def encode_dict(self, dict):
        return list(map(lambda x: x.to_json(), dict.values()))

from json import JSONEncoder
class ComplexJsonEncoder(JSONEncoder):
    def default(self, o):
        print(type(o))

        if isinstance(o, Jsonify):
            model = o.to_json()
            model["_type"] = str(type(o))
            return model

        return "unknown"