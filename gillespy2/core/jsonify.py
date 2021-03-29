from json import JSONEncoder, JSONDecoder

class Jsonify:
    def to_dict(self):
        return self.public_vars()

    def from_json(json_object):
        pass

    def public_vars(self):
        return {k: v for k, v in vars(self).items() if not k.startswith("_")}

    def encode_dict(self, dict):
        return list(map(lambda x: x.to_json(), dict.values()))

class ComplexJsonEncoder(JSONEncoder):
    def default(self, o):
        from numpy import ndarray
        from gillespy2.core.model import Model

        if isinstance(o, ndarray):
            return {
                "start": o[0],
                "end": o[-1],
                "points": o.size
            }

        if not isinstance(o, Jsonify):
            return super().default(o)

        model = o.to_dict()
        if issubclass(o.__class__, Model):
            model["_type"] = f"{Model.__module__}.{Model.__name__}"

        else:
            model["_type"] = f"{o.__class__.__module__}.{o.__class__.__name__}"

        return model

class ComplexTypeDecoder():
    @staticmethod
    def decode_hook(obj):
        from pydoc import locate

        if "_type" not in obj:
            return obj

        obj_type = locate(obj["_type"])

        if obj_type is None:
            raise Exception(f"{obj_type} does not exist.")

        if not issubclass(obj_type, Jsonify):
            raise Exception(f"{obj_type}")

        return obj_type.from_json(obj)