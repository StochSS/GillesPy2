from json import JSONEncoder, JSONDecoder

class Jsonify():
    def to_dict(self):
        return self.public_vars()

    def from_json(json_object):
        pass

    def public_vars(self):
        return {k: v for k, v in vars(self).items() if not k.startswith("_")}

    def encode_dict(self, dict):
        return list(map(lambda x: x.to_json(), dict.values()))

class StaticJsonify():
    @staticmethod
    def to_dict(obj):
        pass

    @staticmethod
    def from_json(json_object):
        pass

    @staticmethod
    def public_vars(obj):
        return {k: v for k, v in vars(obj).items() if not k.startswith("_")}

class ComplexJsonEncoder(JSONEncoder):
    def default(self, o):
        from numpy import ndarray
        from gillespy2.core.model import Model

        if isinstance(o, ndarray):
            return StaticNdarrayCoder.to_dict(o)

        if not isinstance(o, Jsonify):
            return super().default(o)

        model = o.to_dict()
        if issubclass(o.__class__, Model):
            model["_type"] = f"{Model.__module__}.{Model.__name__}"

        else:
            model["_type"] = f"{o.__class__.__module__}.{o.__class__.__name__}"

        return model

class ComplexJsonDecoder():
    @staticmethod
    def decode_hook(obj):
        from pydoc import locate

        if "_type" not in obj:
            return obj

        obj_type = locate(obj["_type"])

        if obj_type is None:
            raise Exception(f"{obj_type} does not exist.")

        if not issubclass(obj_type, Jsonify) and not issubclass(obj_type, StaticJsonify):
            raise Exception(f"{obj_type}")

        return obj_type.from_json(obj)

class StaticNdarrayCoder(StaticJsonify):
    @staticmethod
    def to_dict(obj):
        return {
            "start": obj[0],
            "end": obj[-1],
            "size": obj.size,
            "_type": f"{StaticNdarrayCoder.__module__}.{StaticNdarrayCoder.__name__}"
        }

    @staticmethod
    def from_json(json_object):
        import numpy
        return numpy.linspace(start=json_object["start"], stop=json_object["end"], num=json_object["size"])
