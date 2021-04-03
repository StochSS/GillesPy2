import json, hashlib
from json import JSONEncoder


class Jsonify:
    """
    Interface to allow for instances of arbitrary types to be encoded into json strings and decoded into new objects.
    """

    def to_json(self, translation_table=None):
        import json

        encoder = ComplexJsonEncoder(translation_table)
        return json.dumps(self, indent=4, sort_keys=True, default=encoder.default)

    @classmethod
    def from_json(cls, json_object):
        """
        Convert some json_object into a decoded Python type. This function should return a __new__ instance of the type.

        :param json_object: A json str to be converted into a new type instance.
        """

        # If the json_object is actually a dict, it means we've decoded as much as possible.
        if type(json_object) is dict:
            new = cls.__new__(cls)
            new.__dict__ = json_object

            return new

        return json.loads(json_object, object_hook=ComplexJsonDecoder.decode_hook)

    def to_dict(self):
        """
        Convert the object into a dictionary ready for json encoding.
        Note: Complex types that inherit from Jsonify do not need to be manually encoded.

        By default, this function will return a dictionary of the object's public types.

        :param self: Instance of the object to convert into a dict.
        """
        return self.public_vars()

    def get_translation_table(self):
        """
        Generate a translation table that describes key:value pairs to convert user-defined data into generic equivalents.
        """
        pass

    def public_vars(self):
        """
        Gets a dictionary of public vars that exist on self. Keys starting with '_' are ignored.
        """
        return {k: v for k, v in vars(self).items() if not k.startswith("_")}

    def get_json_hash(self, translation_table):
        """
        Get the hash of the anonymous json representation of self.
        """

        json_str = self.to_json(translation_table)
        return hashlib.md5(str.encode(json_str)).hexdigest()

class ComplexJsonEncoder(JSONEncoder):
    def __init__(self, key_table=None, **kwargs):
        super(ComplexJsonEncoder, self).__init__(**kwargs)
        self.key_table = key_table

    def default(self, o):
        from numpy import ndarray
        from gillespy2.core.model import Model

        if isinstance(o, ndarray):
            return NdArrayCoder.to_dict(o)

        if not isinstance(o, Jsonify):
            return super().default(o)

        model = o.to_dict()

        # If the model is some subclass of gillespy2.core.model.Model, then manually set its type.
        if issubclass(o.__class__, Model):
            model["_type"] = f"{Model.__module__}.{Model.__name__}"

        else:
            model["_type"] = f"{o.__class__.__module__}.{o.__class__.__name__}"

        # If valid, recursively translate keys and values in the current model.
        self.recursive_translate(model)
        return model

    def recursive_translate(self, obj):
        from collections import OrderedDict, Hashable

        if obj is None or self.key_table is None:
            return

        # If the input object is a list, we iterate through it element by element.
        if isinstance(obj, list):
            for i, item in enumerate(list(obj)):
                if isinstance(item, dict) or isinstance(item, list):
                    self.recursive_translate(item)
                    continue

                if not isinstance(item, Hashable):
                    continue

                if item in self.key_table:
                    obj[i] = self.key_table[item]

            return

        # Else, the item is a dictionary, so we iterate through each key/value.
        for k in list(obj.keys()):
            if k in self.key_table:
                obj[self.key_table[k]] = obj.pop(k)
                k = self.key_table[k]

            # If the value is a list, we need to iterate through it.
            if isinstance(obj[k], list):
                self.recursive_translate(obj[k])
                continue

            # OrderedDicts are immutable, so we need to convert it into a dictionary prior to translation.
            if isinstance(obj[k], OrderedDict):
                obj[k] = dict(obj[k])

            # We need to translate all sub-elements in a dictionary, so recurse into it.
            if isinstance(obj[k], dict):
                self.recursive_translate(obj[k])
                continue

            # If the value isn't Hashable, continue.
            if not isinstance(obj[k], Hashable):
                continue

            v = obj[k]
            if v in self.key_table:
                obj[k] = self.key_table[v]


class ComplexJsonDecoder:
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


class NdArrayCoder(Jsonify):
    @staticmethod
    def to_dict(self):
        return {
            "start": self[0],
            "end": self[-1],
            "size": self.size,
            "_type": f"{NdArrayCoder.__module__}.{NdArrayCoder.__name__}"
        }

    @staticmethod
    def from_json(json_object):
        import numpy
        return numpy.linspace(start=json_object["start"], stop=json_object["end"], num=json_object["size"])