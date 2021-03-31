from json import JSONEncoder


class Jsonify:
    """
    Interface to allow for instances of arbitrary types to be encoded into json strings and decoded into new objects.
    """

    def to_dict(self):
        """
        Convert the object into a dictionary ready for json encoding.
        Note: Complex types that inherit from Jsonify do not need to be manually encoded.

        By default, this function will return a dictionary of the object's public types.

        :param self: Instance of the object to convert into a dict.
        """
        return self.public_vars()

    @staticmethod
    def from_json(json_object):
        """
        Convert some json_object into a decoded Python type. This function should return a __new__ instance of the type.

        :param json_object: A json dict to be converted into a new type instance.
        """
        pass

    def public_vars(self):
        """
        Gets a dictionary of public vars that exist on self. Keys starting with '_' are ignored.
        """
        return {k: v for k, v in vars(self).items() if not k.startswith("_")}


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
        def recursive_translate(obj):
            for k in list(obj.keys()):
                from collections import OrderedDict, Hashable

                # OrderedDicts are immutable, so we need to convert it into a dictionary prior to translation.
                if isinstance(obj[k], OrderedDict):
                    obj[k] = dict(obj[k])

                # We need to translate all sub-elements in a dictionary, so recurse into it.
                if isinstance(obj[k], dict):
                    recursive_translate(obj[k])
                    continue

                # If the value isn't Hashable, continue.
                if not isinstance(obj[k], Hashable):
                    continue

                v = obj[k]
                if v in self.key_table:
                    obj[k] = self.key_table[v]

                if k in self.key_table:
                    obj[self.key_table[k]] = obj.pop(k)

        recursive_translate(model)

        return model

    def encode(self, o):
        print(type(o))
        return self.encode(o)


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