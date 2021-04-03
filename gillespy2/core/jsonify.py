import json, hashlib
from json import JSONEncoder


class Jsonify:
    """
    Interface to allow for instances of arbitrary types to be encoded into json strings and decoded into new objects.
    """

    def to_json(self, translation_table=None):
        import json

        encoder = ComplexJsonCoder(translation_table)
        return json.dumps(self, indent=4, sort_keys=True, default=encoder.default)

    @classmethod
    def from_json(cls, json_object, translation_table=None):
        """
        Convert some json_object into a decoded Python type. This function should return a __new__ instance of the type.

        :param json_object: A json str to be converted into a new type instance.
        :param translation_table: A dictionary used to translate anonymous names back into user-defined.
        """

        # If the json_object is actually a dict, it means we've decoded as much as possible.
        if type(json_object) is dict:
            new = cls.__new__(cls)
            new.__dict__ = json_object

            return new

        decoder = ComplexJsonCoder(translation_table)
        return json.loads(json_object, object_hook=decoder.decode)


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


class ComplexJsonCoder(JSONEncoder):
    def __init__(self, translation_table=None, **kwargs):
        super(ComplexJsonCoder, self).__init__(**kwargs)
        self.translation_table = translation_table

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
        if self.translation_table is not None:
            self.recursive_translate(model, self.translation_table.to_anon)

        return model

    def decode(self, obj):
        from pydoc import locate

        # _type is a field embedded by the encoder to indicate which Jsonify instance will be used to decode the json string.
        if "_type" not in obj:
            return obj

        obj_type = locate(obj["_type"])

        if obj_type is None:
            raise Exception(f"{obj_type} does not exist.")

        # If the type is not a subclass of Jsonify, throw an exception. We do this to prevent the execution of arbitrary code.
        if not issubclass(obj_type, Jsonify):
            raise Exception(f"{obj_type}")

        decoded = obj_type.from_json(obj, self.translation_table)

        # If a translation_table is present, we need to convert anon values to named.
        if self.translation_table is not None:
            self.recursive_translate(decoded, self.translation_table.from_anon)

        return decoded

    def recursive_translate(self, obj, translation_table):
        from collections import OrderedDict, Hashable

        if obj is None or translation_table is None:
            return

        # If the input object is a list, we iterate through it element by element.
        if isinstance(obj, list):
            for i, item in enumerate(list(obj)):
                if isinstance(item, dict) or isinstance(item, list):
                    self.recursive_translate(item, translation_table)
                    continue

                if not isinstance(item, Hashable):
                    continue

                if item in translation_table:
                    obj[i] = translation_table[item]

            return

        # Else, the item is a dictionary, so we iterate through each key/value.
        for k in list(obj.keys()):
            if k in translation_table:
                obj[translation_table[k]] = obj.pop(k)
                k = translation_table[k]

            # If the value is a list, we need to iterate through it.
            if isinstance(obj[k], list):
                self.recursive_translate(obj[k], translation_table)
                continue

            # OrderedDicts are immutable, so we need to convert it into a dictionary prior to translation.
            if isinstance(obj[k], OrderedDict):
                obj[k] = dict(obj[k])

            # We need to translate all sub-elements in a dictionary, so recurse into it.
            if isinstance(obj[k], dict):
                self.recursive_translate(obj[k], translation_table)
                continue

            # If the value isn't Hashable, continue.
            if not isinstance(obj[k], Hashable):
                continue

            v = obj[k]
            if v in translation_table:
                obj[k] = translation_table[v]


class TranslationTable(Jsonify):
    def __init__(self, to_anon_table):
        self.to_anon = to_anon_table
        self.from_anon = { (v, k) for k, v in self.to_anon.items() }


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
    def from_json(json_object, translation_table):
        import numpy
        return numpy.linspace(start=json_object["start"], stop=json_object["end"], num=json_object["size"])