import collections, json, hashlib, re, numpy

from json import JSONEncoder
from typing import Hashable


class Jsonify:
    """
    Interface to allow for instances of arbitrary types to be encoded into json strings and decoded into new objects.
    """

    def to_json(self, translation_table=None):
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

        vars = self.__dict__.copy()

        # We remove _hash since it's identifiable to the named status of a model.
        if "_hash" in vars:
            vars.pop("_hash")

        return vars

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

    def get_json_hash(self, translation_table=None):
        """
        Get the hash of the anonymous json representation of self.
        """

        return hashlib.md5(str.encode(self.to_json())).hexdigest()

    def to_anon(self):
        """
        Converts self into an anonymous instance of self.
        """

        jsoned = self.to_json()
        anon = self.get_translation_table().text_to_anon(jsoned)

        return self.from_json(anon)

    def to_named(self):
        """
        Converts self into a named instance of self.
        """

        named_json = self.get_translation_table().text_to_named(self.to_json)
        return self.from_json(named_json)

    def __eq__(self, o):
        """
        Overload to compare the json of two objects that derive from Jsonify. This method will not do any 
        additional translation.
        """

        return self.get_json_hash() == o.get_json_hash()

class ComplexJsonCoder(JSONEncoder):
    def __init__(self, translation_table=None, **kwargs):
        super(ComplexJsonCoder, self).__init__(**kwargs)
        self.translation_table = translation_table

    """
    This function is called when json.dumps() fires. default() is a bad name for the function,
    but anything else makes JSONEncoder freak out.

    :param o: The object that is currently being encoded into JSON.
    """
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
            model = self.recursive_translate(model, self.translation_table.to_anon)

        return model

    def decode(self, obj):
        from pydoc import locate

        # If a translation_table is present, we need to convert anon values to named.
        if self.translation_table is not None:
            obj = self.recursive_translate(obj, self.translation_table.from_anon)

        # _type is a field embedded by the encoder to indicate which Jsonify instance will be used to decode the json string.
        if "_type" not in obj:
            return obj

        obj_type = locate(obj["_type"])

        if obj_type is None:
            raise Exception(f"{obj_type} does not exist.")

        # If the type is not a subclass of Jsonify, throw an exception. We do this to prevent the execution of arbitrary code.
        if not issubclass(obj_type, Jsonify):
            raise Exception(f"{obj_type}")

        return obj_type.from_json(obj, self.translation_table)

    def recursive_translate(self, obj, translation_table):
        if isinstance(obj, list):
            for item in obj:
                item = self.recursive_translate(item, translation_table)

        elif isinstance(obj, dict):
            # Convert the dictionary into a list of tuples. This makes it easier to modify key names.
            obj = list((k, v) for k, v in obj.items())
            new_pairs = [ ]

            for pair in obj:
                new_pairs.append((
                    self.recursive_translate(pair[0], translation_table),
                    self.recursive_translate(pair[1], translation_table)
                ))

            obj = dict((x[0], x[1]) for x in new_pairs)

        elif isinstance(obj, collections.OrderedDict):
            obj = self.recursive_translate(dict(obj), translation_table)

        elif isinstance(obj, Hashable) and obj in translation_table.keys():
            obj = translation_table.get(obj, obj)

        return obj

class TranslationTable(Jsonify):
    def __init__(self, to_anon):
        self.to_anon = to_anon.copy()
        self.to_named = dict((v, k) for k, v in list(self.to_anon.items()))

    def text_to_anon(self, json_texts):
        return self._translate(json_texts, self.to_anon)

    def text_to_named(self, json_texts):
        return self._translate(json_texts, self.to_named)

    def _translate(self, text, table):
        # Grab the indexes of all matching keys via regex.
        # This will grab 1 or more characters comprised of '_' and any other alpha-numeric values.
        matches = list(re.finditer("([\_a-zA-Z0-9])+", text))

        last = 0
        translated = []

        # Iterate through each match, copying last and match boundaries into the list.
        for match in matches:
            translated.append(text[last:match.start()])
            translated.append(table.get(text[match.start():match.end()], text[match.start():match.end()]))

            last = match.end()
        
        # Write the last of the text to the buffer.
        translated.append(text[last:])

        return "".join(translated)

class NdArrayCoder(Jsonify):
    @staticmethod
    def to_dict(self):
        return {
            "data": self.tolist(),
            "_type": f"{NdArrayCoder.__module__}.{NdArrayCoder.__name__}"
        }

    @staticmethod
    def from_json(json_object, translation_table):
        return numpy.array(json_object["data"])