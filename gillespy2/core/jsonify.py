# GillesPy2 is a modeling toolkit for biochemical simulation.
# Copyright (C) 2019-2022 GillesPy2 developers.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import re
import copy
import json
import pydoc
import hashlib

from json import JSONEncoder

import numpy
import gillespy2

class Jsonify:
    """
    Interface to allow for instances of arbitrary types to be encoded into json strings
    and decoded into new objects.
    """

    _translation_table = None

    def to_json(self, encode_private=True) -> str:
        """
        Convert self into a json string.

        :param encode_private: If True all private (prefixed of '_') will be encoded. None if False.
        :type encode_private: bool

        :returns: The JSON representation of self.
        """

        encoder = ComplexJsonCoder(encode_private=encode_private)
        return json.dumps(copy.deepcopy(self), indent=4, sort_keys=True, default=encoder.default)

    @classmethod
    def from_json(cls, json_str: str) -> object:
        """
        Convert some json_str into a decoded Python type. This function should
        return a new instance of the type.

        :param json_str: A json str to be converted into a new type instance.
        :type json_str: str

        :returns: A decoded object.
        """

        # If the json_str is actually a dict, it means we've decoded as much as possible.
        if isinstance(json_str, dict):
            return cls.from_dict(json_str)

        decoder = ComplexJsonCoder()
        return json.loads(json_str, object_hook=decoder.decode)

    def to_dict(self) -> dict:
        """
        Convert the object into a dictionary ready for json encoding.
        Note: Complex types that inherit from Jsonify do not need to be manually encoded.

        By default, this function will return a dictionary of the object's public types.

        :returns: The backing var dictionary of the object.
        """

        return self.__dict__.copy()

    @classmethod
    def from_dict(cls, src_dict: dict) -> object:
        """
        Convert some dict into a new instance of a python type.
        This function will return a __new__ instance of the type.

        :param src_dict: The dictionary to apply onto the new instance.
        :type src_dict: dict

        :returns: A new object with its backing __dict__ set to a copy of src_dict.
        """

        new = cls.__new__(cls)
        new.__dict__ = src_dict.copy()

        return new

    def to_anon(self):
        """
        Converts self into an anonymous instance of self.
        """

        return self.get_translation_table().obj_to_anon(copy.deepcopy(self))

    def to_named(self):
        """
        Converts self into a named instance of self.
        """

        return self.get_translation_table().obj_to_named(copy.deepcopy(self))

    def get_translation_table(self) -> "TranslationTable":
        """
        Make and/or return the translation table.

        :returns: A TranslationTable instance.
        """

        if self._translation_table is None:
            self._translation_table = self.make_translation_table()

        return self._translation_table

    def make_translation_table(self) -> "TranslationTable":
        """
        Make a translation table that describes key:value pairs to convert
        user-define data into generic equivalents.

        :returns: A newly generated TranslationTable instance.
        """
        raise NotImplementedError("make_translation_table() has not been implemented.")

    def public_vars(self) -> dict:
        """
        Gets a dictionary of public vars that exist on self. Keys starting with '_' are ignored.

        :returns: A dictionary containing all public ('_'-prefixed) variables on the object.
        """
        return {k: v for k, v in vars(self).items() if not k.startswith("_")}

    def get_json_hash(self, ignore_whitespace=True, hash_private_vals=False) -> str:
        """
        Get the hash of the json representation of self.

        :param ignore_whitespace: If set to True all whitespace will be stripped from the JSON
            prior to being hashed.
        :type ignore_whitespace: bool

        :param hash_private_vals: If set to True all private and non-private variables will
            be included in the hash.
        :type hash_private_vals: bool

        :returns: An MD5 hash of the object's sorted JSON representation.
        """

        # If _hash_private_vars is set, hash ALL properties on the object.
        model_json = self.to_json(encode_private=hash_private_vals)

        # If ignore_whitespace is set, strip out all whitespace characters.
        if ignore_whitespace:
            model_json = re.sub(r"\s+", "", model_json)

        return hashlib.md5(str.encode(model_json)).hexdigest()

    def __eq__(self, obj: "Jsonify"):
        """
        Overload to compare the json of two objects that derive from Jsonify.
        This method will not do any additional translation.

        :param obj: The Jsonify object to compare against.
        :type obj: Jsonify

        :returns: True if equal, False if not.
        """
        return self.get_json_hash() == obj.get_json_hash()

class ComplexJsonCoder(JSONEncoder):
    """
    This class delegates the encoding and decoding of objects to one or more implementees.

    :param translation_table: A TranslationTable instance that will be used to
        translate objects.
    :type translation_table: TranslationTable

    :param encode_private: If set to True then all private and public variables will
        be converted to JSON. If False, only public.
    :type encode_private: bool
    """

    def __init__(self, translation_table=None, encode_private=True, **kwargs):
        super().__init__(**kwargs)
        self._translation_table = translation_table
        self._encode_private = encode_private

        self._delegation_table = {
            numpy.ndarray: NdArrayCoder,
            numpy.int64: Int64Coder,
            set: SetCoder,
            type: TypeCoder
        }

    def default(self, o: object):
        """
        This function is called when json.dumps() fires. default() is a bad name for the function,
        but anything else makes JSONEncoder freak out.

        :param o: The object that is currently being encoded into JSON.
        """

        # If o is of matching type, use a custom coder.
        for obj_type, coder in self._delegation_table.items():
            if isinstance(o, obj_type):
                return coder.to_dict(o)

        if not isinstance(o, Jsonify):
            return super().default(o)

        if self._encode_private:
            model = o.to_dict()

        else:
            # Strip private variables from the object.
            model = {}

            for key, val in o.to_dict().items():
                if key.startswith("_") and not key.startswith("__"):
                    continue

                model[key] = val

        # If the model is some subclass of gillespy2.core.model.Model, then manually set its type.
        if issubclass(o.__class__, gillespy2.core.Model):
            model["_type"] = f"{gillespy2.core.Model.__module__}.{gillespy2.core.Model.__name__}"

        else:
            model["_type"] = f"{o.__class__.__module__}.{o.__class__.__name__}"

        return model

    def decode(self, json_dict: dict):
        """
        Decode the JSON dictionary into a valid Python object.

        :param json_dict: The JSON dictionary to decode.
        :type json_dict: dict

        :returns: The decoded form of the JSON dictionary.
        """

        # _type is a field embedded by the encoder to indicate which
        # Jsonify instance will be used to decode the json string.
        if "_type" not in json_dict:
            return json_dict

        json_type = pydoc.locate(json_dict["_type"])

        if json_type is None:
            raise Exception(f"{json_type} does not exist.")

        # If the type is not a subclass of Jsonify, throw an exception.
        # We do this to prevent the execution of arbitrary code.
        if not issubclass(json_type, Jsonify):
            raise Exception(f"{json_type}")

        return json_type.from_json(json_dict)

class TranslationTable(Jsonify):
    """
    This class contains functions to enable arbitrary object trees to be "translated" to anonymous
    and named objects. This behavior is defined by a map of 'named' and 'anon' key values.

    :param to_anon: A mapping of 'named' to 'anonymous' strings to be used when converting
        user-defined names to anon.
    :type to_anon: dict[str, str]
    """

    def __init__(self, to_anon: "dict[str, str]"):
        self.to_anon = to_anon.copy()
        self.to_named = dict((v, k) for k, v in list(self.to_anon.items()))

    def obj_to_anon(self, obj: object):
        """
        Recursively anonymise all named properties on the object.

        :param obj: The object to anonymize.
        :type obj: object

        :returns: An anonymized instance of self.
        """

        return self.recursive_translate(obj, self.to_anon)

    def obj_to_named(self, obj):
        """
        Recursively identify all anonymous properties on the object.

        :param obj: The object that will be converted to named.
        :type obj: object

        :returns: A named instance of self.
        """

        return self.recursive_translate(obj, self.to_named)

    def recursive_translate(self, obj: object, translation_table: "dict[str, str]"):
        """
        Recursively search through the object tree searching for property value matches in the
        translation table. If a match is found, substitute.

        :param obj: The object that will be translated.
        :type obj: object

        :param translation_table: The mapping to translate by.
        :type translation_table: TranslationTable
        """

        # If a translation table exists on the object, remove and save it.
        if "_translation_table" in obj.__dict__:
            saved_table = obj.__dict__.pop("_translation_table")

        translated = self._recursive_translate(obj, translation_table)

        # Restore the original translation table, if needed.
        if saved_table is not None:
            obj.__dict__["_translation_table"] = saved_table

        return translated

    def _recursive_translate(self, obj: object, translation_table: "dict[str, str]"):
        # The obj is a class if it's an instance of Jsonify. Class property names *cannot*
        # be changed, so translate just the values.
        if isinstance(obj, Jsonify):
            for key in vars(obj).keys():
                vars(obj)[key] = self._recursive_translate(vars(obj)[key], translation_table)

        elif isinstance(obj, list):
            for item in obj:
                item = self._recursive_translate(item, translation_table)

        elif isinstance(obj, dict):
            # Convert the dictionary into a list of tuples.
            # This makes it easier to modify key names.
            obj = list((k, v) for k, v in obj.items())
            new_pairs = [ ]

            for pair in obj:
                new_pairs.append((
                    self._recursive_translate(pair[0], translation_table),
                    self._recursive_translate(pair[1], translation_table)
                ))

            obj = dict((x[0], x[1]) for x in new_pairs)

        # If the obj is a string, translate it via a regex replace.
        # Note: mathematical functions contain additional characters that should not be translated.
        elif isinstance(obj, str):
            # To handle functions, grab all complete words from the string.
            matches = re.finditer("([0-z])+", obj)

            # For each match, translate the group.
            for match in matches:
                group = match.group()
                obj = obj.replace(group, translation_table.get(group, group))

        return obj

class NdArrayCoder(Jsonify):
    """ This JSON coder enables support for  the `numpy.ndarray` type. """

    @staticmethod
    def to_dict(obj):
        return {
            "data": obj.tolist(),
            "_type": f"{NdArrayCoder.__module__}.{NdArrayCoder.__name__}"
        }

    @staticmethod
    def from_json(obj):
        return numpy.array(obj["data"])

class Int64Coder(Jsonify):
    """ This JSON coder enables support for  the `numpy.int64` type. """

    @staticmethod
    def to_dict(obj):
        return {
            "data": int(obj),
            "_type": f"{Int64Coder.__module__}.{Int64Coder.__name__}"
        }

    @staticmethod
    def from_json(obj):
        return numpy.int64(obj["data"])

class SetCoder(Jsonify):
    """ This JSON coder enables support for the `set` type. """

    @staticmethod
    def to_dict(obj):
        return {
            "data": list(obj),
            "_type": f"{SetCoder.__module__}.{SetCoder.__name__}"
        }

    @staticmethod
    def from_json(obj):
        return set(obj["data"])

class TypeCoder(Jsonify):
    """ This JSON coder enables support for the 'type' type. """

    @staticmethod
    def to_dict(obj):
        return {
            "data": type(obj),
            "_type": f"{TypeCoder.__module__}.{TypeCoder.__name__}"
        }

    @staticmethod
    def from_json(obj):
        return pydoc.locate(obj["data"])
