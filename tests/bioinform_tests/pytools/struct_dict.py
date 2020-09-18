"""
StructDictClass - contains classes StructDict and MatlabStructDict to make it
    possible to access a dictionary with syntax struct.field.
"""

import re
import typing


class StructDict(dict):
    '''Class that adds a structure type syntax to dictionaries,
       i.e. 'dict.field' will invoke dict['field']
    '''

    def __getattr__(self, key):
        '''Implement getattr to support syntax "data.field"'''
        try:
            val = self[key]
        except KeyError:
            val = None
        return val

    def __setattr__(self, key, val):
        '''Implement setattr to support syntax "data.field=x"'''
        self[key] = val

    def __delattr__(self, key):
        '''Implement delattr to support syntax "del data.field"'''
        if key in self:
            del self[key]

    def __getstate__(self):
        '''Needed for pickling, return the underlying dictionary'''
        return dict(self)

    def __setstate__(self, dict_entries):
        '''Needed for pickling, set the underlying dictionary'''
        self.update(dict_entries)

    def copy(self):
        return StructDict(super().copy())


def copy_toplevel(data):
    cptl = StructDict()
    for key, val in data.items():
        if isinstance(val, dict):
            continue
        if type(val) == list:
            if isinstance(val[0], dict):
                continue
        cptl[key] = val
    return cptl


def recurseCreateStructDict(data):
    '''Given a recursive dictionary, i.e. a dictionary that has
            child dictionaries or lists of dictionaries,
            convert each child dictionary to a StructDict.
      '''
    if isinstance(data, dict):
        tmpDict = StructDict()
        for key, value in data.items():
            tmpDict[key] = recurseCreateStructDict(value)
        return tmpDict
    elif isinstance(data, list):
        tmpList = []
        for value in data:
            tmpList.append(recurseCreateStructDict(value))
        return tmpList
    return data
