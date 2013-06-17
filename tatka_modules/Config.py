from configobj import ConfigObj
from validate import Validator
# -*- coding: utf8 -*- 
class Config(dict):
    #The following is to make Config keys accessible as attributes
    def __setitem__(self, key, value):
	    super(Config, self).__setattr__(key, value)
	    super(Config, self).__setitem__(key, value)
    def __getattr__(self, key):
        try:
            return self.__getitem__(key)
        except KeyError, message:
            raise AttributeError, message
    __setattr__ = __setitem__


def Config_Input(ConfigSpec_file,Config_file):
    configspec = ConfigObj(ConfigSpec_file, interpolation=False, list_values=False,
                        _inspec=True)
    _config = ConfigObj(Config_file, configspec=configspec)

    val = Validator()
    test = _config.validate(val)
    if test == True:
        print 'Succeeded'

    return Config(_config.dict().copy())
