import ConfigParser

def updatefromconfig(cat, opt, filename):
	Config = ConfigParser.ConfigParser()
	Config.read(filename)
	if Config.has_option(cat, opt):
		return Config.get(cat, opt)
	else:
		return [opt]
def updatebooleanconfig(cat, opt, filename):
	Config = ConfigParser.ConfigParser()
	Config.read(filename)
	if Config.has_option(cat, opt):
		return Config.getboolean(cat, opt)
	else:
		return globals()[opt]

