class ConditionalList(list):
	# a class that conditionally accepts new elements based on user-defined functions
	def __init__(self, input_list, *checkfuncs):
		if not isinstance(input_list, list):
			input_list = [input_list]
		for checkfunc in checkfuncs:
			if not hasattr(checkfunc, "__call__"):
				raise TypeError("Need to pass a callable function as a parameter")
			for item in input_list:
				checkfunc(item)
		super(ConditionalList, self).__init__(input_list)
		self.__checkfuncs = list(checkfuncs)

		for item in ["__add__", "__iadd__", "append", "extend", "insert", "remove"]:
			setattr(ConditionalList, item, self._check(getattr(list, item)))

	def _check(self, listfunc):
		def decorated(*args):
			print(self)
			if not self:
				return
			items = args[-1]
			if not isinstance(items, list):
				items = [items]
			for checkfunc in self.__checkfuncs:
				for item in items:
					checkfunc(item)
			return listfunc(*args)
		return decorated
