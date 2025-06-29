from pytest import raises
from importlib import import_module

def test_load_pkg():
	import eprsim as es

def test_load_mod():
	modules = ["Convolutions","Direct_conversion_to_Field","EPRsim", "EPRload",
	"FastMotion","Hamiltonian_Eig","Hamiltonian_Point_Group","Interpolation_lib",
	 "Nucdic","Pauli_generators","Presettings","resfield_full","SolidState",
	"spectral_processing","Tools","Validate_input_parameter"] 
	for mod in modules:
		my_module = import_module('eprsim.'+mod)

def test_bad_pkg_struct():
	with raises(ModuleNotFoundError):
		import EPRsim
	with raises(ModuleNotFoundError):	
		import Tools
