from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext
import platform

# MSVC, compiler used in conda uses different flag for the c++ version and does not define the __cplusplus variable by default
# set it so that the breakthrough video chmods don't break (Windows also doesn't have S_IRWXU)
std_version = ["/std:c++17", '/Zc:__cplusplus'] if platform.system() == "Windows" else ["-std=c++17"]
ext_modules = [
    Pybind11Extension(
    "_ruptura",
    sources=glob("src/*.cpp"),
    extra_compile_args= std_version + ["-DPYBUILD=1", "-D_LIBCPP_DISABLE_AVAILABILITY"]
    )
]

setup(
    ext_modules=ext_modules,
    zip_safe=False,
    packages=["ruptura"],
    package_dir={"ruptura":"src"}
)