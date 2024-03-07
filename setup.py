from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext
import platform

std_version = "/std:c++17" if platform.system() == "Windows" else "-std=c++17"
ext_modules = [
    Pybind11Extension(
    "_ruptura",
    sources=glob("src/*.cpp"),
    extra_compile_args= [std_version, "-DPYBUILD=1", "-D_LIBCPP_DISABLE_AVAILABILITY"]
    )
]

setup(
    ext_modules=ext_modules,
    zip_safe=False,
    packages=["ruptura"],
    package_dir={"ruptura":"src"}
)