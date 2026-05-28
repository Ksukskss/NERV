import pybind11
from setuptools import setup, Extension

ext_modules = [
    Extension(
        "ray_tracer_cpp",  # Имя модуля в Python
        ["ray_tracer.cpp"],     # Исходный код C++
        include_dirs=[pybind11.get_include()],
        language='c++',
        extra_compile_args=['-O3', '-std=c++14'],  # -O3 для максимальной скорости!
    ),
]

setup(
    name="seismic_tracer_cpp",
    version="1.0.0",
    ext_modules=ext_modules,
)
