import os
import platform
import re
import subprocess
import sys
import sysconfig
import numpy.distutils
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

_OPT = False
if "--opt" in sys.argv:
    _OPT = True
    sys.argv.remove("--opt")

_DEBUG = False
if "--debug" in sys.argv:
    _DEBUG = True
    sys.argv.remove("--debug")

# Read description from README
with open('README.md') as f:
    long_description = f.read()


class CMakeExtension(Extension):
    def __init__(self, name, cmake_lists_dir='', **kwa):
        Extension.__init__(self, name, **kwa)
        self.cmake_lists_dir = os.path.abspath(cmake_lists_dir)


class CMakeBuild(build_ext):

    def build_extensions(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError('Cannot find CMake executable')

        for ext in self.extensions:

            if self.compiler.compiler_type == 'msvc':
                raise RuntimeError(
                        'You are trying to compile this package with MSVC (Visual Studio). ' +
                        'MSVC does not support standard C99. Please compile this package with mingw. ' +
                        'Refer to the documentation for more details')

            extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
            cfg = 'Debug' if _DEBUG else 'Release'
            build_args = ['--config', cfg]

            cmake_args = [
                '-DDISPTOOLS_DEBUG=%s' % ('ON' if cfg == 'Debug' else 'OFF'),
                '-DDISPTOOLS_OPT=%s' % ('ON' if _OPT else 'OFF'),
                '-DDISPTOOLS_VERBOSE=ON',
                '-DDISPTOOLS_LOW_ORDER_PD=OFF',
                '-DDISPTOOLS_DOUBLE=OFF',
                '-DCMAKE_BUILD_TYPE=%s' % cfg,
                '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
            ]

            if platform.system() == "Windows":
                cmake_args += [
                    '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir),
                    '-G', 'MinGW Makefiles',
                ]

            if not os.path.exists(self.build_temp):
                os.makedirs(self.build_temp)

            subprocess.check_call(['cmake', ext.cmake_lists_dir] + cmake_args,
                                  cwd=self.build_temp)
            subprocess.check_call(['cmake', '--build', '.'] + build_args,
                                  cwd=self.build_temp)

            for e in self.extensions:
                e.extra_compile_args = [
                    '-fopenmp',
                ]
                if _OPT:
                    e.extra_compile_args += ['-O3', '-march=native']
                if _DEBUG:
                    e.extra_compile_args += ['-O0', '-g']

                e.extra_link_args = [
                    '-ldisptools',
                    '-lgomp',
                    '-lm',
                    '-L' + os.path.join(self.build_temp, 'src'),
                ]

            build_ext.build_extensions(self)


disptools_c = CMakeExtension('_disptools',
                    cmake_lists_dir='.',
                    define_macros=[('ORDER_PD', 4)],
                    include_dirs = ['src/headers',
                                    sysconfig.get_path('include')] +
                                    numpy.distutils.misc_util.get_numpy_include_dirs(),
                    libraries = [],
                    library_dirs = [],
                    sources = ['src/_disptools.c'])

version = '0.4.0'

setup(name = 'disptools',
    packages = ['disptools'],
    version = version,
    description = 'Generate displacements fields with known volume changes',
    author = 'Martino Pilia',
    author_email = 'martino.pilia@gmail.com',
    url = 'https://github.com/m-pilia/disptools',
    download_url = 'https://github.com/m-pilia/disptools/archive/v%s.tar.gz' % version,
    keywords = ['jacobian', 'displacement', 'image processing'],
    long_description = long_description,
    long_description_content_type='text/markdown',
    install_requires=[
        'numpy',
        'scipy',
        'SimpleITK',
    ],
    ext_modules = [disptools_c],
    cmdclass = {'build_ext': CMakeBuild},
    zip_safe=False,
    )

