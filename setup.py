import re
import sys
import sysconfig
import numpy.distutils
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

_PEDANTIC = False
if "--pedantic" in sys.argv:
    _PEDANTIC = True
    sys.argv.remove("--pedantic")

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

# C preprocessor macros
define_macros = [
    ('DISPTOOLS_VERBOSE',  '1'),
    ('ORDER_PD', '4'),
    ('FLOAT_SIZE', '32'),
]
if _DEBUG:
    define_macros.append(('DISPTOOLS_DEBUG', '1'))
else:
    define_macros.append(('DISPTOOLS_DEBUG', '0'))


# C compiler and linker flags
cflags = {
    'msvc'   : ['/openmp', '/Ox', '/fp:fast'],
    'mingw32': ['-Wall', '-std=c99', '-fPIC', '-fopenmp'],
    'unix'   : ['-Wall', '-std=c99', '-fPIC', '-fopenmp'],
    'cygwin' : ['-Wall', '-std=c99', '-fPIC', '-fopenmp'],
}
cflags_opt = {
    'msvc'   : [],
    'mingw32': ['-O3', '-march=native'],
    'unix'   : ['-O3', '-march=native'],
    'cygwin' : ['-O3', '-march=native'],
}
cflags_debug = {
    'msvc'   : [],
    'mingw32': ['-O0', '-g'],
    'unix'   : ['-O0', '-g'],
    'cygwin' : ['-O0', '-g'],
}
ldflags = {
    'msvc'    : [],
    'mingw32' : ['-lm', '-lgomp'],
    'unix'    : ['-lm', '-lgomp'],
    'cygwin'  : ['-lm', '-lgomp'],
}


class build_ext_xplatform(build_ext):
    """ Class to build the extension with a cross platform setup.
    """

    def build_extensions(self):
        compiler = self.compiler.compiler_type

        if compiler == 'msvc':
            raise Exception(
                    'You are trying to compile this package with msvc (Visual Studio). ' +
                    'msvc does not support C99. Please compile this package with mingw. ' +
                    'Refer to the documentation for more details')

        if compiler in cflags.keys():
            for e in self.extensions:
                e.extra_compile_args = cflags[compiler]
                if _PEDANTIC:
                    e.extra_compile_args.append('--pedantic')
                if _DEBUG:
                    e.extra_compile_args += cflags_debug[compiler]
                elif _OPT:
                    e.extra_compile_args += cflags_opt[compiler]
                print('extra_compile_args: ' + str(e.extra_compile_args))

        if compiler in ldflags.keys():
            for e in self.extensions:
                e.extra_link_args = ldflags[compiler]
                print('extra_link_args: ' + str(e.extra_link_args))

        build_ext.build_extensions(self)


disptools_c = Extension('_disptools',
                    define_macros = define_macros,
                    include_dirs = ['src/headers',
                                    sysconfig.get_path('include')] +
                                    numpy.distutils.misc_util.get_numpy_include_dirs(),
                    libraries = [],
                    library_dirs = [],
                    sources = ['src/displacement_field_gradient.c',
                               'src/displacement_field_greedy.c',
                               'src/field.c',
                               'src/jacobian.c',
                               'src/shape_descriptors.c',
                               'src/VolumeMatching3D.c',
                               'src/_disptools.c'
                               ])

version = '0.2.2'

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
    cmdclass = {'build_ext': build_ext_xplatform}
    )

