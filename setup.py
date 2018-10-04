import os
import platform
import re
import subprocess
import sys
import sysconfig
from pprint import pprint
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

# Filename for the C extension module library
c_module_name = '_disptools'

# Parse command line flags
options = {k: 'OFF' for k in ['--opt', '--debug', '--cuda']}
for flag in options.keys():
    if flag in sys.argv:
        options[flag] = 'ON'
        sys.argv.remove(flag)

# Command line flags forwarded to CMake
cmake_cmd_args = []
for f in sys.argv:
    if f.startswith('-D'):
        cmake_cmd_args.append(f)
        sys.argv.remove(f)


class CMakeExtension(Extension):
    def __init__(self, name, cmake_lists_dir='.', sources=[], **kwa):
        Extension.__init__(self, name, sources=sources, **kwa)
        self.cmake_lists_dir = os.path.abspath(cmake_lists_dir)


class CMakeBuild(build_ext):

    def build_extensions(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError('Cannot find CMake executable')

        for ext in self.extensions:

            extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
            cfg = 'Debug' if options['--debug'] == 'ON' else 'Release'

            cmake_args = [
                '-DDISPTOOLS_DEBUG=%s' % ('ON' if cfg == 'Debug' else 'OFF'),
                '-DDISPTOOLS_OPT=%s' % options['--opt'],
                '-DDISPTOOLS_VERBOSE=ON',
                '-DDISPTOOLS_LOW_ORDER_PD=OFF',
                '-DDISPTOOLS_DOUBLE=OFF',
                '-DDISPTOOLS_CUDA_SUPPORT=%s' % options['--cuda'],
                '-DDISPTOOLS_CUDA_ERROR_CHECK=ON',
                '-DDISPTOOLS_CUDA_ERROR_CHECK_SYNC=ON',
                '-DDISPTOOLS_PYTHON_SUPPORT=ON',
                '-DDISPTOOLS_PYTHON_C_MODULE_NAME=%s' % c_module_name,
                '-DCMAKE_BUILD_TYPE=%s' % cfg,
                '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir),
                '-DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), self.build_temp),
                '-DPYTHON_EXECUTABLE={}'.format(sys.executable),
            ]

            if platform.system() == 'Windows':
                plat = ('x64' if platform.architecture()[0] == '64bit' else 'Win32')
                cmake_args += [
                    '-DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=TRUE',
                    '-DCMAKE_RUNTIME_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir),
                ]
                if self.compiler.compiler_type == 'msvc':
                    cmake_args += [
                        '-DCMAKE_GENERATOR_PLATFORM=%s' % plat,
                    ]
                else:
                    cmake_args += [
                        '-G', 'MinGW Makefiles',
                    ]

            cmake_args += cmake_cmd_args

            pprint(cmake_args)

            if not os.path.exists(self.build_temp):
                os.makedirs(self.build_temp)

            # Config and build the extension
            subprocess.check_call(['cmake', ext.cmake_lists_dir] + cmake_args,
                                  cwd=self.build_temp)
            subprocess.check_call(['cmake', '--build', '.', '--config', cfg],
                                  cwd=self.build_temp)


# The following line is parsed by Sphinx
version = '0.4.0'

setup(name = 'disptools',
      packages = ['disptools'],
      version = version,
      description = 'Generate displacement fields with known volume changes',
      author = 'Martino Pilia',
      author_email = 'martino.pilia@gmail.com',
      url = 'https://github.com/m-pilia/disptools',
      download_url = 'https://github.com/m-pilia/disptools/archive/v%s.tar.gz' % version,
      keywords = ['jacobian', 'displacement', 'image processing'],
      long_description = open('README.md').read(),
      long_description_content_type='text/markdown',
      install_requires=['numpy', 'scipy', 'SimpleITK'],
      extras_require = {
          'ITK': ['itk>=4.12'],
          'VTK': ['vtk>=7.0'],
          'elastix': ['ply'],
          'testing': ['scikit-image'],
      },
      ext_modules = [CMakeExtension(c_module_name)],
      cmdclass = {'build_ext': CMakeBuild},
      zip_safe=False,
      )

