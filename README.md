disptools
=========
Generate displacement fields with known volume changes
------------------------------------------------------
[![GitHub release](https://img.shields.io/github/release/m-pilia/disptools.svg)](https://github.com/m-pilia/disptools/releases/latest)
[![PyPI](https://img.shields.io/pypi/v/disptools.svg)](https://pypi.python.org/pypi/disptools)
[![Wheels](https://img.shields.io/pypi/wheel/disptools.svg)](https://pypi.org/project/disptools)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/m-pilia/disptools/blob/master/LICENSE)
[![GitHub Actions](https://github.com/m-pilia/disptools/workflows/ChecksLinux/badge.svg)](https://github.com/m-pilia/disptools/actions/workflows/checks_linux.yml)
[![codecov](https://codecov.io/gh/m-pilia/disptools/branch/master/graph/badge.svg)](https://codecov.io/gh/m-pilia/disptools/branch/master)
[![Downloads](https://pepy.tech/badge/disptools)](https://pepy.tech/project/disptools)

This library provides utilities to generate and manipulate displacement fields with known volume changes. It implements three search-based algorithms for the generation of deformation fields, along with a small collection of utility functions, and provides optional GPU acceleration through a CUDA implementation.

The three algorithms implemented are referred as:
+ <tt>gradient</tt>: a gradient descent method (default).
+ <tt>greedy</tt>: a greedy search method proposed in [[1]](#1).
+ <tt>matching</tt>: a volume matching method proposed in [[2]](#2) and [[3]](#3). The implementation comes from the [PREDICT atrophysim tool](https://www.nitrc.org/projects/atrophysim).

The library is built on top of SimpleITK, in order to provide a simple yet powerful set of functionalities for image analysis. Images stored as numpy arrays can be easily converted from and to [SimpleITK](http://simpleitk.github.io/SimpleITK-Notebooks/01_Image_Basics.html) and [ITK](https://blog.kitware.com/convert-itk-data-structures-to-numpy-arrays/) image objects.

### Documentation

The complete documentation for this package is available on https://martinopilia.com/disptools.

### Quick example
```python
import SimpleITK as sitk
import disptools.displacements as dsp
import disptools.drawing as drw

# Create an example Jacobian map
# A spherical ROI with a Jacobian of 1.1 (expansion)
jacobian = drw.create_sphere(10, 40, fg_val=1.1, bg_val=1.0)

# Create a binary mask for the ROI
mask = drw.create_sphere(10, 40) > 0

# Generate the displacement
displacement = dsp.displacement(jacobian, mask=mask)

# Check the correctness of the result within the ROI
error = jacobian - dsp.jacobian(displacement)
error = sitk.Mask(error, mask)
```

A 3D rendering of the resulting displacement field with [ParaView](https://www.paraview.org/), and a visualisation of the the error on the Jacobian within the ROI:

<img src="https://github.com/m-pilia/disptools/blob/master/sphinx/img/example_2.png?raw=true" alt="displacement" width=45% /> <img src="https://github.com/m-pilia/disptools/blob/master/sphinx/img/example_1.png?raw=true" alt="error" width=45% />

### Architecture

The project is structured in three layers:
+ a pure standard [C99](https://en.wikipedia.org/wiki/C99) library (whose headers are in <tt>src/headers</tt>), with no external dependencies, implementing the core algorithms for the generation of deformation fields. It is a standalone library that can be directly included in a C or C++ project. It is paired with an optional [CUDA](https://en.wikipedia.org/wiki/CUDA) library, whose headers are in <tt>cuda/headers</tt>, that depends on <tt>src/headers</tt> and provides a GPGPU implementation of the key routines.
+ a [Python C extension](https://docs.python.org/3.6/extending/extending.html) package <tt>_disptools</tt> (whose source is in the file <tt>python_c_extension/_disptools.c</tt>), providing a bare Python wrapper to the aforementioned library, using the [NumPy C API](https://docs.scipy.org/doc/numpy-1.14.0/reference/c-api.html) to pass arrays. This can be directly included in a Python project with no dependencies other than NumPy.
+ a Python package (<tt>disptools</tt>), that wraps the <tt>_disptools</tt> package providing file IO (through SimpleITK) and implementing high-level features (such as the multi-resolution framework) and auxiliary utilities and functions.
    - <tt>disptools.displacements</tt>: module providing the main functions for the generation and manipulation of displacement fields.
    - <tt>disptools.drawing</tt>: collection of utilities to create test images.
    - <tt>disptools.io</tt>: collection of utilities to read and write to file.
    - <tt>disptools.measure</tt>: collection of utilities to measure some image features.
    - <tt>disptools.simulatrophy</tt>: routines to interface with the [Simul@atrophy](https://github.com/Inria-Asclepios/simul-atrophy) tool.
    - <tt>disptools.predict</tt>: routines to interface with the [PREDICT](https://www.nitrc.org/projects/atrophysim) tool.

### Install

This package is available on [PyPI](https://pypi.org/project/disptools) both as source distribution and as a Windows pre-compiled binary wheel. You can install it with <tt>pip</tt>:
```bash
 python3 -m pip install disptools
```

As always, it is recommended to use the package inside a [virtual environment](https://docs.python.org/3/tutorial/venv.html).

### Build from source

#### Requirements

Requirements are specified by the <tt>requirements.txt</tt> file and can be installed with <tt>pip</tt>.
```bash
python3 -m pip install -r requirements.txt
```

The library is a cross-platform Python 3.5+ package, with a compiled C extension. The Python dependencies are:
+ [numpy](https://github.com/numpy/numpy) ([pypi package](https://pypi.python.org/pypi/numpy))
+ [scipy](https://github.com/scipy/scipy) ([pypi package](https://pypi.org/pypi/scipy))
+ [SimpleITK](https://github.com/SimpleITK/SimpleITK) ([pypi package](https://pypi.org/pypi/SimpleITK))

Build dependencies are a standard C compiler (tested with gcc 8.2 on Linux, mingw-w64 7.2 and MSVC 19 on Windows 10), [CMake](https://cmake.org/), the [numpy](https://pypi.python.org/pypi/numpy) and the [setuptools](https://pypi.python.org/pypi/setuptools) packages. [scikit-build](https://pypi.python.org/pypi/scikit-build) may be required to build the other Python dependencies.

Some optional dependencies are required only for a limited set of features, and the package should build and run without them:
+ [itk](https://github.com/InsightSoftwareConsortium/ITK) ([pypi package](https://pypi.org/project/itk)): for <tt>disptools.drawing.sitk_to_itk</tt>
+ [vtk](https://github.com/Kitware/VTK) ([pypi package](https://pypi.org/project/vtk)): for <tt>disptools.io.write_vtk_points</tt>
+ [ply](https://github.com/dabeaz/ply) ([pypi package](https://pypi.org/project/ply)): for <tt>disptools.io.read_elastix_parameters</tt>
+ [scikit-image](https://github.com/scikit-image/scikit-image) ([pypi package](https://pypi.org/project/scikit-image)): for some features of <tt>disptools.drawing.extract_slice</tt>, and to run the unit tests

#### Build options

The following environment variables affect the <tt>setup.py</tt>:
+ <tt>DISPTOOLS_OPT=ON</tt>: enable non-portable optimisations.
+ <tt>DISPTOOLS_DEBUG=ON</tt>: disable optimisations, compile with debug symbols.
+ <tt>DISPTOOLS_CUDA_SUPPORT=ON</tt>: enable CUDA support.

#### Windows (Visual Studio) and Linux

Install the dependencies with your favourite package manager. For example, with <tt>pip</tt>:
```bash
python3 -m pip install -r requirements.txt
```

The package provides a <tt>setuptools</tt> based install script. To install the library, run from the project root folder
```bash
python3 setup.py install
```
which should build the C extension and install the Python package.

#### Windows (MinGW)

1. First, be sure that [mingw](https://mingw-w64.org), CMake and Python are installed and their executables [are in your PATH](https://docs.python.org/3/using/windows.html#excursus-setting-environment-variables).

2. Ensure that <tt>gcc</tt> is working correctly:
```none
> gcc --version
gcc (x86_64-posix-seh-rev1, Built by MinGW-W64 project) 7.2.0
Copyright (C) 2017 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```

3. Ensure that <tt>distutils</tt> correctly recognises your version of Visual Studio (even when using <tt>mingw</tt>). Open the file <tt>C:\Users\yourname\AppData\Local\Programs\Python\Python3x\Lib\distutils\cygwinccompiler.py</tt> (the exact location may vary according to your setup) and check that your version of Visual Studio is present in the function <tt>get_msvcr()</tt>; if not, adjust it according to the following:
```python
def get_msvcr():
    """Include the appropriate MSVC runtime library if Python was built
    with MSVC 7.0 or later.
    """
    msc_pos = sys.version.find('MSC v.')
    if msc_pos != -1:
        msc_ver = sys.version[msc_pos+6:msc_pos+10]
        if msc_ver == '1300':
            # MSVC 7.0
            return ['msvcr70']
        elif msc_ver == '1310':
            # MSVC 7.1
            return ['msvcr71']
        elif msc_ver == '1400':
            # VS2005 / MSVC 8.0
            return ['msvcr80']
        elif msc_ver == '1500':
            # VS2008 / MSVC 9.0
            return ['msvcr90']
        elif msc_ver == '1600':
            # VS2010 / MSVC 10.0
            return ['msvcr100']
        elif msc_ver == '1700':
            # Visual Studio 2012 / Visual C++ 11.0
            return ['msvcr110']
        elif msc_ver == '1800':
            # Visual Studio 2013 / Visual C++ 12.0
            return ['msvcr120']
        elif msc_ver == '1900':
            # Visual Studio 2015 / Visual C++ 14.0
            # "msvcr140.dll no longer exists" http://blogs.msdn.com/b/vcblog/archive/2014/06/03/visual-studio-14-ctp.aspx
            return ['vcruntime140']
        else:
            raise ValueError("Unknown MS Compiler version %s " % msc_ver)
```

4. Ensure that the library <tt>vcruntime140.dll</tt> is present in your library path. Otherwise, download it and place it in <tt>C:\Users\yourname\AppData\Local\Programs\Python\Python3x\libs</tt> (the exact path may vary according to your setup).

5. Clone the sources of this package with <tt>git</tt>, or download and extract them as a <tt>zip</tt> archive. Move to the root folder of the sources (<tt>C:\Users\yourname\disptools</tt> in this example), specify the right compiler, and launch the setup script to build and install the package.
```cmd
> cd C:\Users\yourname\disptools
> python setup.py setopt --command=build --option=compiler --set-value=mingw32
> python -m pip install -r requirements.txt
> python setup.py install
```

### References
+ <a id="1"></a>[1] van Eede, M. C., Scholz, J., Chakravarty, M. M., Henkelman, R. M., and Lerch, J. P. *Mapping registration sensitivity in MR mouse brain images.* Neuroimage 82 (2013), 226–236.
+ <a id="2"></a>[2] Karaçali, B., and Davatzikos, C. *Estimating topology preserving and smooth displacement fields.* IEEE Transactions on Medical Imaging 23, 7 (2004), 868–880.
+ <a id="3"></a>[3] Karaçali, B., and Davatzikos, C. *Simulation of tissue atrophy using a topology preserving transformation model.* IEEE transactions on medical imaging 25, 5 (2006), 649–652.

### License

The software is distributed under the MIT license.
