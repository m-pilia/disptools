disptools
=========
Generate displacement fields with known volume changes
------------------------------------------------------

This library provides utilities to generate and manipulate displacement fields with known volume changes. It implements three search-based algorithms for the generation of deformation fields, along with a small collection of utility functions.

The three algorithms implemented are referred as:
+ `gradient`: a gradient descent method (default).
+ `greedy`: a greedy search method based on the method proposed in [[1]](#1).
+ `matching`: a volume matching routine based on gradient descent,
  published in [[2]](#2) and [[3]](#3). The implementation comes from the [atrophysim tool](https://www.nitrc.org/projects/atrophysim).

The library is built on top of SimpleITK, in order to provide a simple yet powerful set of functionalities for image analysis. Images stored as numpy arrays can be easily converted from and to [SimpleITK](http://simpleitk.github.io/SimpleITK-Notebooks/01_Image_Basics.html) and [ITK](https://blog.kitware.com/convert-itk-data-structures-to-numpy-arrays/) image objects.

### Documentation

The complete documentation for this package is available on https://martinopilia.com/disptools.

### Architecture

The project is structured in three layers:
+ a pure standard [C99](https://en.wikipedia.org/wiki/C99) library (whose headers are in `src/headers/`), with no external dependencies, implementing the core algorithms for the generation of deformation fields. It is a standalone library that can be directly included in a C or C++ project.
+ a [Python C extension](https://docs.python.org/3.6/extending/extending.html) package `_disptools` (whose source is in the file `src/_disptools.c`), providing a bare Python wrapper to the aforementioned library, using the [NumPy C API](https://docs.scipy.org/doc/numpy-1.14.0/reference/c-api.html) to pass arrays. This can be directly included in a Python project with no dependencies other than NumPy.
+ a Python package (`disptools`), that wraps the `_disptools` package providing file IO (through SimpleITK) and implementing high-level features (such as the multi-resolution framework) and auxiliary utilities and functions.

### Requirements

The library is a cross-platform Python 3.5+ package, with a compiled C extension. The Python dependencies are:
+ [numpy](https://github.com/numpy/numpy) ([pypi package](https://pypi.python.org/pypi/numpy))
+ [scipy](https://github.com/scipy/scipy) ([pypi package](https://pypi.org/pypi/scipy))
+ [SimpleITK](https://github.com/SimpleITK/SimpleITK) ([pypi package](https://pypi.org/pypi/SimpleITK))

Build dependencies are a standard C99 compiler (tested with gcc 7.3 on Linux, mingw-w64 7.2 on Windows 10), the [numpy](https://pypi.python.org/pypi/numpy) and the [setuptools](https://pypi.python.org/pypi/setuptools) packages. [scikit-build](https://pypi.python.org/pypi/scikit-build) may be required to build the other Python dependencies.

Some optional dependencies are required only for a limited set of features, and the package should build and run without them:
+ [itk](https://github.com/InsightSoftwareConsortium/ITK) ([pypi package](https://pypi.org/project/itk)): for `disptools.drawing.sitk_to_itk`
+ [vtk](https://github.com/Kitware/VTK) ([pypi package](https://pypi.org/project/vtk)): for `disptools.io.write_vtk_points`
+ [ply](https://github.com/dabeaz/ply) ([pypi package](https://pypi.org/project/ply)): for `disptools.io.read_elastix_parameters`
+ [scikit-image](https://github.com/scikit-image/scikit-image) ([pypi package](https://pypi.org/project/scikit-image)): to run the unit tests

### Install

This package is available on [PyPI](https://pypi.org/project/disptools) both as source distribution and as a Windows pre-compiled binary wheel. You can install it with `pip`:
```bash
 python3 -m pip install disptools
```

### Build from source

#### Linux

Install the dependencies with your favourite package manager. For example, with `pip`:
```bash
python3 -m pip install scikit-build numpy scipy SimpleITK
```

The package provides a `setuptools` based install script. To install the library, run from the project root folder
```bash
python3 setup.py install
```
which should build the C extension and install the Python package.

#### Windows

1. First, be sure that the Python executables [are in your PATH](https://docs.python.org/3/using/windows.html#excursus-setting-environment-variables).

2. Since `msvc` (Visual Studio) does not support the C99 standard, this package should be compiled with [mingw](https://mingw-w64.org) (or any other C99-compliant compiler of your choice; in the following, instructions for `mingw` are provided). Install `mingw` and add its binary folder to your path. Ensure that `gcc` is working correctly:
```none
> gcc --version
gcc (x86_64-posix-seh-rev1, Built by MinGW-W64 project) 7.2.0
Copyright (C) 2017 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```

3. Ensure that `distutils` correctly recognises your version of Visual Studio (even when using `mingw`). Open the file `C:\Users\yourname\AppData\Local\Programs\Python\Python3x\Lib\distutils\cygwinccompiler.py` (the exact location may vary according to your setup) and check that your version of Visual Studio is present in the function `get_msvcr()`; if not, adjust it according to the following:
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

4. Ensure that the library `vcruntime140.dll` is present in your library path. Otherwise, download it and place it in `C:\Users\yourname\AppData\Local\Programs\Python\Python3x\libs` (the exact path may vary according to your setup).

5. Install the dependencies:
```cmd
> python -m pip install scikit-build numpy scipy SimpleITK
```

6. Clone the sources of this package with `git`, or download and extract them as a `zip` archive. Move to the root folder of the sources (`C:\Users\yourname\disptools` in this example), specify the right compiler, and launch the setup script to build and install the package.
```cmd
> cd C:\Users\yourname\disptools
> python setup.py setopt --command=build --option=compiler --set-value=mingw32
> python setup.py install
```

### Content

+ `disptools.displacements`: module providing the main functions for the generation and manipulation of displacement fields.
+ `disptools.drawing`: collection of utilities to create test images.
+ `disptools.io`: collection of utilities to read and write to file.
+ `disptools.measure`: collection of utilities to measure some image features.
+ `disptools.predict`: routines to interface with the PREDICT tool.
+ `disptools.simulatrophy`: routines to interface with the Simulatrophy tool.

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
error = jacobian - sitk.DisplacementFieldJacobianDeterminant(displacement)
error = sitk.Mask(error, mask)
```

A 3D rendering of the resulting displacement field:
![error](https://github.com/m-pilia/disptools/blob/master/sphinx/img/example_2.png?raw=true)

And a visualisation of the the error on the Jacobian:
![error](https://github.com/m-pilia/disptools/blob/master/sphinx/img/example_1.png?raw=true)

### References
+ <a id="1"></a>[1] van Eede, M. C., Scholz, J., Chakravarty, M. M., Henkelman, R. M., and Lerch, J. P. "Mapping registration sensitivity in MR mouse brain images." Neuroimage 82 (2013), 226–236.
+ <a id="2"></a>[2] Karaçali, B., and Davatzikos, C. "Estimating topology preserving and smooth displacement fields." IEEE Transactions on Medical Imaging 23, 7 (2004), 868–880.
+ <a id="3"></a>[3] Karaçali, B., and Davatzikos, C. "Simulation of tissue atrophy using a topology preserving transformation model." IEEE transactions on medical imaging 25, 5 (2006), 649–652.

### License

The software is distributed under the MIT license.
