import numpy as np
import os
import re
import SimpleITK as sitk
from typing import *
from disptools import *

try:
    import vtk
except ImportError as e:
    print("Warning: cannot import 'vtk' module. " +
          "Some functionalities depending upon it may be unavailable.")

try:
    import ply
    import ply.lex as lex
    import ply.yacc as yacc
except ImportError as e:
    print("Warning: cannot import 'ply' module. " +
          "Some functionalities depending upon it may be unavailable.")


class _ElastixParametersLexer():
    r""" Lexer for Elastix parameter files.
    """

    def __init__(self, **kwargs):
        self.lexer = lex.lex(module=self, **kwargs)

    def test(self, data):
        self.lexer.input(data)
        while True:
             tok = self.lexer.token()
             if not tok:
                 break
             print(tok)

    t_ignore = " \t"

    tokens = (
        'LPAREN',
        'RPAREN',
        'IDENTIFIER',
        'STRING',
        'NUMBER',
        'COMMENT'
    )

    t_LPAREN     = r'\('
    t_RPAREN     = r'\)'
    t_IDENTIFIER = r'[a-zA-Z_][a-zA-Z0-9_]*'

    def t_NUMBER(self, t):
        r'[+-]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?'
        t.value = float(t.value)
        if t.value.is_integer():
            t.value = int(t.value)
        return t

    def t_STRING(self, t):
        r'"[^"]*"'
        t.value = t.value[1:-1]
        if t.value == 'true':
            t.value = True
        elif t.value == 'false':
            t.value = False
        return t

    def t_error(self, t):
        raise Exception("Lexer error while reading parameter file at line %d, token '%s'"
                        % (t.lexer.lineno, t.value[0]))

    def t_COMMENT(self, t):
        r'//.*'
        pass

    def t_newline(self, t):
        r'\n+'
        t.lexer.lineno += t.value.count("\n")


class _ElastixParametersParser():
    r""" Parser object for Elastix parameter files.

    The `parse(text)` method accepts as input the content of an Elastix
    parameter file, and returns a dictionary containing couples of
    parameters with their respective values.
    """

    def __init__(self, lexer_opts={}, **kwargs):
        self.names = {}
        self.lexer = _ElastixParametersLexer(**lexer_opts)
        self.tokens = self.lexer.tokens
        self.parser = yacc.yacc(module=self, **kwargs)

    def parse(self, text):
        self.names = {}
        self.parser.parse(text, lexer=self.lexer.lexer)
        return self.names

    precedence = ()

    def p_program(self, p):
        '''program : program statement
                   | statement'''
        pass

    def p_statement_assign(self, p):
        '''statement : LPAREN IDENTIFIER rvaluelist RPAREN'''
        self.names[p[2]] = p[3]

    def p_statement_comment(self, p):
        '''statement : COMMENT'''
        pass

    def p_expression_rvalue(self, p):
        '''rvalue : NUMBER
                  | STRING'''
        p[0] = p[1]

    def p_expression_rvaluelist_2(self, p):
        '''rvaluelist : rvalue
                      | rvaluelist rvalue'''
        if len(p) == 2:
            p[0] = [p[1]]
        else:
            p[0] = p[1] + [p[2]]

    def p_error(self, p):
        raise Exception("Syntax error while reading parameter file at '%s'" % p.value)


def make_unique_directory(path: str) -> str:
    r""" Create a unique directory.

    If a directory with the given name already exists,
    a suffix in the form `_x` with `x` integer will be
    added at the end of its name.

    Parameters
    ----------
    path : str
        Relative or absolute path for the new directory.

    Returns
    -------
    str
        A string containing the path of the new directory.
    """

    if not os.path.exists(path):
        os.makedirs(path)
        return path

    base_path, name = os.path.split(path)
    name, extension = os.path.splitext(name)

    i = 1
    while True:
        path = os.path.join(base_path, '%s_%d%s' % (name, i, extension))
        if not os.path.exists(path):
            os.makedirs(path)
            return path
        i += 1


def read_rvf(filename: str) -> sitk.Image:
    r""" Read an RVF file.

    Read image data from an RVF file.

    Parameters
    ----------
    filename : str
        Filename for the input vector field.

    Returns
    -------
    result : sitk.Image
        Vector field.
    """

    with open(filename, 'rb') as f:
        size = [int(i) for i in f.readline().split()]
        spacing = [float(i) for i in f.readline().split()]
        data = np.fromfile(f, dtype=np_float_type)

    image = sitk.GetImageFromArray(data.reshape((*size, 3), order='C'))
    image.SetSpacing(tuple(spacing))

    return image


def write_vtk_points(points: np.ndarray, filename: str) -> None:
    r""" Write a set of points to a VTK PolyData file.

    The points are given as a numpy bidimensional array, with a
    row for each point, and three columns, one per component.

    Parameters
    ----------
    points : np.ndarray
        A `n × m` array containing `n` points with `m` components.
    filename : str
        Output file.
    """

    if 'vtk' not in sys.modules:
        raise Exception('write_vtk_points: vtk module is required to use this feature.')

    vtk_points = vtk.vtkPoints()

    for i in range(0, points.shape[0]):
        vtk_points.InsertNextPoint(points[i,0], points[i,1], points[i,2])

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(vtk_points)

    writer = vtk.vtkPolyDataWriter()
    writer.SetInputData(polydata)
    writer.SetFileName(filename)
    writer.Update()


def read_elastix_points(filename: str) -> np.ndarray:
    r""" Read a set of points from a file in Elastix format.

    The points are returned as a two-dimensional numpy array, with a
    row for each point and three columns, one per coordinate.

    Parameters
    ----------
    filename : str
        Output file.

    Returns
    -------
    np.ndarray
        A `n × m` array containing `n` points with `m` components.
    """

    with open(filename, 'r') as f:

        line = f.readline()
        if not re.match(r'\s*point|index\s*', line):
            raise Exception('Invalid point file format')

        line = f.readline()
        try:
            n = int(re.match(r'([0-9]+)', line).group(1))
        except:
            raise Exception('Invalid point file format')

        points = np.empty([n, 3])

        i = 0
        for line in f:
            try:
                m = re.match(r'([0-9.e+-]+)\s+([0-9.e+-]+)\s+([0-9.e+-]+)', line)
                points[i,:] = [float(x) for x in m.groups()]
            except:
                raise Exception('Invalid point file format')

            i += 1

        return points


def write_elastix_points(
        points: np.ndarray,
        filename: str,
        point_format: str = 'point'
        ) -> None:
    r""" Write a set of points to a file in Elastix format.

    The points are passed as a two-dimensional numpy array, with a row
    for each point and three columns, one per coordinate.

    Parameters
    ----------
    points : np.ndarray
        A `n × m` array containing `n` points with `m` components.
    filename : str
        Output file.
    point_format : str
        One of 'point' (default) or 'index'.
    """

    if point_format not in ['point', 'index']:
        raise Exception('Unsupported point format %s' % point_format)

    if len(points.shape) != 2 and points.shape[1] != 3:
        raise Exception('Invalid point array')

    with open(filename, 'w') as f:
        f.write('%s\n' % point_format)
        f.write('%d\n' % points.shape[0])
        for i in range(0, points.shape[0]):
            f.write('%.10e\t%.10e\t%.10e\n' % (points[i,0], points[i,1], points[i,2]))


def read_elastix_parameters(filename: str) -> dict:
    r""" Read an Elastix parameter file.

    Read an Elastix parameter file, returning a dictionary mapping each
    parameter name to a list of values. Booleans are converted to Python
    `bool`s, same for `int` and `float` values. Strings do not need to be
    quoted. Comments are discarded.

    .. warning::
        This function will not work if Python is allowed to discard docstrings
        (e.g. due to the option -OO).


    Parameters
    ----------
    filename : str
        Name of the parameter file to be read.

    Returns
    -------
    dict
        Dictionary of parameters.
    """

    if 'ply' not in sys.modules:
        raise Exception('read_elastix_parameters: ply module is required to use this feature.')

    if not hasattr(read_elastix_parameters, "parser"):
        try:
            read_elastix_parameters.parser = _ElastixParametersParser(debug=False, write_tables=False)
        except SyntaxError:
            raise Exception('read_elastix_parameters: unable to build the lexer. '
                            'Please make sure you are not running Python with '
                            'docstring discarding enabled (-OO).')

    with open(filename, 'r') as f:
        return read_elastix_parameters.parser.parse(f.read())


def write_elastix_parameters(parameters: dict, filename: str) -> None:
    r""" Write Elastix parameters to file.

    Write an Elastix parameter file. Parameters are passed as a dictionary
    mapping each parameter name to a list of values. Python `bool`s, `int`s,
    and `float`s are converted automatically.

    Parameters
    ----------
    parameters : dict
        Dictionary of parameters.

    filename : str
        Name of the parameter file to be read.
    """
    def _format(p):
        if isinstance(p, bool):
            return '"true"' if p else '"false"'
        elif isinstance(p, str):
            return '"%s"' % p
        else:
            return '%s' % p

    with open(filename, 'w') as f:
        for k, v in parameters.items():
            f.write('(%s %s)\n' % (k, ' '.join([_format(p) for p in v])))

