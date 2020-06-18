import math
import numpy
numpy.sin

def __piecewise(*args):
    # Eval entry for piecewise functions
    args = list(args)
    sol = None
    if len(args) % 2: args.append(True)
    for i, arg in enumerate(args):
        if not i % 2: continue
        if arg:
            sol = args[i-1]
            break
    return int(sol)


def __xor(*args):
    # Eval entry for MathML xor function
    args = list(args)
    if not len(args):
        return False
    if len(args) == 1:
        return args
    from operator import ixor
    from functools import reduce
    return reduce(ixor, args)


# Next 4 functions, see test case 01112
def __plus(*args):
    # Eval entry for EMPTY plus
    args = list(args)
    if not len(args):
        return 0
    if len(args) == 1:
        return args[0]

def __times(*args):
    # Eval entry for EMPTY times
    args = list(args)
    if not len(args):
        return 1
    if len(args) == 1:
        return args[0]

def __and(*args):
    # Eval entry for EMPTY and
    args = list(args)
    if not len(args):
        return True
    if len(args) == 1:
        if args[0] == 'true':
            return int(True)
        if args[0] == 'false':
            return int(False)
        return args[0]

def __or(*args):
    args = list(args)
    if len(args) == 1:
        return args[0]


def __quotient(a, b):
    # Eval entry for quotient

    return a // b


def __rem(a, b):
    # Eval entry for remainder

    return a % b


def __implies(a, b):
    # Eval entry for implies

    return (not(a) or b)


def __max(*args):
    # Eval entry for max

    args = list(args)
    return max(args)


def __min(*args):
    # Eval entry for min

    args = list(args)
    return min(args)


def __sec(cos):
    # Eval entry for sec

    return numpy.arccos(cos)


def __ln(a):
    # Eval entry for ln

    return numpy.log(a)


def __csc(x):
    # Eval entry for cosecant, input in radians.
    return 1/math.sin(x)


def __cot(x):
    # Eval entry for cotangent, input in radians.
    return 1/math.tan(x)


def __arcsec(x):
    # Eval entry for arcsecant, input in radians.
    return math.acos(1/x)


def __arccsc(x):
    # Eval entry for arccosecant, input in radians
    return __csc(1/x)

def __arccot(x):
    # Eval entry for arccotangent, input in radians
    if x == 0:
        return numpy.pi/2
    return __cot(1/x)


def __arcsinh(x):
    return numpy.arcsinh(x)


def __arccosh(x):
    return numpy.arccosh(x)


def __arctanh(x):
    return numpy.arctanh(x)


def __arcsech(x):
    if x == 0:
        return numpy.inf
    return numpy.arccosh(1/x)

def __arccsch(x):
   return __arccsc(1/x)

def __arccoth(x):
    if x == 0:
        return numpy.inf

eval_globals = math.__dict__
eval_globals['false'] = False
eval_globals['true'] = True
eval_globals['piecewise'] = __piecewise
eval_globals['xor'] = __xor
eval_globals['quotient'] = __quotient
eval_globals['rem'] = __rem
eval_globals['implies'] = __implies
eval_globals['max'] = __max
eval_globals['min'] = __min
eval_globals['sec'] = __sec
eval_globals['__and'] = __and
eval_globals['plus'] = __plus
eval_globals['times'] = __times
eval_globals['or'] = __or
eval_globals['exponentiale'] = numpy.e
eval_globals['avogadro'] = 6.02214179*10**23
eval_globals['ln'] = __ln
eval_globals['csc'] = __csc
eval_globals['cot'] = __cot
eval_globals['arcsec'] = __arcsec
eval_globals['arccsc'] = __arccsc
eval_globals['arccot'] = __arccot
eval_globals['arcsinh'] = __arcsinh
eval_globals['arccosh'] = __arccosh
eval_globals['arctanh'] = __arctanh
eval_globals['arcsech'] = __arcsech
eval_globals['arccsch'] = __arccsch
eval_globals['arccoth'] = __arccoth
#eval_globals['time'] = 0
eval_globals['true'] = True
eval_globals['false'] = False

