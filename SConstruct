import sysconfig

sys_info_include = sysconfig.get_path('include')

env = Environment(
    CPPPATH=[
        '#/lib/cgillespy/src',
        '#/lib/cgillespy/src/template',
        sys_info_include,
    ],
    CXXFLAGS=['-std=c++14', '-fPIC'],
)

SConscript([
    'lib/cgillespy/src/ssa_cpp_solver/SConscript',
    'lib/cgillespy/src/SConscript',
], exports='env')

env.SharedLibrary('cgillespy', [
    'lib/cgillespy/src/extension.cpp',
    'lib/cgillespy/src/model_context.cpp',
    'lib/cgillespy/src/template/template.cpp',
])
