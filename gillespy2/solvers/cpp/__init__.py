from gillespy2.solvers.cpp.ssa_c_solver import SSACSolver

def check_cpp_support():
    from gillespy2.example_models import Example
    try:
        model = Example()
        results = model.run(solver=SSACSolver)
        return True
    except:
        return False

can_use_cpp = check_cpp_support()

__all__ = ['SSACSolver'] if can_use_cpp else []
