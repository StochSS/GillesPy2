from gillespy2.solvers.cpp.ssa_c_solver import SSACSolver

def check_cpp_support():
    import example_models
    try:
        model = Example()
        results = model.run(solver=SSACSolver)
        return True
    except:
        return False

can_use_cpp = check_cpp_support()

__all__ = ['SSACSolver'] if can_use_cpp else []
