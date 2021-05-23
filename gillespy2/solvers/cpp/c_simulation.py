import subprocess, os, signal, threading

from gillespy2.core import Model
from gillespy2.solvers.cpp.build.build_engine import BuildEngine

class CSimulation():
    def __init__(self, model: Model = None, output_directory: str = None, delete_directory: bool = False, variable: bool = True, debug: bool = False):
        self.build_engine = BuildEngine(debug, output_directory)
        self.build_engine.prepare(model, variable)

        self.executable = self.build_engine.build_simulation(self.type)

    def __exec_simulation(self, timeout: int = 0, *args):
        sim_args = [self.executable] + args

        if os.name == "nt":
            proc_kill = lambda sim: sim.send_signal(signal.CTRL_BREAK_EVENT)
            platform_args = {
                "creationargs": subprocess.CREATE_NEW_PROCESS_GROUP,
                "start_new_session": True
            }

        else:
            proc_kill = lambda sim: os.killpg(sim.pid, signal.SIGINT)
            platform_args = { "start_new_session": True }

        buffer = []
        timeout_event = False

        with subprocess.Popen(args, stdout=subprocess.PIPE, **platform_args) as simulation:
            def timeout_kill():
                timeout_event = True
                proc_kill(simulation)

            timeout_thread = threading.Timer(timeout, timeout_kill)

            try:
                read_thread = threading.Thread(name = "SimulationHandlerThread", target=self.__sim_delegate, args=(simulation, buffer))
                read_thread.start()

                if timeout > 0:
                    timeout_thread.start()

            except KeyboardInterrupt:
                proc_kill(simulation)
                pause = True

            finally:
                return_code = simulation.wait()
                read_thread.join()

                if timeout_thread.is_alive():
                    timeout_thread.cancel()

                simulation_output = "".join(buffer).split(",")

                if timeout_event:
                    pause = True
                    return_code = 33

    def __sim_delegate(self, sim, sim_buffer):
        def read_next():
            line = sim.stdout.read().decode("utf-8")
            length = len(line)

            if length > 0:
                sim_buffer.append(line)

            return length
        
        block_size = read_next()
        while block_size > 0 and sim.poll() is None:
            block_size = read_next()

    def clean(self):
        self.build_engine.clean()