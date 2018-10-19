from .gillespyError import *
class GillesPySolver():
    """ 
    Abstract class for a solver. This is generally called from within a
    gillespy Model through the Model.run function. Returns simulation 
    trajectories.
        
    Attributes
    ----------
    model : gillespy.Model
        The model on which the solver will operate.
    t : float
        The end time of the solver.
    number_of_trajectories : int
        The number of times to sample the chemical master equation. Each
        trajectory will be returned at the end of the simulation.
    increment : float
        The time step of the solution.
    seed : int
        The random seed for the simulation. Defaults to None.
    stochkit_home : str
        Path to stochkit. This is set automatically upon installation, but 
        may be overwritten if desired.
    algorithm : str
        The solver by which to simulate the model. 'ssa' or 'tau_leaping'
        are the available options. If 'ssa' is chosen, StochKit will choose
        from the available ssa options.
    job_id : str
        If given, this will be the name of the solver run. Usually not set.
    extra_args : str
        Any extra arguments for the stochkit solver. See StochKit2 
        documentation for details.
    debug : bool (False)
        Set to True to provide additional debug information about the     
        simulation.
    show_labels : bool (True)
        Use names of species as index of result object rather than position numbers.
    """

    def run(self, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, stochkit_home=None, algorithm=None,
            job_id=None, extra_args='', debug=False, show_labels=False):
        """ 
        Call out and run the solver. Collect the results.
        """
        
        if algorithm is None:
            raise SimuliationError("No algorithm selected")
        
        # We write all StochKit input and output files to a temporary folder
        prefix_basedir = tempfile.mkdtemp()
        prefix_outdir = os.path.join(prefix_basedir, 'output')
        os.mkdir(os.path.join(prefix_basedir, 'output'))
        
        if job_id is None:
            job_id = str(uuid.uuid4())
        
        # Write a temporary StochKit2 input file.
        if isinstance(model, Model):
            outfile =  os.path.join(prefix_basedir, 
                                        "temp_input_"+job_id+".xml")
            mfhandle = open(outfile, 'w')
            #document = StochMLDocument.from_model(model)

        # If the model is a Model instance, we serialize it to XML,
        # and if it is an XML file, we just make a copy.
        if isinstance(model, Model):
            document = model.serialize()
            mfhandle.write(document)
            mfhandle.close()
        elif isinstance(model, str):
            outfile = model

        # Assemble argument list for StochKit
        ensemblename = job_id
    
        directories = os.listdir(prefix_outdir)
        
        
        outdir = prefix_outdir+'/'+ensemblename
        

        # Algorithm, SSA or Tau-leaping?
        executable = None
        if stochkit_home is not None:
            if os.path.isfile(os.path.join(stochkit_home, algorithm)):
                executable = os.path.join(stochkit_home, algorithm)
            else:
                raise SimuliationError("stochkit executable '{0}' not found \
                stochkit_home={1}".format(algorithm, stochkit_home))
        elif os.environ.get('STOCHKIT_HOME') is not None:
            if os.path.isfile(os.path.join(os.environ.get('STOCHKIT_HOME'), 
                                           algorithm)):
                executable = os.path.join(os.environ.get('STOCHKIT_HOME'), 
                                          algorithm)
        if executable is None:
            # try to find the executable in the path
            if os.environ.get('PATH') is not None:
                for dir in os.environ.get('PATH').split(':'):
                    if os.path.isfile(os.path.join(dir, algorithm)):
                        executable = os.path.join(dir, algorithm)
                        break
        if executable is None:
            raise SimulationError("stochkit executable '{0}' not found. \
                Make sure it is your path, or set STOCHKIT_HOME envronment \
                variable'".format(algorithm))



        # Assemble the argument list
        args = ''
        args += '--model '
        args += outfile
        args += ' --out-dir '+outdir
        args += ' -t '
        args += str(t)
        if increment == None:
            increment = t/20.0
        num_output_points = str(int(float(t/increment)))
        args += ' -i ' + num_output_points
        if ensemblename in directories:
            print('Ensemble '+ensemblename+' already existed, using --force.')
            args+=' --force'
        

        # If we are using local mode, shell out and run StochKit 
        # (SSA or Tau-leaping or ODE)
        cmd = executable+' '+args+' '+extra_args
        if debug:
            print("cmd: {0}".format(cmd))

        # Execute
        try:
            #print "CMD: {0}".format(cmd)
            handle = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            return_code = handle.wait()
        except OSError as e:
            raise SimuliationError("Solver execution failed: \
            {0}\n{1}".format(cmd, e))
        
        try:
            stderr = handle.stderr.read()
        except Exception as e:
            stderr = 'Error reading stderr: {0}'.format(e)
        try:
            stdout = handle.stdout.read()
        except Exception as e:
            stdout = 'Error reading stdout: {0}'.format(e)

        if return_code != 0:
            #print stdout
            #print stderr
            raise SimuliationError("Solver execution failed: \
            '{0}' output: {1}{2}".format(cmd,stdout,stderr))

        # Get data using solver specific function
        try:
            if show_labels:
                labels, trajectories = self.get_trajectories(outdir, debug=debug, show_labels=True)
            else:
                trajectories = self.get_trajectories(outdir, debug=debug, show_labels=False)
        except Exception as e:
            fname = os.path.join(prefix_basedir,'temp_input_{0}_generated_code'.format(ensemblename),'compile-log.txt')
            if os.path.isfile(fname):
                with open(fname) as f:
                    cerr = f.read()
                raise SimulationError("Error compiling custom propensities: {0}\n{1}\n".format(fname,cerr))

            fname = os.path.join(prefix_outdir,ensemblename,'log.txt')
            if os.path.isfile(fname):
                with open(fname) as f:
                    cerr = f.read()
                raise SimulationError("Error running simulation: {0}\n{1}\n".format(fname,cerr))
            
            raise SimulationError("Error using solver.get_trajectories('{0}'): {1}".format(outdir, e))

        if len(trajectories) == 0:
            #print stdout
            #print stderr
            raise SimuliationError("Solver execution failed: \
            '{0}' output: {1}{2}".format(cmd,stdout,stderr))

        # Clean up
        if debug:
            print("prefix_basedir={0}".format(prefix_basedir))
            print("STDOUT: {0}".format(stdout))
            print("STDERR: {0}".format(stderr))
        else:
            shutil.rmtree(prefix_basedir)
        # Return data
        if show_labels:
            results2 = []
            for r in trajectories:
                ret = {}
                for n,l in enumerate(labels):
                    ret[l] = r[:,n]
                results2.append(ret)
            return results2
        else:
            return trajectories

class StochKitSolver(GillesPySolver):
    """ 
    Abstract class for StochKit solver derived from the GillesPySolver class.
    This is generally used to set up the solver.
        
    Attributes
    ----------
    model : gillespy.Model
        The model on which the solver will operate.
    t : float
        The end time of the solver.
    number_of_trajectories : int
        The number of times to sample the chemical master equation. Each
        trajectory will be returned at the end of the simulation.
    increment : float
        The time step of the solution.
    seed : int
        The random seed for the simulation. Defaults to None.
    stochkit_home : str
        Path to stochkit. This is set automatically upon installation, but 
        may be overwritten if desired.
    algorithm : str
        The solver by which to simulate the model. 'ssa' or 'tau_leaping'
        are the available options. If 'ssa' is chosen, StochKit will choose
        from the available ssa options.
    job_id : str
        If given, this will be the name of the solver run. Usually not set.
    method : str
        The specific SSA to call. NOT YET FUNCTIONAL.
    debug : bool (False)
        Set to True to provide additional debug information about the     
        simulation.
    """
    
    @classmethod
    def run(cls, model, t=20, number_of_trajectories=1,
            increment=0.05, seed=None, stochkit_home=None, algorithm='ssa',
            job_id=None, method=None,debug=False, show_labels=False):
    
        # all this is specific to StochKit
        if model.units == "concentration":
            raise SimuliationError("StochKit can only simulate population "+
                "models, please convert to population-based model for "+
                "stochastic simulation. Use solver = StochKitODESolver "+
                "instead to simulate a concentration model deterministically.")

        if seed is None:
            seed = random.randint(0, 2147483647)
        # StochKit breaks for long ints
        if seed.bit_length()>=32:
            seed = seed & ((1<<32)-1)
            if seed > (1 << 31) -1:
                seed -= 1 << 32

        # Only use on processor per StochKit job.
        args = ' -p 1'
      
        # We keep all the trajectories by default.
        args += ' --keep-trajectories'
        args += ' --label'

        args += ' --seed '
        args += str(seed)
        
        realizations = number_of_trajectories
        args += ' --realizations '
        args += str(realizations)

        if method is not None:  #This only works for StochKit 2.1
            args += ' --method ' + str(method)

        
        self = StochKitSolver()
        return GillesPySolver.run(self, model,t, number_of_trajectories, 
                                  increment, seed, stochkit_home,
                                  algorithm, 
                                  job_id, extra_args=args, debug=debug,
                                  show_labels=show_labels)


    def get_trajectories(self, outdir, debug=False, show_labels=False):
        # Collect all the output data
        files = os.listdir(outdir + '/stats')
        trajectories = []
        files = os.listdir(outdir + '/trajectories')
        labels = []
        if show_labels:
            with open(outdir + '/trajectories/trajectory0.txt', 'r') as f:
                first_line= f.readline()
                labels = first_line.split()
        for filename in files:
            if 'trajectory' in filename:
                trajectories.append(numpy.loadtxt(outdir + '/trajectories/' + 
                                        filename, skiprows=1))
            else:
                raise SimuliationError("Couldn't identify file '{0}' found in \
                                        output folder".format(filename))
        if show_labels:
            return (labels, trajectories)
        else:
            return trajectories


class StochKitODESolver(GillesPySolver):
    """ 
    Abstract class for StochKit solver derived from the GillesPySolver class.
    This is generally used to set up the solver.
        
    Attributes
    ----------
    model : gillespy.Model
        The model on which the solver will operate.
    t : float
        The end time of the solver.
    number_of_trajectories : int
        The number of times to sample the chemical master equation. Each
        trajectory will be returned at the end of the simulation.
    increment : float
        The time step of the solution.
    seed : int
        The random seed for the simulation. Defaults to None.
    stochkit_home : str
        Path to stochkit. This is set automatically upon installation, but 
        may be overwritten if desired.
    algorithm : str
        Already set to 'stochkit_ode.py'
    job_id : str
        If given, this will be the name of the solver run. Usually not set.
    debug : bool (False)
        Set to True to provide additional debug information about the     
        simulation.
    """
    
    @classmethod
    def run(cls, model, t=20, number_of_trajectories=1,
                increment=0.05, seed=None, stochkit_home=None, 
                algorithm='stochkit_ode.py',
                job_id=None, debug=False, show_labels=False):
        self = StochKitODESolver()
        return GillesPySolver.run(self,model,t, number_of_trajectories, 
                                  increment, seed, stochkit_home,
                                  algorithm, 
                                  job_id, debug=debug,
                                  show_labels=show_labels)

    def get_trajectories(self, outdir, debug=False, show_labels=False):
        if debug:
            print("StochKitODESolver.get_trajectories(outdir={0}".format(outdir))
        # Collect all the output data
        trajectories = []
        with open(outdir + '/output.txt') as fd:
            fd.readline()
            headers = fd.readline()
            fd.readline()
            data = []
            data.append([float(x) for x in fd.readline().split()])
            fd.readline()
            for line in fd:
                data.append([float(x) for x in line.split()])
        trajectories.append(numpy.array(data))
        if show_labels:
            return (headers.split(), trajectories)
        else:
            return trajectories
