# ======================================================================
#         Import modules
# ======================================================================
import os
import argparse
import shutil
from mpi4py import MPI
from adflow import ADFLOW
import inspect
import time
from multipoint import multiPointSparse
from baseclasses.utils import redirectIO
from pprint import pprint as pp

from SETUP import setup_problem

parser = argparse.ArgumentParser()
parser.add_argument("--output", type=str, default="debug", help="This is the solution directory.")
parser.add_argument("--outDir", type=str, default="debug", help="This is the base output directory. Run case directories are then created here.")
parser.add_argument("--geometry", choices=["G1", "G2"], default="G2")
parser.add_argument("--gridFamily", choices=["A", "B"], default="A")
parser.add_argument("--gridLevel", choices=["L0", "L1", "L2"], default="L0")
parser.add_argument("--mach", type=float, default=0.05, help="[default: %(default)s]")
parser.add_argument("--alpha", type=float, default=0.0, help="[default: %(default)s]")
args = parser.parse_args()

comm = MPI.COMM_WORLD

# Reference all files in a dictionary so we can auto-save them
files = {}
files["self"] = inspect.getfile(inspect.currentframe())
files["gridFile"] = "circle_rans_129x65.cgns"
#"cyl.cgns"#f"./INPUT/mesh/OAT15A_{args.geometry}_{args.gridFamily}_{args.gridLevel}.cgns"

# Set output base directory
outDir = os.path.join(args.outDir, f"m{args.mach}_a{args.alpha}", args.output)
solDir = os.path.join(outDir, "output")

# Make backup of what was actually run
if comm.rank == 0:
    # Create a folder called input in the output directory to save all input files.
    copyDir = os.path.join(outDir, "INPUT")
    os.system(f"mkdir -p {copyDir}")
    for key in files:
        shutil.copy(files[key], copyDir)

    # Create a folder of the output files
    os.system(f"mkdir -p {solDir}")

comm.barrier()
'''
#Redirect STDOUT
if comm.rank == 0:
    #fName = os.path.join(outDir, "%s_%d.out"%(setName, ptID))
    fName = os.path.join(outDir, "std.out")
    outFile = open(fName, "w")
    redirectIO(outFile)
    #sys.stdout = outFile
'''
if comm.rank == 0:
    print("This script was executed: ", time.strftime("%Y-%m-%d %H:%M"))

    print("+------------------------------------------------+")
    print("+            Command Line Options                +")
    print("+------------------------------------------------+")
    pp(vars(args))

# ======================================================================
#         Create multipoint communication object
# ======================================================================
MP = multiPointSparse(MPI.COMM_WORLD)
MP.addProcessorSet("cruise", nMembers=1, memberSizes=MPI.COMM_WORLD.size)
comm, setComm, setFlags, groupFlags, ptID = MP.createCommunicators()


# ======================================================================
#         Setup problem
# ======================================================================
ap = setup_problem.setup(args)

# Echo the various options:
if comm.rank == 0:
    print("+------------------------------------------------+")
    print("|            AeroProblem Options                 |")
    print("+------------------------------------------------+")
    pp(vars(ap))


# ======================================================================
#         ADflow Set-up
# ======================================================================


aeroOptions = {
    # Common Parameters
    "gridFile": files["gridFile"],
    "outputDirectory": solDir,
    # Physics Parameters
    "nCycles": 22000,
    "monitorvariables": ["resrho", "cl", "resTurb"],
    "useNKSolver": True,
    "useANKsolver": True,

    "writeTecplotSurfaceSolution": False,
    "writeSurfaceSolution": True,
    "writeVolumeSolution": True,
    "L2ConvergenceCoarse": 1e-04,
    # Adjoint Parameters
    "adjointL2Convergence": 1e-11,
    "L2Convergence": 1e-14, #6E-11 use it for Re 60, Ma 0.05 for cylinder to get match. fine mesh 385 by 385
    "nsubiterturb": 15,
    "equationType": "laminar NS",#"RANS",
    #"turbulenceModel": "SA",
    "ANKSwitchTol" : 1e+2,
    "NKSwitchTol" : 1e-9,
    #"turbResScale" : 1.0,
    "CFL" : 2.0,
    #"mgCycle" : "sg",
    "discretization": "central plus matrix dissipation",
    #"restartFile": ".cgns",
}


# Update options based on grid config
if args.geometry == "G2" and args.gridFamily == "B":
    specialAeroOptions = {
        "nsubiterturb": 15,
        #"ANKUseTurbDADI": False,
        #"ANKNSubiterTurb": 1,
        # "ANKTurbKSPDebug": False,
        #"ANKCoupledSwitchTol": 1e-5,

    }
    aeroOptions.update(specialAeroOptions)



# Create solver
CFDSolver = ADFLOW(options=aeroOptions, comm=comm)

# Add slices
CFDSolver.addSlices("z",[0.5])
CFDSolver.addLiftDistribution(2,"z")
CFDSolver(ap)

# Evaluate functions
funcs = {}
CFDSolver.evalFunctions(ap, funcs)
CFDSolver.checkSolutionFailure(ap, funcs)

# Print the evaluated functions
if comm.rank == 0:
    print(funcs)


