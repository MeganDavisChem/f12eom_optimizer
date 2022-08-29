"""Main module."""
import qcengine as qcng
import qcelemental as qcel
import psi4
import optking
from optking import molsys, stepAlgorithms, hessian, displace, intcosMisc, history, optwrapper, convCheck
from optking import optparams as op
import subprocess
import numpy as np


#TODO may have to feed in a Zmat to molpro. <-- can do this for efficiency but i'd rather not
#TODO use logger and consolidate all outputs to one file, cleanup molpro cruft. Probably don't even need to write
#molpro input files now that I think about it, just use stdin or python or something
#rewrite these to use a python dictionary


#TODO remove default values? Could lead to confusing bugs
optparams = {
        'maxiter': 50,
        'natoms': 3,
        'charge': 0,
        'spin': 0,
        'psispin': 1,
        'memory': 8,
        'mem': 2,
        'nproc': 4,
        'roots': [0, 1],
        'theory': 'dz',
        'units': 'bohr',
        'freeze_core': 'false',
        'xyz': '',
        'mol_core': '',
        'mol_basis': '',
        'psi_basis': ''
        }
input_file = 'opt.inp'

maxiter = 50
natoms = 3
charge = 0
spin = 0
memory = 8
nproc = 4
mem = memory*1000/nproc
roots = [0, 1]

initial_geom="""
  H      0.000000000    1.442306288   -2.754574215
  C      0.000000000   -0.186314094   -1.423042607
  F      0.000000000    0.041269805    1.045804233
"""

#this should be tz or dz and allow switching between the two
theory = "dz"
#bohr or angstrom
units = 'bohr'

#TODO move this back to the run molpro function and/or fix global hacky stuff
molpro_input = """*** part of geometry opt run
memory, {mem}, m;
gthresh,energy=1.d-12,zero=1.d-22,oneint=1.d-22,twoint=1.d-22;
gthresh,optgrad=1.d-8,optstep=1.d-8;
nocompress;
ang;
geometry={{
     {natoms}
     {xyz}
     }}

set,charge={charge}
set,spin={spin}
basis=sto-3g
hf,maxit=500;accu,20;
basis={mol_basis}
hf,maxit=500;accu,20;
ccsd(t)-f12,maxit=250;{mol_core}orbital,IGNORE_ERROR;
ee(1)=energy(2)
forces,varsav
SHOW[f18.12],ee(1)
SHOW[f18.12],GRADX
SHOW[f18.12],GRADY
SHOW[f18.12],GRADZ
"""

def read_input(_input):
    def strip_line(_line):
        """helper function"""
        stripped_line = line.split('=')
        stripped_arg = stripped_line[-1].split('\n')[0].strip()
        return stripped_arg
    print('hello world')

    #this is stupid

    with open(_input, 'r') as f:
        input_lines = f.readlines()
    print(input_lines)

    for i,line in enumerate(input_lines):
        for key in optparams.keys():
            if key in line:
                optparams[key] = strip_line(line)
        if "STARTGEOM" in line:
            startgeom = i
        elif "ENDGEOM" in line:
            endgeom = i


    #convert root from a string to an actual list there is surely a better way to do this
    optparams['natoms'] = input_lines[startgeom+1].strip()
    optparams['roots'] = [int(root.strip()) for root in optparams['roots'].strip('[').strip(']').split(',')]
    optparams['psispin'] = str(int(optparams['spin'])+1)
    initial_geom = ''.join(input_lines[startgeom+3:endgeom])
    optparams['xyz'] = initial_geom
    optparams['mem'] = int(optparams['memory'])*1000-50
    #shouldn't actually do this cause it needs to be converted first
    #actually yeah do that in this function

    #converts geometry to right unit I know this is dumb
    mol = psi4.geometry("""
    {charge} {psispin}
    {xyz}
    units {units}
         """.format_map(optparams))
    molsys_obj = molsys.Molsys.from_psi4_molecule(mol)
    optparams['xyz'] = molsys_obj.show_geom()

    if optparams['freeze_core'] == 'true':
        optparams['mol_basis'] = 'v{theory}-f12'.format_map(optparams)
        optparams['psi_theory'] = 'aug-cc-pv{theory}'.format_map(optparams)
    elif optparams['freeze_core'] == 'false':
        optparams['mol_basis'] = 'cc-pcv{theory}-f12'.format_map(optparams)
        optparams['mol_core'] = 'core;'
        optparams['psi_theory'] = 'aug-cc-pcv{theory}'.format_map(optparams)

    return optparams, molsys_obj
                




def run_psi4(_mol):
    """runs psi4 using qcengine, for troubleshooting but should be able to delete"""
    inpE = qcel.models.AtomicInput(
            molecule=_mol,
            driver="energy",
            model={"method": "SCF", "basis": "aug-cc-pvdz"},
            keywords={'reference': 'rohf', 
        'max_force_g_convergence': '1.0e-8', 
        'max_energy_g_convergence': '1.0e-12',
        'max_disp_g_convergence': '1.0e-8',
        'basis': f'aug-cc-pv{theory}',
        'freeze_core': 'true',
        'cachelevel': '0',
        'maxiter': '500',
        'dertype': 'none',
        'ints_tolerance': '20',
        'roots_per_irrep': roots,
        'scf__d_convergence': '13',
        'ccenergy__r_convergence': '15',
        'cceom__r_convergence': '8'
        })
    inpG = qcel.models.AtomicInput(
            molecule=_mol,
            driver="gradient",
            model={"method": "SCF", "basis": "aug-cc-pvdz"},
            keywords={'reference': 'rohf', 
        'max_force_g_convergence': '1.0e-8', 
        'max_energy_g_convergence': '1.0e-12',
        'max_disp_g_convergence': '1.0e-8',
        'basis': f'aug-cc-pv{theory}',
        'freeze_core': 'true',
        'cachelevel': '0',
        'maxiter': '500',
        'dertype': 'none',
        'ints_tolerance': '20',
        'roots_per_irrep': roots,
        'scf__d_convergence': '13',
        'ccenergy__r_convergence': '15',
        'cceom__r_convergence': '8'
        })
    
    #compute energy and gradient
    retE = qcng.compute(inpE, "psi4")
    retG = qcng.compute(inpG, "psi4")
    
    E = retE.return_result
    gradient = retG.return_result.flatten()
    print("Here is the qcengine gradient")
    print(gradient)
    return gradient, E

def real_run_psi4(_molsys_obj, _optparams):
    """Runs an EOM-CCSD/aug-cc-pv(d,t)z calculation in PSI4
   
    Returns an EOM-CCSD excitation energy gradient and EOM-CCSD XS energy
    """
    psi4.set_memory('{mem} MB'.format_map(_optparams))
    #might need to fix this
    psi4.core.set_num_threads(int(_optparams['nproc']))
    #formats this for a psi4 geomtetry

    mol = psi4.geometry(_molsys_obj.show_geom().replace("\tFragment 1\n\t", "") + f"\nunits ang")
    

    
    #compute energy and gradient
    psi4.set_options({'reference': 'rohf', 
        'max_force_g_convergence': '1.0e-8', 
        'max_energy_g_convergence': '1.0e-12',
        'max_disp_g_convergence': '1.0e-8',
        'basis': _optparams['psi_theory'],
        'freeze_core': _optparams['freeze_core'],
        'cachelevel': '0',
        'maxiter': '500',
        'dertype': 'none',
        'ints_tolerance': '20',
        'roots_per_irrep': _optparams['roots'],
        'scf__d_convergence': '13',
        'ccenergy__r_convergence': '15',
        'cceom__r_convergence': '8'
        })
    ccsd_gradient = psi4.gradient('ccsd', dertype='energy').np.flatten()
    ##    gradient,wfn = psi4.gradient('eom-ccsd', dertype='energy', return_wfn='true').np.flatten()
    eom_gradient,wfn = psi4.gradient('eom-ccsd', dertype='energy', return_wfn='true')
    #
    eom_gradient = eom_gradient.np.flatten()
    #E = wfn.variables()['CURRENT ENERGY']
    gradient = np.subtract(eom_gradient, ccsd_gradient)
    E = wfn.variables()['CURRENT ENERGY'] - wfn.variables()['CCSD TOTAL ENERGY']

#    gradient = psi4.gradient('scf').np.flatten()
#    E = psi4.energy('scf')
#    E = psi4.energy('eom-ccsd')
    print("Here is the gradient")
    print(gradient)
    return gradient, E

def run_molpro_v2(_molpro_input, _molpro_keys):
    """Runs molpro with a given input. Will make this function run without writing files"""
    with open('molpro_step.inp', 'w') as molpro_file:
        molpro_file.write(_molpro_input.format_map(_molpro_keys))
    subprocess.Popen(["/home/qc/bin/molpro2020.sh", "1", 
        "{nproc}".format_map(_molpro_keys), "molpro_step.inp"]).wait()


def run_molpro(_molpro_input):
    """Runs molpro with input given above (i.e. a CCSD(T)-F12b/v(d,t)z-f12 calculation).
    This structure should be completely redone (when I have time) to
    not generate useless molpro files"""
    with open('molpro_step.inp', 'w') as molpro_file:
        molpro_file.write(_molpro_input.format_map(globals()))
    #TODO let this run parallel with psi4 job? Or remove reliance on the tsch script
    subprocess.Popen(["/home/qc/bin/molpro2020.sh", f"{nproc}", "1", "molpro_step.inp"]).wait()

def read_molpro():
    """parses molpro output for gradients and energy"""
    def grad_to_array(_line):
        """helper function"""
        char1 = '['
        char2 = ']'
        str_array = _line[_line.find(char1)+1 : line.find(char2)].split()
        float_array = [float(grad) for grad in str_array]
        return float_array

    with open('molpro_step.out', 'r') as molpro_output:
        molpro_output = molpro_output.readlines()
    
    for line in molpro_output:
        if 'GRADX(' in line:
            xline = grad_to_array(line)
        elif 'GRADY(' in line:
            yline = grad_to_array(line)
        elif 'GRADZ(' in line:
            zline = grad_to_array(line)
        elif 'EE(1)' in line:
            energy = line
    #hacky, do here to grab the variable from end of file
    energy = float(energy.split()[2])
    molpro_gradients = np.empty(len(xline)*3)
    #convert to flattened array in correct order for optking
    for i in range(len(xline)):
        molpro_gradients[i*3] = xline[i]
        molpro_gradients[i*3+1] = yline[i]
        molpro_gradients[i*3+2] = zline[i]
    return molpro_gradients, energy

def main():
    optwrapper.initialize_options({"g_convergence": "gau_verytight"})
    converged = False

    #TODO change this to use the click interface or something
    optparams, molsys_obj = read_input(input_file)
    optking.make_internal_coords(molsys_obj)
    energies = []
    for i in range(int(maxiter)):
        print("iteration: " + str(i))

        psi_gradient, psi_energy = real_run_psi4(molsys_obj, optparams)

        run_molpro_v2(molpro_input, optparams)
        mol_gradient, mol_energy = read_molpro()

        gradient = np.add(psi_gradient,mol_gradient)
        E = mol_energy + psi_energy

        energies.append(E)


        #optking stuff
        H = hessian.guess(molsys_obj, guessType=op.Params.intrafrag_hess)

        #this converts the cartesian gradients to internal coordinate forces or something
        fq = intcosMisc.q_forces(molsys_obj.intcos, molsys_obj.geom, gradient)
        
        #I don't really know what these do but they make optimization work good
        intcosMisc.apply_fixed_forces(molsys_obj, fq, H, i)
        intcosMisc.project_redundancies_and_constraints(molsys_obj, fq, H)
        molsys_obj.q_show()

        #append to optimization history
        history.oHistory.append(molsys_obj.geom, E, fq)

        #simplified backstep logic TODO see if I need to add anything to this
        lastStepOK = history.oHistory.current_step_report()
        if len(history.oHistory.steps) < 5:
            pass
        #take a backstep if we made a bad step. TODO soft-code the 5 here and see if it's a reasonable value
        elif history.History.consecutiveBacksteps < 5 and not lastStepOK:
            history.History.consecutiveBacksteps += 1
            stepAlgorithms.take_step(molsys_obj, None, fq, H, stepType="BACKSTEP")
            continue
        # None is energy, not needed for any algorithms except for the linesearch
        #TODO consider allowing algorithm control uptop
        step = stepAlgorithms.take_step(molsys_obj, None, fq, H, stepType='RFO')
        #update hessian (idk if i actually need this)
        history.oHistory.hessian_update(H, molsys_obj)

        #TODO again, switch to using logger
        print("Here is the updated geometry")
        print(molsys_obj.show_geom()) 
        optparams['xyz'] = molsys_obj.show_geom()
    
        #check convergence
        converged = convCheck.conv_check(i, molsys_obj, step, fq, energies)
        if converged:
            break
    

def main_old():
    #set convergence options TODO add more of these, possibly up top? (yes up top)
    optwrapper.initialize_options({"g_convergence": "gau_verytight"})
    converged = False

    #define geo as psi4 ZMAT input
    h2o = psi4.geometry(f"""
    {charge} {spin+1}
  {initial_geom}
  units {units}
         """)
    #convert to molsys object
    molsys_obj = molsys.Molsys.from_psi4_molecule(h2o)
    optking.make_internal_coords(molsys_obj)
    
    energies = []
    for i in range(int(maxiter)):
        #TODO redundant but yeah change this to use the python logger
        print("iteration: " + str(i))
        #calculate psi4 stuff
        #TODO change these function names around
        psi_gradient, psi_energy = real_run_psi4(molsys_obj)
        #calculate psi4 stuff in qcengine
        #psi_gradient, psi_energy = run_psi4(molsys_obj.molsys_to_qc_molecule())
    
        #calculate molpro stuff
        global xyz 
        xyz = molsys_obj.show_geom()
        run_molpro(molpro_input)
        mol_gradient, mol_energy = read_molpro()

        #molpro version
        #gradient, E = mol_gradient, mol_energy

        #psi4 version
#        gradient, E = psi_gradient, psi_energy

        #(T)+EOM version
        gradient = np.add(psi_gradient,mol_gradient)
        E = mol_energy + psi_energy

        energies.append(E)

        
        #all the optking stuff now that we have our gradient
        #form hessian
        # Same guess types as c++ optking
        # intcosMisc has a convert_hessian_to_internals method if needed
        H = hessian.guess(molsys_obj, guessType=op.Params.intrafrag_hess)

        #this converts the cartesian gradients to internal coordinate forces or something
        fq = intcosMisc.q_forces(molsys_obj.intcos, molsys_obj.geom, gradient)
        
        #I don't really know what these do but they make optimization work good
        intcosMisc.apply_fixed_forces(molsys_obj, fq, H, i)
        intcosMisc.project_redundancies_and_constraints(molsys_obj, fq, H)
        molsys_obj.q_show()

        #append to optimization history
        history.oHistory.append(molsys_obj.geom, E, fq)

        #simplified backstep logic TODO see if I need to add anything to this
        lastStepOK = history.oHistory.current_step_report()
        if len(history.oHistory.steps) < 5:
            pass
        #take a backstep if we made a bad step. TODO soft-code the 5 here and see if it's a reasonable value
        elif history.History.consecutiveBacksteps < 5 and not lastStepOK:
            history.History.consecutiveBacksteps += 1
            stepAlgorithms.take_step(molsys_obj, None, fq, H, stepType="BACKSTEP")
            continue
        # None is energy, not needed for any algorithms except for the linesearch
        #TODO consider allowing algorithm control uptop
        step = stepAlgorithms.take_step(molsys_obj, None, fq, H, stepType='RFO')
        #update hessian (idk if i actually need this)
        history.oHistory.hessian_update(H, molsys_obj)

        #TODO again, switch to using logger
        print("Here is the updated geometry")
        print(molsys_obj.show_geom()) 
    
        #check convergence
        converged = convCheck.conv_check(i, molsys_obj, step, fq, energies)
        if converged:
            break
if __name__=="__main__":
    main()

