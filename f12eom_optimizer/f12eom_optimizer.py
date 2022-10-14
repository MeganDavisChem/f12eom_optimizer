"""Main module."""
import qcengine as qcng
import qcelemental as qcel
import psi4
import optking
from optking import molsys, stepAlgorithms, hessian, displace, intcosMisc, history, optwrapper, convCheck
from optking import optparams as op
import subprocess
import numpy as np
import os


#TODO may have to feed in a Zmat to molpro. <-- can do this for efficiency but i'd rather not
#TODO use logger and consolidate all outputs to one file, cleanup molpro cruft. Probably don't even need to write
#molpro input files now that I think about it, just use stdin or python or something
#rewrite these to use a python dictionary


#TODO remove default values? Could lead to confusing bugs
optparams = {
        'maxiter': 50,
        'natoms': 2,
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
        'psi_basis': '',
        'canonical': 'false',
        'IP': 'false',
        'cfour_states': '0/0/1/0',
        'cfour_core': 'ON',
        'cfour_basis': '',
        'cfour_grad': '',
        'cfour_ref': 'RHF',
        'cfour_excite': 'EOMIP',
        'cfouryz': ''
        }
input_file = 'opt.inp'

maxiter = 50
natoms = 2
charge = 0
spin = 0
memory = 8
nproc = 4
mem = memory*1000/nproc
roots = [0, 1]
molpro_executable = "/home/qc/bin/molpro2020.sh"
#TODO get molpro working on maple
#molpro_executable = "molpro"

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

molpro_input_canonical = """*** part of geometry opt run
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
basis={mol_basis}
hf,maxit=500;accu,20;
ccsd(t),maxit=250;{mol_core}orbital,IGNORE_ERROR;
ee(1)=energy(1)
forces,varsav
SHOW[f18.12],ee(1)
SHOW[f18.12],GRADX
SHOW[f18.12],GRADY
SHOW[f18.12],GRADZ
"""

cfour_ccsd_input = """%GRD={cfour_grad}
CO AUG-PVTZ CCSD
{cfouryz}



*CFOUR(CALC=CCSD,BASIS={cfour_basis},COORDINATES=CARTESIAN,UNITS=ANGSTROM
CC_CONV=10
LINEQ_CONV=12
SCF_CONV=10
ESTATE_CONV=8
REFERENCE={cfour_ref}
CHARGE={charge}
MULTIPLICITY={psispin}
FROZEN_CORE={cfour_core}
VIB=FINDIF
FD_CALCTYPE=GRADONLY
MEM_UNIT=GB,MEMORY_SIZE={memory})
"""

cfour_eomip_input = """%GRD={cfour_grad}
CO AUG-PVTZ CCSD
{cfouryz}



*CFOUR(CALC=CCSD,BASIS={cfour_basis},COORDINATES=CARTESIAN,UNITS=ANGSTROM
CC_CONV=10
LINEQ_CONV=12
SCF_CONV=10
ESTATE_CONV=8
REFERENCE={cfour_ref}
CHARGE={charge}
MULTIPLICITY={psispin}
FROZEN_CORE={cfour_core}
VIB=FINDIF
FD_CALCTYPE=GRADONLY
EXCITE={cfour_excite}
ESTATE_SYM={cfour_states}
MEM_UNIT=GB,MEMORY_SIZE={memory})
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
    print("""
    {charge} {psispin}
    {xyz}
    units {units}
         """.format_map(optparams))
    mol = psi4.geometry("""
    {charge} {psispin}
    {xyz}
    units {units}
         """.format_map(optparams))
    molsys_obj = molsys.Molsys.from_psi4_molecule(mol)
    optparams['xyz'] = molsys_obj.show_geom()

    #TODO this program needs to be rewritten lol
    if optparams['theory'] == 'mt':
        optparams['psi_theory'] = 'mt'
    elif optparams['freeze_core'] == 'true':
        if optparams['canonical'] == 'true':
            optparams['mol_basis'] = 'aug-cc-pv{theory}'.format_map(optparams)
        else:
            optparams['mol_basis'] = 'v{theory}-f12'.format_map(optparams)
        optparams['cfour_basis'] = 'aug-pv{theory}'.format_map(optparams)
        optparams['psi_theory'] = 'aug-cc-pv{theory}'.format_map(optparams)
    elif optparams['freeze_core'] == 'false':
        if optparams['canonical'] == 'true':
            optparams['mol_basis'] = 'aug-cc-pcv{theory}'.format_map(optparams)
        else:
            optparams['mol_basis'] = 'cc-pcv{theory}-f12'.format_map(optparams)
        optparams['mol_core'] = 'core;'
        optparams['cfour_core'] = 'OFF'
        optparams['cfour_basis'] = 'aug-pcv{theory}'.format_map(optparams)
        optparams['psi_theory'] = 'aug-cc-pcv{theory}'.format_map(optparams)
    optparams['cfour_basis'] = optparams['cfour_basis'].upper()

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

    mol.set_molecular_charge(int(_optparams['charge']))
    mol.set_multiplicity(int(_optparams['psispin']))
    

    
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

def run_psi4_canonical(_molsys_obj, _optparams):
    """
    Runs CCSD(T) in Psi4
    """
    psi4.set_memory('{mem} MB'.format_map(_optparams))
    #might need to fix this
    psi4.core.set_num_threads(int(_optparams['nproc']))
    #formats this for a psi4 geomtetry

    mol = psi4.geometry(_molsys_obj.show_geom().replace("\tFragment 1\n\t", "") + f"\nunits ang")
    
    mol.set_molecular_charge(int(_optparams['charge']))
    mol.set_multiplicity(int(_optparams['psispin']))

    
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
    gradient = psi4.gradient('ccsd(t)', dertype='energy').np.flatten()
    E = psi4.energy('ccsd(t)')
    ##    gradient,wfn = psi4.gradient('eom-ccsd', dertype='energy', return_wfn='true').np.flatten()
#    eom_gradient,wfn = psi4.gradient('eom-ccsd', dertype='energy', return_wfn='true')
    #
#    eom_gradient = eom_gradient.np.flatten()
    #E = wfn.variables()['CURRENT ENERGY']

 #   gradient = np.subtract(eom_gradient, ccsd_gradient)
#    E = wfn.variables()['CURRENT ENERGY'] - wfn.variables()['CCSD TOTAL ENERGY']

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
    subprocess.Popen([molpro_executable, "1", 
        "{nproc}".format_map(_molpro_keys), "molpro_step.inp"]).wait()

def run_cfour(_cfour_ccsd_input, _cfour_eomip_input, _count):
    """This will run cfour, similar to the molpro command"""
    #run cfour for CCSD gradient
    cwd = os.getcwd()
    ccsd_dir = f'CCSD_{_count}'
    optparams['cfour_grad'] = cwd + '/' + ccsd_dir + '/' + 'GRD'
    os.mkdir(ccsd_dir)
    os.chdir(ccsd_dir)
    cfourgeom = optparams['xyz'].split('\n')[1:]

    for i in range(len(cfourgeom)):
        cfourgeom[i] = cfourgeom[i].strip('\t').split()
        cfourgeom[i].append('\n')
        print(cfourgeom[i])
        cfourgeom[i] = " ".join(cfourgeom[i])

    cfourgeom = "".join(cfourgeom)



    optparams['cfouryz'] = cfourgeom

#    optparams['cfouryz'] = optparams['xyz']

    with open('ZMAT', 'w') as ZMAT:
        ZMAT.write(_cfour_ccsd_input.format_map(optparams))
    subprocess.Popen(['/home/qc/bin/c4_new.sh', '1']).wait()
    #run cfour for EOMIP gradient
    #make directory
    os.chdir(cwd)
    eomip_dir = f'EOMIP_{_count}'
    optparams['cfour_grad'] = cwd + '/' + eomip_dir + '/' + 'GRD'
    os.mkdir(eomip_dir)
    os.chdir(eomip_dir)
    with open('ZMAT', 'w') as ZMAT:
        ZMAT.write(_cfour_eomip_input.format_map(optparams))
    subprocess.Popen(['/home/qc/bin/c4_new.sh', '1']).wait()
    os.chdir(cwd)


    #format input file

    #run cfour
    pass

def read_cfour(_count):
    """This will parse the cfour output for a gradient and put it into the proper format which I'm sure will be a nightmare"""
    cwd = os.getcwd()
    ccsd_dir = f'CCSD_{_count}'
    os.chdir(ccsd_dir)

    with open('GRD') as GRD:
        grdlines_ccsd = GRD.readlines()
        #parse from here into right format

    with open('output.dat') as output:
        ccsd_output = output.readlines()

    os.chdir(cwd)
    eomip_dir = f'EOMIP_{_count}'
    os.chdir(eomip_dir)
    with open('GRD') as GRD:
        grdlines_eomip = GRD.readlines()
    with open('output.dat') as output:
        eomip_output = output.readlines()

    os.chdir(cwd)

    na = int(optparams['natoms'])
    grdlines_ccsd = grdlines_ccsd[na+1:]
    grdlines_eomip = grdlines_eomip[na+1:]

    grd_ccsd = []
    for line in grdlines_ccsd:
        grd = line.strip('\n').split()[1:]
        for g in grd:
            grd_ccsd.append(float(g))

    grd_eomip = []
    for line in grdlines_eomip:
        grd = line.strip('\n').split()[1:]
        for g in grd:
            grd_eomip.append(float(g))

    grd_fin = np.empty(na*3)

    print('printing this stupid stuff')
    print(grd_eomip)
    print(grd_ccsd)
    for i,grd in enumerate(grd_eomip):
        grd_comb = grd - grd_ccsd[i]
        grd_fin[i] = grd_comb
    print(grd_fin)





#    grd_ccsd = [float(grd) for i,grd in enumerate(line.strip('\n').split()) if i > 1 for line in grdlines_ccsd]
    
#    grd_ccsd = []
#    for line in grdlines_ccsd:
#        grds = [float(grd) for i,grd in enumerate(line.strip('\n').split()) if i > 0] 
#        for grd in grds:
#            grd_ccsd.append(grd)
#
#    
#    print("I am printing grd_ccsd")
    #this is all fucking stupid
 #   grd_eomip = [float(grd) for i,grd in enumerate(line.strip('\n').split()) if i > 1 for line in grdlines_eomip]
#    grd_eomip = [float(grd) for grd in line.strip('\n').split() for line in grdlines_eomip]
#    grd_eomip = []
#    for line in grdlines_eomip:
#        grds = [float(grd) for i,grd in enumerate(line.strip('\n').split()) if i > 0] 
#        for grd in grds:
#            grd_eomip.append(grd)
#
#    grd_fin = np.empty(na*3)
#
#    for i,grd in enumerate(grd_eomip):
#        grd_comb = grd - grd_ccsd[i]
#        grd_fin[i] = grd_comb
    
#    print(grd_fin)
#    print(grdlines_eomip)

    _gradient = grd_fin


    for line in eomip_output:
        if "Total EOMIP-CCSD electronic energy" in line:
            _energy = float(line.split()[-2])

    #parse for CCSD energy

    #parse for CCSD gradient

    #parse for EOMIP energy

    #parse for EOMIP gradient

    #form pure EOM-IP gradient



    return _gradient, _energy


def run_molpro(_molpro_input):
    """Runs molpro with input given above (i.e. a CCSD(T)-F12b/v(d,t)z-f12 calculation).
    This structure should be completely redone (when I have time) to
    not generate useless molpro files"""
    with open('molpro_step.inp', 'w') as molpro_file:
        molpro_file.write(_molpro_input.format_map(globals()))
    #TODO let this run parallel with psi4 job? Or remove reliance on the tsch script
    subprocess.Popen([molpro_executable, f"{nproc}", "1", "molpro_step.inp"]).wait()

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

        if optparams['IP'] == 'true':
            run_cfour(cfour_ccsd_input, cfour_eomip_input, i)
            xs_gradient, xs_energy = read_cfour(i)
        else:
            xs_gradient, xs_energy = real_run_psi4(molsys_obj, optparams)


        if optparams['canonical'] == 'true':
        #    canon_gradient, canon_energy = run_psi4_canonical(molsys_obj, optparams)
            run_molpro_v2(molpro_input_canonical, optparams)
            mol_gradient, mol_energy = read_molpro()
            gradient = np.add(xs_gradient,mol_gradient)
            E = mol_energy + xs_energy
            #gradient = np.add(psi_gradient,canon_gradient)
            #E = psi_energy + canon_energy
        else:
            run_molpro_v2(molpro_input, optparams)
            mol_gradient, mol_energy = read_molpro()

            gradient = np.add(xs_gradient,mol_gradient)
            E = mol_energy + xs_energy

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

