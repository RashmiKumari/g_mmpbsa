#!/usr/bin/env python
#
# This file is part of g_mmpbsa.
#
# Author: Rajendra Kumar
#
# Copyright (C) 2013, 2014 Rashmi Kumari and Andrew Lynn
#
# g_mmpbsa is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# g_mmpbsa is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with g_mmpbsa.  If not, see <http://www.gnu.org/licenses/>.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#
#

import os, sys, shlex, subprocess
from subprocess import Popen, PIPE
import re
import argparse

from modeller import *
from modeller.optimizers import molecular_dynamics, conjugate_gradients
from modeller.automodel import autosched


def main():
	parser = ParseOptions()
	args = CheckInput(parser)
	residues2mutate = get_residue_list_for_scan(args)
		
	wildtype_path = os.path.abspath(args.dirWildType)
	mutation_path = os.path.abspath(args.dirMutations)
	file_list = [ os.path.join(wildtype_path, f) for f in os.listdir(wildtype_path) if os.path.isfile(os.path.join(wildtype_path, f)) ]
	file_list.sort()

	mt = Mutations(file_list)
	for residue in residues2mutate:
		sys.stdout.write("\n----------------------------------------------------")
		sys.stdout.write("\nMutating residue no. %s of chain %s to %s " % (residue.pos, residue.chain, residue.name))
		sys.stdout.write("\n----------------------------------------------------\n")
		sys.stdout.flush()
		
		# Add residue into object
		mt.add_mutation(residue)

		# Mutate all frames using modeller
		mt.mutate_frames(mutation_path, args)

		# Generate trajectory, topology and tpr file of mutated protein complex	
		mt.generate_traj(mutation_path, args)

		# Remove residue from mutations
		mt.del_mutations()

		sys.stdout.write("+++++++++++++++++++ FINISHED +++++++++++++++++++++++\n")
		sys.stdout.flush()


class resid2mutate:
	def __init__(self, name, pos, chain):
		self.name = name
		self.pos = pos
		self.chain = chain

class Mutations:
	def __init__(self, file_list):
		self.file_list = file_list
		self.residues=[]
	
	def add_mutation(self, residues):
		self.residues.append(residues)

	def del_mutations(self):
		for i in range(len(self.residues)):
			self.residues.pop()

	def dirname(self):
		dirname=""
		for residue in self.residues:
			dirname = dirname + residue.name + residue.pos + residue.chain + "-"
		return dirname[:-1]

	def mutate_frames(self, mutation_path, args):
		try:
			os.makedirs(mutation_path)
		except OSError:
			if not os.path.isdir(mutation_path):
				raise

		dest_path = os.path.join(mutation_path, self.dirname())	
		try:
			os.makedirs(dest_path)
		except OSError:
			if not os.path.isdir(dest_path):
				raise

		cd = ChDir(dest_path)
		n = 1
		for a_file in self.file_list:
			sys.stdout.write("\r    ...Mutating Frame: %d/%d" % (n, len(self.file_list)))
			sys.stdout.flush()
		
			# Change of atomname CD to CD1 of ILE, modeller does not detect CD atomname
			CD_2_CD1_for_ILE(a_file, "modeller_in.pdb")
			
			# Mutations of selected residue by modeller
			mutate_a_frame("modeller_in.pdb", os.path.basename(a_file), self.residues)
							
			os.remove("modeller_in.pdb")
			n=n+1
		sys.stdout.write("\n    ...Finished\n")
		sys.stdout.flush()

		del cd

	def generate_traj(self, mutation_path, args):
		try:
			os.path.exists(mutation_path)
		except OSError:
			if os.path.isdir(mutation_path):
				raise

		dest_path = os.path.join(mutation_path, self.dirname())	
		try:
			os.path.exists(dest_path)
		except OSError:
			if os.path.isdir(dest_path):
				raise

		cd = ChDir(dest_path)
		
		#Generate gro and toplogy file for alanine scanning
		gen_top_gro_for_alascan(args, self.residues, self.file_list)	
	
		#Energy minimization
		n = 1
		for a_file in self.file_list:
			sys.stdout.write("\r    ...Energy Minimization Frame: %d/%d" % (n, len(self.file_list)))
			sys.stdout.flush()
			energy_minimization(args, a_file)
			n=n+1
		
		sys.stdout.write("\n    ...Finished\n")
		sys.stdout.flush()
		
		sys.stdout.write("\r    ...Generating trajectory and tpr file")
		sys.stdout.flush()

		# Generate trajectory and tpr from energy minimized mutated frames
		gen_traj_tpr(args)
		
		sys.stdout.write("\n    ...Finished\n")
		sys.stdout.flush()

		del cd
	
def gen_traj_tpr(args):
	files = os.listdir(os.getcwd())

	frames = []
	
	# Deleting top and itp files
	for a_file in files:
		if a_file.endswith("em.gro"):
			frames.append(a_file)
	frames.sort()

	# Running trjcat
	command = 'trjcat -f '
	for frame in frames:
		command = command + frame + ' '
	command = '{0} -cat -o {1}' .format(command, 'traj_temp.xtc')
	p = Popen(shlex.split(command), stdin=PIPE, stdout=PIPE, stderr=PIPE)
	stdout, stderr = p.communicate()
	if p.returncode != 0:
		print stdout, stderr
		sys.exit(1)

	# Generating tpr file
	gen_tpr(args, files[0])
	
	# Running trjconv to set time in trajectory
	command = 'trjcat -f '
	for frame in frames:
		command = command + frame + ' '
	command = 'trjconv -f traj_temp.xtc -o trajout.xtc -t0 {0} -timestep {1}' .format(args.time0, args.timestep)
	p = Popen(shlex.split(command), stdin=PIPE, stdout=PIPE, stderr=PIPE)
	stdout, stderr = p.communicate()
	if p.returncode != 0:
		print stdout, stderr
		sys.exit(1)

	for i in range(1, len(frames)):
		os.remove(frames[i])
	os.remove('traj_temp.xtc')
	
def gen_top_gro_for_alascan(args, residues, file_list):

	# Update protonation states of the selected residues
	n = 0
	for aFile in file_list:
		n=n+1
		sys.stdout.write("\r    ...Updating protonation state for frame %d/%d" % (n, len(file_list)) )
		sys.stdout.flush()
		InFile = os.path.basename(aFile)
		update_prot_state(args, InFile)
	
	sys.stdout.write("\n    ...Finished\n")
	sys.stdout.flush()

	# Generating pdb file using pdb2gmx, recovering hydrogen atoms and generating periodic box
	n=0
	for aFile in file_list:
		n=n+1
		sys.stdout.write("\r    ...Generating pdb file using pdb2gmx, recovering hydrogen atoms and generating periodic box for frame %d/%d" % (n, len(file_list)) )
		sys.stdout.flush()
		InFile = os.path.basename(aFile)
		OutFile1 = '{0}_nobox.pdb' .format(os.path.splitext(InFile)[0])
		OutFile2 = '{0}.gro' .format(os.path.splitext(InFile)[0])
	
		# Running pdb2gmx	
		command = 'pdb2gmx -f {0} -water tip3p -o {1}' .format(InFile, OutFile1)
		p = Popen(shlex.split(command), stdin=PIPE, stdout=PIPE, stderr=PIPE)
		p.stdin.write('{0}\n' .format(args.force_field))
		stdout, stderr = p.communicate()
		if p.returncode != 0:
			print stdout, stderr
			sys.exit(1)

		# Recover coordinates of hydrogen atoms
		if(not args.no_orig_h_pos):
			recover_hydrogen_coords(residues, aFile, OutFile1)
	
		# Generating dodecahedron box
		command = 'editconf -f {0} -bt dodecahedron -d 1 -c -o {1}' .format(OutFile1, OutFile2)
		p = Popen(shlex.split(command), stdin=PIPE, stdout=PIPE, stderr=PIPE)
		p.stdin.write('0\n0\n')
		stdout, stderr = p.communicate()
		if p.returncode != 0:
			print stdout, stderr
			sys.exit(1)
	
		files = os.listdir(os.getcwd())

		#Deleting top and itp files
		for top in files:
			if top.endswith(".top") or top.endswith(".itp"):
				os.remove(top)
		os.remove(OutFile1)
		
	sys.stdout.write("\n    ...Finished\n")
	sys.stdout.flush()

	#Generate topolgy using pdb2gmx	
	sys.stdout.write("\r    ...Generating topology files using pdb2gmx...")
	sys.stdout.flush()

	InFile = os.path.basename(file_list[0])
	
	command = ['pdb2gmx', '-f', InFile, '-water', 'tip3p', '-o', 'temp.gro']
	p = Popen(command, stdin=PIPE, stdout=PIPE, stderr=PIPE)
	p.stdin.write('{0}\n' .format(args.force_field))
	stdout, stderr = p.communicate()
	if p.returncode != 0:
		print stdout, stderr
		sys.exit(1)
		 
	os.remove('temp.gro')

	for aFile in file_list:
		InFile = os.path.basename(aFile)
		os.remove(InFile)
	
	sys.stdout.write("\n    ...Finished")
	sys.stdout.flush()

def energy_minimization(args, filename):
	pdbfile = os.path.basename(filename)
	InFile = '{0}.gro' .format(os.path.splitext(pdbfile)[0])
	OutFile = '{0}_em.gro' .format(os.path.splitext(pdbfile)[0])

	# Generating tpr file
	gen_tpr(args, InFile)	

	# Running mdrun for minimization
	command = 'mdrun -s topol.tpr -c {0}' .format(OutFile)
	p = Popen(shlex.split(command), stdin=PIPE, stdout=PIPE, stderr=PIPE)
	stdout, stderr = p.communicate()
	if p.returncode != 0:
		print stdout, stderr
		sys.exit(1)
	
	#Deleting files generated during minimization
	files = os.listdir(os.getcwd())
	for tfile in files:
		if tfile.endswith(".trr") or tfile.endswith(".tpr") or tfile.endswith(".edr") or tfile.endswith(".log") or tfile.endswith(".mdp"):
			os.remove(tfile)
	os.remove(InFile)

def gen_tpr(args, infile):
	try:
		fmdp = open("em.mdp", 'w')
	except IOError:
		print "\nCould not write em.mdp\n"
		raise

	fmdp.write("constraints         =  none\n")
	fmdp.write("integrator          =  steep\n")
	fmdp.write("nsteps              =  {0}\n" .format(args.em_nsteps))
	fmdp.write("emtol               =  {0}\n" .format(args.emtol))
	fmdp.write("emstep              =  {0}\n" .format(args.emstep))
	fmdp.write("nstxout		=  1\n")
	fmdp.write("nstcgsteep          =  10\n")
	fmdp.write("energygrps	        =  system\n")
	fmdp.write("nstcomm             =  1\n")
	fmdp.write("ns_type             =  {0}\n" .format(args.em_ns_type))
	fmdp.write("rlist               =  {0}\n" .format(args.em_cutoff))
	fmdp.write("rcoulomb            =  {0}\n" .format(args.em_cutoff))
	fmdp.write("rvdw                =  {0}\n" .format(args.em_cutoff))
	fmdp.write("Tcoupl              =  no\n")
	fmdp.write("Pcoupl              =  no\n")
	fmdp.write("gen_vel             =  no\n")
	fmdp.close()	

	command = 'grompp -f em.mdp -c {0} -p topol.top -o topol.tpr' .format(infile) 
	p = Popen(shlex.split(command), stdin=PIPE, stdout=PIPE, stderr=PIPE)
	stdout, stderr = p.communicate()
	if p.returncode != 0:
		print stdout, stderr
		sys.exit(1)


def CD_2_CD1_for_ILE(infile, outfile):
	#Opening input pdb file
	try:
		finPDB = open(infile, 'r')
	except IOError:
		print '\nCould not open File: {0}\n' .format(infile)
		raise
	
	PDB = []
	for line in finPDB:
		PDB.append(line)

	finPDB.close()

	for i in range(len(PDB)):
		tmp = re.split('\s+', PDB[i])
		if (tmp[0]=='ATOM') and (tmp[3]=='ILE') and (tmp[2] == 'CD'):
			PDB[i] = PDB[i].replace('CD ', 'CD1')
	
	#Writing output PDB file	
	try:
		fout = open(outfile, 'w')
	except IOError:
		print '\nCould not open File: {0}\n' .format(filename)
		raise
	for line in PDB:
		fout.write(line)
	fout.close()

def recover_hydrogen_coords(residues, infile, outfile):
	#Opening input pdb file
	try:
		finPDB = open(infile, 'r')
	except IOError:
		print '\nCould not open File: {0}\n' .format(infile)
		raise
	
	inPDB = []
	for line in finPDB:
		inPDB.append(line)

	finPDB.close()

	#Opening output pdb file
	try:
		foutPDB = open(outfile, 'r')
	except IOError:
		print '\nCould not open File: {0}\n' .format(outfile)
		raise
	
	outPDB = []
	for line in foutPDB:
		outPDB.append(line)
	foutPDB.close()
	os.remove(outfile)

	mutateIN = False
	mutateOUT = False
	for i in range(len(inPDB)):
		tmpIN = re.split('\s+', inPDB[i])
		if (tmpIN[0] != 'ATOM'):
			continue

		#Checking for mutated residues
		for residue in residues:
			if (residue.pos == tmpIN[5]) and (residue.chain == tmpIN[4]):
				mutateIN = True
				break
		if mutateIN:
			mutateIN = False
			continue

		## Maximum atom difference between mutated and wild-type
		natoms=len(residues)*20
		if (i<=natoms):
			back = 0
		else:
			back = i-natoms
		####
		for j in range(back, len(outPDB)):
			tmpOUT = re.split('\s+', outPDB[j])
			if (tmpOUT[0] !='ATOM'):
				continue
			
			#Checking for mutated residues
			for residue in residues:
				if (residue.pos == tmpOUT[5]) and (residue.chain == tmpOUT[4]):
					mutateOUT = True
					break
			if mutateOUT:
				mutateOUT = False
				continue
		
			#Replacing coordinates of hydrogen atoms
			if ( (re.search(r"^H", tmpOUT[2]) != None) or (re.search(r"^\dH", tmpOUT[2]) != None)) and (tmpOUT[2] == tmpIN[2]) and (tmpOUT[4] == tmpIN[4]) and (tmpOUT[5] == tmpIN[5]):
				outPDB[j] = inPDB[i]
				break
	
	#Writing output PDB file	
	try:
		fout = open(outfile, 'w')
	except IOError:
		print '\nCould not open File: {0}\n' .format(outfile)
		raise
	for line in outPDB:
		fout.write(line)
	fout.close()

def CheckInput(parser):
	args = parser.parse_args()
	if args.residue_range==None:
		print "ERROR: Enter -rr or --residue_range!!!\n\n"
		parser.print_help()
		sys.exit(1)
	return args

def get_residue_list_for_scan(args):
	line = args.residue_range
	residues2mutate = []
	chains = line.split(';')
	for value in chains:
		rlist = []
		chain, all_resids = value.split(':')
		resids = all_resids.split(',')
		for resid in resids:
			if re.search(r"-", resid) != None:
				rstart, rend = resid.split('-')
				rlist = rlist + range(int(rstart), int(rend)+1)
			else:
				rlist.append(resid)
		for i in range(len(rlist)):
			residues2mutate.append(resid2mutate(args.residue_name, str(rlist[i]), chain))
	#for i in range(len(residues2mutate)):	
		#print residues2mutate[i].name, residues2mutate[i].pos, residues2mutate[i].chain
	#sys.exit(0)
	return residues2mutate

def update_prot_state(args, filename):
	#Get residue list from command prompt
	def get_resid_list(line):
		resid, chain = [], []
		csv = line.split(';')
		for csv_value in csv:
			ch, rr = csv_value.split(':')
			rlist = rr.split(',')
			for i in range(len(rlist)):
				chain.append(ch)
				resid.append(rlist[i])
		return resid, chain
	
	#Replace specific residue
	def replace_residue(PDB, resid, chain, old, new):
		for r in range(len(resid)):
			for i in range(len(PDB)):
				tmp = re.split('\s+', PDB[i])
				if (tmp[0]=='ATOM') and (resid[r]==tmp[5]) and (chain[r]==tmp[4]):
					PDB[i] = PDB[i].replace(old, new)
		return PDB
	
	#Replace HIS with HIE
	def replace_HIS_to_HIE(PDB):
		for i in range(len(PDB)):
			tmp = re.split('\s+', PDB[i])
			if (tmp[0]=='ATOM') and (tmp[3]=='HIS'):
				PDB[i] = PDB[i].replace('HIS', 'HIE')
		return PDB
		
	#Opening input pdb file
	try:
		finPDB = open(filename, 'r')
	except IOError:
		print '\nCould not open File: {0}\n' .format(infile)
		raise
	
	PDB = []
	for line in finPDB:
		PDB.append(line)

	finPDB.close()
	os.remove(filename)
	
	#Replacing HIS with HIP
	if (args.hip != None):
		resid, chain = get_resid_list(args.hip)
		PDB = replace_residue(PDB, resid, chain, 'HIS', 'HIP')	
			
	#Replacing HIS with HID
	if (args.hid != None):
		resid, chain = get_resid_list(args.hid)
		PDB = replace_residue(PDB, resid, chain, 'HIS', 'HID')

	#Replacing HIS with HIE
	replace_HIS_to_HIE(PDB)
	
	#Replacing GLU with GLH
	if (args.glh != None):
		resid, chain = get_resid_list(args.glh)
		PDB = replace_residue(PDB, resid, chain, 'GLU', 'GLH')
	
	#Replacing ASP with ASH
	if (args.ash != None):
		resid, chain = get_resid_list(args.ash)
		PDB = replace_residue(PDB, resid, chain, 'ASP', 'ASH')
	
	#Replacing LYS with LYN
	if (args.lyn != None):
		resid, chain = get_resid_list(args.lyn)
		PDB = replace_residue(PDB, resid, chain, 'LYS', 'LYN')
	
	#Replacing CYS with CYM
	if (args.cym != None):
		resid, chain = get_resid_list(args.cym)
		PDB = replace_residue(PDB, resid, chain, 'CYS', 'CYM')

	#Writing output PDB file	
	try:
		fout = open(filename, 'w')
	except IOError:
		print '\nCould not open File: {0}\n' .format(filename)
		raise
	for line in PDB:
		fout.write(line)
	fout.close()


def ParseOptions():
	parser = argparse.ArgumentParser()
	
	parser.add_argument("-drWT", "--dirWildType", help='Directory containing wild-type frames from MD trajectory', action="store", default='wildtype', metavar='wildtype')
	parser.add_argument("-drMT", "--dirMutations", help='Directory containing mutated frames from MD trajectory', action="store", default='mutations', metavar='mutations')
	parser.add_argument("-ff", "--force_field", help='Force-field number to choose in pdb2gmx. Default value is 6, which corrosponds to AMBER99SB-ILDN force-field in standard GROMACS package', action="store", default='6', metavar='6')
	parser.add_argument("-rr", "--residue_range", help='Input residue range and/or list with respective chain for mutating residues to ALA. e.g. -rr \'A:27,45,67-70;B:40,43,78-80\'', action="store")
	parser.add_argument("-rnm", "--residue_name", help='Target residue name for mutation. e.g. for Alanine scanning: -rnm \'ALA\'', action="store", default="ALA", metavar="ALA")
	parser.add_argument("-nrch", "--no_orig_h_pos", help='Do not recover original coordinates of hydrogen atoms.', action="store_true", default=False)
	parser.add_argument("-hip", "--hip", help='HIS residues need to be protonated as HIP. Usage: e.g. -hip \'A:23,34,56;B:32,56\'', action="store", default=None, metavar=None)
	parser.add_argument("-hid", "--hid", help='HIS residues need to be as HID. Usage: e.g. -hid \'A:23,34,56;B:32,56\'', action="store", default=None, metavar=None)
	parser.add_argument("-glh", "--glh", help='GLU residues need to be protonated as GLH. Usage: e.g. -glh \'A:23,34,56;B:32,56\'', action="store", default=None, metavar=None)
	parser.add_argument("-ash", "--ash", help='ASP residues need to be protonated as ASH. Usage: e.g. -ash \'A:23,34,56;B:32,56\'', action="store", default=None, metavar=None)
	parser.add_argument("-lyn", "--lyn", help='LYS residues need to be neutral as LYN. Usage: e.g. -lyn \'A:23,34,56;B:32,56\'', action="store", default=None, metavar=None)
	parser.add_argument("-cym", "--cym", help='CYS residues need to be deprotonated as CYM. Usage: e.g. -cym \'A:23,34,56;B:32,56\'', action="store", default=None, metavar=None)
	parser.add_argument("-emnsteps", "--em_nsteps", help='Number of minimization steps (nsteps)', action="store", default=10000, type=int, metavar=10000)
	parser.add_argument("-emct", "--em_cutoff", help='Cut-off in nm during energy minimization', action="store", default=1.4, type=float, metavar=1.4)
	parser.add_argument("-emtol", "--emtol", help='Tolerance value (emtol) for energy minimization', action="store", default=10, type=float, metavar=10)
	parser.add_argument("-emstep", "--emstep", help='Minimization step (emstep) for energy minimization', action="store", default=0.01, type=float, metavar=0.01)
	parser.add_argument("-emnstyp", "--em_ns_type", help='Cut-off scheme during energy minimization', action="store", default="grid", metavar="grid")
	parser.add_argument("-t0", "--time0", help='Starting time (ps) of mutated trajectory', action="store", default=0, type=float, metavar=0)
	parser.add_argument("-tstep", "--timestep", help='Time step (ps) between frames of mutated trajectory', action="store", default=500, type=float, metavar=500)
	
	return parser
	

#def mutate(modelname, respos, restyp, chain):
def mutate_a_frame(modelname, outfile, residues):


	def optimize(atmsel, sched):
		#conjugate gradient
		for step in sched:
			step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
		#md
		refine(atmsel)
		cg = conjugate_gradients()
		cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)


	#molecular dynamics
	def refine(atmsel):
		# at T=1000, max_atom_shift for 4fs is cca 0.15 A.
		md = molecular_dynamics(cap_atom_shift=0.39, md_time_step=4.0, md_return='FINAL')
		init_vel = True
		for (its, equil, temps) in ((200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)), (200, 600, (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0))):
			for temp in temps:
				md.optimize(atmsel, init_velocities=init_vel, temperature=temp, max_iterations=its, equilibrate=equil)
		init_vel = False


	#use homologs and dihedral library for dihedral angle restraints
	def make_restraints(mdl1, aln):
   		rsr = mdl1.restraints
		rsr.clear()
		s = selection(mdl1)
   		for typ in ('stereo', 'phi-psi_binormal'):
			rsr.make(s, restraint_type=typ, aln=aln, spline_on_site=True)
		for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'):
			rsr.make(s, restraint_type=typ+'_dihedral', spline_range=4.0,spline_dx=0.3, spline_min_points = 5, aln=aln, spline_on_site=True)

	log.level(output=0, warnings=0, errors=1)	
	#log.verbose()
	# Set a different value for rand_seed to get a different final model
	env = environ(rand_seed=-49837)

	env.io.hetatm = True
	#soft sphere potential
	env.edat.dynamic_sphere=False
	#lennard-jones potential (more accurate)
	env.edat.dynamic_lennard=True
	env.edat.contact_shell = 4.0
	env.edat.update_dynamic = 0.39

	# Read customized topology file with phosphoserines (or standard one)
	env.libs.topology.read(file='$(LIB)/top_heav.lib')

	# Read customized CHARMM parameter library with phosphoserines (or standard one)
	env.libs.parameters.read(file='$(LIB)/par.lib')


	# Read the original PDB file and copy its sequence to the alignment array:
	mdl1 = model(env, file=modelname)
	ali = alignment(env)
	ali.append_model(mdl1, atom_files=modelname, align_codes=modelname)

	for residue in residues:
		#set up the mutate residue selection segment
		#s = selection(mdl1.chains[chain].residues[respos])
		s = selection(mdl1.chains[residue.chain].residues[residue.pos])

		#perform the mutate residue operation
		#s.mutate(residue_type=restyp)
		s.mutate(residue_type=residue.name)
	
	#get two copies of the sequence.  A modeller trick to get things set up
	ali.append_model(mdl1, align_codes=modelname)

	# Generate molecular topology for mutant
	mdl1.clear_topology()
	mdl1.generate_topology(ali[-1])


	# Transfer all the coordinates you can from the template native structure
	# to the mutant (this works even if the order of atoms in the native PDB
	# file is not standard):
	#here we are generating the model by reading the template coordinates
	mdl1.transfer_xyz(ali)

	# Build the remaining unknown coordinates
	mdl1.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')

	#yes model2 is the same file as model1.  It's a modeller trick.
	mdl2 = model(env, file=modelname)

	#required to do a transfer_res_numb
	#ali.append_model(mdl2, atom_files=modelname, align_codes=modelname)
	#transfers from "model 2" to "model 1"
	mdl1.res_num_from(mdl2,ali)

	#It is usually necessary to write the mutated sequence out and read it in
	#before proceeding, because not all sequence related information about MODEL
	#is changed by this command (e.g., internal coordinates, charges, and atom
	#types and radii are not updated).

	mdl1.write(file=outfile+'.tmp')
	mdl1.read(file=outfile+'.tmp')

	#set up restraints before computing energy
	#we do this a second time because the model has been written out and read in,
	#clearing the previously set restraints
	make_restraints(mdl1, ali)

	#a non-bonded pair has to have at least as many selected atoms
	mdl1.env.edat.nonbonded_sel_atoms=1

	sched = autosched.loop.make_for_model(mdl1)

	#only optimize the selected residue (in first pass, just atoms in selected
	#residue, in second pass, include nonbonded neighboring atoms)
	#set up the mutate residue selection segment
	for residue in residues:
		#s = selection(mdl1.chains[chain].residues[respos])
		s = selection(mdl1.chains[residue.chain].residues[residue.pos])

	mdl1.restraints.unpick_all()
	mdl1.restraints.pick(s)

	s.energy()

	s.randomize_xyz(deviation=4.0)

	mdl1.env.edat.nonbonded_sel_atoms=2
	optimize(s, sched)

	#feels environment (energy computed on pairs that have at least one member
	#in the selected)
	mdl1.env.edat.nonbonded_sel_atoms=1
	optimize(s, sched)

	s.energy()

	#give a proper name
	mdl1.write(file=outfile)

	#delete the temporary file
	os.remove(outfile+'.tmp')

class ChDir:
	def __init__(self, newPath):
		self.savedPath = os.getcwd()
		os.chdir(newPath)
	def __del__(self):
		os.chdir(self.savedPath)


if __name__=="__main__":
	main()
