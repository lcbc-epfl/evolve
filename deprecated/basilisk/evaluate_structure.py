#  BASILISK: A generative, probabilistic model of side chains in proteins
#
#  Copyright (C) 2010 	Tim Harder and Jes Frellsen 
#
#  BASILISK is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  BASILISK is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with BASILISK.  If not, see <http://www.gnu.org/licenses/>.
#
############################################################################
#
#	evaluate_structure.py
#
#	this script is an example script showing how to make use of the 
#	BASILISK dbn. It is a quite simple program just evaluating the 
#	likelihood of the side chain conformations in a given pdb structure.
#
############################################################################
#
# first a bunch of standard python modules
import sys
import os
import optparse
from math import log, exp

#
# lets see if we can find biopython which helps us 
# reading those pdb files
try:
    import Bio.PDB
except ImportError :                    
    sys.stderr.write("This program requires biopython to run. Please make sure you have biopython installed. \nSee http://www.biopython.org for download and help.\n\n")
    sys.exit() 

#
#import basilisk_lib 
from basilisk_lib import *
import  basilisk_lib.CalcAngles as CalcAngles

#############################################################################
#  get a little help intepreting the command line arguments.
#
def parse_commandline () :
	"""
	Just parse the commandline options and see what we got there. 
	This functions uses the optparse module.
	@return options object, containing all the options set.
	"""
	
	#
	# setting up the parser
	parser = optparse.OptionParser(usage = "%prog [options]",version="%prog 0.1 alpha")
	#
	# the different files we need
	parser.add_option("-p", "--pdb-file", dest="pdb", help="Specify the PDB file to be used [%default]", default='')
	#
	# some general switches 
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="be verbose with the output.", default=False)
	#
	# and go for it then
	options, args = parser.parse_args()
		
	#
	# be done ..
	return options


#############################################################################
if __name__ == "__main__" :
	
	#
	# lets get all the input together.
	options = parse_commandline() 
	
	if options.pdb == '' and len(sys.argv) < 2 :
		sys.stderr.write("Usage %s <pdb file> \n" % sys.argv[0])
		sys.exit()
	elif options.pdb == "" :
		options.pdb = sys.argv[1]
	
	#
	#  create a new BASILISK dbn object and load the given structure
	dbn = basilisk_dbn()
	pdb = basilisk_utils.load_pdb_file(options.pdb)

	# lets first calculate all the angles in the 
	# structure .. 
	for model in pdb :
		for chain in model :
			pd = 0
			for residue in chain: 	
				if residue.get_id()[0] != ' ':
					#don't care for the HET Atoms
					continue
				CalcAngles.calc_phi_psi(residue, pd)
				CalcAngles.calc_chi_angles(residue)
				pd = residue

	# 
	# print a column description header
	if options.verbose :
		print "Model\tChain\tResID\tResType\tphi\tpsi\tx1\tx2\tx3\tx4\tlog_likelihood"	
		
	#
	# and loop over the residues again .. this time we calculate the likelihood for the 
	# angles we found.
	for model in pdb :
		for chain in model :
			ll_sum = 0
			for residue in chain: 	
				if residue.get_id()[0] != ' ':
					#don't care for the HET Atoms
					continue
				# generate the output right away
				out = ""
				out+= "%s\t%s\t%s\t%s\t" % (model.id, chain.id , residue.get_id()[1], residue.resname)
				#
				phi = -5.
				psi = -5.
				chis = []
				aa = basilisk_utils.three_to_index(residue.resname)
				#
				# PHI angle
				if residue.xtra.has_key("phi") and residue.xtra["phi"] > -5.:
					out+= "% .3f\t" % (residue.xtra["phi"])
					phi = residue.xtra["phi"]
				else: 
					out+= " --- \t" 
				#
				# PSI angle
				if residue.xtra.has_key("psi") and residue.xtra["psi"] > -5.:
					out+= "% .3f\t" % (residue.xtra["psi"])
					psi = residue.xtra["psi"]
				else: 
					out+= " --- \t" 
				#
				#chis 
				if residue.xtra.has_key("chi") :
					for i in range (1,5) :
						x = "x%d" % i
						if  residue.xtra["chi"].has_key(x) and residue.xtra["chi"][x] > -5.:
							out+= "% .3f\t" % (residue.xtra["chi"][x])
							chis.append(residue.xtra["chi"][x])
						else: 
							out+= " --- \t" 
				else :
					out+= " --- \t --- \t --- \t --- \t" 
						
				ll =  dbn.get_log_likelihood(aa, chis, phi, psi) 

				out+= "% .3f\t" % (ll)
				ll_sum += exp(ll)
				if options.verbose:
					print out			
			ll_sum = log(ll_sum)
			if options.verbose :
				print "-"*90
				print "\t\t\t\t\t\t\t\t\t\t% .3f" % ll_sum
				print "\n"			
			else :
				print ll_sum
			
		


