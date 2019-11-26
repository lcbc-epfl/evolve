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
#	get_angles.py
#
#	this script is an example script showing how to make use of the 
#	BASILISK dbn. It is a quite simple program just drawing a single
#	set of chi angles for a given amino acid typ. 
#	The script allows to specify a pair of backbone angles, which will
#	lead to a backbone dependent sample being drawn. Otherwise the 
#	backbone information will be integrated out. 
#
############################################################################
#
# first a bunch of standard python modules
import sys
import os
#
#import basilisk_lib 
from basilisk_lib import *

#
if __name__ == "__main__" :
   
   	# make sure that we acutally got some input
   	# otherwise print the usage message
	if len(sys.argv) < 2 :
		sys.stderr.write("Usage %s <residue type> [<phi> <psi>] \n" % sys.argv[0])
		sys.stderr.write("   <residue type>\tamino acid type in either one or three letter notation or as index (int 0-19)\n" )
		sys.stderr.write("   <phi>, <psi>\t\toptional, backbone angle in radians in the interval [-pi, pi) \n" )
		sys.stderr.write(" \n" )
		sys.exit()
	
	#
	# this is the amino acid to sample for.
	aa = sys.argv[1]
	#
	# lets find out what we got here 
	if basilisk_utils.is_numeric(aa) :
		aa = int(aa)
	elif len(aa) == 1 :
		aa = basilisk_utils.one_to_index(aa)
	elif len(aa) == 3 :
		aa = basilisk_utils.three_to_index(aa)
	else :
		sys.stderr.write("Sorry, could not match \"%s\" to any standard amino acid.\n\n" )
		sys.exit()
	#
	# keep the somewhat longer name for output later
	aa_3 = basilisk_utils.index_to_three(aa)
	
	#
	# alanine and glycine are not welcome .. nothing to sample 
	# for those
	if aa in [0,5] :
		sys.stderr.write("%s does not have any rotameric degrees of freedom .. nothing to sample here.\n\n" % aa_3 )
		sys.exit()

	#
	# did we also get some backbone angles?
	phi = -5.
	psi = -5.
	if len(sys.argv) >=4 :
		phi = sys.argv[2]
		if basilisk_utils.is_numeric(phi) :
			phi = float(phi)
		
		psi = sys.argv[3]
		if basilisk_utils.is_numeric(psi) :
			psi = float(psi)
	
	############################################################################
	# now for the actual calculations. 
	#
	# load the dbn 
	dbn = basilisk_dbn()
	#
	# sample a new set of angles
	# what we get back from the dbn is a tuple consisting of a
	# list of chi angles, a list with the sampled or given 
	# backbone angles and the log likelihood that was calculated 
	# for that sample 
	chis, bb, ll = dbn.get_sample(aa, phi, psi)
	phi, psi = bb
	
	#
	# and that was it already .. very quick and easy
	#
	############################################################################
	#
	# now be a bit verbose about what we found out. 
	out = ""
	out_short = ""
	
	out+= "# Sampling Sidechain angles for %s \n" % aa_3
	out+= "#  log likelihood : % .3f  \n" % (ll)
	out_short+= "% .3f\t" % (ll)
	
	out+= "#  phi : % .3f\t psi : % .3f  \n" % (phi,psi)
	out_short+= "% .3f\t% .3f\t" % (phi,psi)

	out+= "#  "
	for i in range (len(chis)) :
		out+= "x%s : % .3f\t" % (i+1, chis[i])
		out_short+= "% .3f\t" % (chis[i])

	out+= "\n"
	
	print out
	print out_short
	


