#!/usr/bin/python
#This script is written by Jingwei Xu, ETH Zurich

try:
        from optparse import OptionParser
except:
        from optik import OptionParser

import sys, os, math, time
import numpy as np
from EMAN2 import *

def main():
	usage = """python subtomo_helixtool_mapback.py [starfile] --options\n
	For example:
	python ~/JX_data/scripts/subtomo_helixtool_mapback.py Refine3D/actin_2Dproj_run1_data.star --ini_file actin_bin2_ctf_test/crop.tbl --column20 actin_bin2_ctf_test/indices_column20.doc --helix --crop_folder actin_bin2_ctf_test --polarity_voting"""
	parser = OptionParser(usage=usage)
	parser.add_option("--dynamo", dest="dynamo", action="store_true", help="To convert to a dynamo table? The default is False.", default=False) 
	parser.add_option("--ini_file", dest="ini_file", type="string", help="The initial orientation file for mapping back. It could be Relion star file or dynamo tbl file. The default is none.", default="")
	parser.add_option("--column20", dest="column20", type="string", help="The column20 file for convention to a dynamo table. The default is none.", default="")
	parser.add_option("--helix", dest="helix", action="store_true", help="Was the orientation of each segments generated based on the neighboring points? The default is False.", default=False)
	parser.add_option("--crop_folder", dest="crop_folder", type="string", help="The relative directory of crop folder from dynamo. The default is current directory.", default=".")
	parser.add_option("--polarity_voting", dest="polarity_voting", action="store_true", help="To vote the polarity of filament segments based on tube ID? The default is False.", default=False)
	parser.add_option("--match_bycor", dest="match_bycor", action="store_false", help="To find the corresponding paired orientation based on the coordinate X/Y? The default is True. If not, it will be based on the ImageName.", default=True)
	(options, args) = parser.parse_args()

	if len(args) < 1:
                print "ERROR: please provide the input file for processing. Exit!"
                print usage
                sys.exit(-1)

	ptclfile = args[0]
	dynamo = options.dynamo
		
	ini_file = options.ini_file
	ini_dynamo = False
	if ini_file.split('.')[-1] == "tbl":
		ini_dynamo = True
	column20 = options.column20
	if ini_dynamo and len(column20) == 0:
		print "ERROR: Please provide the corresponding column20 file to convert to a dynamo table. Exit!"
		sys.exit(-1)
	column20_lst = []
	helix = options.helix
	crop_folder = options.crop_folder
	ptcl_info_ini = ptcl_header_ini = optic_info_ini = []
	ptcl_ini_info_lst = []
	rot_index_ini = tilt_index_ini = psi_index_ini = ptcl_index_ini = tomo_index_ini = corx_index_ini = cory_index_ini = corz_index_ini = helixID_index_ini = helix_track_index_ini = psi_prior_flip_index_ini = randomset_index_ini = prior_tilt_index_ini = prior_psi_index_ini = -1
	if ini_dynamo:
		column20_info = open(column20, "r").readlines()
		for i in column20_info:
			index = int(i.split()[0])
			tomo_name = i.split()[1].split('/')[-1]
			column20_lst.append((index, tomo_name))

#Process the initial file and generate a list (because the orientation from dynamo need to be recalculated if it is previously generated based on the neighboring points.	
	ini_txt = open(ini_file, "r").readlines()
	if ini_dynamo:
		line_num = 0
		column20_start = -1
		tomo_ptcl_lst = []
		ptcl_dynamo_lst = []
		randomset = 0
		for i in ini_txt:
			ptcl_index = i.split()[0]
			dx, dy, dz = i.split()[3:6]
			posx, posy, posz = i.split()[23:26]
			posx_star = float(posx) + float(dx)
			posy_star = float(posy) + float(dy)
			posz_star = float(posz) + float(dz)
			tdrot = i.split()[6]
			tilt = i.split()[7]
			narot = i.split()[8]
			column20 = int(i.split()[19])
			column21 = i.split()[20]
			tomo_name = ""
			ptcl = "%s/particle_%s"%(crop_folder, ptcl_index.zfill(6))
			for line in column20_lst:
				if column20 == line[0]:
					tomo_name = line[1]
			if line_num %2 == 1:
				randomset = 1
			else:
				randomset = 2
			info = [ptcl, float(tdrot), float(tilt), float(narot), posx_star, posy_star, posz_star, randomset, tomo_name, int(column21), line_num]
			line_num += 1
			if int(column20) != column20_start:
				column20_start = int(column20)
				if len(tomo_ptcl_lst) != 0:
					ptcl_dynamo_lst.append(tomo_ptcl_lst)
				tomo_ptcl_lst = []
				tomo_ptcl_lst.append(info)
			else:
				tomo_ptcl_lst.append(info)
		ptcl_dynamo_lst.append(tomo_ptcl_lst)

		for tomo_lst in ptcl_dynamo_lst:
			ptcl_sort_lst = []
			tube_ID = initial_id = 0
			seg_ptcl_lst = []
			for i in tomo_lst:
				column21 = i[-2]
				if column21 != initial_id:
					tube_ID += 1
					initial_id = column21
					i.append(tube_ID)
					if len(seg_ptcl_lst) != 0:
						ptcl_sort_lst.append(seg_ptcl_lst)
					seg_ptcl_lst = []
					seg_ptcl_lst.append(i)
				else:
					i.append(tube_ID)
					seg_ptcl_lst.append(i)
			ptcl_sort_lst.append(seg_ptcl_lst)
			
			for subtube_lst in ptcl_sort_lst:
#Calculate prior tilt/psi based on start/end points
				xstart_px = subtube_lst[0][4]
				ystart_px = subtube_lst[0][5]
				zstart_px = subtube_lst[0][6]
				dist = 0.0
				xend_px = subtube_lst[-1][4]
				yend_px = subtube_lst[-1][5]
				zend_px = subtube_lst[-1][6]
				diff_x = xend_px - xstart_px
				diff_y = yend_px - ystart_px
				diff_z = zend_px - zstart_px
				length = math.sqrt(diff_x**2 + diff_y**2 + diff_z**2)
				psi_prior = tilt_prior = 0.0
				if helix:
					psi_cal = math.atan2(diff_y, diff_x) / math.pi * 180
					if psi_cal > 0:
						psi_prior = 180 - psi_cal
					else:
						psi_prior = -1 * (180 + psi_cal)
					tilt_prior = math.acos(diff_z/length) / math.pi * 180

				for num, ptcl in enumerate(subtube_lst):
					ptcl_info = ptcl[0]
					tdrot = ptcl[1]
					tilt = ptcl[2]
					narot = ptcl[3]
					posx_star = ptcl[4]
					posy_star = ptcl[5]
					posz_star = ptcl[6]
					randomset = ptcl[7]
					tomo_name = ptcl[8]
					line_num = ptcl[10]
					helical_tube_ID = ptcl[-1]
					rot_ptcl = tilt_ptcl = psi_ptcl = 0.0
					if not helix:
						rotation_matrix = euler2matrix(tdrot, tilt, narot, True)
						rot_ptcl, tilt_ptcl, psi_ptcl = matrix2euler(rotation_matrix, False)
					else:
						particle_helical_track_length = math.sqrt((posx_star - xstart_px)**2 + (posy_star - ystart_px)**2 + (posz_star - zstart_px)**2)
						rot_ptcl = tdrot
						if num == 0:
							tilt_ptcl = tilt_prior
							psi_ptcl = psi_prior
						else:
							posx_last = subtube_lst[num-1][4]
							posy_last = subtube_lst[num-1][5]
							posz_last = subtube_lst[num-1][6]
							diffx = posx_star - posx_last
							diffy = posy_star - posy_last
							diffz = posz_star - posz_last
							tilt_ptcl = math.acos(diffz/math.sqrt((diffx)**2 + (diffy)**2 + (diffz)**2)) / math.pi * 180
							psi_cal = math.atan2(diffy, diffx) / math.pi * 180
							if psi_cal > 0:
								psi_ptcl = 180 - psi_cal
							else:
								psi_ptcl = -1 * (180 + psi_cal)
					info = [ptcl_info, rot_ptcl, tilt_ptcl, psi_ptcl, helical_tube_ID, particle_helical_track_length, posx_star, posy_star, posz_star, randomset, tomo_name, line_num] 
					ptcl_ini_info_lst.append(info)
	else:
#Assume that the star file is generated from relion3.1 or higher, which contains the optic information.
#Make sure that the star file contains the following tags.
		optic_info_ini, ptcl_info_ini = fetch_optic(ini_txt)
		ptcl_index_ini = fetch_index(ptcl_info_ini, "_rlnImageName")
		rot_index_ini = fetch_index(ptcl_info_ini, "_rlnAngleRot")
		tilt_index_ini = fetch_index(ptcl_info_ini, "_rlnAngleTilt")
		psi_index_ini = fetch_index(ptcl_info_ini, "_rlnAnglePsi")
		prior_tilt_index_ini = fetch_index(ptcl_info_ini, "_rlnAngleTiltPrior")
		prior_psi_index_ini = fetch_index(ptcl_info_ini, "_rlnAnglePsiPrior")
		tomo_index_ini = fetch_index(ptcl_info_ini, "_rlnTomoName")
		corx_index_ini = fetch_index(ptcl_info_ini, "_rlnCoordinateX")
		cory_index_ini = fetch_index(ptcl_info_ini, "_rlnCoordinateY")
		corz_index_ini = fetch_index(ptcl_info_ini, "_rlnCoordinateZ") 
		offx_index_ini = fetch_index(ptcl_info_ini, "_rlnOriginXAngst")
		offy_index_ini = fetch_index(ptcl_info_ini, "_rlnOriginYAngst")
		offz_index_ini = fetch_index(ptcl_info_ini, "_rlnOriginZAngst")
		helixID_index_ini = fetch_index(ptcl_info_ini, "_rlnHelicalTubeID")
		randomset_index_ini = fetch_index(ptcl_info_ini, "_rlnRandomSubset")
		helix_track_index_ini = fetch_index(ptcl_info_ini, "_rlnHelicalTrackLength")
		psi_prior_flip_index_ini = fetch_index(ptcl_info_ini, "_rlnAnglePsiFlipRatio")
		line_num = 0

		for i in ptcl_info_ini:
			if len(i.split()) < 3 or i.startswith('#'):
				line_num += 1
				ptcl_header_ini.append(i)
				continue
			ptcl = i.split()[ptcl_index_ini].split('.')[0].split('/')[-1]
			tomo_name = i.split()[tomo_index_ini].split('/')[-1]
			corx = float(i.split()[corx_index_ini])
			cory = float(i.split()[cory_index_ini])
			corz = float(i.split()[corz_index_ini])
			rot = float(i.split()[rot_index_ini])
			tilt = float(i.split()[tilt_index_ini])
			psi = float(i.split()[psi_index_ini])
			helixID = int(i.split()[helixID_index_ini])
			helix_track = float(i.split()[helix_track_index_ini])
			psi_prior_flip = float(i.split()[psi_prior_flip_index_ini])
			randomset = int(i.split()[randomset_index_ini])
			offx = float(i.split()[offx_index_ini])
			offy = float(i.split()[offy_index_ini])
			offz = float(i.split()[offz_index_ini])
			info = [ptcl, rot, tilt, psi, helixID, helix_track, corx, cory, corz, randomset, offx, offy, offz, tomo_name, line_num]
			line_num += 1
			ptcl_ini_info_lst.append(info)	

	ptcltxt = open(ptclfile, "r").readlines()
	optic_info, ptcl_info = fetch_optic(ptcltxt)
	angpix_index = fetch_index(optic_info, "_rlnImagePixelSize")
	rot_index = fetch_index(ptcl_info, "_rlnAngleRot")
	tilt_index = fetch_index(ptcl_info, "_rlnAngleTilt")
	psi_index = fetch_index(ptcl_info, "_rlnAnglePsi")
	corx_index = fetch_index(ptcl_info, "_rlnCoordinateX")
	cory_index = fetch_index(ptcl_info, "_rlnCoordinateY")
#	offx_index = fetch_index(ptcl_info, "_rlnOriginXAngst")
#	offy_index = fetch_index(ptcl_info, "_rlnOriginYAngst")
	micro_index = fetch_index(ptcl_info, "_rlnMicrographName")
	ptcl_index = fetch_index(ptcl_info, "_rlnImageName")
	helixID_index = fetch_index(ptcl_info, "_rlnHelicalTubeID")

#Voting filament polarity
	polarity_voting = options.polarity_voting
        ptcl_voted = ptcl_info
	if polarity_voting:
                ptcl_voted = helcal_voting(ptcl_info, helixID_index, micro_index, psi_index)

	angpix = 0.0
	for i in optic_info:
		if len(i.split()) < 3 or i.startswith('#'):
			continue
		angpix = float(i.split()[angpix_index])

	outprefix = ptclfile.split('.star')[0]
	outfile = "%s_mod.star"%(outprefix)
	if dynamo:
		outfile = "%s_mod.tbl"%(outprefix)
	ptcl_lst = []
	ptcl_update_lst = []
	index_num = 0
	for i in ptcl_voted:
		if len(i.split()) < 3 or i.startswith('#'):
			continue
		micro = i.split()[micro_index]
		rot = float(i.split()[rot_index])
		tilt = float(i.split()[tilt_index])
		psi = float(i.split()[psi_index])
		corx = float(i.split()[corx_index])
		cory = float(i.split()[cory_index])
#		offx = float(i.split()[offx_index])
#		offx_pixel = offx / angpix
#		offy = float(i.split()[offy_index])
#		offy_pixel = offy / angpix
		ptcl = i.split()[ptcl_index].split('/')[-1].split('_proj')[0]
		rot_ini, tilt_ini, psi_ini, ini_line, ptcl_old_info = match_info(ptcl, corx, cory, ptcl_ini_info_lst, options.match_bycor)
		rot_ref = rot
		tilt_ref = tilt_ini + (tilt - 90)
#need to check whether it needs to add 180
		psi_ref = psi_ini + psi + 180
		ptcl_update = ""
		if dynamo:
			rotation_matrix = euler2matrix(rot_ref, tilt_ref, psi_ref, False)
			tdrot, tilt, narot = matrix2euler(rotation_matrix, True)
			for num, item in enumerate(ini_txt[ini_line].split()):
#				if num == 3:
#					dx_new = float(item) + offx
#					ptcl_update += "%s "%(dx_new)
#				elif num == 4:
#					dy_new = float(item) + offy
#					ptcl_update += "%s "%(dy_new)
				if num == 6:
					ptcl_update += "%s "%(tdrot)
				elif num == 7:
					ptcl_update += "%s "%(tilt)
				elif num == 8:
					ptcl_update += "%s "%(narot)
				else:
					ptcl_update += "%s "%(item)
		else:
			if not(ini_dynamo):
				for num, item in enumerate(ptcl_info_ini[ini_line].split()):
					if num == rot_index_ini:
						ptcl_update += "%s\t"%(rot_ref)
					elif num == psi_index_ini:
						ptcl_update += "%s\t"%(psi_ref)
					elif num == tilt_index_ini:
						ptcl_update += "%s\t"%(tilt_ref)
					elif num == prior_tilt_index_ini:
						ptcl_update += "%s\t"%(tilt_ref)
					elif num == prior_psi_index_ini:
						ptcl_update += "%s\t"%(psi_ref)
					else:
						ptcl_update += "%s\t"%(item)
#				ptcl_update += "%s\t%s\t"%(tilt_ref, psi_ref)
				index_num = num
#Convert dynamo to relion
#Using the info from lst [ptcl_info, rot_ptcl, tilt_ptcl, psi_ptcl, helical_tube_ID, particle_helical_track_length, posx_star, posy_star, posz_star, randomset, tomo_name, line_num]
			else:
				ptcl_name = "%s.mrc"%(ptcl_old_info[0])
				helical_tube_ID = ptcl_old_info[4]
				particle_helical_track_length = ptcl_old_info[5]
				posx_star, posy_star, posz_star = ptcl_old_info[6:9]
				randomset = ptcl_old_info[9]
				tomo_name = ptcl_old_info[10]
				ptcl_update = "%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t0.000000\t0.000000\t0.000000\t%s\t%s\t0.500000\t%s\t1\t%.6f\t%.6f\t%s\t"%(ptcl_name, rot_ref, tilt_ref, psi_ref, posx_star, posy_star, posz_star, helical_tube_ID, particle_helical_track_length, tomo_name, tilt_ref, psi_ref, randomset)
		ptcl_update += "\n"
		ptcl_update_lst.append(ptcl_update)

	outtxt = open(outfile, "w")
	if not dynamo:
		if not(ini_dynamo):
			for i in optic_info_ini:
				outtxt.write(i)
			for i in ptcl_header_ini:
				outtxt.write(i)
#			outtxt.write("_rlnAngleTiltPrior #%s\n_rlnAnglePsiPrior #%s\n"%(index_num + 1, index_num + 2))
		else:
			outtxt.write("# version 30001\n\ndata_optics\n\nloop_\n_rlnOpticsGroup #1\n_rlnOpticsGroupName #2\n_rlnSphericalAberration #3\n_rlnVoltage #4\n_rlnTomoTiltSeriesPixelSize #5\n_rlnCtfDataAreCtfPremultiplied #6\n_rlnImageDimensionality #7\n_rlnTomoSubtomogramBinning #8\n_rlnImagePixelSize #9\n_rlnImageSize #10\n")
			outtxt.write("1\topticsGroup1\t2.700000\t300.000000\t4.510000\t1\t3\t2.000000\t9.020000\t48\n\n")
			outtxt.write("# version 30001\n\ndata_particles\n\nloop_\n_rlnImageName #1\n_rlnAngleRot #2\n_rlnAngleTiltPrior #3\n_rlnAnglePsiPrior #4\n_rlnCoordinateX #5\n_rlnCoordinateY #6\n_rlnCoordinateZ #7\n_rlnOriginXAngst #8\n_rlnOriginYAngst #9\n_rlnOriginZAngst #10\n_rlnHelicalTubeID #11\n_rlnHelicalTrackLength #12\n_rlnAnglePsiFlipRatio #13\n_rlnMicrographName #14\n_rlnOpticsGroup #15\n_rlnAngleTilt #16\n_rlnAnglePsi #17\n_rlnRandomSubset #18\n")

	for i in ptcl_update_lst:
		outtxt.write(i)
	outtxt.close()	

def helcal_voting(ptcl_info_lst, helixID_index, micro_index, psi_index):
        helixID_lst = []
        helix_individual_lst = []
        helix_total_lst = []
        ptcl_voted = []
        for i in ptcl_info_lst:
                if len(i.split()) < 3 or i.startswith('#'):
                        ptcl_voted.append(i)
                        continue
                micro_name = i.split()[micro_index]
                helixID = i.split()[helixID_index]
                helixID_info = "%s@%s"%(helixID, micro_name)
                if helixID_info not in helixID_lst:
			helixID_lst.append(helixID_info)
                        helix_individual_lst.append(i)
			helix_total_lst.append(helix_individual_lst)
			helix_individual_lst = []
                else:
			helix_lst_index = helixID_lst.index(helixID_info)
                        helix_total_lst[helix_lst_index].append(i)

        for index, helix_individual in enumerate(helix_total_lst):
                invert = False
                invert_num = no_invert_num = 0
                for i in helix_individual:
                        psi = float(i.split()[psi_index])
                        if abs(psi) > 90:
                                invert_num += 1
                        else:
                                no_invert_num += 1
                invert_percentage = float(invert_num)/(float(invert_num) + float(no_invert_num))*100
                print "There are %.2f percentage of inverted polairty in filament %s. The number of total segment is %d."%(invert_percentage, helixID_lst[index], float(invert_num) + float(no_invert_num))
                if invert_percentage > 50:
                        invert = True
		for i in helix_individual:
			psi = float(i.split()[psi_index])
                        if invert:
                                if abs(psi) < 90:
                                        psi = psi + 180
                        else:
                                if abs(psi) > 90:
                                        psi = psi - 180
                        ptcl_voted_info = ""
                        for num, item in enumerate(i.split()):
                                if num == psi_index:
                                        ptcl_voted_info += "%.6f\t"%(psi)
                                else:
                                        ptcl_voted_info += "%s\t"%(item)
                        ptcl_voted_info += "\n"
                        ptcl_voted.append(ptcl_voted_info)
	return ptcl_voted
										
def match_info(ptcl, corx, cory, ptcl_ini_lst, bycor):
        rot_ini = tilt_ini = psi_ini = line = 0
	ptcl_old_info = []
        for i in ptcl_ini_lst:
		if bycor:
			corx_ini = i[6]
			cory_ini = i[7]
			if abs(corx_ini - corx) < 0.5 and abs(cory_ini - cory) < 0.5:
				rot_ini, tilt_ini, psi_ini = i[1:4]
				line = i[-1]
				ptcl_old_info = i
		else:
	                if ptcl == i[0]:
	                        rot_ini, tilt_ini, psi_ini = i[1:4]
	                        line = i[-1]
				ptcl_old_info = i

        return rot_ini, tilt_ini, psi_ini, line, ptcl_old_info

def fetch_optic(star_info):
        optic_info = []
        ptcl_info = []
        line = 0
        line_start = line_end = 0
        for i in star_info:
                if i.startswith('data_optics'):
                        line_start = line
                if i.startswith('data_particles'):
                        line_end = line
                line += 1
        optic_info = star_info[:line_end]
        ptcl_info = star_info[line_end:]
        return optic_info, ptcl_info

def fetch_index(txt, item):
        imgindex = ""
        for i in txt:
                if len(i.split()) < 3:
                        if item in i:
                                imgindex = i.split('#')[-1]
        if imgindex == "":
                imgindex = 0
        return int(imgindex) - 1

def euler2matrix(rot, tilt, psi, zxz):
        deg2rad = 180/math.pi
        alpha = rot / deg2rad
        beta = tilt / deg2rad
        gamma = psi / deg2rad
        ca = math.cos(alpha)
        cb = math.cos(beta)
        cg = math.cos(gamma)
        sa = math.sin(alpha)
        sb = math.sin(beta)
        sg = math.sin(gamma)
        a00 = a01 = a02 = a10 = a11 = a12 = a20 = a21 = a22 = 0
        if zxz:
                cc = cg
                sc = sg
                a00 = ca * cc - cb * sa * sc
                a01 = ca * sc + cb * cc * sa
                a02 = sb * sa
                a10 = cc * sa + cb * ca * sc
                a10 *= -1
                a11 = cb * ca * cc - sa * sc
                a12 = ca * sb
                a20 = sb * sc
                a21 = -1 * cc * sb
                a22 = cb
        else:
                cc = cb * ca
                cs = cb * sa
                sc = sb * ca
                ss = sb * sa
                a00 = cg * cc - sg * sa
                a01 = cg * cs + sg * ca
                a02 = -cg * sb
                a10 = -sg * cc - cg * sa
                a11 = -sg * cs + cg * ca
                a12 = sg * sb
                a20 =  sc
                a21 =  ss
                a22 = cb
        t = np.array([[a00,a01,a02],[a10,a11,a12],[a20,a21,a22]])
        return t

def matrix2euler(matrix, zxz):
        t10 = matrix[1,0]
        t00 = matrix[0,0]
        t21 = matrix[2,1]
        t20 = matrix[2,0]
        t12 = matrix[1,2]
        t02 = matrix[0,2]
        t22 = matrix[2,2]
        abs_sinbeta = math.sqrt(t02**2 + t12**2)
        epsilon = 10 ** -5
        if (abs_sinbeta > 16 * epsilon):
                alpha = math.atan2(t21,t20)
		if zxz:
			alpha = math.atan2(t20, -t21)
                gamma = math.atan2(t12, -t02)
		if zxz:
			gamma = math.atan2(t02, t12)
                if abs(math.sin(gamma)) < epsilon:
                        sign_sb = sgn(-t02 / cos(gamma))
                else:
                        if math.sin(gamma) > 0:
                                sign_sb = sgn(t12)
                        else:
                                sign_sb = -sgn(t12)
                beta = math.atan2(sign_sb * abs_sinbeta, t22)
        else:
                if sgn(t22) > 0:
                        alpha = 0
                        beta = 0
                        gamma = math.atan2(-t10,t00)
                else:
                        alpha = 0
                        beta = math.pi
                        gamma = math.atan2(t10,-t00)
        rot2 = "%.5f"%(alpha * 180/math.pi)
        tilt2 = "%.5f"%(beta * 180/math.pi)
        psi2 = "%.5f"%(gamma *180/math.pi)
        return rot2, tilt2, psi2

def sgn(x):
        if x > 0:
                t = 1
        elif x == 0:
                t = 0
        else:
                t = -1
        return t

if __name__ == "__main__":
	main()
