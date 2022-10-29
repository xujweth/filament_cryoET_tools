#!/usr/bin/python
#This script is written by Jingwei Xu, ETH Zurich
#Modify		2022/03/05	fix the issue to orient sub-volume
#				update the script on dynamo table part
#		2022/03/14	update the script for the central projection from the refined particles (not only filament)
#		2022/04/01	fix the bug to calculate the psi angle (180 - psi_cal or -1*(180 + psi_cal))
#		2022/04/04	fix the bug about the psi definition difference between EMAN and Relion
#		2022/04/12	add the function to align with reference
#		2022/04/14	add the function to mask out the density out of filament.
#		2022/06/13	fix the bug for non_helix when using dynamo tbl
try:
        from optparse import OptionParser
except:
        from optik import OptionParser

import sys, os, math, time
import numpy as np
from EMAN2 import *

def main():
        usage = """python subtomo_helixtool.py [starfile] --options\n
		For example:
		python subtomo_helixtool_v2.py exampleData/crop.tbl --dynamo --only_center --percentage 0.4 --prior --crop_folder exampleData/ --column20 exampleData/indices_column20.doc
		python subtomo_helixtool_v2.py crop.tbl --dynamo --only_center --prior --crop_folder . --column20 indices_column20.doc --non_helix --percentage 0.2
		python ~/JX_data/scripts/subtomo_helixtool_v2.py actin_bin2_ctf_test/crop.tbl --dynamo --only_center --percentage 0.2 --prior --crop_folder actin_bin2_ctf_test/ --column20 actin_bin2_ctf_test/indices_column20.doc
		python ~/JX_data/scripts/subtomo_helixtool_v2.py particles_bin1/crop.tbl --dynamo --only_center --prior --crop_folder particles_bin1 --column20 particles_bin1/indices_column20.doc --angpix 4.51 --align --ref actin_mask.mrc --binning 2\n"""
	parser = OptionParser(usage=usage)
	parser.add_option("--dynamo", dest="dynamo", action="store_true", help="Is the input file from dynamo table? The default is False.", default=False)
	parser.add_option("--only_center", dest="only_center", action="store_true", help="To only use center slice to do the projection?", default=False)
	parser.add_option("--percentage", dest="percentage", type="float", help="The percentage of center slice to do the projection. The default is 40%.", default=0.4)
	parser.add_option("--prior", dest="prior", action="store_true", help="To use prior knowledge of tilt/psi? The default is False.", default=False)
	parser.add_option("--crop_folder", dest="crop_folder", type="string", help="The name of crop folder from dynamo. The default is none.", default="")
	parser.add_option("--column20", dest="column20_file", type="string", help="The name of column20 file from dynamo cropping. The default is none.", default="")
	parser.add_option("--non_helix", dest="non_helix", action="store_true", help="To process the orientation file from non-filament data? Or use the orientation determined from dynamo? The default is False.", default=False)
	parser.add_option("--no_invert", dest="no_invert", action="store_true", help="Do not invert contrast? The default is False.", default=False)
	parser.add_option("--angpix", dest="angpix", type="float", help="The pixel value of subvolume. The default is 1.0.", default=1.0)
	parser.add_option("--align", dest="align", action="store_true", help="To align with the reference before projection? The default is False.", default=False)
	parser.add_option("--ref", dest="ref", type="string", help="The name of the refernece. The default is none.", default="")
	parser.add_option("--binning", dest="binning", type="int", help="The binning level to align with reference. The default is 2.", default=2)
	parser.add_option("--mask", dest="mask", type="string", help="The mask file to mask out the density out of filament. The default is none.", default="")
	(options, args) = parser.parse_args()

	if len(args) < 1:
		print "ERROR: please provide the input file for processing. Exit!"
		print usage
		sys.exit(-1)

	ptclfile = args[0]
	only_center = options.only_center
	z_per = options.percentage
	dynamo = options.dynamo
	prior = options.prior

	ptcltxt = open(ptclfile, "r").readlines()
	optic_info = []
	ptcl_info_lst = []
	current_directory = os.getcwd()
	crop_folder = options.crop_folder
	non_helix = options.non_helix
	angpix = options.angpix
	mask = options.mask

	if not(dynamo):
		optic_info, ptcl_info = fetch_optic(ptcltxt)
		angpix_index = fetch_index(optic_info, "_rlnImagePixelSize")
		ptcl_index = fetch_index(ptcl_info, "_rlnImageName")
		rot_index = fetch_index(ptcl_info, "_rlnAngleRot")
		tilt_index = fetch_index(ptcl_info, "_rlnAngleTilt")
		psi_index = fetch_index(ptcl_info, "_rlnAnglePsi")
#Please check the label
		tomo_index = fetch_index(ptcl_info, "_rlnMicrographName")
		corx_index = fetch_index(ptcl_info, "_rlnCoordinateX")
		cory_index = fetch_index(ptcl_info, "_rlnCoordinateY")
		shiftx_index = fetch_index(ptcl_info, "_rlnOriginXAngst")
		shifty_index = fetch_index(ptcl_info, "_rlnOriginYAngst")
		shiftz_index = fetch_index(ptcl_info, "_rlnOriginZAngst")
		randomset_index = fetch_index(ptcl_info, "_rlnRandomSubset")
		optic_group_index = fetch_index(ptcl_info, "_rlnOpticsGroup")
		helixID_index = helix_track_index = psi_prior_flip_index = 0
		if prior:
			tilt_index = fetch_index(ptcl_info, "_rlnAngleTiltPrior")
			psi_index = fetch_index(ptcl_info, "_rlnAnglePsiPrior")
			helixID_index = fetch_index(ptcl_info, "_rlnHelicalTubeID")
			helix_track_index = fetch_index(ptcl_info, "_rlnHelicalTrackLength")
			psi_prior_flip_index = fetch_index(ptcl_info, "_rlnAnglePsiFlipRatio")
	
		for i in optic_info:
			if len(i.split()) < 3 or i.startswith('#'):
				continue
			angpix = float(i.split()[angpix_index])
		for i in ptcl_info:
			if len(i.split()) < 3:
				continue
#			ptcl = "%s/%s"%(current_directory, i.split()[ptcl_index])
			ptcl = i.split()[ptcl_index]
			tomo_name = i.split()[tomo_index]
			corx = i.split()[corx_index]
			cory = i.split()[cory_index]
			randomset = i.split()[randomset_index]
			optic_group = i.split()[optic_group_index]
			rot = float(i.split()[rot_index])
			tilt = float(i.split()[tilt_index])
			psi = float(i.split()[psi_index])
			helixID = helix_track = psi_prior_flip = 0
			shiftx = float(i.split()[shiftx_index]) / options.angpix
			shifty = float(i.split()[shifty_index]) / options.angpix
			shiftz = float(i.split()[shiftz_index]) / options.angpix
			if prior:
				helixID = int(i.split()[helixID_index])
				helix_track = float(i.split()[helix_track_index])
				psi_prior_flip = float(i.split()[psi_prior_flip_index])
			if not(os.path.exists(ptcl)):
				print "ERROR: the pathway of particle %s seems not correct. Please check it. Exit!"%(ptcl)
			info = [ptcl, rot, tilt, psi, helixID, helix_track, psi_prior_flip, corx, cory, randomset, tomo_name, shiftx, shifty, shiftz, optic_group]
			ptcl_info_lst.append(info)
	else:
#Modified 2022/03/05
#Rewrite this part
#fix the psi/tilt angle based on the start/end points in the filemant
		ptcl_dynamo_lst = []
		column20_start = -1
		tomo_ptcl_lst = []
		randomset = 0
		num = 1
		column20_file = options.column20_file
		column20_info = open(column20_file, "r").readlines()
		column20_lst = []
		for i in column20_info:
			index = int(i.split()[0])
			tomo_name = i.split()[1].split('/')[-1]
			column20_lst.append((index, tomo_name))

#Sort particles first based on the column20 (tomograms)
		for i in ptcltxt:
			ptcl_index = i.split()[0]
			dx, dy, dz = i.split()[3:6]
			posx, posy, posz = i.split()[23:26]
#Currently do not consider binning
			posx_star = float(posx) + float(dx)
			posy_star = float(posy) + float(dy)
			posz_star = float(posz) + float(dz)
			tdrot = i.split()[6]
			tilt = i.split()[7]
			narot = i.split()[8]
			column20 = int(i.split()[19])
			column21 = i.split()[20]
			tomo_name = ""
			ptcl = "%s/particle_%s.em"%(crop_folder, ptcl_index.zfill(6))
			for line in column20_lst:
				if column20 == line[0]:
					tomo_name = line[1]
			if not(os.path.exists(ptcl)):
				print "ERROR: the pathway of particle %s seems not correct. Please check it. Exit!"%(ptcl)
			rotation_matrix = euler2matrix(float(tdrot), float(tilt), float(narot), True)
			rot, tilt, psi = matrix2euler(rotation_matrix)
			if num%2 == 1:
				randomset = 1
			else:
				randomset = 2
			num += 1
			info = [ptcl, float(tdrot), float(posx_star), float(posy_star), float(posz_star), randomset, tomo_name, float(dx), float(dy), float(dz), int(column21)]
			if non_helix:
				info = [ptcl, float(rot), float(posx_star), float(posy_star), float(posz_star), randomset, tomo_name, float(tilt), float(psi), float(dx), float(dy), float(dz), int(column21)]
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
				column21 = i[-1]
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
				xstart_px = subtube_lst[0][2]
				ystart_px = subtube_lst[0][3]
				zstart_px = subtube_lst[0][4]
				dist = 0.0
				xend_px = subtube_lst[-1][2]
				yend_px = subtube_lst[-1][3]
				zend_px = subtube_lst[-1][4]
				diff_x = xend_px - xstart_px
				diff_y = yend_px - ystart_px
				diff_z = zend_px - zstart_px
				length = math.sqrt(diff_x**2 + diff_y**2 + diff_z**2)
				psi_prior = tilt_prior = 0.0
				if not(non_helix):
					psi_cal = math.atan2(diff_y, diff_x) / math.pi * 180
#Modified: 2022/04/01
#Fix the bug on psi angle
					if psi_cal > 0:
						psi_prior = 180 - psi_cal
					else:
						psi_prior = -1 * (180 + psi_cal)
				
				if not(non_helix):
					tilt_prior = math.acos(diff_z/length) / math.pi * 180
				for num, ptcl in enumerate(subtube_lst):
					posx_star = ptcl[2]
					posy_star = ptcl[3]
					posz_star = ptcl[4]
					rot = ptcl[1]
					randomset = ptcl[5]
					tomo_name = ptcl[6]
					helical_tube_ID = ptcl[-1]
					particle_helical_track_length = math.sqrt((posx_star - xstart_px)**2 + (posy_star - ystart_px)**2 + (posz_star - zstart_px)**2)
					tilt_ptcl = psi_ptcl = 0.000000
					if non_helix:
						tilt_ptcl = ptcl[7]
						psi_ptcl = ptcl[8]
					else:
						if num == 0:
							tilt_ptcl = tilt_prior
							psi_ptcl = psi_prior
						else:
							posx_last = subtube_lst[num-1][2]
							posy_last = subtube_lst[num-1][3]
							posz_last = subtube_lst[num-1][4]
							diffx = posx_star - posx_last
							diffy = posy_star - posy_last
							diffz = posz_star - posz_last
#							print diffz, particle_helical_track_length, ptcl
#Modified: 2022/04/08
#Fix the bug for tilt_ptcl: it should be length between two neighboring points, not particle_helical_trach_length!
							tilt_ptcl = math.acos(diffz/math.sqrt((diffx)**2 + (diffy)**2 + (diffz)**2)) / math.pi * 180
							psi_cal = math.atan2(diffy, diffx) / math.pi * 180
							if psi_cal > 0:
								psi_ptcl = 180 - psi_cal
							else:
								psi_ptcl = -1 * (180 + psi_cal)

					psi_prior_flip_ratio = 0.500000
					ptcl_name = ptcl[0]
					optic_group = 1
					shiftx, shifty, shiftz = ptcl[7:10]
					if non_helix:
						shiftx, shifty, shiftz = ptcl[9:12]
					info = [ptcl_name, rot, tilt_ptcl, psi_ptcl, helical_tube_ID, particle_helical_track_length, psi_prior_flip_ratio, posx_star, posy_star, randomset, tomo_name, shiftx, shifty, shiftz, optic_group]
					ptcl_info_lst.append(info)
#					print info[2], info[3]

	print "There are %d particles in file %s."%(len(ptcl_info_lst), ptclfile)

	d1 = EMData()
	d1.read_image(ptcl_info_lst[0][0])
	mapsize = d1.get_zsize()
	z_limit = int(mapsize * z_per)
	if z_limit %2 != 0:
		z_limit += 1
	
	half_z_limit = z_limit / 2
	z_start = mapsize/2 - half_z_limit
	z_end = mapsize/2 + half_z_limit
	if only_center:
		print "Only %d center slices of volume: %d-%d will be used for projection."%(z_limit, z_start, z_end)
	write_num = 0	
	for i in ptcl_info_lst:
		ptcl = i[0]
		ptcl_prefix = ptcl.split('/')[-1].split('.')[0]
		ptcl_path = ptcl.split(ptcl_prefix)[0]
		rot, tilt, psi = i[1:4]
		shiftx, shifty, shiftz = i[11:14]
		d = EMData()
		d.read_image(ptcl)
#Dynamo: shift + rotation; relion/EMAN: rotation + shift
#Here for Dynamo
		if dynamo:
			t0 = Transform()
			t0.set_trans((-1 * shiftx, -1 * shifty, -1 *shiftz))
			d.process_inplace("xform", {"transform":t0})
#		d.write_image("temp.mrc")
		t = Transform()
		t.set_params({"type":"spider", "phi":rot, "theta":tilt, "psi":psi})
		t2 = t.get_params("eman")
		invert_t = Transform()
#to orient particle along z axis first and then make it along Y? axis by rotating 90 degree
		invert_t.set_rotation({"type":"eman", "az": 0, "alt":0, "phi":-1 * t2["phi"]})
		d.transform(invert_t)
		invert_t.set_rotation({"type":"eman", "az": 0, "alt": -1 * t2["alt"], "phi":0})
		d.transform(invert_t)
#Here for Relion/EMAN shift
		if not(dynamo):
			t0 = Transform()
			t0.set_trans((-1 * shiftx, -1 * shifty, -1 *shiftz))
			d.process_inplace("xform", {"transform":t0})
#Modified 2022/03/09    To invert the contrast
		if not(options.no_invert):
			d.mult(-1)
#Modified 2022/04/11 Add align function
		if options.align:
			binning = options.binning
			ref = EMData(options.ref)
			d2 = d.copy()
#Lowpass filter 50 Angstrom for rough alignment (only shift!)
			d2.process_inplace("filter.lowpass.gauss", {"cutoff_freq":1.0/50, "apix": options.angpix*binning})
			if binning > 1:
				d2.process_inplace("math.meanshrink", {"n":binning})
			ref2 = ref.copy()
			ref2.process_inplace("filter.lowpass.gauss", {"cutoff_freq":1.0/50, "apix": options.angpix*binning})
			if binning > 1:
				ref2.process_inplace("math.meanshrink", {"n":binning})
#                       t_align = d2.align("refine_3d_grid", ref2, {'delta':0, 'range':2, 'search': 1, 'verbose':1}, "ccc.tomo")
#Only consider translation on x/y
			t_align = d2.align("rotate_translate_3d_grid", ref2, {'phi0':0, 'phi1':1, 'alt0':0, 'alt1':1, 'az0':0, 'az1':1, 'dphi':2, 'daz':2, 'dalt':2, 'searchx':1, 'searchy':1, 'searchz':0, 'verbose':0})
			t_aln = t_align["xform.align3d"].get_params("eman")
			t0 = Transform()
			t0.set_trans((-1 *float(t_aln["tx"])*binning, -1 * float(t_aln["ty"])*binning, -1 * float(t_aln["tz"]) * binning))
			d.process_inplace("xform", {"transform":t0})
#Modified 2022/04/14 add function of mask
		if len(mask) != 0:
			m = EMData(mask)
			d2 = d.copy()
			d = d2 * m		
	
#Modified 2022/04/04 Relion psi 0 is equal to 90 in EMAN?
		rot_t = Transform({"type":"eman", "az":0, "alt":90, "phi":0})
		d.transform(rot_t)
		rot_psi = Transform({"type":"eman", "az":0, "alt":0, "phi":90})
		d.transform(rot_psi)
		d.process_inplace("normalize.edgemean")
#The projection orientation will be only alt=90 compared to the reference.

		transformed_tmp = "%s_tmp_transformed.mrcs"%(ptcl_prefix)
		d.write_image(transformed_tmp)
		
		if only_center:
			part_tmp = "%s_tmp_center.mrcs"%(ptcl_prefix)
			for z in range(z_limit):
				z_num = z_start + z - 1
				transformed_map = EMData()
				transformed_map.read_image(transformed_tmp, z_num)
				transformed_map.write_image(part_tmp, z)
			transformed_tmp = part_tmp
		target_map = "%s.mrc"%(transformed_tmp.split('.mrcs')[0])
		os.rename(transformed_tmp, target_map)
		d2 = EMData()
		d2.read_image(target_map)
		proj = Transform()
		proj.set_params({"type":"eman", "az":0, "alt":0, "phi":0})
		proj_image = d2.project("standard", proj)
#		out_proj = "%s%s_proj.mrc"%(ptcl_path, ptcl_prefix)
		out_proj = "projections.mrcs"
		proj_image.write_image(out_proj, write_num)
		if only_center:
			temp = "%s_tmp_transformed.mrcs"%(ptcl_prefix)
			os.remove(temp)
		os.remove(target_map)
		i[0] = "%s@%s"%(str(write_num + 1).zfill(6), out_proj)
#The tilt and psi is changed to 0
#Modified tilt should be 90
		i[2] = 90.000000
		i[3] = 0.000000
		write_num += 1
	
#Write out the star file
	out_star = "%s_helix.star"%(ptclfile.split('.')[0])
	out_txt = open(out_star, "w")
	if not dynamo:
		if len(optic_info) != 0:
			for i in optic_info:
				out_txt.write(i)
		else:
			out_txt.write("\n# version 30001\n\ndata_optics\n\nloop_\n_rlnOpticsGroup #1\n_rlnOpticsGroupName #2\n_rlnSphericalAberration #3\n_rlnVoltage #4\n_rlnMicrographOriginalPixelSize #5\n_rlnAmplitudeContrast #6\n_rlnImageDimensionality #7\n_rlnImagePixelSize #8\n_rlnImageSize #9\n")
			out_txt.write("1\topticsGroup1\t2.700000\t300.000000\t4.510000\t0.100000\t2\t9.020000\t48\n")
	out_txt.write("\n# version 30001\n\ndata_particles\n\nloop_\n_rlnImageName #1\n_rlnAngleRot #2\n_rlnAngleTiltPrior #3\n_rlnAnglePsiPrior #4\n_rlnHelicalTubeID #5\n_rlnHelicalTrackLength #6\n_rlnAnglePsiFlipRatio #7\n_rlnCoordinateX #8\n_rlnCoordinateY #9\n_rlnRandomSubset #10\n_rlnMicrographName #11\n_rlnOriginXAngst #12\n_rlnOriginYAngst #13\n_rlnOriginZAngst #14\n_rlnOpticsGroup #15\n")
	for i in ptcl_info_lst:
		ptcl_update = ""
		for num, item in enumerate(i):
			if num == 11:
				ptcl_update += "0.000000\t"
			elif num == 12:
				ptcl_update += "0.000000\t"
			elif num == 13:
				ptcl_update += "0.000000\t"
			else:
				ptcl_update += "%s\t"%(item)
		ptcl_update += "\n"
		out_txt.write(ptcl_update)
		

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

def matrix2euler(matrix):
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
                gamma = math.atan2(t12, -t02)
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

		

