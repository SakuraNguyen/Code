#!/usr/bin/env python

import numpy as np
import sys
import math
import MDAnalysis as mda 

u = mda.Universe(sys.argv[1], sys.argv[2])

donors = u.select_atoms(sys.argv[3])
immediates = u.select_atoms(sys.argv[4])
acceptors = u.select_atoms(sys.argv[5])

cutoff = float(sys.argv[6])
angle = float(sys.argv[7]))
def abc_to_hmatrix(a, b, c, alpha, beta, gamma, degrees=True):
    if degrees:
        alpha, beta, gamma = map(math.radians, (alpha, beta, gamma))
    result = np.zeros((3, 3))
    a = np.array((a, 0, 0))
    b = b*np.array((math.cos(gamma), math.sin(gamma),0))
    bracket = (math.cos(alpha)-math.cos(beta)*math.cos(gamma))/math.sin(gamma)
    c = c*np.array((math.cos(beta), bracket, math.sin(beta)**2-bracket**2))
    result[:, 0] = a
    result[:, 1] = b
    result[:, 2] = c
    return result

def smallest_distance_to(u, pos, ag):
	hmat = abc_to_hmatrix(*u.dimensions)
	result = np.ones(len(ag))*1000
	resultpos = np.zeros((len(ag), 3))
	rpos = []
	for x in range(-1, 2):
		for y in range(-1, 2):
			for z in range(-1, 2):
				rpos.append(pos + x*hmat[:, 0] + y*hmat[:, 1] + z*hmat[:, 2])
	for idx in range(len(ag)):
		ds = np.linalg.norm(rpos-ag.get_positions()[idx], axis=1)
		f, = np.where(ds > 0)
		m = min(ds[f])
		result[idx] = min(result[idx], m)
		f, = np.where(ds == m)
		resultpos[idx] = rpos[f[0]]
	return result, resultpos

def _angle_between(a, b):
    a=np.array(a, dtype=np.float32)
    b=np.array(b, dtype=np.float32)
    try:
        a /= np.linalg.norm(a)
        b /= np.linalg.norm(b)
    except:
        raise ValueError('Got zero vector.')
    angle = np.arccos(np.dot(a, b))
    if np.isnan(angle):
        if (a == b).all():
            return 0.0
        else:
            return np.pi
    return angle

tmatches = dict()
for ts in u.trajectory[:]:
	count=0
	#print ts.frame
	matches = dict()
	for acc in acceptors:
		# distance criterion
		adist, rpos = smallest_distance_to(u, acc.pos, donors)
		#print donor.name, adist
		acceptable, = np.where(adist < cutoff)
		
		# angle criterion
		for candidate in acceptable:
			for hydrogen in immediates:
				a1 = donors[candidate].pos - hydrogen.pos
				a2 = rpos[candidate] - hydrogen.pos
				ta = _angle_between(a1, a2)
				#print acceptors[candidate].name,hydrogen.name, np.rad2deg(_angle_between(a1, a2))
				if np.rad2deg(ta) > 180-angle:	
					#print 'match', acc.name, hydrogen.name, donors[candidate].name, np.rad2deg(_angle_between(a1, a2))
					key = '-'.join([u.atoms[_].name for _ in ([donors[candidate].number, hydrogen.number, acc.number])])
					count+=1
					if key in matches:
						matches[key] += 1
					else:
						#print ta, np.rad2deg(ta)
						matches[key] = 1#(acc.number, hydrogen.number, donors[candidate].number)
	for e in matches:
		if e in tmatches:
			tmatches[e] += matches[e]
		else:
			tmatches[e] = matches[e]
	print ts.frame, count
for e in tmatches:
	print e, float(tmatches[e])#/len(u.trajectory)*100
