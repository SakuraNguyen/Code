# hb.py
"""Calculates X-H...Y hydrogen-bonding occupancy."""

import glob
import math
import os
import time


class HB:
    def __init__(self, models):
        self.version = "2018-02-13"
        self.cwd = os.getcwd()
        self.models = models
    def method1(self):
        """Creates timestamp."""
        year = time.localtime()[0]
        month = time.localtime()[1]
        day = time.localtime()[2]
        self.timestamp = "{0}-{1}{2}-{3}{4}".format(year, "0" * (2 - len(str(month))), month, "0" * (2 - len(str(day))), day)
    def method2(self):
        """Defines hydrogen bond(s)."""
        self.D1 = {}
        self.D1["01"] = {} # HB1: A-Lys30-NZ-HZ1 ... O-29-Leu-B
        self.D1["01"]["A"] = {}
        self.D1["01"]["A"]["resid"] = 30
        self.D1["01"]["A"]["atid1"] = 399 # NZ
        self.D1["01"]["A"]["atid2"] = 400 # HZ1
        self.D1["01"]["B"] = {}
        self.D1["01"]["B"]["resid"] = 29
        self.D1["01"]["B"]["atid3"] = 831 # O
        self.D1["02"] = {} # HB2: A-Lys30-NZ-HZ2 ... O-29-Leu-B
        self.D1["02"]["A"] = {}
        self.D1["02"]["A"]["resid"] = 30
        self.D1["02"]["A"]["atid1"] = 399 # NZ
        self.D1["02"]["A"]["atid2"] = 401 # HZ2
        self.D1["02"]["B"] = {}
        self.D1["02"]["B"]["resid"] = 29
        self.D1["02"]["B"]["atid3"] = 831 # O
        self.D1["03"] = {} # HB3: A-Lys30-NZ-HZ2 ... O-29-Leu-B
        self.D1["03"]["A"] = {}
        self.D1["03"]["A"]["resid"] = 30
        self.D1["03"]["A"]["atid1"] = 399 # NZ
        self.D1["03"]["A"]["atid2"] = 402 # HZ3
        self.D1["03"]["B"] = {}
        self.D1["03"]["B"]["resid"] = 29
        self.D1["03"]["B"]["atid3"] = 831 # O
        self.hbs = sorted([x for x in self.D1])
    def method3(self):
        """Calculates hydrogen-bonding distances (H...Y)."""
        self.D2 = {}
        for hb in self.hbs:
            self.D2[hb] = []
            for model in glob.glob(self.models):
                for line in open(model):
                    if line.split()[4] == "A":
                        if int(line.split()[1]) == self.D1[hb]["A"]["atid2"]:
                            atx2 = float(line.split()[-3])
                            aty2 = float(line.split()[-2])
                            atz2 = float(line.split()[-1])
                    if line.split()[4] == "B":
                        if int(line.split()[1]) == self.D1[hb]["B"]["atid3"]:
                            atx3 = float(line.split()[-3])
                            aty3 = float(line.split()[-2])
                            atz3 = float(line.split()[-1])
                distance = ((atx3 - atx2)**2 + (aty3 - aty2)**2 + (atz3 - atz2)**2)**0.5
                self.D2[hb].append(distance)
    def method4(self):
        """Calculates hydrogen-bonding angles (X-H...Y)."""
        self.D3 = {}
        for hb in self.hbs:
            self.D3[hb] = []
            for model in glob.glob(self.models):
                for line in open(model):
                    if line.split()[4] == "A":
                        if int(line.split()[1]) == self.D1[hb]["A"]["atid1"]:
                            atx1 = float(line.split()[-3])
                            aty1 = float(line.split()[-2])
                            atz1 = float(line.split()[-1])
                    if line.split()[4] == "A":
                        if int(line.split()[1]) == self.D1[hb]["A"]["atid2"]:
                            atx2 = float(line.split()[-3])
                            aty2 = float(line.split()[-2])
                            atz2 = float(line.split()[-1])
                    if line.split()[4] == "B":
                        if int(line.split()[1]) == self.D1[hb]["B"]["atid3"]:
                            atx3 = float(line.split()[-3])
                            aty3 = float(line.split()[-2])
                            atz3 = float(line.split()[-1])
                dist1x2 = atx1 - atx2
                dist1y2 = aty1 - aty2
                dist1z2 = atz1 - atz2
                dist3x2 = atx3 - atx2
                dist3y2 = aty3 - aty2
                dist3z2 = atz3 - atz2
                dist1xyz2 = (dist1x2**2 + dist1y2**2 + dist1z2**2)**0.5
                dist3xyz2 = (dist3x2**2 + dist3y2**2 + dist3z2**2)**0.5
                dist1xyz2i = 1 / dist1xyz2
                dist3xyz2i = 1 / dist3xyz2
                costheta = ((dist1x2 * dist3x2) + (dist1y2 * dist3y2) + (dist1z2 * dist3z2)) * dist1xyz2i * dist3xyz2i
                angle_rad = math.acos(costheta) # [angle] = rad
                angle_deg = math.degrees(angle_rad) # [angle] = deg
                self.D3[hb].append(angle_deg)
    def method5(self, cutoffs):
        """Calculates hydrogen-bonding presence (1) and absence (0) for distances and angles."""
        self.D4 = {}
        for hb in self.hbs:
            self.D4[hb] = {}
            self.D4[hb]["distance"] = [1 if x < cutoffs["distance"] else 0 for x in self.D2[hb]]
            self.D4[hb]["angle"] = [1 if x > cutoffs["angle"] else 0 for x in self.D3[hb]]
    def method6(self):
        """Calculates the presence of individual hydrogen bonds (1 = hydrogen bond is present; 0 = hydrogen bond is absent)."""
        self.D5 = {}
        for hb in self.hbs:
            L1 = list(zip(self.D4[hb]["distance"], self.D4[hb]["angle"]))
            L2 = [sum(x) for x in L1]
            self.D5[hb] = [1 if x == 2 else 0 for x in L2]
        print("method6", self.D5)
    def method7(self):
        """Calculates total hydrogen-bonding presence (1) and absence (0)."""
        self.D6 = {}
        L1 = list(zip(self.D5[self.hbs[0]], self.D5[self.hbs[1]], self.D5[self.hbs[2]]))
        self.D6["/".join(self.hbs)] = [max(x) for x in L1]
        print("method7", self.D6)
    def method8(self):
        """Calculates hydrogen-bonding occupancies (1.0 = 100.0 %; 0.0 = 0.0 %)."""
        self.D7 = {}
        for hbonds in self.D6:
            self.D7[hbonds] = {}
            self.D7[hbonds] = sum(self.D6[hbonds]) / float(len(self.D6[hbonds]))
        print("method8", self.D7)


i1 = HB(models="*.pdb")
i1.method1()
i1.method2()
i1.method3()
i1.method4()
i1.method5(cutoffs={"distance": 2.5, "angle": 120.0})
i1.method6()
i1.method7()
i1.method8()
