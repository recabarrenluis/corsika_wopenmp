import ROOT
import numpy as np
import argparse
from corsikaio import CorsikaParticleFile, as_dict
from numba import njit

@njit(nopython=True)
def sec_energy(px, py, pz):
    return np.sqrt(px**2 + py**2 + pz**2)

parser = argparse.ArgumentParser()
parser.add_argument('--corsika', type=str, help='DAT file')
parser.add_argument('--output', type=str, help='path_to_write')
args = parser.parse_args()

path_to_data = '/home/lure_wsl/5_HPC_exam/' #edit yours
simfile = path_to_data + args.corsika
nev = 100
output_file = args.output + '/' + args.corsika + ".root"

# Open output ROOT file
f_out = ROOT.TFile(output_file, "RECREATE")
tree = ROOT.TTree("events", "corsika c++ readable")

# Define C++ vectors
from array import array
import ctypes

from ROOT import std

energy = std.vector('float')()
particle_type = std.vector('short')()
pos_x = std.vector('float')()
pos_y = std.vector('float')()
time = std.vector('float')()

# Set branches
tree.Branch("energy", energy)
tree.Branch("particle_type", particle_type)
tree.Branch("pos_x", pos_x)
tree.Branch("pos_y", pos_y)
tree.Branch("time", time)

# Fill tree
with CorsikaParticleFile(simfile) as f:
    for ievt in range(nev):
        event = next(f)
        part = as_dict(event.particles)

        # Clear vectors
        energy.clear()
        particle_type.clear()
        pos_x.clear()
        pos_y.clear()
        time.clear()

        # Fill vectors
        e = sec_energy(part['px'], part['py'], part['pz'])
        for i in range(len(e)):
            energy.push_back(float(e[i]))
            particle_type.push_back(int(part['particle_description'][i] * 1e-3))
            pos_x.push_back(float(part['y'][i] * 1e-2))
            pos_y.push_back(float(part['x'][i] * 1e-2))
            time.push_back(float(part['t'][i]))

        tree.Fill()

# Write file
f_out.Write()
f_out.Close()
print(f"ROOT file written: {output_file}")
