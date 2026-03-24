import ROOT
import numpy as np
import argparse
import ctypes
from array import array
from ROOT import std
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
    # First event header (primary info)
    event0 = next(f)
    hdr0 = as_dict(event0.header)

    primary_energy = float(hdr0['total_energy'])
    primary_zenith = float(hdr0['zenith'])
    primary_azimuth = float(hdr0['azimuth'])
    primary_type = int(hdr0['particle_id'])

    print("=== PRIMARY INFO (first event) ===")
    print(f"Energy  = {primary_energy}")
    print(f"Zenith  = {primary_zenith}")
    print(f"Azimuth = {primary_azimuth}")
    print(f"Type    = {primary_type}")

    part = as_dict(event0.particles)

    energy.clear()
    particle_type.clear()
    pos_x.clear()
    pos_y.clear()
    time.clear()

    e = sec_energy(part['px'], part['py'], part['pz'])
    for i in range(len(e)):
        energy.push_back(float(e[i]))
        particle_type.push_back(int(part['particle_description'][i] * 1e-3))
        pos_x.push_back(float(part['y'][i] * 1e-2))
        pos_y.push_back(float(part['x'][i] * 1e-2))
        time.push_back(float(part['t'][i]))

    tree.Fill()
    
    for ievt in range(1, nev):
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

    # Save primary info
    ROOT.TParameter(float)("primary_energy", primary_energy).Write()
    ROOT.TParameter(float)("primary_zenith", primary_zenith).Write()
    ROOT.TParameter(float)("primary_azimuth", primary_azimuth).Write()
    ROOT.TParameter(int)("primary_type", primary_type).Write()

# Write file
f_out.Write()
f_out.Close()
print(f"ROOT file written: {output_file}")
