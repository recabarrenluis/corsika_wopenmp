import ROOT
import numpy as np
import argparse
from corsikaio import CorsikaParticleFile, as_dict

parser = argparse.ArgumentParser()
parser.add_argument('--corsika', type=str, help='DAT file')
parser.add_argument('--output', type=str, help='path_to_write')
args = parser.parse_args()

path_to_data = '/mnt/d/' #'/home/lure_wsl/5_HPC_exam/'
simfile = path_to_data + args.corsika
nev = 100
output_file = args.output + '/' + args.corsika + ".root"

# Open output ROOT file
f_out = ROOT.TFile(output_file, "RECREATE")
f_out.SetCompressionLevel(0)

tree = ROOT.TTree("events", "corsika c++ readable")

from ROOT import std

px = std.vector('float')()
py = std.vector('float')()
pz = std.vector('float')()
particle_type = std.vector('short')()
pos_x = std.vector('float')()
pos_y = std.vector('float')()
time = std.vector('float')()

tree.Branch("px", px)
tree.Branch("py", py)
tree.Branch("pz", pz)
tree.Branch("particle_type", particle_type)
tree.Branch("pos_x", pos_x)
tree.Branch("pos_y", pos_y)
tree.Branch("time", time)

with CorsikaParticleFile(simfile) as f:
    event0 = next(f)
    hdr0 = as_dict(event0.header)

    primary_energy = float(hdr0['total_energy'])
    primary_zenith = float(hdr0['zenith'])
    primary_azimuth = float(hdr0['azimuth'])
    primary_type = int(hdr0['particle_id'])

    print("=== PRIMARY PARTICLE INFO ===")
    print(f"Energy  = {primary_energy}")
    print(f"Zenith  = {primary_zenith}")
    print(f"Azimuth = {primary_azimuth}")
    print(f"Type    = {primary_type}")

    def fill_event(event):
        part = as_dict(event.particles)

        px_arr = np.asarray(part['px'], dtype=np.float32)
        py_arr = np.asarray(part['py'], dtype=np.float32)
        pz_arr = np.asarray(part['pz'], dtype=np.float32)
        pid_arr = np.asarray(part['particle_description'] * 1e-3, dtype=np.int16)
        pos_x_arr = np.asarray(part['y'] * 1e-2, dtype=np.float32)
        pos_y_arr = np.asarray(part['x'] * 1e-2, dtype=np.float32)
        time_arr = np.asarray(part['t'], dtype=np.float32)

        n = len(px_arr)

        px.clear()
        py.clear()
        pz.clear()
        particle_type.clear()
        pos_x.clear()
        pos_y.clear()
        time.clear()

        if px.capacity() < n:
            px.reserve(n)
            py.reserve(n)
            pz.reserve(n)
            particle_type.reserve(n)
            pos_x.reserve(n)
            pos_y.reserve(n)
            time.reserve(n)

        for i in range(n):
            px.push_back(float(px_arr[i]))
            py.push_back(float(py_arr[i]))
            pz.push_back(float(pz_arr[i]))
            particle_type.push_back(int(pid_arr[i]))
            pos_x.push_back(float(pos_x_arr[i]))
            pos_y.push_back(float(pos_y_arr[i]))
            time.push_back(float(time_arr[i]))

        tree.Fill()

    fill_event(event0)

    for _ in range(1, nev):
        fill_event(next(f))

    ROOT.TParameter(float)("primary_energy", primary_energy).Write()
    ROOT.TParameter(float)("primary_zenith", primary_zenith).Write()
    ROOT.TParameter(float)("primary_azimuth", primary_azimuth).Write()
    ROOT.TParameter(int)("primary_type", primary_type).Write()

f_out.Write()
f_out.Close()
print(f"ROOT file written: {output_file}")
