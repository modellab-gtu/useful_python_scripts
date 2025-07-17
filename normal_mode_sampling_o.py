#!/usr/bin/env python3
import numpy as np
import argparse
import os
from math import exp, sqrt

# Global constants
kB = 0.69503476  # Boltzmann constant in cm-1/K

atomic_number_to_symbol = {
    1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 16: "S", 17: "Cl", 35: "Br", 53: "I"
}

def parse_arguments():
    parser = argparse.ArgumentParser(description="Normal mode sampling from Gaussian log file.")
    parser.add_argument("log_file", help="Gaussian log file containing frequencies and normal modes.")
    parser.add_argument("--num_total_conformers", type=int, default=100, help="Total number of conformers to generate.")
    parser.add_argument("--scale", type=float, default=0.2, help="Scaling factor for displacement amplitude.")
    parser.add_argument("--T", type=float, default=298.15, help="Temperature in Kelvin.")
    parser.add_argument("--output_dir", type=str, default="conformers", help="Directory to save conformers.")
    parser.add_argument("--include_imaginary", action="store_true", help="Include imaginary mode in sampling.")
    parser.add_argument("--use_rms_amplitudes", action="store_true", help="Use RMS amplitude instead of Boltzmann weights.")
    parser.add_argument("--uniform_sampling", action="store_true", help="Uniformly distribute conformers across modes.")
    parser.add_argument("-fo", "--file_output", choices=["unified", "individual", "mode_by_mode"], default="unified", help="Output mode for XYZ files: 'unified' (default), 'individual' (one file per conformer), or 'mode_by_mode' (grouped by vibrational mode).")
    return parser.parse_args()

def extract_coordinates(log_file):
    coords = []
    atomic_numbers = []
    with open(log_file) as f:
        lines = f.readlines()
    start = False
    for line in lines:
        if "Standard orientation:" in line:
            coords = []
            atomic_numbers = []
            start = True
        elif start and "---" in line:
            if len(coords) > 0:
                break
            continue
        elif start:
            tokens = line.split()
            if len(tokens) == 6 and tokens[0].isdigit():
                coords.append([float(x) for x in tokens[-3:]])
                atomic_numbers.append(int(tokens[1]))
    coords = np.array(coords)
    print(f"\n[INFO] Extracted {len(coords)} atoms from '{log_file}'")
    for i, (Z, xyz) in enumerate(zip(atomic_numbers, coords)):
        print(f"  Atom {i+1:3d}: Z={Z:<2}  x={xyz[0]: .4f}  y={xyz[1]: .4f}  z={xyz[2]: .4f}")
    return coords, atomic_numbers
#    return np.array(coords), atomic_numbers

def extract_frequencies_and_modes(log_file, include_imaginary=False):
    with open(log_file) as f:
        lines = f.readlines()

    freqs = []
    displacements = []
    mode_blocks = []
    mode_block = []

    reading_block = False
    for line in lines:
        if "Frequencies --" in line:
            current_freqs = [float(x) for x in line.split()[2:]]
            freqs.extend(current_freqs)
            mode_blocks.append([])
        elif "Atom  AN" in line:
            reading_block = True
            continue
        elif reading_block:
            tokens = line.split()
            if not tokens or not tokens[0].isdigit():
                if mode_block:
                    mode_blocks[-1] = mode_block
                    mode_block = []
                reading_block = False
                continue
            floats = [float(x) for x in tokens[2:]]
            mode_block.append(floats)

    num_atoms = len(mode_blocks[0]) - 1
    num_modes = len(freqs)
    modes = []
    print(f"\n[INFO] Found {num_modes} modes (including imaginary if allowed)")
    for m in range(num_modes):
        displacement = []
        for a in range(num_atoms):
            vec = [
                mode_blocks[m // 3][a][3 * (m % 3) + 0],
                mode_blocks[m // 3][a][3 * (m % 3) + 1],
                mode_blocks[m // 3][a][3 * (m % 3) + 2],
            ]
            displacement.append(vec)
        modes.append(np.array(displacement))
        print(f"  Mode {m+1:3d}: Frequency = {freqs[m]: .2f} cm⁻¹  | Displacement (atom 1) = {displacement[0]}")
    if not include_imaginary:
        real_modes = [(f, d) for f, d in zip(freqs, modes) if f > 0]
        print(f"[INFO] Retained {len(real_modes)} real modes (> 0 cm⁻¹)")
    else:
        real_modes = list(zip(freqs, modes))
    return real_modes

def compute_rms_amplitudes(frequencies, T):
    amplitudes = []
    for freq in frequencies:
        if freq == 0.0:
            amplitudes.append(0.0)
            continue
        zpe = 0.5
        try:
            occupation = 1.0 / (exp(freq / (kB * T)) - 1.0)
        except OverflowError:
            occupation = 0.0
        E_vib = zpe + occupation
        amp = sqrt(2.0 * E_vib / freq)
        amplitudes.append(amp)
    return amplitudes

def boltzmann_weights(frequencies, T):
    weights = [exp(-f / (kB * T)) if f > 0 else 1.0 for f in frequencies]
    total = sum(weights)
    return [w / total for w in weights]

def generate_conformers(coords, modes, num_total, scale, T, use_rms=False, use_uniform=False):
    frequencies = [f for f, _ in modes]
    displacements = [d for _, d in modes]
    num_modes = len(modes)

    if use_uniform:
        print("Using uniform mode-based sampling.")
        base_count = num_total // num_modes
        counts = np.full(num_modes, base_count, dtype=int)
        remainder = num_total - base_count * num_modes
        if remainder > 0:
            counts[:remainder] += 1
        amplitudes = [scale] * num_modes
    elif use_rms:
        print("Using RMS amplitude-based sampling.")
        amplitudes = compute_rms_amplitudes(frequencies, T)
        total = sum(amplitudes)
        raw_counts = np.array(amplitudes) / total * num_total
        counts = np.floor(raw_counts).astype(int)
        remainder = num_total - np.sum(counts)
        if remainder > 0:
            fractional = raw_counts - counts
            for idx in np.argsort(fractional)[-remainder:]:
                counts[idx] += 1
    else:
        print("Using Boltzmann weight-based sampling.")
        weights = boltzmann_weights(frequencies, T)
        raw_counts = np.array(weights) * num_total
        counts = np.floor(raw_counts).astype(int)
        remainder = num_total - np.sum(counts)
        if remainder > 0:
            fractional = raw_counts - counts
            for idx in np.argsort(fractional)[-remainder:]:
                counts[idx] += 1
        amplitudes = [scale] * len(modes)

    conformers = []
    metadata = []
    for mode_idx, (disp, amp, count) in enumerate(zip(displacements, amplitudes, counts), start=1):
        if count == 1:
            displacements_linspace = [0.0]
        else:
            displacements_linspace = np.linspace(-1, 1, count)
        for step_idx, a in enumerate(displacements_linspace, start=1):
            conformers.append(coords + a * scale * amp * disp)
            metadata.append((mode_idx, step_idx, count))
    return conformers, metadata

#def save_xyz(conformers, atomic_numbers, output_dir, metadata):
#    if not os.path.exists(output_dir):
#        os.makedirs(output_dir)
#    combined_path = os.path.join(output_dir, "combined.xyz")
#    with open(combined_path, "w") as f:
#        for i, (conf, (mode_idx, step_idx, total_steps)) in enumerate(zip(conformers, metadata)):
#            f.write(f"{len(conf)}\nConformer {i+1} (mode {mode_idx}, step {step_idx} of {total_steps})\n")
#            for Z, xyz in zip(atomic_numbers, conf):
#                f.write(f"{atomic_number_to_symbol.get(Z, 'X')} {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f}\n")


def save_xyz(conformers, atomic_numbers, output_dir, metadata, file_output="unified"):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if file_output == "unified":
        combined_path = os.path.join(output_dir, "combined.xyz")
        with open(combined_path, "w") as f:
            for i, (conf, (mode_idx, step_idx, total_steps)) in enumerate(zip(conformers, metadata)):
                f.write(f"{len(conf)}\nConformer {i+1} (mode {mode_idx}, step {step_idx} of {total_steps})\n")
                for Z, xyz in zip(atomic_numbers, conf):
                    f.write(f"{atomic_number_to_symbol.get(Z, 'X')} {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f}\n")

    elif file_output == "individual":
        for i, (conf, (mode_idx, step_idx, total_steps)) in enumerate(zip(conformers, metadata)):
            fname = os.path.join(output_dir, f"conformer_{i+1:03d}.xyz")
            with open(fname, "w") as f:
                f.write(f"{len(conf)}\nConformer {i+1} (mode {mode_idx}, step {step_idx} of {total_steps})\n")
                for Z, xyz in zip(atomic_numbers, conf):
                    f.write(f"{atomic_number_to_symbol.get(Z, 'X')} {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f}\n")

    elif file_output == "mode_by_mode":
        # Group by mode
        mode_dict = {}
        for conf, (mode_idx, step_idx, total_steps) in zip(conformers, metadata):
            mode_dict.setdefault(mode_idx, []).append((conf, step_idx, total_steps))

        for mode_idx, conf_list in mode_dict.items():
            fname = os.path.join(output_dir, f"mode_{mode_idx:02d}.xyz")
            with open(fname, "w") as f:
                for i, (conf, step_idx, total_steps) in enumerate(conf_list, start=1):
                    f.write(f"{len(conf)}\nMode {mode_idx}, Conformer {i} (step {step_idx} of {total_steps})\n")
                    for Z, xyz in zip(atomic_numbers, conf):
                        f.write(f"{atomic_number_to_symbol.get(Z, 'X')} {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f}\n")


def main():
    args = parse_arguments()
    coords, atomic_numbers = extract_coordinates(args.log_file)
    modes = extract_frequencies_and_modes(args.log_file, include_imaginary=args.include_imaginary)
    conformers, metadata = generate_conformers(
        coords, modes, args.num_total_conformers, args.scale, args.T,
        use_rms=args.use_rms_amplitudes, use_uniform=args.uniform_sampling
    )
    save_xyz(conformers, atomic_numbers, args.output_dir, metadata,file_output=args.file_output)
    print(f"Generated {len(conformers)} conformers in '{args.output_dir}/combined.xyz'.")

if __name__ == "__main__":
    main()

