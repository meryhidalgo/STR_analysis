#!/usr/bin/env python3
import pysam
import regex as re
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
import random
import glob
import os

def calculate_mode(values):
    counts = {}
    for v in values:
        if v in counts:
            counts[v] += 1
        else:
            counts[v] = 1
    # Encontrar el valor con la frecuencia máxima
    max_freq = max(counts.values())
    moda = [k for k, v in counts.items() if v == max_freq]
    return moda

def reading_bam(bam_file, motif, qscore_threshold=20, minimum=9):
	"""
	Reads a BAM file and returns the counts of consecutive motifs in the specified region, only those found more than 3 times (CAGCAGCAG)
	"""
	# === LECTURA DEL BAM ===
	bam = pysam.AlignmentFile(bam_file, "rb")

	repeat_counts_per_motif = defaultdict(list)
	total_read_count = 0
	q_read_count = 0

	for read in bam.fetch():
		total_read_count += 1
		#if read.is_unmapped or read.mapping_quality < min_mapq:
		if read.mapping_quality < qscore_threshold:
			q_read_count += 1
			continue
		seq = read.query_sequence
		if seq:
			matches = re.findall(f"(?:{motif})+", seq, overlapped=False)
			if matches:
				total_repeats = sum(len(m) // len(motif) for m in matches if len(m) >= minimum)
				if total_repeats > 0:
					repeat_counts_per_motif[motif].append(total_repeats)

	sample_name = os.path.splitext(os.path.basename(bam_file))[0]
	output_txt = os.path.join("stats", f"{sample_name}_stats.txt")
	with open(output_txt, "w") as f:
		print(f"Estadísticas para {sample_name}", file=f)
		print(f"Total reads: {total_read_count}", file=f)
		print(f"Q reads: {q_read_count} with Q score < {qscore_threshold}", file=f)
		print("\nRepeticiones encontradas (motivo -> lista de cuentas):", file=f)
		if repeat_counts_per_motif[motif]:
			mean_value = sum(repeat_counts_per_motif[motif]) / len(repeat_counts_per_motif[motif])
			median_value = sorted(repeat_counts_per_motif[motif])[len(repeat_counts_per_motif[motif]) // 2]
			mode_value = calculate_mode(repeat_counts_per_motif[motif])[0]
			print(f"\nMedia de repeticiones para {motif}: {mean_value:.2f}", file=f)
			print(f"Mediana de repeticiones para {motif}: {median_value}", file=f)
			print(f"Moda de repeticiones para {motif}: {mode_value}\n", file=f)

		for motif, counts in repeat_counts_per_motif.items():
			print(f"  {motif} -> {counts}", file=f)
	
	bam.close()
	return repeat_counts_per_motif


def plot_relative_histogram(values, sample):
	plt.figure(figsize=(8, 5))
	bins = range(min(values), max(values) - 10)
	color = "#%06x" % random.randint(0, 0xFFFFFF)
	counts, bins, patches = plt.hist(values, bins=bins, edgecolor='black', align='left', density=True, color=color)
	for p in patches:
		p.set_height(p.get_height() * 100)
	plt.ylim(0, max(p.get_height() for p in patches) * 1.1)
	plt.xlim(min(bins), max(bins))
	plt.title(f'{sample}')
	plt.xlabel('Valor')
	plt.ylabel('Frecuencia relativa (%)')
	plt.grid(axis='y', linestyle='--', alpha=0.7)
	plt.tight_layout()
	plt.savefig(os.path.join("plots", f"{sample}_hist.png"))
	plt.close()


def main():
	os.makedirs("stats", exist_ok=True)
	os.makedirs("plots", exist_ok=True)
	print("Análisis de repeticiones en archivos BAM\n")
	print("Los archivos bam deben encontrarse en la misma carpeta desde la que se lanza el script.\n")

	bam_path = "."
	motif = input("Motivo a buscar (ej: CAG): ").strip().upper() or "CAG"

		
	qscore = input("Umbral de calidad de mapeo [por defecto: 20]: ").strip()
	qscore = int(qscore) if qscore else 20
	minimum = input("Longitud mínima de la repetición [por defecto: 9]: ").strip()
	minimum = int(minimum) if minimum else 9

	bam_files = glob.glob(os.path.join(bam_path, "*.bam"))
	if not bam_files:
		print("No se encontraron archivos BAM en el directorio especificado.")
		return

	for sample in bam_files:
		sample_name = os.path.splitext(os.path.basename(sample))[0]
		repeat_counts = reading_bam(sample, motif, qscore, minimum)
		if repeat_counts[motif]:
			plot_relative_histogram(repeat_counts[motif], sample_name)
		else:
			print(f"No se encontraron repeticiones '{motif}' en {sample_name}")

if __name__ == "__main__":
	main()
