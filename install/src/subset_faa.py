import os
from Bio import SeqIO

# Parámetros
input_dir = "/axolote/diana/algas/Data/Articulo/Proteoma/"
output_dir = "test/faa"
num_sequences = 20                  # Número de secuencias por archivo


# Crear directorio de salida si no existe
os.makedirs(output_dir, exist_ok=True)

# Procesar cada archivo .faa en el directorio
for filename in os.listdir(input_dir):
    if filename.endswith(".faa"):
        input_path = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, filename)

        try:
            # Leer y tomar las primeras N secuencias con codificación latin-1
            with open(input_path, encoding="latin-1") as infile:
                records = list(SeqIO.parse(infile, "fasta"))[:num_sequences]

            # Guardar las secuencias en el nuevo archivo
            with open(output_path, "w") as outfile:
                SeqIO.write(records, outfile, "fasta")

            print(f"{filename}: {len(records)} secuencias escritas.")
        except Exception as e:
            print(f"Error procesando {filename}: {e}")

print("¡Listo! Se generaron los archivos en el directorio:", output_dir)
