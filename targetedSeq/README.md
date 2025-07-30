# 🔬 Análisis de repeticiones CAG en archivos BAM

Este script en Python permite analizar repeticiones de motivos específicos (como `CAG`) en archivos BAM y generar gráficos de distribución, además de estadísticas básicas (media, mediana, moda). Está pensado para facilitar su uso incluso sin experiencia en terminal.

---

## Requisitos

### 1. Tener instalado Python (>= 3.8)

### 2. Instalar las dependencias:
Puedes usar `conda` o `pip`.

#### Opción A: conda (recomendado)
```bash
conda create -n repeats python=3.10
conda activate repeats
conda install -c bioconda pysam
pip install matplotlib regex
```




📁 Estructura de carpetas de salida
Al ejecutarlo, el script crea automáticamente:
  - stats/ → contiene archivos .txt con estadísticas por muestra.
  - plots/ → contiene histogramas .png por muestra.


🚀 Cómo usarlo (modo interactivo)
1. Ejecuta el script:
```bash
python analyze_repeats.py
```


2. El script pedirá:

  - Motivo a buscar (por ejemplo CAG)
  - Umbral de calidad (Q score)
  - Longitud mínima de la expansión. Es recomendable emplear aquí, por ejemplo para CAG, 6 o 9. Esto asegurará no cuenta aquellas expansiones que no tengan más de 2/3 CAG seguidos. 

3. Resultado:
Por cada archivo .bam, el script:
  - Calcula las repeticiones del motivo en la región indicada
  - Genera un histograma de frecuencias relativas
  - Guarda estadísticas (media, mediana, moda, total lecturas, etc.) en un .txt

🧪 Ejemplo de salida
Archivo de estadísticas (sample1_stats.txt):

```yaml
📊 Estadísticas para sample1
Total reads: 12345
Q reads: 678 with Q score < 20

Repeticiones encontradas (motivo -> lista de cuentas):
CAG: [11, 12, 13, 15, 12, 13, 13, 12, 13]

Moda de repeticiones para CAG: [13]
Media de repeticiones para CAG: 12.67
Mediana de repeticiones para CAG: 13
```
Archivo de gráfico (plots/sample1_hist.png):
→ Histograma de frecuencias relativas de tamaños de repetición.


Los colores o formato de las figuras


