# Proyecto_BBDD_Macromoleculas

Trabajo final de la asignatura BBDD de macromoleculas de la Universidad Nacional de Quilmes.

Consiste en una aplicación de consola (terminal de comandos) que permite buscar proteínas usando el identificador de Uniprot o el número de accession de NCBI protein para:

1. Reconocer el tipo de identificador (input) del usuario : Incluye también un input de texto para Uniprot.
2. Realizar la búsqueda en la BBDD correspondiente y descargar de la proteína sus Features (Tabla de Características) utilizando la API REST.
3. Buscar los identificador PDB de la proteína.
4. Permitir al usuario seleccionar dos PDB de la proteína.
5. Realizar un análisis de alineamiento estructural computando el RMSD en una ventana de 5 residuos (Esto puede ser un parámetro libre seteable por el usuario) y devolver como output un gráfico de RMSD vs. posición centroide-residuo de cada ventana.
6. Emplear un visualizador para mostrar las Features y el alineamiento residual de las proteínas.

## Instalación

```bash
# Clonar el repositorio
git clone <https://github.com/Javyne/Proyecto_BBDD_Macromoleculas.git>
cd protein_cli

# Crear entorno virtual (recomendado)
python -m venv venv
venv\Scripts\activate  # En Windows
# source venv/bin/activate  # En Linux/Mac

# Instalar dependencias
pip install -r requirements.txt
```

## Uso

### Comandos Disponibles

```bash
# Ver todos los comandos disponibles
python main.py --help
```

### 1. Buscar proteínas

```bash
# Buscar por Numero de Accession o ID de UniProt o NCBI
python main.py buscar "Accession / ID"

# Buscar por texto descriptivo
python main.py buscar "Texto"
```

### 2. Buscar estructuras PDB

```bash
# Buscar estructuras PDB asociadas a un numero de Accession
python main.py buscar-pdb "Accession"
```

### 3. Descargar features de proteínas

```bash
# Descargar características de una proteína desde UniProt (formato JSON por defecto)
python main.py features P01308

# OPCIONALES
# Descargar en formato específico
python main.py features P01308 --formato json
python main.py features P01308 --formato txt
python main.py features P01308 --formato xml
python main.py features P01308 --formato gff

# Usando la opción corta
python main.py features P01308 -f xml
```

### 4. Análisis RMSD

```bash
# Análisis RMSD local entre dos estructuras PDB
python main.py rmsd-pdb "PDB1" "PDB2"

# OPCIONALES
# Especificar cadena
python main.py rmsd-pdb "PDB1" "PDB2" --cadena X

# Cambiar tamaño de ventana (default: 5)
python main.py rmsd-pdb "PDB1" "PDB2" --ventana N

# Combinar opciones
python main.py rmsd-pdb "PDB1" "PDB2" --cadena X --ventana N
```

## Ejemplos de Uso

```bash
# Buscar Proteina por Accession
python main.py buscar P01308

# Buscar información por Texto
python main.py buscar "albumin human"

# Descargar features
python main.py features P01308 --formato json

# Buscar estructuras PDB
python main.py buscar-pdb P01308

# Analizar RMSD entre dos estructuras encontradas
python main.py rmsd-pdb 1HHO 2HHB 
```

## Funcionalidades Detalladas

### Análisis RMSD Local

El comando `rmsd-local` realiza el cálculo de RMSD entre dos PDB:

1. **Superposición Global**: Alinea las estructuras completas usando átomos CA de aminoácidos estándar
2. **Ventana Deslizante**: Calcula RMSD en ventanas de tamaño configurable (default: 5 residuos)
3. **Fórmula Estándar**: RMSD = √(Σ(coord₁ - coord₂)² / n)
4. **Visualización**: Genera gráficos con estadísticas completas

**Características:**

- Maneja estructuras de diferentes longitudes
- Filtra solo aminoácidos estándar
- Gráficos guardados automáticamente en carpeta `graficos/`
- Soporte para diferentes cadenas y tamaños de ventana

### Búsqueda de PDB con Pandas

El comando `buscar-pdb` utiliza pandas para presentar los resultados en formato tabular.

### Descarga de Features

Soporte para múltiples formatos de salida:

- **JSON**: Estructurado y fácil de procesar
- **TXT**: Formato de texto plano
- **XML**: Formato extensible
- **GFF**: Formato genómico

## Estructura del Proyecto

```bash
protein_cli/
├── main.py                 # Punto de entrada principal
├── requirements.txt        # Dependencias del proyecto
├── README.md              # Este archivo
├── data/                  # Módulo para APIs
│   ├── __init__.py
│   ├── fetch_ncbi.py      # Funciones para NCBI
│   └── fetch_uniprot.py   # Funciones para UniProt
└── utils/                 # Utilidades
    ├── __init__.py
    ├── prote_search.py    # Lógica principal de búsqueda de proteínas
    ├── pdb_search.py      # Lógica para búsqueda de PDB (con pandas)
    ├── features_search.py # Lógica para búsqueda y descarga de features
    └── rmsd_analysis.py   # Análisis RMSD local entre estructuras PDB
```

## Dependencias

```bash
click>=8.0.0          # Interfaz de línea de comandos
requests>=2.25.0      # Peticiones HTTP
biopython>=1.79       # Análisis de estructuras biológicas
pandas>=1.3.0         # Manipulación de datos
matplotlib>=3.5.0     # Generación de gráficos
numpy>=1.21.0         # Cálculos numéricos
seaborn>=0.11.0       # Estilos de gráficos
```
