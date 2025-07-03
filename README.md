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

### Requisitos Previos

- Python 3.8 o superior
- Conexión a internet (para descargar estructuras PDB y acceder a APIs)
- Navegador web (para visualización 3D)

### Instalación Paso a Paso

```bash
# 1. Clonar el repositorio
git clone <https://github.com/Javyne/Proyecto_BBDD_Macromoleculas.git>
cd protein_cli

# 2. Crear entorno virtual (recomendado)
python -m venv venv

# 3. Activar entorno virtual
# En Windows:
venv\Scripts\activate
# En Linux/Mac:
# source venv/bin/activate

# 4. Instalar dependencias
pip install -r requirements.txt

# 5. Verificar instalación
python main.py --help
```

### Instalación Rápida (Sin Entorno Virtual)

```bash
# Solo para desarrollo/testing
pip install -r requirements.txt
python main.py --help
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
# Especificar cadena común para ambos PDBs
python main.py rmsd-pdb "PDB1" "PDB2" --cadena1 X

# Especificar cadenas independientes para cada PDB
python main.py rmsd-pdb "PDB1" "PDB2" --cadena1 X --cadena2 Y

# Cambiar tamaño de ventana (default: 5)
python main.py rmsd-pdb "PDB1" "PDB2" --ventana N

# Combinar opciones
python main.py rmsd-pdb "PDB1" "PDB2" --cadena1 X --ventana N
python main.py rmsd-pdb "PDB1" "PDB2" --cadena1 X --cadena2 Y --ventana N

# NOTA: No se pueden usar --cadena y --cadenas simultáneamente
```

### 5. Visualización de estructura terciaria de proteínas

```bash
# Visualización de una proteína
python main.py mostrar-PDB-simple "PDB"

# Mostrar Features
python main.py mostrar-PDB-features "PDB" 

# Visualización de alineamiento estructural
python main.py mostrar-alineamiento "PDB1" "PDB2"

# OPCIONALES

# Mostrar solo una cadena
python main.py mostrar-PDB-simple "PDB" --cadena X

# Mostrar Features y asignar un feature personalizado
python main.py mostrar-PDB-features "PDB" --feature "Color HEX"(str) inicio(int) fin(int) "NombreFeature"(str)

# Especificar cadena en alineamiento
python main.py mostrar-alineamiento "PDB1" "PDB2" --cadena X

# Especificar tamaño de ventana (default: 50)
python main.py mostrar-alineamiento "PDB1" "PDB2" --ventana N

# Cambiar colores (default: blue red)
python main.py mostrar-alineamiento "PDB1" "PDB2" --colores "col1" "col2"

# Combinar opciones en alineamiento
python main.py mostrar-alineamiento "PDB1" "PDB2" --cadena X --colores "col1" "col2" --ventana 25
```

## Ejemplos de Uso

### Búsqueda y Descarga de Datos

```bash
# Buscar Proteina por Accession
python main.py buscar P01308

# Buscar información por Texto
python main.py buscar "albumin human"

# Descargar features
python main.py features P01308 --formato json

# Buscar estructuras PDB
python main.py buscar-pdb P01308
```

### Análisis RMSD

#### 1. Misma Proteína, Diferentes Estados Conformacionales

```bash
# Hemoglobina en diferentes estados (misma molécula, diferentes conformaciones)
python main.py rmsd-pdb 1HHO 2HHB --cadena1 A
# Resultado: RMSD bajo (< 2 Å) - estructuras muy similares
# UniProt IDs compartidos: P01922 (Hemoglobina alfa humana)
```

#### 2. Misma Proteína, Diferentes Cadenas

```bash
# Hemoglobina - comparar cadenas A y B de la misma estructura
python main.py rmsd-pdb 1HHO 1HHO --cadena1 A --cadena2 B
# Resultado: ADVERTENCIA + RMSD bajo entre cadenas homólogas
```

#### 3. Proteínas Diferentes

```bash
# Hemoglobina vs Albúmina (moléculas completamente diferentes)
python main.py rmsd-pdb 1A00 1AO6
# Resultado: ADVERTENCIA - proteínas diferentes
# UniProt IDs: P01966 (Hemoglobina) vs P02768 (Albúmina)
# El análisis continúa solo si el usuario confirma
```

#### 4. Proteínas Relacionadas - Familia de Proteínas

```bash
# Diferentes isoformas de hemoglobina
python main.py rmsd-pdb 1HHO 1A00 --cadena1 A
# Resultado: RMSD moderado - proteínas relacionadas pero no idénticas
```

#### 5. Análisis con Diferentes Tamaños de Ventana

```bash
# Ventana pequeña (3 residuos) - más sensible a cambios locales
python main.py rmsd-pdb 1HHO 2HHB --ventana 3

# Ventana grande (10 residuos) - más suavizado
python main.py rmsd-pdb 1HHO 2HHB --ventana 10
```

#### 6. Cadenas Independientes

```bash
# Comparar cadenas A y B de la misma estructura
python main.py rmsd-pdb 1HHO 1HHO --cadena1 A --cadena2 B

# Comparar cadenas específicas de diferentes estructuras
python main.py rmsd-pdb 1HHO 2HHB --cadena1 A --cadena2 B

# Comparar cadenas diferentes de proteínas relacionadas
python main.py rmsd-pdb 1HHO 1A00 --cadena1 A --cadena2 B

```

### Interpretación de Resultados RMSD

#### RMSD Bajo (< 2 Å)

- Estructuras muy similares
- Misma proteína en diferentes condiciones
- Diferentes conformaciones de la misma molécula

#### RMSD Moderado (2-5 Å)

- Proteínas relacionadas (familia de proteínas)
- Diferentes isoformas
- Mutaciones puntuales

#### RMSD Alto (> 5 Å)

- Proteínas diferentes
- Estructuras no relacionadas
- Comparación biológicamente sin sentido

### Visualización de Estructuras

```bash
# Visualizar una estructura PDB (simple sin features):
python main.py mostrar_PDB_simple 1HHO

# Visualizar una estructura PDB (Con features):
python main.py mostrar_PDB_features 1HHO

# Visualizar alineamiento entre dos estructuras
python main.py mostrar-alineamiento 1HHO 2HHB --cadena B --colores blue green --ventana 10
```

## Flujo de Trabajo Típico

### **Análisis Completo de una Proteína**

```bash
# 1. Buscar información de la proteína
python main.py buscar P01308

# 2. Descargar características detalladas
python main.py features P01308 --formato gff

# 3. Encontrar estructuras PDB disponibles
python main.py buscar-pdb P01308

# 4. Analizar RMSD entre dos estructuras seleccionadas
python main.py rmsd-pdb 1HHO 2HHB --cadena1 A --ventana 5

# 5. Visualizar las estructuras
python main.py mostrar-alineamiento 1HHO 2HHB --cadena1 A
```

## Funcionalidades Detalladas

### Análisis RMSD Local

El comando `rmsd-pdb` realiza el cálculo de RMSD entre dos PDB:

1. **Verificación UniProt**: Verifica automáticamente si las estructuras pertenecen a la misma proteína
2. **Superposición Global**: Alinea las estructuras completas usando átomos CA de aminoácidos estándar
3. **Ventana Deslizante**: Calcula RMSD en ventanas de tamaño configurable (default: 5 residuos)
4. **Fórmula Estándar**: RMSD = √(Σ(coord₁ - coord₂)² / n)
5. **Gráficos RMSD**: Genera gráficos con estadísticas completas
6. **Visualización de estructura**: Genera una instancia y un archivo HTML/CSS para mostrar la estructura 3D de las proteínas interactiva.

**Características:**

- **Verificación automática de compatibilidad**: Detecta si las estructuras pertenecen a la misma proteína
- **Manejo de estructuras de diferentes longitudes**: Corta automáticamente a la longitud mínima
- **Filtra solo aminoácidos estándar**: Excluye residuos no estándar del análisis
- **Gráficos guardados automáticamente**: En carpeta `graficos/` con nombres descriptivos
- **Soporte para diferentes cadenas y tamaños de ventana**: Flexibilidad total en el análisis
- **Advertencias inteligentes**: Informa al usuario sobre comparaciones biológicamente cuestionables

### Búsqueda de PDB con Pandas

El comando `buscar-pdb` utiliza pandas para presentar los resultados en formato tabular.

### Descarga de Features

Soporte para múltiples formatos de salida:

- **JSON**: Estructurado y fácil de procesar
- **TXT**: Formato de texto plano
- **XML**: Formato extensible
- **GFF**: Formato genómico

### Búsqueda de Proteínas

El comando `buscar` permite buscar proteínas por diferentes criterios:

- **Accession/ID**: Búsqueda directa por identificador UniProt o NCBI
- **Texto descriptivo**: Búsqueda por palabras clave o descripción
- **Reconocimiento automático**: Detecta el tipo de identificador automáticamente

### Visualización de Estructuras PDB

#### Visualización Simple

- Comando: `mostrar-PDB-simple`
- Funcionalidad: Visualización 3D básica de estructuras PDB

#### Visualización con Features

- Comando: `mostrar-PDB-features`
- Funcionalidad: Visualización con características biológicas destacadas

#### Alineamiento Estructural

- Comando: `mostrar-alineamiento`
- Funcionalidad: Comparación visual de dos estructuras PDB

## Estructura del Proyecto

```bash
protein_cli/
├── main.py                 # Punto de entrada principal
├── requirements.txt        # Dependencias del proyecto
├── README.md              # Este archivo
├── data/                  # Módulo para APIs
│   ├── __init__.py
│   ├── fetch_ncbi.py      # Funciones para NCBI
│   └── fetch_pdb.py       # Funciones para PDB
│   └── fetch_uniprot.py   # Funciones para UniProt
└── utils/                 # Utilidades
    ├── __init__.py
    ├── prote_search.py    # Lógica principal de búsqueda de proteínas
    ├── pdb_search.py      # Lógica para búsqueda de PDB (con pandas)
    ├── features_search.py # Lógica para búsqueda y descarga de features
    ├── pdb_viewer.py      # Visualización de archivos PDB
    └── rmsd_analysis.py   # Análisis RMSD local con verificación UniProt
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
py3Dmol>=2.5.1        # Manipulación y visualización 3D de archivos PDB.
```

## Licencia

Este proyecto es parte del trabajo final de la asignatura BBDD de macromoléculas de la Universidad Nacional de Quilmes.

### **Uso Académico**

- Libre para uso educativo y de investigación
- Citar el proyecto en publicaciones académicas
- Respetar las licencias de las dependencias utilizadas

### **Dependencias Externas**

- **RCSB PDB**: Datos de estructuras proteicas
- **UniProt**: Información de proteínas
- **NCBI**: Datos de secuencias

## Contacto

- **Universidad**: Universidad Nacional de Quilmes
- **Asignatura**: BBDD de Macromoléculas
- **Año**: 2025

---

**Nota**: Este proyecto está diseñado para análisis académico y de investigación. Para uso comercial, consultar las licencias de las bases de datos externas utilizadas.
