import os
import tempfile
import warnings
import io

import matplotlib.pyplot as plt
import numpy as np
import requests
import seaborn as sns
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB import PDBIO, Select

warnings.filterwarnings("ignore")

# Configurar estilo de gráficos
plt.style.use("seaborn-v0_8")
sns.set_palette("husl")


# Descarga un archivo PDB desde la base de datos PDB
# Entrada = ID de PDB (string)
# Salida = ruta del archivo temporal descargado
def descargar_pdb(pdb_id):
    pdb_id = pdb_id.upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    try:
        response = requests.get(url, timeout=30)

        # Manejar específicamente el error 404
        if response.status_code == 404:
            raise Exception(
                f"El código PDB '{pdb_id}' no existe en la base de datos RCSB. "
                f"Verifica que el código sea correcto o consulta https://www.rcsb.org/"
            )

        response.raise_for_status()

        # Verificar que el contenido no esté vacío
        if not response.text.strip():
            raise Exception(
                f"El archivo PDB '{pdb_id}' está vacío o no contiene datos válidos"
            )

        # Crear archivo temporal
        temp_file = tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False)
        temp_file.write(response.text)
        temp_file.close()

        return temp_file.name

    except requests.exceptions.Timeout:
        raise Exception(
            f"Timeout al descargar PDB {pdb_id}. Verifica tu conexión a internet."
        )
    except requests.exceptions.ConnectionError:
        raise Exception(
            f"Error de conexión al descargar PDB {pdb_id}. Verifica tu conexión a internet."
        )
    except requests.RequestException as e:
        raise Exception(f"Error al descargar PDB {pdb_id}: {e}")
    except Exception as e:
        raise Exception(f"Error inesperado al descargar PDB {pdb_id}: {e}")


# Carga una estructura PDB usando BioPython
# Entrada = ruta del archivo PDB
# Salida = objeto estructura de BioPython
def cargar_estructura(archivo_pdb):
    parser = PDBParser(QUIET=True)
    try:
        estructura = parser.get_structure("protein", archivo_pdb)
        return estructura
    except Exception as e:
        raise Exception(f"Error al cargar estructura: {e}")


# Encuentra las cadenas que están presentes en ambas estructuras
# Entrada = dos objetos estructura de BioPython
# Salida = lista de IDs de cadenas comunes
def obtener_cadenas_comunes(estructura1, estructura2):
    cadenas1 = set([chain.id for chain in estructura1.get_chains()])
    cadenas2 = set([chain.id for chain in estructura2.get_chains()])
    cadenas_comunes = cadenas1.intersection(cadenas2)

    if not cadenas_comunes:
        raise Exception("No se encontraron cadenas comunes entre las estructuras")

    return list(cadenas_comunes)


# Calcula RMSD local usando una ventana
# Entrada = dos estructuras, ID de cadena, tamaño de ventana (5 por default)
# Salida = listas de posiciones y valores RMSD locales
def calcular_rmsd_local(estructura1, estructura2, cadena_id, ventana=5):
    """
    Algoritmo científico estándar para RMSD local:
    1. Superponer globalmente las estructuras completas
    2. Calcular RMSD local en ventanas sin superponer nuevamente
    """
    # Obtener cadenas
    cadena1 = estructura1[0][cadena_id]
    cadena2 = estructura2[0][cadena_id]

    # Extraer coordenadas CA de aminoácidos estándar
    coords1 = []
    coords2 = []
    residuos1 = []
    residuos2 = []

    for residuo in cadena1:
        if is_aa(residuo, standard=True) and "CA" in residuo:
            coords1.append(residuo["CA"].get_coord())
            residuos1.append(residuo)

    for residuo in cadena2:
        if is_aa(residuo, standard=True) and "CA" in residuo:
            coords2.append(residuo["CA"].get_coord())
            residuos2.append(residuo)

    # Convertir a arrays numpy
    coords1 = np.array(coords1)
    coords2 = np.array(coords2)

    # Verificar que tenemos suficientes residuos
    if len(coords1) < ventana or len(coords2) < ventana:
        raise Exception(f"Se necesitan al menos {ventana} residuos para el análisis")

    # Cortar para que tengan igual longitud (algoritmo estándar)
    min_len = min(len(coords1), len(coords2))
    coords1 = coords1[:min_len]
    coords2 = coords2[:min_len]
    residuos1 = residuos1[:min_len]
    residuos2 = residuos2[:min_len]

    # PASO 1: Superposición global (estándar científico)
    print("Realizando superposición global...")

    # Crear átomos para la superposición global
    atomos_global1 = []
    atomos_global2 = []

    for residuo in residuos1:
        atomos_global1.append(residuo["CA"])

    for residuo in residuos2:
        atomos_global2.append(residuo["CA"])

    superimposer_global = Superimposer()
    superimposer_global.set_atoms(atomos_global1, atomos_global2)
    superimposer_global.apply(atomos_global2)

    # Actualizar coordenadas después de la superposición
    coords2 = np.array([atomo.get_coord() for atomo in atomos_global2])

    # PASO 2: Calcular RMSD local en ventanas (sin superponer nuevamente)
    rmsd_local = []
    posiciones = []

    for i in range(min_len - ventana + 1):
        # Tomar ventana de coordenadas (ya superpuestas globalmente)
        ventana_coords1 = coords1[i : i + ventana]
        ventana_coords2 = coords2[i : i + ventana]

        # Calcular RMSD directamente (fórmula estándar)
        # RMSD = sqrt(sum((coord1 - coord2)^2) / n)
        diff_squared = np.sum((ventana_coords1 - ventana_coords2) ** 2, axis=1)
        rmsd = np.sqrt(np.mean(diff_squared))

        rmsd_local.append(rmsd)

        # Posición central de la ventana (usar índice de residuo)
        pos_central = residuos1[i + ventana // 2].id[1]
        posiciones.append(pos_central)

    return posiciones, rmsd_local


# Genera un gráfico de RMSD local
# Entrada = posiciones, valores RMSD, IDs de PDB, cadena, ventana
# Salida = objeto figura de matplotlib
def generar_grafico_rmsd(
    posiciones, rmsd_values, pdb1_id, pdb2_id, cadena_id, ventana=5
):

    fig, ax = plt.subplots(figsize=(12, 6))

    # Crear gráfico
    ax.plot(posiciones, rmsd_values, "o-", linewidth=2, markersize=4, alpha=0.7)

    # Configurar gráfico
    ax.set_xlabel("Posición del residuo", fontsize=12)
    ax.set_ylabel("RMSD Local (Å)", fontsize=12)
    ax.set_title(
        f"RMSD Local entre {pdb1_id} y {pdb2_id} (Cadena {cadena_id}, Ventana={ventana})",
        fontsize=14,
        fontweight="bold",
    )

    # Agregar grid
    ax.grid(True, alpha=0.3)

    # Agregar estadísticas
    rmsd_mean = np.mean(rmsd_values)
    rmsd_std = np.std(rmsd_values)
    rmsd_max = np.max(rmsd_values)

    stats_text = f"Promedio: {rmsd_mean:.3f} Å\nDesv. Est.: {rmsd_std:.3f} Å\nMáximo: {rmsd_max:.3f} Å"
    ax.text(
        0.02,
        0.98,
        stats_text,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    # Ajustar layout
    plt.tight_layout()

    return fig


# Función principal para analizar RMSD local entre dos estructuras PDB
# Entrada = IDs de PDB, cadena opcional, tamaño de ventana
# Salida = ruta del archivo guardado, posiciones, valores RMSD
def analizar_rmsd_local(pdb1_id, pdb2_id, cadena_id=None, ventana=5):

    print(f"Analizando RMSD local entre {pdb1_id} y {pdb2_id}...")

    try:
        # Descargar estructuras en archivos temporales
        print("Descargando estructuras PDB...")
        archivo1 = descargar_pdb(pdb1_id)
        archivo2 = descargar_pdb(pdb2_id)

        # Cargar estructuras utilizando PDBParser de BioPython
        print("Cargando estructuras...")
        estructura1 = cargar_estructura(archivo1)
        estructura2 = cargar_estructura(archivo2)

        # Encontrar cadenas comunes
        cadenas_comunes = obtener_cadenas_comunes(estructura1, estructura2)

        if cadena_id is None:
            cadena_id = cadenas_comunes[0]
            print(f"Usando cadena {cadena_id} (primera disponible)")
        elif cadena_id not in cadenas_comunes:
            raise Exception(
                f"Cadena {cadena_id} no encontrada. Cadenas disponibles: {cadenas_comunes}"
            )

        # Calcular RMSD local
        print(f"Calculando RMSD local para cadena {cadena_id}...")
        posiciones, rmsd_values = calcular_rmsd_local(
            estructura1, estructura2, cadena_id, ventana
        )

        # Generar gráfico
        print("Generando gráfico...")
        fig = generar_grafico_rmsd(
            posiciones, rmsd_values, pdb1_id, pdb2_id, cadena_id, ventana
        )

        # Guardar gráfico
        nombre_archivo = f"rmsd_local_{pdb1_id}_{pdb2_id}_{cadena_id}.png"
        carpeta='graficos'; os.makedirs(carpeta,exist_ok=True)
        ruta_completa = os.path.join("graficos", nombre_archivo)
        fig.savefig(ruta_completa, dpi=300, bbox_inches="tight")
        print(f"Gráfico guardado como: {ruta_completa}")

        # Mostrar estadísticas
        print(f"\nEstadísticas RMSD Local:")
        print(f"Promedio: {np.mean(rmsd_values):.3f} Å")
        print(f"Desviación estándar: {np.std(rmsd_values):.3f} Å")
        print(f"Máximo: {np.max(rmsd_values):.3f} Å")
        print(f"Mínimo: {np.min(rmsd_values):.3f} Å")

        # Limpiar archivos temporales
        os.unlink(archivo1)
        os.unlink(archivo2)

        return ruta_completa, posiciones, rmsd_values

    except Exception as e:
        print(f"Error: {e}")
        return None, None, None


# Para alineamiento estructural gráfico => Con alineamiento global.
# Extraer los C-alfa de la cadena X
def conseguir_atomos_CA(estructura, cadenaID=None):
    """Devuelve un diccionario {resid: atom} para los C-alfa de una cadena."""
    modelo = estructura[0]
    if cadenaID is None:
        cadenaID = 'A'
    cadena = modelo[cadenaID]
    return {residuo.id: residuo['CA'] for residuo in cadena if 'CA' in residuo}

# Alinear estructuras usando Superimposer
def alinear_estructuras(estructuraReferencia, estructuraOtra, cadenaID, tolerancia=3):
    ca_ref = conseguir_atomos_CA(estructuraReferencia, cadenaID)
    ca_otro = conseguir_atomos_CA(estructuraOtra, cadenaID)

    # Encontrar residuos comunes -> Para alineamiento con proteinas de diferente tamaño
    residuos_comunes = set(ca_ref.keys()) & set(ca_otro.keys())
    
    if len(residuos_comunes) < tolerancia:
        raise ValueError("Muy pocos residuos comunes para alinear.")
    
    atomos_referencia = [ca_ref[residuo] for residuo in sorted(residuos_comunes)]
    atomos_comunes = [ca_otro[residuo] for residuo in sorted(residuos_comunes)]
    
    if len(atomos_referencia) != len(atomos_comunes):
        raise ValueError("Las listas de C-Alfa no tienen el mismo largo.")
    
    si = Superimposer()
    si.set_atoms(atomos_referencia, atomos_comunes)
    si.apply([atom for atom in estructuraOtra.get_atoms()]) # Aplica la transformación sobre todos los átomos
    return si.rms

#Convertir estructura Biopython a string PDB para py3Dmol
class SeleccionarCadena(Select):
    def __init__(self, cadenaID):
        self.cadenaID = cadenaID
    def accept_chain(self, cadena):
        return cadena.id == self.cadenaID

def estructura_PDB_a_str(estructura, cadenaID=None):
    io_pdb = PDBIO()
    io_pdb.set_structure(estructura)
    string_io = io.StringIO()
    if cadenaID is None:
        cadenaID = 'A'
    io_pdb.save(string_io, select=SeleccionarCadena(cadenaID))
    return string_io.getvalue()