import io
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from Bio.PDB.PDBIO import PDBIO, Select
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.Superimposer import Superimposer

from data.fetch_pdb import descargar_pdb
from data.fetch_uniprot import buscar_pdb_accessions

warnings.filterwarnings("ignore")

# Configurar estilo de gráficos
plt.style.use("seaborn-v0_8")
sns.set_palette("husl")


# =============================================================================
# RESUMEN DE FUNCIONES DEL MÓDULO RMSD_ANALYSIS
# =============================================================================
#
# FUNCIONES PRINCIPALES:
# ---------------------
# analizar_rmsd_local() - Función principal que coordina todo el análisis RMSD
#
# FUNCIONES DE VERIFICACIÓN Y VALIDACIÓN:
# ---------------------------------------
# verificar_compatibilidad_uniprot() - Verifica si las estructuras tienen IDs UniProt compatibles
# determinar_cadena_analisis() - Determina qué cadena usar para el análisis
#
# FUNCIONES DE MANEJO DE DATOS:
# -----------------------------
# cargar_estructura() - Carga una estructura PDB usando BioPython
# cargar_estructuras_pdb() - Descarga y carga dos estructuras PDB
# obtener_cadenas_comunes() - Encuentra cadenas comunes entre dos estructuras
# extraer_coordenadas_ca() - Extrae coordenadas CA de aminoácidos estándar
# preparar_coordenadas_para_analisis() - Prepara coordenadas para el análisis RMSD
#
# FUNCIONES DE CÁLCULO RMSD:
# --------------------------
# superponer_estructuras_globalmente() - Realiza superposición global de estructuras
# calcular_rmsd_ventana() - Calcula RMSD para una ventana específica
# calcular_rmsd_local() - Algoritmo principal para RMSD local con ventanas
#
# FUNCIONES DE VISUALIZACIÓN:
# ---------------------------
# generar_grafico_rmsd() - Genera gráfico de RMSD local
# generar_y_guardar_grafico() - Genera y guarda el gráfico en archivo
# mostrar_estadisticas_rmsd() - Muestra estadísticas descriptivas del análisis
#
# FUNCIONES AUXILIARES:
# ---------------------
# conseguir_atomos_CA() - Extrae átomos CA de una cadena
# alinear_estructuras() - Alinea estructuras usando Superimposer
# SeleccionarCadena() - Clase para seleccionar cadenas específicas
# estructura_PDB_a_str() - Convierte estructura BioPython a string PDB
#
# =============================================================================


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
def obtener_cadenas_comunes(estructura1, estructura2, cadena1_id=None, cadena2_id=None):

    if cadena1_id is None:
        cadenas1 = set([chain.id for chain in estructura1.get_chains()])
    else:
        cadenas1 = set([cadena1_id])

    if cadena2_id is None:
        cadenas2 = set([chain.id for chain in estructura2.get_chains()])
    else:
        cadenas2 = set([cadena2_id])

    cadenas_comunes = cadenas1.intersection(cadenas2)

    if not cadenas_comunes:
        raise Exception("No se encontraron cadenas comunes entre las estructuras")

    return list(cadenas_comunes)


# Calcula RMSD local usando una ventana
# Entrada = Estructura, ID de cadena, tamaño de ventana (5 por default)
# Salida = listas de posiciones y valores RMSD locales
def extraer_coordenadas_ca(estructura1, cadena_id):

    # Extrae las coordenadas de los átomos C-alfa (CA) de ambas estructuras
    try:
        cadena1 = estructura1[0][cadena_id]
    except Exception as e:
        print(
            f"No existe cadena {cadena_id} en los PDBs. Se utilizará por defecto el valor de cadena A."
        )
        cadena_id = "A"
        cadena1 = estructura1[0][cadena_id]

    coords1, residuos1 = [], []

    # Extraer coordenadas CA de la estructura
    for residuo in cadena1:
        if is_aa(residuo, standard=True) and "CA" in residuo:
            coords1.append(residuo["CA"].get_coord())
            residuos1.append(residuo)

    return np.array(coords1), residuos1


# Prepara las coordenadas para el análisis RMSD
# Entrada = arrays de coordenadas, listas de residuos, tamaño de ventana
# Salida = coordenadas preparadas y longitud mínima
def preparar_coordenadas_para_analisis(coords1, coords2, residuos1, residuos2, ventana):

    # Verificar que tenemos suficientes residuos para el análisis
    if len(coords1) < ventana or len(coords2) < ventana:
        raise Exception(f"Se necesitan al menos {ventana} residuos para el análisis")

    # Cortar para que tengan igual longitud (algoritmo estándar)
    min_len = min(len(coords1), len(coords2))
    coords1 = coords1[:min_len]
    coords2 = coords2[:min_len]
    residuos1 = residuos1[:min_len]
    residuos2 = residuos2[:min_len]

    return coords1, coords2, residuos1, residuos2, min_len


# Realiza la superposición global de las estructuras usando los átomos CA
# Entrada = listas de residuos de ambas estructuras
# Salida = coordenadas actualizadas de la segunda estructura
def superponer_estructuras_globalmente(residuos1, residuos2):

    print("Realizando superposición global...")

    # Crear listas de átomos CA para la superposición
    atomos_global1 = [residuo["CA"] for residuo in residuos1]
    atomos_global2 = [residuo["CA"] for residuo in residuos2]

    # Aplicar superposición usando BioPython
    superimposer_global = Superimposer()
    superimposer_global.set_atoms(atomos_global1, atomos_global2)
    superimposer_global.apply(atomos_global2)

    # Retornar las coordenadas actualizadas de la segunda estructura
    return np.array([atomo.get_coord() for atomo in atomos_global2])


# Calcula el RMSD para una ventana específica de residuos
# Entrada = coordenadas de dos ventanas de residuos
# Salida = valor RMSD calculado
def calcular_rmsd_ventana(ventana_coords1, ventana_coords2):

    diff_squared = np.sum((ventana_coords1 - ventana_coords2) ** 2, axis=1)
    return np.sqrt(np.mean(diff_squared))


# Algoritmo científico estándar para RMSD local
# Entrada = dos estructuras, ID de cadena, tamaño de ventana (5 por default)
# Salida = listas de posiciones y valores RMSD locales
def calcular_rmsd_local(estructura1, estructura2, cadena1_id, cadena2_id, ventana=5):

    # Extraer coordenadas CA de ambas estructuras
    coords1, residuos1 = extraer_coordenadas_ca(estructura1, cadena1_id)
    coords2, residuos2 = extraer_coordenadas_ca(estructura2, cadena2_id)

    # Preparar coordenadas para el análisis
    coords1, coords2, residuos1, residuos2, min_len = (
        preparar_coordenadas_para_analisis(
            coords1, coords2, residuos1, residuos2, ventana
        )
    )

    # PASO 1: Superposición global (estándar científico)
    coords2 = superponer_estructuras_globalmente(residuos1, residuos2)

    # PASO 2: Calcular RMSD local en ventanas (sin superponer nuevamente)
    rmsd_local, posiciones = [], []

    for i in range(min_len - ventana + 1):
        # Tomar ventana de coordenadas (ya superpuestas globalmente)
        ventana_coords1 = coords1[i : i + ventana]
        ventana_coords2 = coords2[i : i + ventana]

        # Calcular RMSD para esta ventana
        rmsd = calcular_rmsd_ventana(ventana_coords1, ventana_coords2)
        rmsd_local.append(rmsd)

        # Posición central de la ventana (usar índice de residuo)
        pos_central = residuos1[i + ventana // 2].id[1]
        posiciones.append(pos_central)

    return posiciones, rmsd_local


# Genera un gráfico de RMSD local
# Entrada = posiciones, valores RMSD, IDs de PDB, cadena, ventana
# Salida = objeto figura de matplotlib
def generar_grafico_rmsd(
    posiciones, rmsd_values, pdb1_id, pdb2_id, cadena1_id, cadena2_id, ventana=5
):

    fig, ax = plt.subplots(figsize=(12, 6))

    # Crear gráfico
    ax.plot(posiciones, rmsd_values, "o-", linewidth=2, markersize=4, alpha=0.7)

    # Configurar gráfico
    ax.set_xlabel("Posición del residuo", fontsize=12)
    ax.set_ylabel("RMSD Local (Å)", fontsize=12)
    ax.set_title(
        f"RMSD Local entre {pdb1_id} y {pdb2_id} (Cadena {cadena1_id} y {cadena2_id}, Ventana={ventana})",
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
def verificar_compatibilidad_uniprot(pdb1_id, pdb2_id, cadena1_id, cadena2_id):

    # Verifica si las estructuras PDB tienen IDs de UniProt compatibles
    uni1 = buscar_pdb_accessions(pdb1_id, cadena1_id)
    uni2 = buscar_pdb_accessions(pdb2_id, cadena2_id)

    # Caso 1: Al menos una cadena no tiene mapeo UniProt
    if not uni1 or not uni2:
        print(
            "⚠️  Advertencia: al menos una de las cadenas no tiene mapeo UniProt "
            "en PDBe-SIFTS; no se puede verificar si son la misma molécula."
        )
        return True

    # Caso 2: Las cadenas no comparten ningún accession UniProt
    if uni1.isdisjoint(uni2):
        print(
            "⚠️  Advertencia: las cadenas seleccionadas no comparten ningún "
            "accession UniProt ({}  vs  {}).  Es probable que pertenezcan a "
            "moléculas distintas, de modo que el RMSD resultante puede carecer "
            "de sentido biológico.".format(
                ", ".join(sorted(uni1)), ", ".join(sorted(uni2))
            )
        )
        respuesta = (
            input("¿Desea continuar con el análisis RMSD? (s/n): ").lower().strip()
        )
        return respuesta in ["s", "si", "sí", "y", "yes"]

    # Caso 3: Las cadenas comparten al menos un accession UniProt
    return True


# Descarga y carga las estructuras PDB desde la base de datos
# Entrada = dos IDs de PDB
# Salida = archivos temporales y estructuras cargadas
def cargar_estructuras_pdb(pdb1_id, pdb2_id):

    print("Descargando estructuras PDB...")
    archivo1 = descargar_pdb(pdb1_id)
    archivo2 = descargar_pdb(pdb2_id)

    print("Cargando estructuras...")
    estructura1 = cargar_estructura(archivo1)
    estructura2 = cargar_estructura(archivo2)

    return archivo1, archivo2, estructura1, estructura2


# Genera el gráfico de RMSD local y lo guarda en la carpeta 'graficos'
# Entrada = posiciones, valores RMSD, IDs de PDB, cadena, ventana
# Salida = ruta completa del archivo guardado
def generar_y_guardar_grafico(
    posiciones, rmsd_values, pdb1_id, pdb2_id, cadena1_id, cadena2_id, ventana
):

    print("Generando gráfico...")
    fig = generar_grafico_rmsd(
        posiciones, rmsd_values, pdb1_id, pdb2_id, cadena1_id, cadena2_id, ventana
    )

    # Crear carpeta si no existe y guardar gráfico
    nombre_archivo = f"rmsd_local_{pdb1_id}_{pdb2_id}_{cadena1_id}_{cadena2_id}.png"
    carpeta = "graficos"
    os.makedirs(carpeta, exist_ok=True)
    ruta_completa = os.path.join("graficos", nombre_archivo)
    fig.savefig(ruta_completa, dpi=300, bbox_inches="tight")
    print(f"Gráfico guardado como: {ruta_completa}")

    return ruta_completa


# Muestra las estadísticas descriptivas del análisis RMSD
# Entrada = lista de valores RMSD
# Salida = impresión de estadísticas en consola
def mostrar_estadisticas_rmsd(rmsd_values):

    print(f"\nEstadísticas RMSD Local:")
    print(f"Promedio: {np.mean(rmsd_values):.3f} Å")
    print(f"Desviación estándar: {np.std(rmsd_values):.3f} Å")
    print(f"Máximo: {np.max(rmsd_values):.3f} Å")
    print(f"Mínimo: {np.min(rmsd_values):.3f} Å")


# Función principal para analizar RMSD local entre dos estructuras PDB
# Entrada = IDs de PDB, cadena opcional, tamaño de ventana
# Salida = ruta del archivo guardado, posiciones, valores RMSD
def analizar_rmsd_local(pdb1_id, pdb2_id, cadena1_id=None, cadena2_id=None, ventana=5):

    print(f"Analizando RMSD local entre {pdb1_id} y {pdb2_id}...")

    try:
        # PASO 1: Descargar y cargar estructuras PDB
        archivo1, archivo2, estructura1, estructura2 = cargar_estructuras_pdb(
            pdb1_id, pdb2_id
        )

        # PASO 2: si cadena1_id o cadena2_id es None, obtener cadenas en común
        if cadena1_id is None or cadena2_id is None:
            cadenas_comunes = obtener_cadenas_comunes(
                estructura1, estructura2, cadena1_id, cadena2_id
            )
            if cadena1_id is None:
                cadena1_id = cadenas_comunes[0]
            if cadena2_id is None:
                cadena2_id = cadenas_comunes[0]

        # PASO 3: Verificar compatibilidad UniProt antes de descargar estructuras
        print("Verificando anotaciones UniProt...")

        if not verificar_compatibilidad_uniprot(
            pdb1_id, pdb2_id, cadena1_id, cadena2_id
        ):
            print("Análisis cancelado por el usuario.")
            return None, None, None

        # PASO 5: Calcular RMSD local
        print(f"Calculando RMSD local...")
        posiciones, rmsd_values = calcular_rmsd_local(
            estructura1, estructura2, cadena1_id, cadena2_id, ventana
        )

        # PASO 6: Generar gráfico y mostrar estadísticas
        ruta_completa = generar_y_guardar_grafico(
            posiciones, rmsd_values, pdb1_id, pdb2_id, cadena1_id, cadena2_id, ventana
        )
        mostrar_estadisticas_rmsd(rmsd_values)

        # PASO 7: Limpiar archivos temporales
        os.unlink(archivo1)
        os.unlink(archivo2)

        return ruta_completa, posiciones, rmsd_values

    except Exception as e:
        print(f"Error: {e}")
        return None, None, None


# Para alineamiento estructural gráfico => Con alineamiento global.
# Extraer los C-alfa de la cadena X
def conseguir_atomos_CA(estructura, cadenaID=None):
    modelo = estructura[0]
    cadena = modelo[cadenaID]

    return {residuo.id: residuo["CA"] for residuo in cadena if "CA" in residuo}


# Alinear estructuras usando Superimposer
def alinear_estructuras(estructuraReferencia, estructuraOtra, cadenaID, tolerancia=3):
    ca_ref = conseguir_atomos_CA(estructuraReferencia, cadenaID)
    ca_otro = conseguir_atomos_CA(estructuraOtra, cadenaID)

    # Encontrar residuos comunes -> Para alineamiento con proteinas de diferente tamaño
    residuos_comunes = set(ca_ref.keys()) & set(ca_otro.keys())
    print(cadenaID, residuos_comunes)
    if len(residuos_comunes) < tolerancia:
        raise ValueError("Muy pocos residuos comunes para alinear.")

    atomos_referencia = [ca_ref[residuo] for residuo in sorted(residuos_comunes)]
    atomos_comunes = [ca_otro[residuo] for residuo in sorted(residuos_comunes)]

    if len(atomos_referencia) != len(atomos_comunes):
        raise ValueError("Las listas de C-Alfa no tienen el mismo largo.")

    si = Superimposer()
    si.set_atoms(atomos_referencia, atomos_comunes)
    si.apply(
        [atom for atom in estructuraOtra.get_atoms()]
    )  # Aplica la transformación sobre todos los átomos
    return si.rms


# Convertir estructura Biopython a string PDB para py3Dmol
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
        cadenaID = "A"
    io_pdb.save(string_io, select=SeleccionarCadena(cadenaID))
    return string_io.getvalue()
