import re

import pandas as pd

from data.fetch_uniprot import buscar_pdb_uniprot


# Verifica si el accession es de UniProt
def es_accession_uniprot(accession):
    return bool(
        re.match(
            r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^[A-Z][0-9][A-Z0-9]{3}[0-9]$|^A0A[A-Z0-9]{7}$",
            accession,
        )
    )


# Formatea los resultados de PDB en una tabla usando pandas
# Entrada = resultados de PDB
# Salida = tabla formateada con información de las estructuras PDB
def formatear_resultados_pdb(resultados):

    if isinstance(resultados, str):
        return resultados  # Es un mensaje de error

    if not resultados:
        return "No se encontraron estructuras PDB"

    # Crear DataFrame con la nueva estructura
    df = pd.DataFrame(resultados)

    # Renombrar columnas para mejor presentación
    df = df.rename(
        columns={
            "identifier": "Identifier",
            "method": "Method",
            "resolution": "Resolution",
            "chain": "Chain/Position",
        }
    )

    # Configurar pandas para mejor visualización
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", None)
    pd.set_option("display.max_colwidth", 50)

    return df.to_string(index=False, justify="left")


# Busca estructuras PDB asociadas a un accession
# Entrada = accession de UniProt
# Salida = tabla con información detallada de PDBs encontrados
def lista_pdb(accession):

    print(f"Buscando estructuras PDB para: {accession}")

    # Verificar si es un accession de UniProt
    if not es_accession_uniprot(accession):
        return f"Error: '{accession}' no es un accession válido de UniProt"

    # Buscar PDBs en UniProt
    pdb_info = buscar_pdb_uniprot(accession)

    if not pdb_info:
        return f"No se encontraron estructuras PDB para el accession {accession}"

    # Usar pandas para formatear y mostrar los resultados
    tabla_formateada = formatear_resultados_pdb(pdb_info)

    print(f"\nEncontradas {len(pdb_info)} estructuras PDB:")
    print("=" * 120)
    print(tabla_formateada)
    print("=" * 120)
