import re  # re es una libreria para evaluar expreciones regulares

import pandas as pd

from data.fetch_ncbi import buscar_acn_ncbi
from data.fetch_uniprot import buscar_id_uniprot, buscar_uniprot


# identifica si el texto es un id de uniprot utilizando re y las expresiones de los id de uniprot
# Entrada = texto
# Salida = booleano
def es_id_uniprot(texto):
    return bool(
        re.fullmatch(
            r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^A0A[A-Z0-9]{7}$",
            texto,
        )
    )


# identifica si el texto es un id de ncbi utilizando re y las expresiones de los id de proteinas de ncbi (NP_, XP_, YP_, etc...)
# Entrada = texto
# Salida = booleano
def es_id_ncbi(texto):
    return bool(re.fullmatch(r"^(NP|XP|YP|WP|ZP|AP)_\d{6,9}(\.\d+)?$", texto))


# Formatea los resultados de UniProt en una tabla usando pandas
# Entrada = datos de UniProt
# Salida = tabla con información de la proteína
def formatear_resultados_uniprot(data):

    if "results" not in data or not data["results"]:
        return "No se encontraron resultados"

    # Preparar datos
    datos = []

    for i, result in enumerate(data["results"][:10], 1):
        accession = result.get("primaryAccession", "N/A")
        protein_id = result.get("uniProtkbId", "N/A")

        # Extraer nombre de proteína
        protein_name = "N/A"
        protein_desc = result.get("proteinDescription", {})
        if "recommendedName" in protein_desc:
            protein_name = (
                protein_desc["recommendedName"].get("fullName", {}).get("value", "N/A")
            )
        elif "submissionNames" in protein_desc and protein_desc["submissionNames"]:
            protein_name = (
                protein_desc["submissionNames"][0]
                .get("fullName", {})
                .get("value", "N/A")
            )

        # Extraer organismo
        organism = result.get("organism", {}).get("scientificName", "N/A")

        # Extraer longitud
        length = result.get("sequence", {}).get("length", "N/A")

        datos.append(
            {
                "#": i,
                "Accession": accession,
                "ID": protein_id,
                "Nombre de Proteína": protein_name,
                "Organismo": organism,
                "Longitud": length,
            }
        )

    # Crear DataFrame y formatear
    df = pd.DataFrame(datos)

    # Configurar pandas para mejor visualización
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", None)
    pd.set_option("display.max_colwidth", 50)

    return df.to_string(index=False, justify="left")


# Formatea los resultados de NCBI en una tabla usando pandas
# Entrada = datos de NCBI
# Salida = tabla con información de la proteína
def formatear_resultados_ncbi(data):

    if not isinstance(data, dict):
        return f"Error: Formato de datos inesperado: {type(data)}"

    # Verificar si hay error
    if "error" in data:
        return data["error"]

    # Buscar la estructura correcta de NCBI
    if "result" not in data:
        return "Formato de respuesta de NCBI no reconocido"

    result = data["result"]
    if not isinstance(result, dict):
        return "Estructura de resultado de NCBI inválida"

    # Preparar datos
    datos = []

    # Buscar el UID de la proteína (excluir 'uids' que es una lista)
    protein_uid = None
    for key in result.keys():
        if key != "uids" and key.isdigit():
            protein_uid = key
            break

    if not protein_uid:
        return "No se encontró información de proteína válida"

    protein_data = result[protein_uid]
    if not isinstance(protein_data, dict):
        return "Datos de proteína inválidos"

    # Extraer información de la proteína
    accession = protein_data.get("accessionversion", "N/A")
    title = protein_data.get("title", "N/A")
    organism = protein_data.get("organism", "N/A")
    length = protein_data.get("slen", "N/A")

    datos.append(
        {
            "#": 1,
            "UID": protein_uid,
            "Accession": accession,
            "Título": title,
            "Organismo": organism,
            "Longitud": length,
        }
    )

    # Crear DataFrame y formatear
    df = pd.DataFrame(datos)

    # Configurar pandas para mejor visualización
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", None)
    pd.set_option("display.max_colwidth", 50)

    return df.to_string(index=False, justify="left")


# Busca información de una proteína por su ID
# Entrada = ID de UniProt, NCBI o texto descriptivo
# Salida = tabla con información de la proteína
def buscar(prompt):
    if es_id_uniprot(prompt):
        resultado = buscar_id_uniprot(prompt)
        if isinstance(resultado, dict) and "error" not in resultado:
            # Para IDs específicos de UniProt, crear un formato similar a búsquedas múltiples
            data = {"results": [resultado]}
            return formatear_resultados_uniprot(data)
        else:
            return resultado
    elif es_id_ncbi(prompt):
        resultado = buscar_acn_ncbi(prompt)
        if isinstance(resultado, dict) and "error" not in resultado:
            return formatear_resultados_ncbi(resultado)
        else:
            return resultado
    else:
        resultado = buscar_uniprot(prompt)
        if isinstance(resultado, dict) and "error" not in resultado:
            return formatear_resultados_uniprot(resultado)
        else:
            return resultado
