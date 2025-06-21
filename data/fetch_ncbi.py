import json

from Bio import Entrez, SeqIO


# Busca una proteina en NCBI por ID
# Entrada = texto
# Salida = data json
def buscar_acn_ncbi(accession, email="tucorreo@example.com"):

    Entrez.email = email
    print("Buscando el ID en NCBI")

    try:
        # Paso 1: Buscar el UID del accession
        handle = Entrez.esearch(db="protein", term=accession, retmode="json")
        raw_search = handle.read()
        search_data = json.loads(raw_search)
        handle.close()

        id_list = search_data["esearchresult"]["idlist"]
        if not id_list:
            return {"error": f"No se encontró el accession '{accession}' en NCBI"}

        uid = id_list[0]
        print(f"UID encontrado: {uid}")

        # Paso 2: Obtener resumen en JSON
        handle = Entrez.esummary(db="protein", id=uid, retmode="json")
        raw_summary = handle.read()
        summary_data = json.loads(raw_summary)
        handle.close()

        return summary_data

    except Exception as error:
        print(f"Error en la búsqueda de NCBI: {str(error)}")
        return {"error": f"Error en la búsqueda de NCBI: {str(error)}"}
