import requests


# Busca una proteina en UniProt por ID
# Entrada = texto
# Salida = data json
def buscar_id_uniprot(id):

    print("Buscando el ID en UniProt")

    base_url = f"https://rest.uniprot.org/uniprotkb/{id}"

    params = {
        "format": "json",
        "fields": "accession,id,sequence,protein_name,organism_name",
    }

    try:

        response = requests.get(base_url, params=params, timeout=30)
        response.raise_for_status()

        return response.json()
    except Exception as error:
        print("Error en la busqueda:", str(error))


# Busca una proteina en UniProt por query
# Entrada = texto
# Salida = data
def buscar_uniprot(query):

    print("Buscando en UniProt")

    base_url = f"https://rest.uniprot.org/uniprotkb/search"

    params = {
        "query": query,
        "format": "json",
        "fields": "accession,id,sequence,protein_name,organism_name",
        "size": 10,
    }

    try:

        response = requests.get(base_url, params=params, timeout=30)
        response.raise_for_status()

        return response.json()

    except Exception as error:
        print("Error en la busqueda:", str(error))
