import requests


# Busca una proteina en UniProt por ID
# Entrada = texto
# Salida = data json
def buscar_id_uniprot(id):

    print("Buscando el ID en UniProt")

    base_url = f"https://rest.uniprot.org/uniprotkb/{id}"

    params = {
        "format": "json",
        "fields": "accession,id,sequence,protein_name,organism_name,organism_id,length",
    }

    try:

        response = requests.get(base_url, params=params, timeout=30)
        response.raise_for_status()

        return response.json()
    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 404:
            return {"error": f"No se encontró el ID '{id}' en UniProt"}
        else:
            return {"error": f"Error HTTP {e.response.status_code}: {e.response.text}"}
    except requests.exceptions.RequestException as error:
        return {"error": f"Error de conexión: {str(error)}"}
    except Exception as error:
        return {"error": f"Error inesperado: {str(error)}"}


# Busca una proteina en UniProt por query
# Entrada = texto
# Salida = data
def buscar_uniprot(query):

    print("Buscando en UniProt")

    base_url = f"https://rest.uniprot.org/uniprotkb/search"

    params = {
        "query": query,
        "format": "json",
        "fields": "accession,id,sequence,protein_name,organism_name,organism_id,length",
        "size": 10,
    }

    try:

        response = requests.get(base_url, params=params, timeout=30)
        response.raise_for_status()

        return response.json()

    except requests.exceptions.HTTPError as e:
        return {"error": f"Error HTTP {e.response.status_code}: {e.response.text}"}
    except requests.exceptions.RequestException as error:
        return {"error": f"Error de conexión: {str(error)}"}
    except Exception as error:
        return {"error": f"Error inesperado: {str(error)}"}


# Busca estructuras PDB asociadas a un accession de UniProt
# Entrada = accession de UniProt
# Salida = lista de diccionarios con información de PDB
def buscar_pdb_uniprot(accession):

    print(f"Buscando estructuras PDB para accession: {accession}")

    # URL de la API de UniProt con extensión .json
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"

    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        data = response.json()

        # Extraer información detallada de PDB
        pdb_info = []
        for xref in data.get("uniProtKBCrossReferences", []):
            if xref.get("database") == "PDB":
                pdb_id = xref.get("id")
                if pdb_id:
                    # Extraer información de las propiedades
                    method = "N/A"
                    resolution = "N/A"
                    chains = "N/A"

                    for prop in xref.get("properties", []):
                        if prop.get("key") == "Method":
                            method = prop.get("value", "N/A")
                        elif prop.get("key") == "Resolution":
                            resolution = prop.get("value", "N/A")
                        elif prop.get("key") == "Chains":
                            chains = prop.get("value", "N/A")

                    pdb_info.append(
                        {
                            "identifier": pdb_id,
                            "method": method,
                            "resolution": resolution,
                            "chain": chains,
                        }
                    )

        return pdb_info

    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 404:
            print(f"Error: No se encontró el accession '{accession}' en UniProt")
        else:
            print(f"Error HTTP {e.response.status_code}: {e.response.text}")
        return []
    except requests.exceptions.RequestException as error:
        print(f"Error al buscar PDB: {str(error)}")
        return []


# Descarga features de una proteína desde UniProt usando accession y formato
# Entrada = accession de UniProt, formato deseado
# Salida = archivo descargado o mensaje de error
def buscar_features_uniprot(accession, formato="json"):

    print(f"Iniciando búsqueda de features para accession: {accession}")

    # Validar formato
    formatos_permitidos = ["json", "txt", "xml", "gff"]
    if formato not in formatos_permitidos:
        return f"Error: Formato '{formato}' no válido. Formatos permitidos: {formatos_permitidos}"

    # URL de la API de UniProt
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.{formato}"

    try:
        print(f"Realizando consulta a: {url}")
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        return response.content

    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 404:
            return f"Error: No se encontró el accession '{accession}' en UniProt"
        else:
            return f"Error HTTP {e.response.status_code}: {e.response.text}"
    except requests.exceptions.RequestException as error:
        return f"Error de conexión: {str(error)}"
    except Exception as error:
        return f"Error inesperado: {str(error)}"


# Busca los accession asociados a un PDB en UniProt
# Entrada = PDB ID, cadena ID
# Salida = lista de accessions
def buscar_pdb_accessions(pdb_id: str, chain_id: str, timeout: int = 20) -> set[str]:

    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
    try:
        resp = requests.get(url, timeout=timeout)
        resp.raise_for_status()
        data = resp.json()
    except Exception as exc:
        print(f"Advertencia: no se pudo contactar con PDBe-SIFTS ({exc}).")
        return set()

    acceso = set()
    for record in data.get(pdb_id.lower(), {}).get("UniProt", {}).values():
        for seg in record.get("mappings", []):
            if seg.get("chain_id").strip() == chain_id.strip():
                acceso.add(record["identifier"])
    return acceso
