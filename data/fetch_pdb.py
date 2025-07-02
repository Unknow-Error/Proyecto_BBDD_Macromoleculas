import tempfile

import requests


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
