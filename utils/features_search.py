import json
import os

from data.fetch_uniprot import buscar_features_uniprot
from utils.pdb_search import es_accession_uniprot, map_ncbi_to_uni

# Guarda features en un archivo
# Entrada = contenido, accession, formato
# Salida = archivo guardado o mensaje de error
def guardar_features_archivo(contenido, accession, formato):
    # Crear directorio si no existe
    directorio = "features_files"
    if not os.path.exists(directorio):
        os.makedirs(directorio)

    # Nombre del archivo
    nombre_archivo = f"{accession}_features.{formato}"
    archivo_salida = os.path.join(directorio, nombre_archivo)

    try:
        # Para JSON, formatear con indentación
        if formato == "json":
            contenido_texto = contenido.decode("utf-8")
            datos_json = json.loads(contenido_texto)
            with open(archivo_salida, "w", encoding="utf-8") as f:
                json.dump(datos_json, f, indent=2, ensure_ascii=False)
        else:
            with open(archivo_salida, "wb") as f:
                f.write(contenido)

        # Verificación y retorno
        if os.path.exists(archivo_salida):
            if os.path.getsize(archivo_salida) > 0:
                print(f"Archivo guardado exitosamente: {archivo_salida}")
                return archivo_salida  
            else:
                print(f"Error: El archivo se creó pero está vacío: {archivo_salida}")
        else:
            print(f"Error: No se pudo crear el archivo: {archivo_salida}")

    except Exception as error:
        print(f"Error al guardar archivo: {str(error)}")

    return None  # Solo en caso de error o fallo

# Descarga features de una proteína por accession
# Entrada = accession de UniProt, formato
# Salida = mensaje de estado de la descarga
def descargar_features(accession, formato):
    archivos_guardados = []

    # Accession es de UniProt directamente
    if es_accession_uniprot(accession):
        try:
            contenido = buscar_features_uniprot(accession, formato)
            ruta = guardar_features_archivo(contenido, accession, formato)
            if ruta:
                archivos_guardados.append(ruta)
        except Exception as e:
            print(f"Error con '{accession}': {e}")
    
    # Es un ID de NCBI, mapear a UniProt
    else:
        try:
            uniprot_ids = map_ncbi_to_uni(accession)
            if not uniprot_ids:
                return f"El ID '{accession}' no se pudo mapear a UniProt."
            print(f"IDs UniProt mapeados desde '{accession}': {uniprot_ids}")
        except Exception as e:
            return f"Error durante el mapeo: {e}"

        for uid in uniprot_ids:
            try:
                print(f"Descargando features para UniProt ID: {uid}")
                contenido = buscar_features_uniprot(uid, formato)
                ruta = guardar_features_archivo(contenido, uid, formato)
                if ruta:
                    archivos_guardados.append(ruta)
            except Exception as e:
                print(f"No se pudieron obtener features para '{uid}': {e}")
    
    if not archivos_guardados:
        return "No se guardó ningún archivo."
    
    print("\nArchivos guardados:")
    for archivo in archivos_guardados:
        print(f"  - {archivo}")
    
    return None  # No es necesario devolver texto si ya imprime todo