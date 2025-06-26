import json
import os

from data.fetch_uniprot import buscar_features_uniprot


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
            # Decodificar contenido y parsear JSON
            contenido_texto = contenido.decode("utf-8")
            datos_json = json.loads(contenido_texto)

            # Guardar con formato legible
            with open(archivo_salida, "w", encoding="utf-8") as f:
                json.dump(datos_json, f, indent=2, ensure_ascii=False)
        else:
            # Para otros formatos, guardar como bytes
            with open(archivo_salida, "wb") as f:
                f.write(contenido)

        # Verificar que el archivo se guardó correctamente
        if os.path.exists(archivo_salida):
            # Verificar que el archivo no esté vacío
            if os.path.getsize(archivo_salida) > 0:
                print(f"Archivo guardado exitosamente: {archivo_salida}")
            else:
                print(f"Error: El archivo se creó pero está vacío: {archivo_salida}")
        else:
            print(f"Error: No se pudo crear el archivo: {archivo_salida}")

        return None
    except Exception as error:
        print(f"Error al guardar archivo: {str(error)}")
        return None


# Descarga features de una proteína por accession
# Entrada = accession de UniProt, formato
# Salida = mensaje de estado de la descarga
def descargar_features(accession, formato):
    contenido = buscar_features_uniprot(accession, formato)
    return guardar_features_archivo(contenido, accession, formato)
