# Proyecto_BBDD_Macromoleculas

Trabajo final de la asignatura BBDD de macromoleculas de la Universidad Nacional de Quilmes.

Consiste en una aplicación de consola (terminal de comandos) que permite buscar proteínas usando el identificador de Uniprot o el número de accession de NCBI protein para:

1- Reconocer el tipo de indetificador (input) del usuario : Incluye también un input de texto para Uniprot.
2- Realizar la búsqueda en la BBDD correpondiente y descargar de la proteína sus Features (Tabla de Características) utilizando la API REST.
3- Buscar los identificador PDB de la proteína. 
4- Permitir al usuario seleccionar dos PDB de la proteína.
5- Realizar un análisis de alineamiento estructural computando el RMSD en una ventana de 5 residuos (Esto puede ser un parámetro libre seteable por el usuario) y devolver como output un gráfico de RMSD vs. posición centroide-residuo de cada ventana.
6- Emplear un visualizador para mostrar las Features y el alineamiento residual de las proteínas.
