# Análisis de Expresión Diferencial y Propagación en Redes de Genes de *Arabidopsis thaliana*

Este repositorio contiene los pasos para el análisis de expresión diferencial y propagación en redes de genes de *Arabidopsis thaliana* utilizando datos del proyecto [PRJNA156363](https://www.ebi.ac.uk/ena/browser/view/PRJNA156363?show=reads). 

Los datos fueron [preprocesados y mapeados con *Kallisto*](https://usegalaxy.eu/published/history?id=1e4ffa727de55dbc), y aquí se lleva a cabo un análisis funcional a partir de los genes sobreexpresados, seguido de la descarga y análisis de una red de interacciones de proteínas utilizando el algoritmo DIAMOnD para la propagación en red.

> [!note]  
> Este proyecto está orientado al análisis funcional y la propagación de información en redes biológicas, especificamente para proyectos de biología en sistemas (asignatura de la Universidad De Málaga).

## Datos utilizados

- Se realiza un análisis funcional de los genes sobreexpresados, filtrando aquellos con `logFC > 0` y `adj.P.Val < 0.05`. Dichos genes puedes encontrarlo en la pestaña dataset con mayor información de su obtención.

> [!tip]  
> Los valores de `logFC` y `adj.P.Val` se utilizan en análisis de expresión diferencial para priorizar genes relevantes.

## Proceso General

1. **Análisis de expresión diferencial**:
    - Se leen los datos de expresión diferencial desde archivos `.tsv`.
    - Se filtran los genes sobreexpresados.
    
2. **Análisis funcional**:
    - Se utiliza la API de STRINGdb para obtener un enriquecimiento funcional de los genes sobreexpresados.
    - Se filtran los procesos significativos con un `FDR < 0.001`.
    - Los resultados se guardan en archivos `.csv`.

> [!important]  
> La utilización de APIs como STRINGdb permite un análisis funcional pero recurrimos a la API lo cual puede disminuir la velocidad del script.

3. **Descarga de red de interacción proteica**:
    - Se descargan las interacciones de proteínas de STRINGdb para *Arabidopsis thaliana*.
    - Se filtra la red utilizando un `combined_score > 400` para reducir el ruido y estudiar relaciones más directas (para un trabajo científico se recomienda el uso de por lo menos 400, si queremos simplificar los datos 700 o 900 pueden ser buenas opciones.

4. **Propagación en red con DIAMOnD**:
    - Se implementa el algoritmo [DIAMOnD](http://diamond.barabasilab.com/) para identificar proteínas adicionales relevantes en la red a partir de los genes semilla sobreexpresados.

> [!warning]  
> El algoritmo DIAMOnD no está optimizado para tiempos de ejecución cortos, y puede tardar en compilar y ejecutarse, especialmente en redes grandes. Considera usar versiones más parecidas al [código original.](http://diamond.barabasilab.com/)

5. **Nuevo Análisis funcional**:
   - Se repite el paso 2 pero usando los datos obtenidos con la propagación de red.

> [!caution]  
> Recuerda que la propagación en redes puede añadir proteínas que no necesariamente están relacionadas de manera directa con los genes semilla originales, por lo que es crucial interpretar los resultados con cuidado.
