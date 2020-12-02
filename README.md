# seqBio: An R package

El paquete seqBio contiene dos objetos: las funciones seq_alignment() y generateSeqsWithMultinomialModel(). La primera función es útil para estudiar la significación del alineamiento de secuencias y para ello, reliza el cálculo de la probabilidad de que la puntuación del alineamiento ha sido obtenida por azar (aplicando la estrategia de shuffling). La función generateSeqsWithMultinomialModel() es una función auxiliar dentro de seq_alignment(), necesaria para generar un numero determinado de secuencias aleatorias.

## Getting Started

Estas instrucciones proporcionan los pasos a seguir para instalarse el paquete en su propia máquina local de manera que podrá ejecutar las funciones, ir a la ayuda de R para ver su utilización y visitar un informe dinamico (vignette) para obtener una explicación más larga y detallada sobre el paquete y sus funciones. 

### Prerequisitos

El paquete necesita una instalación báscia del software R. Además, será necesaria también la instalación de varios paquetes: 

```
install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("Biostrings")
install.packages("seqinr")
install.packages("extRemes")
install.packages("gplots")
```

### Instalación

Para instalar este paquete, utiliza la siguiente linia de codigo:

```
devtools::install_github("laurajuliamelis/Bioinformatics_Task3", subdir="seqBio", build_vignettes=TRUE)
```

Luego, necesitarás cargar el paquete:

```
library(seqBio)
```


### Ayuda 

Para obtener más información acerca del paquete:
```
?seqBio
```

Para buscar en la ayuda de la función seq_alignment(), utilizar el siguiente código:
```
?seq_alignment
```
En la ayuda se podrá encontrar una breve descripción del paquete, los argumentos y sus definiciones, los valores de salida, las referencias, y varios ejemplos.  

Para obtener más información de la función generateSeqsWithMultinomialModel(), visitar la ayuda utilizando el siguiente código:

```
?generateSeqsWithMultinomialModel
```
Esta es una función necesaria para la ejecucion de seq_alignment().

### Vignette (informe dinámico)

Adicionalmente, se puede acceder al informe dinámico del paquete para obtener una explicación más detallada y extensiva de la función seq_align().

Para acceder a esta viñeta, se debe utilizar el siguiente código:

```
browseVignettes('seqBio')
```

A continuación, pulsar el boton "**HTML**".

## Licencia

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

